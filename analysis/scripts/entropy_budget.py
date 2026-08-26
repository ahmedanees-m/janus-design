"""Gate 5: how much sequence freedom each backbone carries, and where.

Secondary structure comes from biotite's P-SEA implementation rather than DSSP,
which needs a binary that is not available in this container. The substitution is
recorded rather than left implicit. Relative accessibility uses the Shrake and
Rupley construction with the theoretical maxima of Tien et al. 2013.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np

from janus.lattice import amino_acid_shell
from janus.objectives.mpnn import load_unconditional

# Tien et al. 2013, theoretical maximum accessible surface area, square angstrom.
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
TOPOLOGY = re.compile(r"^(HHH|EHEE|EEHEE|HEEH)(?:_[A-Z]{2})?_rd\d+")
DELTAS = (0.5, 1.0, 2.0, 3.0)


def structure_features(path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    sse = struc.annotate_sse(atoms)
    area = struc.sasa(atoms, vdw_radii="Single")
    per_residue = struc.apply_residue_wise(atoms, area, np.nansum)
    names = struc.get_residues(atoms)[1]
    residues = [THREE_TO_ONE.get(n, "X") for n in names]
    relative = np.array([
        per_residue[i] / MAX_ASA[r] if r in MAX_ASA else np.nan
        for i, r in enumerate(residues)
    ])
    return sse, relative, residues


def family(name):
    match = TOPOLOGY.match(name)
    if match:
        return match.group(1)
    if name.startswith(("EA:", "GG:")) or "TrROS" in name:
        return "hallucination"
    if PDBLIKE.match(name):
        return "natural"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    skipped = []
    files = sorted(Path(args.marginals).glob("*.npz"))
    for number, path in enumerate(files, start=1):
        name = path.stem
        marginals = load_unconditional(path)
        entropy = -(np.exp(marginals) * marginals).sum(axis=1)
        shells = {d: [len(s) for s in amino_acid_shell(marginals, d)] for d in DELTAS}

        try:
            sse, relative, residues = structure_features(Path(args.backbones) / f"{name}.pdb")
        except Exception as exc:
            skipped.append((name, type(exc).__name__))
            continue
        if len(sse) != len(marginals):
            skipped.append((name, f"length {len(sse)} against {len(marginals)}"))
            continue

        length = len(marginals)
        rows.append({
            "backbone": name,
            "origin": "natural" if PDBLIKE.match(name) else "designed",
            "family": family(name),
            "length": length,
            "entropy": entropy.tolist(),
            "sse": "".join(sse),
            "relative_sasa": [None if np.isnan(v) else round(float(v), 4) for v in relative],
            "terminal_distance": [min(i, length - 1 - i) for i in range(length)],
            **{f"shell_{d}": shells[d] for d in DELTAS},
        })
        if number % 100 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    print(f"\nwrote {args.out}: {len(rows)} backbones, {len(skipped)} skipped")
    for name, why in skipped[:10]:
        print(f"    skipped {name}: {why}")

    summarise(rows)


def summarise(rows):
    for origin in ("designed", "natural"):
        subset = [r for r in rows if r["origin"] == origin]
        if not subset:
            continue
        entropy = np.concatenate([r["entropy"] for r in subset])
        print(f"\n=== {origin} ({len(subset)} backbones) ===")
        lengths = [r["length"] for r in subset]
        print(f"  length: min {min(lengths)} median {int(np.median(lengths))} max {max(lengths)}")
        print(f"  per-position entropy: mean {entropy.mean():.3f} median {np.median(entropy):.3f}")
        print(f"  fixed (H below 0.1): {100 * (entropy < 0.1).mean():.1f}%   "
              f"free (H above 2.0): {100 * (entropy > 2.0).mean():.1f}%")
        for d in DELTAS:
            sizes = np.concatenate([r[f"shell_{d}"] for r in subset])
            print(f"  mean shell at delta={d}: {sizes.mean():.2f}")
        total = [np.sum(r["entropy"]) for r in subset]
        print(f"  total freedom per design: median {np.median(total):.0f} nats")

    print("\n=== freedom against structure, all backbones ===")
    entropy = np.concatenate([r["entropy"] for r in rows])
    sse = np.concatenate([list(r["sse"]) for r in rows])
    sasa = np.concatenate([[np.nan if v is None else v for v in r["relative_sasa"]] for r in rows])
    terminal = np.concatenate([r["terminal_distance"] for r in rows])

    for code, label in [("a", "helix"), ("b", "strand"), ("c", "coil")]:
        mask = sse == code
        if mask.sum():
            print(f"  {label:<7} n={mask.sum():>6}  mean entropy {entropy[mask].mean():.3f}")

    ok = ~np.isnan(sasa)
    buried = ok & (sasa < 0.15)
    exposed = ok & (sasa > 0.40)
    print(f"  buried  (rel SASA below 0.15) n={buried.sum():>6}  mean entropy {entropy[buried].mean():.3f}")
    print(f"  exposed (rel SASA above 0.40) n={exposed.sum():>6}  mean entropy {entropy[exposed].mean():.3f}")

    from scipy.stats import spearmanr
    print(f"  Spearman entropy against relative SASA: {spearmanr(entropy[ok], sasa[ok]).statistic:.3f}")
    print(f"  Spearman entropy against distance from a terminus: "
          f"{spearmanr(entropy, terminal).statistic:.3f}")

    near = terminal < 3
    print(f"  within 3 of a terminus: mean entropy {entropy[near].mean():.3f}   "
          f"interior: {entropy[~near].mean():.3f}")


if __name__ == "__main__":
    main()
