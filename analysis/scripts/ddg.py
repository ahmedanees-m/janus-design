"""Do the substitutions JANUS proposes cost measured stability?

MegaScale reports folding free energy for single mutants of each parent, so a
substitution's cost is the difference between the mutant and the wild type,
negative when destabilising. Mutant names carry a position, and the numbering is
verified against the translated construct before any of it is used rather than
assumed.

Comparators are drawn from the same shell so the contrast isolates which residue
was chosen and not how many substitutions were made.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from janus import Weights, design, hosts
from janus.genetic_code import AA_ALPHABET, translate
from janus.lattice import amino_acid_shell
from janus.objectives.mpnn import load_unconditional

SINGLE = re.compile(r"^(?P<parent>.+?)\.pdb_(?P<wt>[A-Z])(?P<pos>\d+)(?P<mut>[A-Z])$")
BARE = re.compile(r"^(?P<parent>[A-Za-z0-9_:.\-]+?)_(?P<wt>[A-Z])(?P<pos>\d+)(?P<mut>[A-Z])$")
AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}


def parse_single(name):
    for pattern in (SINGLE, BARE):
        match = pattern.match(name)
        if match:
            return (match.group("parent"), match.group("wt"),
                    int(match.group("pos")), match.group("mut"))
    return None


def load_measurements(parquet_dir):
    """Wild-type free energies and single-mutant effects, keyed by parent."""
    wild = {}
    constructs = {}
    singles = defaultdict(dict)

    for path in sorted(Path(parquet_dir).glob("Lib*_K50dG.parquet")):
        table = pq.read_table(path, columns=["name", "dna_seq", "deltaG"])
        names = table.column("name").to_pylist()
        dna = table.column("dna_seq").to_pylist()
        energies = table.column("deltaG").to_pylist()

        for name, seq, energy in zip(names, dna, energies, strict=True):
            if not name or energy is None or not np.isfinite(energy):
                continue
            base = name.partition(".pdb")[0] if ".pdb" in name else name
            if name in (base, base + ".pdb"):
                wild[base] = float(energy)
                if seq:
                    constructs[base] = translate(seq)
                continue
            parsed = parse_single(name)
            if parsed:
                parent, wt, position, mut = parsed
                singles[parent][position, wt, mut] = float(energy)
    return wild, constructs, singles


def infer_offset(construct, mutations):
    """Mutant-name positions are offset against the construct by a per-parent
    amount. Recover it by agreement rather than assuming a convention."""
    best, score = None, 0
    for offset in range(-5, 30):
        agree = sum(
            1 for position, wt, _ in mutations
            if 0 <= position - 1 + offset < len(construct)
            and construct[position - 1 + offset] == wt
        )
        if agree > score:
            best, score = offset, agree
    return best, (score / len(mutations) if mutations else 0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--agreement", type=float, default=0.95)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb

    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
        "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V",
    }

    def backbone_sequence(path):
        atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
        atoms = atoms[struc.filter_amino_acids(atoms)]
        return "".join(three_to_one.get(n, "X") for n in struc.get_residues(atoms)[1])

    host = hosts.load(args.host)
    rng = np.random.default_rng(args.seed)
    wild, constructs, singles = load_measurements(args.parquet)
    print(f"{len(wild)} parents with a wild-type measurement, "
          f"{len(singles)} with single mutants")

    rows = []
    examined = aligned = 0
    reasons = defaultdict(int)

    for path in sorted(Path(args.marginals).glob("*.npz")):
        name = path.stem
        if name not in singles or name not in wild or name not in constructs:
            reasons["no measurement"] += 1
            continue
        examined += 1
        construct = constructs[name]
        mutations = list(singles[name])

        offset, agreement = infer_offset(construct, mutations)
        if offset is None or agreement < args.agreement:
            reasons["numbering unresolved"] += 1
            continue

        native = backbone_sequence(Path(args.backbones) / f"{name}.pdb")
        start = construct.find(native)
        if start < 0:
            reasons["backbone not in construct"] += 1
            continue

        marginals = load_unconditional(path)
        if len(marginals) != len(native):
            reasons["length mismatch"] += 1
            continue
        aligned += 1

        shells = amino_acid_shell(marginals, args.delta)
        proposal = design(marginals, host, weights=Weights(mpnn=1.0, cai=0.5, cpb=0.3),
                          delta=args.delta, k=1)[0].protein

        for index in range(len(native)):
            if native[index] == "X":
                continue
            position = start + index - offset + 1
            arms = {
                "janus": proposal[index],
                "shell_random": shells[index][int(rng.integers(len(shells[index])))],
                "mpnn_sample": AA_ALPHABET[int(rng.choice(
                    len(AA_ALPHABET), p=_tempered(marginals[index], 0.1)))],
            }
            for arm, mutant in arms.items():
                if mutant == native[index]:
                    continue
                energy = singles[name].get((position, native[index], mutant))
                if energy is None:
                    continue
                rows.append({
                    "backbone": name, "position": index, "arm": arm,
                    "wild_type": native[index], "mutant": mutant,
                    "ddg": energy - wild[name],
                })

    print(f"\n{examined} backbones had measurements, {aligned} aligned cleanly")
    for reason, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"    {reason}: {n}")
    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows, aligned)


def _tempered(row, temperature):
    logits = row / temperature
    logits -= logits.max()
    weights = np.exp(logits)
    return weights / weights.sum()


def summarise(rows, parents):
    from scipy.stats import mannwhitneyu

    print(f"\n=== measured stability cost, {len(rows)} covered substitutions ===")
    arms = ["janus", "shell_random", "mpnn_sample"]
    values = {a: np.array([r["ddg"] for r in rows if r["arm"] == a]) for a in arms}
    covered = {a: len({(r["backbone"], r["position"]) for r in rows if r["arm"] == a})
               for a in arms}

    print(f"{'arm':<14}{'n':>7}{'positions':>11}{'median':>9}{'mean':>9}"
          f"{'below -0.5':>12}{'below -1.0':>12}")
    for arm in arms:
        v = values[arm]
        if len(v) == 0:
            continue
        print(f"{arm:<14}{len(v):>7}{covered[arm]:>11}{np.median(v):>9.3f}{v.mean():>9.3f}"
              f"{100 * np.mean(v < -0.5):>11.1f}%{100 * np.mean(v < -1.0):>11.1f}%")

    if len(values["janus"]) and len(values["shell_random"]):
        stat = mannwhitneyu(values["janus"], values["shell_random"], alternative="greater")
        print(f"\n  one-sided test that JANUS substitutions are not more destabilising")
        print(f"  than random draws from the same shell: p = {stat.pvalue:.4g}")
    if len(values["janus"]) and len(values["mpnn_sample"]):
        stat = mannwhitneyu(values["janus"], values["mpnn_sample"], alternative="greater")
        print(f"  against ProteinMPNN sampling at T=0.1: p = {stat.pvalue:.4g}")

    print(f"\n  backbones contributing: "
          f"{len({r['backbone'] for r in rows})} of {parents} with verified numbering")


if __name__ == "__main__":
    main()
