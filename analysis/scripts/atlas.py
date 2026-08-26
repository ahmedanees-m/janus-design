"""The liability atlas: designs against four control sets.

Scrambled and random controls keep the accessibility profile of the design they
replace, so accessibility-weighted features are comparable across every group:
same structure, different sequence. Natural controls carry their own structures.

There are two natural arms. The first is the natural domains in the same
MegaScale release, which is the like-for-like comparison but is mostly excised:
only 15 of 104 entries are whole UniProt proteins. The second is whole natural
proteins of the same size assembled separately, which is the comparison that can
carry a claim about proteostatic liability, since a motif outside an excision is
not absent, only unseen.

Gene-layer features are computed on coding sequences generated the same way for
every group, by taking the highest relative-adaptiveness codon for each residue.
Differences therefore reflect the residue sequence rather than the codon policy.
The MegaScale as-ordered sequences are not used; per the data audit they are
oligo-pool constrained and cannot carry a descriptive claim.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np

from features import gene_features, protein_features
from janus import hosts
from janus.genetic_code import SYNONYMOUS
from janus.objectives.proteostasis import load_classes

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


def read_backbone(path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    sse = "".join(struc.annotate_sse(atoms))
    area = struc.apply_residue_wise(atoms, struc.sasa(atoms, vdw_radii="Single"), np.nansum)
    residues = [THREE_TO_ONE.get(n, "X") for n in struc.get_residues(atoms)[1]]
    rsa = [
        None if r not in MAX_ASA or np.isnan(area[i]) else float(area[i] / MAX_ASA[r])
        for i, r in enumerate(residues)
    ]
    return "".join(residues), rsa, sse


def family(name):
    match = TOPOLOGY.match(name)
    if match:
        return match.group(1)
    if name.startswith(("EA:", "GG:")) or "TrROS" in name:
        return "hallucination"
    return "other"


def best_codon_cds(protein, host):
    adaptiveness = host.relative_adaptiveness
    best = {
        aa: max(family_codons, key=lambda c: adaptiveness.get(c, 0.0))
        for aa, family_codons in SYNONYMOUS.items()
    }
    return "".join(best[r] for r in protein if r in best)


COUNT_FEATURES = (
    "c_degron", "n_degron", "internal_degron", "all_degron",
    "protease_site", "targeting", "modification",
)


def annotate(name, group, protein, rsa, sse, classes, host, topology):
    if "X" in protein:
        return None
    row = {"name": name, "group": group, "topology": topology}
    row.update(protein_features(protein, classes, rsa, sse))
    row.update(gene_features(best_codon_cds(protein, host), host))

    # Motif counts scale with chain length, and the design and natural sets
    # differ in length, so every count also gets a per-100-residue density.
    for stem in COUNT_FEATURES:
        for suffix in ("_raw", "_weighted"):
            key = f"{stem}{suffix}"
            if key in row:
                row[f"{stem}_density{suffix}"] = 100.0 * row[key] / len(protein)
    for key in ("exposed_hydrophobic_area", "exposed_hydrophobic_patches",
                "free_cysteines", "longest_exposed_hydrophobic_run"):
        if key in row:
            row[f"{key}_density"] = 100.0 * row[key] / len(protein)
    return row


def composition(sequence):
    counts = np.zeros(20)
    for index, residue in enumerate("ACDEFGHIKLMNPQRSTVWY"):
        counts[index] = sequence.count(residue)
    return counts / max(len(sequence), 1)


def read_manifest(manifest, backbones):
    """Records for the whole-protein arm, taken from the manifest not the directory.

    The fetcher leaves models in place when a later run narrows the selection, so
    the manifest is the authority on which of them belong to this arm.
    """
    records = []
    lines = Path(manifest).read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        accession = line.split("\t")[0]
        path = Path(backbones) / f"{accession}.pdb"
        if not path.exists():
            continue
        protein, rsa, sse = read_backbone(path)
        records.append((accession, protein, rsa, sse))
    return records


def matched_arm(designs, records, label, classes, host):
    """One control row per design, matched on length then on composition.

    Within a 10 percent length band, take the nearest neighbour in composition
    space, which is D8's second stratum. Without it a difference in low-complexity
    content could be a pure composition artefact.
    """
    if not records:
        return []
    lengths = np.array([len(protein) for _, protein, _, _ in records])
    compositions = np.array([composition(protein) for _, protein, _, _ in records])

    rows, unmatched = [], 0
    for path in designs:
        protein, _, _ = read_backbone(path)
        target = len(protein)
        close = np.flatnonzero(np.abs(lengths - target) <= 0.10 * target)
        if len(close) == 0:
            close = np.array([int(np.argmin(np.abs(lengths - target)))])
            unmatched += 1
        distances = np.linalg.norm(compositions[close] - composition(protein), axis=1)
        name, sequence, rsa, sse = records[int(close[int(np.argmin(distances))])]
        # Carry the design's topology, as the scrambled and random arms do, so
        # the per-topology comparison has matched pairs to work with.
        row = annotate(f"{name}_for_{path.stem}", label, sequence, rsa, sse,
                       classes, host, family(path.stem))
        if row:
            rows.append(row)
    print(f"\n{label}: {len(rows)} drawn from {len(records)} candidates, "
          f"{unmatched} outside the 10 percent length band")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--whole", help="manifest of whole natural proteins")
    parser.add_argument("--whole-backbones", help="directory of their models")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=2)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    print(f"{len(classes)} ELM classes loaded")
    rng = np.random.default_rng(args.seed)

    paths = sorted(Path(args.backbones).glob("*.pdb"))
    designs = [p for p in paths if not PDBLIKE.match(p.stem)]
    naturals = [p for p in paths if PDBLIKE.match(p.stem)]
    if args.limit:
        designs = designs[: args.limit]
        naturals = naturals[: args.limit]
    print(f"{len(designs)} designs, {len(naturals)} natural domains")

    rows = []
    pool = []
    for path in designs:
        protein, rsa, sse = read_backbone(path)
        pool.extend(protein)

    aggregate = np.array(pool)
    for number, path in enumerate(designs, start=1):
        protein, rsa, sse = read_backbone(path)
        topology = family(path.stem)

        row = annotate(path.stem, "design", protein, rsa, sse, classes, host, topology)
        if row:
            rows.append(row)

        shuffled = list(protein)
        rng.shuffle(shuffled)
        row = annotate(path.stem, "scrambled", "".join(shuffled), rsa, sse,
                       classes, host, topology)
        if row:
            rows.append(row)

        drawn = "".join(rng.choice(aggregate, size=len(protein)))
        row = annotate(path.stem, "random", drawn, rsa, sse, classes, host, topology)
        if row:
            rows.append(row)

        if number % 50 == 0:
            print(f"  {number}/{len(designs)} designs", flush=True)

    # Length-match the natural controls to the designs, per the analysis plan.
    natural_records = [(p.stem, *read_backbone(p)) for p in naturals]
    rows.extend(matched_arm(designs, natural_records, "natural", classes, host))

    if args.whole:
        whole_records = read_manifest(args.whole, args.whole_backbones)
        print(f"{len(whole_records)} whole natural proteins")
        rows.extend(matched_arm(designs, whole_records, "whole_natural", classes, host))

    groups = ("design", "natural", "whole_natural")
    lengths = {g: [row["length"] for row in rows if row["group"] == g] for g in groups}
    for group, values in lengths.items():
        if values:
            print(f"  {group:<14} median {int(np.median(values))}, "
                  f"range {min(values)} to {max(values)}")

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    counts = {}
    for row in rows:
        counts[row["group"]] = counts.get(row["group"], 0) + 1
    print(f"\nwrote {args.out}: {len(rows)} records")
    for group, n in sorted(counts.items()):
        print(f"  {group:<12} {n}")


if __name__ == "__main__":
    main()
