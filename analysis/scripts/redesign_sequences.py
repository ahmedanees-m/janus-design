"""Emit the residue sequence every arm returns, for structural scoring.

The method comparison records scores but not sequences, and the structural
measurement needs the sequences themselves. This runs the two search arms over
the same backbones, operating point and objective as the comparison, reads the
proposal arms from the files their own containers wrote, and writes one record
per arm and backbone.

The output feeds `fold_sequences.py` and then `structure_recovery.py`.
"""

from __future__ import annotations

import argparse
import csv
import json
import zlib
from pathlib import Path

import numpy as np
from baselines import PDBLIKE, TIER1, coordinate_descent, hill_climb, make_objective, read_backbone

from janus import design, hosts
from janus.genetic_code import translate
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights
from janus.sample import ShellSearch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--codonmpnn")
    parser.add_argument("--proposals", nargs="*", default=[])
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--max-steps", type=int, default=40)
    parser.add_argument("--folding", type=float, default=0.25)
    parser.add_argument("--penalty", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--sequences", required=True, help="JSONL for the folder")
    parser.add_argument("--index", required=True, help="CSV of arm, backbone, sequence")
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    weights = FoldingWeights(initiation=args.folding, liability=args.penalty)

    inverse = {}
    if args.codonmpnn:
        inverse = json.loads(Path(args.codonmpnn).read_text(encoding="utf-8"))
    proposals = {}
    for item in args.proposals:
        arm, _, path = item.partition("=")
        proposals[arm] = json.loads(Path(path).read_text(encoding="utf-8"))

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    if args.shards > 1:
        files = files[args.shard::args.shards]
    print(f"{len(files)} designed backbones", flush=True)

    records = []
    for number, path in enumerate(files, start=1):
        name = path.stem
        marginals = load_unconditional(path)
        try:
            native, rsa = read_backbone(Path(args.backbones) / f"{name}.pdb")
        except Exception:
            continue
        if len(rsa) != len(marginals) or "X" in native:
            continue

        rng = np.random.default_rng(args.seed + zlib.crc32(name.encode()))
        total, _ = make_objective(marginals, host, rsa, classes, args.delta,
                                  native, args.pool, weights, rng)
        search = ShellSearch(marginals, host, TIER1, args.delta, anchor=native)
        pinned = design(marginals, host, weights=TIER1, delta=args.delta, k=1,
                        fixed=dict(enumerate(native)), anchor=native)[0]

        found = {"design": (native, pinned.cds)}

        descended, _, _ = coordinate_descent(search, pinned, total, marginals,
                                             args.rounds)
        found["coordinate_descent"] = (descended.protein, descended.cds)

        joint, _, _, _ = hill_climb(search, pinned, total, args.max_steps)
        found["joint"] = (joint.protein, joint.cds)

        for arm, table in proposals.items():
            entry = table.get(name)
            if entry is not None and len(entry["protein"]) == len(native):
                found[arm] = (entry["protein"], None)

        entry = inverse.get(name)
        if entry is not None:
            protein = translate(entry["cds"])
            if len(protein) == len(native):
                found["codonmpnn"] = (protein, entry["cds"])

        for arm, (protein, cds) in found.items():
            records.append({
                "arm": arm, "backbone": name, "length": len(protein),
                "protein": protein, "cds": cds or "",
                "identity_to_design": sum(a == b for a, b in zip(protein, native, strict=True))
                / len(native),
            })

        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    index = Path(args.index)
    index.parent.mkdir(parents=True, exist_ok=True)
    with index.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["arm", "backbone", "length",
                                                "protein", "cds",
                                                "identity_to_design"])
        writer.writeheader()
        writer.writerows(records)

    sequences = Path(args.sequences)
    sequences.parent.mkdir(parents=True, exist_ok=True)
    with sequences.open("w", encoding="utf-8") as fh:
        for row in records:
            fh.write(json.dumps({"name": f"{row['arm']}__{row['backbone']}",
                                 "seq": row["protein"]}) + "\n")

    arms = sorted({row["arm"] for row in records})
    print(f"\n{len(records)} sequences across {len(arms)} arms: {', '.join(arms)}")
    print(f"wrote {index} and {sequences}")


if __name__ == "__main__":
    main()
