"""Fold the benchmark designs so the liability score can use accessibility.

The Garcia benchmark supplies sequences, labels and a pLDDT, but no coordinates.
Without a structure the atlas's accessibility weighting cannot be applied, so
every motif feature in the H5 score is a raw count and the score is the weaker
sequence-only form. Weighting roughly doubled the degron effect in the atlas, so
whether it helps here is a question the sequence-only run cannot answer.

This produces the coordinates. pLDDT is recorded alongside so the structures can
be filtered on model confidence rather than trusted uniformly: a liability read
off a 40-pLDDT model is a liability of the prediction, not of the design.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import torch

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def read_benchmark(path):
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for row in rows:
        sequence = (row.get("sequence") or "").strip().upper()
        if not sequence or set(sequence) - STANDARD:
            continue
        out.append((row["Name"], sequence))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--out", required=True, help="directory for the models")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--recycles", type=int, default=3)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    designs = read_benchmark(args.benchmark)
    pending = [(name, seq) for name, seq in designs
               if not (out / f"{name}.pdb").exists()]
    print(f"{len(designs)} designs, {len(pending)} still to fold", flush=True)

    if pending:
        import esm

        model = esm.pretrained.esmfold_v1().eval().cuda()
        model.set_chunk_size(args.chunk_size)
        print("model loaded", flush=True)

        for number, (name, sequence) in enumerate(pending, start=1):
            with torch.no_grad():
                pdb = model.infer_pdb(sequence, num_recycles=args.recycles)
            (out / f"{name}.pdb").write_text(pdb, encoding="utf-8")
            if number % 25 == 0:
                print(f"  {number}/{len(pending)}", flush=True)

    # pLDDT is written into the B-factor column, one value per atom.
    records = []
    for name, sequence in designs:
        path = out / f"{name}.pdb"
        if not path.exists():
            continue
        values = [float(line[60:66]) for line in path.read_text(encoding="utf-8").splitlines()
                  if line.startswith("ATOM")]
        records.append({"name": name, "length": len(sequence),
                        "mean_plddt": sum(values) / len(values) if values else None})

    summary = Path(args.summary)
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps(records), encoding="utf-8")
    folded = [r["mean_plddt"] for r in records if r["mean_plddt"] is not None]
    print(f"\n{len(records)} models, mean pLDDT {sum(folded) / len(folded):.1f}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
