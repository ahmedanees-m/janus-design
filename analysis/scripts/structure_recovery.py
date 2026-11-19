"""Structural agreement between a predicted model and the target backbone.

Every arm that moves residues returns a different protein, so fold retention is
measured directly. Each returned sequence is folded independently and the model
is superposed on the backbone the design was made for, which gives a TM-score and
an RMSD per arm and backbone on a scale set outside this work.

TM-align is called as an external binary and its output parsed, so the alignment
is the reference implementation rather than a reimplementation of it.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path

ALIGNED = re.compile(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+),"
                     r"\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)")
SCORE = re.compile(r"TM-score=\s*([\d.]+)\s*\(if normalized by length of "
                   r"Chain_(\d)")


def mean_plddt(text):
    values = [float(line[60:66]) for line in text.splitlines()
              if line.startswith("ATOM")]
    if not values:
        return None
    mean = sum(values) / len(values)
    # The writer emits confidence on a zero to one scale; the rest of this work
    # reports pLDDT out of a hundred.
    return mean * 100 if mean <= 1.0 else mean


def align(binary, model, target):
    out = subprocess.run([binary, str(model), str(target)],
                         capture_output=True, text=True, check=True).stdout
    match = ALIGNED.search(out)
    scores = {int(chain): float(value) for value, chain in SCORE.findall(out)}
    if match is None or 2 not in scores:
        return None
    return {"aligned": int(match.group(1)), "rmsd": float(match.group(2)),
            "seq_identity": float(match.group(3)),
            "tm_score": scores[2], "tm_score_model": scores.get(1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", required=True,
                        help="directory of <arm>__<backbone>.pdb predictions")
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--tmalign", required=True)
    parser.add_argument("--index", required=True,
                        help="CSV from redesign_sequences.py, for identity and length")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    identity = {}
    with Path(args.index).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            identity[(row["arm"], row["backbone"])] = row

    rows, missing = [], 0
    models = sorted(Path(args.models).glob("*__*.pdb"))
    print(f"{len(models)} models", flush=True)
    for number, path in enumerate(models, start=1):
        arm, _, backbone = path.stem.partition("__")
        target = Path(args.backbones) / f"{backbone}.pdb"
        if not target.exists():
            missing += 1
            continue
        result = align(args.tmalign, path, target)
        if result is None:
            missing += 1
            continue
        source = identity.get((arm, backbone), {})
        rows.append({
            "arm": arm, "backbone": backbone,
            "length": source.get("length", ""),
            "identity_to_design": source.get("identity_to_design", ""),
            "mean_plddt": mean_plddt(path.read_text(encoding="utf-8")),
            **result,
        })
        if number % 200 == 0:
            print(f"  {number}/{len(models)}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {out} ({len(rows)} rows, {missing} without an alignment)")

    import statistics
    by_arm = {}
    for row in rows:
        by_arm.setdefault(row["arm"], []).append(row)
    print(f"\n{'arm':<20}{'n':>5}{'TM-score':>12}{'RMSD':>10}{'pLDDT':>9}")
    for arm in sorted(by_arm):
        group = by_arm[arm]
        print(f"{arm:<20}{len(group):>5}"
              f"{statistics.median(r['tm_score'] for r in group):>12.3f}"
              f"{statistics.median(r['rmsd'] for r in group):>10.2f}"
              f"{statistics.median(r['mean_plddt'] for r in group):>9.1f}")


if __name__ == "__main__":
    main()
