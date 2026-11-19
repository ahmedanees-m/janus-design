"""Degradation signals read under each host's own proteostasis machinery.

The motif panel is built from Eukaryotic Linear Motif classes, so the degron
terms it carries describe a eukaryotic degradation system and transfer directly
to the mammalian host only. The bacterial host has a different and much smaller
set of sequence-encoded degradation signals, of which the N-end rule is the one
determined by the designed sequence: methionine aminopeptidase removes the
initiator methionine when the second residue is small, and ClpS then recognises
a mature N-terminal leucine, phenylalanine, tryptophan or tyrosine, with
arginine and lysine reaching the same pathway after leucylation.

This computes that signal for every arm and reports it beside the eukaryotic
panel, so each host is read under its own machinery.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

SMALL = set("ACGPSTV")
PRIMARY = set("LFWY")
SECONDARY = set("RK")


def mature_terminus(protein):
    if protein.startswith("M") and len(protein) > 1 and protein[1] in SMALL:
        return protein[1], True
    return protein[0], False


def n_end_risk(protein):
    residue, cleaved = mature_terminus(protein)
    if residue in PRIMARY:
        return 1.0, residue, cleaved
    if residue in SECONDARY:
        return 0.5, residue, cleaved
    return 0.0, residue, cleaved


def read_sequences(path):
    rows = []
    with Path(path).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", required=True,
                        help="CSV from redesign_sequences.py")
    parser.add_argument("--benchmark",
                        help="outcome benchmark CSV, for the success split")
    parser.add_argument("--label-column", default="Expressed",
                        help="benchmark column to group the summary by")
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    rows = []
    for row in read_sequences(args.index):
        risk, residue, cleaved = n_end_risk(row["protein"])
        rows.append({"source": "design set", "arm": row["arm"],
                     "name": row["backbone"], "length": int(row["length"]),
                     "n_end_residue": residue, "initiator_cleaved": cleaved,
                     "n_end_risk": risk})

    if args.benchmark:
        with Path(args.benchmark).open(newline="", encoding="utf-8",
                                       errors="replace") as fh:
            for record in csv.DictReader(fh):
                sequence = (record.get("sequence") or "").strip().upper()
                if not sequence:
                    continue
                risk, residue, cleaved = n_end_risk(sequence)
                label = (record.get(args.label_column) or "").strip()
                rows.append({"source": "outcome benchmark",
                             "arm": f"{args.label_column} {label or 'unknown'}",
                             "name": record.get("Name", ""),
                             "length": len(sequence),
                             "n_end_residue": residue,
                             "initiator_cleaved": cleaved,
                             "n_end_risk": risk})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {}
    for row in rows:
        key = (row["source"], row["arm"])
        bucket = summary.setdefault(key, {"n": 0, "primary": 0, "secondary": 0,
                                          "cleaved": 0})
        bucket["n"] += 1
        bucket["cleaved"] += int(bool(row["initiator_cleaved"]))
        if row["n_end_risk"] == 1.0:
            bucket["primary"] += 1
        elif row["n_end_risk"] == 0.5:
            bucket["secondary"] += 1

    print(f"{'source':<20}{'arm':<22}{'n':>6}{'primary':>10}{'secondary':>11}"
          f"{'cleaved':>10}")
    table = []
    for (source, arm), bucket in sorted(summary.items()):
        share = bucket["primary"] / bucket["n"]
        table.append({"source": source, "arm": arm, **bucket,
                      "primary_fraction": round(share, 4),
                      "secondary_fraction": round(bucket["secondary"] / bucket["n"], 4),
                      "cleaved_fraction": round(bucket["cleaved"] / bucket["n"], 4)})
        print(f"{source:<20}{arm:<22}{bucket['n']:>6}"
              f"{bucket['primary'] / bucket['n']:>10.3f}"
              f"{bucket['secondary'] / bucket['n']:>11.3f}"
              f"{bucket['cleaved'] / bucket['n']:>10.3f}")

    path = Path(args.summary)
    path.write_text(json.dumps(table, indent=1), encoding="utf-8")
    print(f"\nwrote {out} and {path}")


if __name__ == "__main__":
    main()
