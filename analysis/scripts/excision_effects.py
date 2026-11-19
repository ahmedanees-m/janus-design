"""Effect sizes on each side of the excised arm's parent-coverage split.

The excised control arm is drawn from natural MegaScale domains, most of which
are fragments of a larger parent. A fragment presents a terminus to the solvent
that was internal in the parent, so terminal features and the degron panel that
reads them are not comparable across the two kinds of entry.

This splits the arm by how much of its UniProt parent each entry covers and
reports Cliff's delta against the design set separately on each side, with
cluster bootstrap intervals over source entries, the same resampling unit the
pooled analysis uses.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from cluster_bootstrap import SKIP, cliffs_delta, cluster_of, grouped, interval

FEATURES = ["all_degron_density_weighted", "n_degron_density_weighted",
            "c_degron_density_weighted", "internal_degron_density_weighted",
            "n_term_rsa", "c_term_rsa", "mean_rsa",
            "low_complexity_fraction", "longest_repeat",
            "exposed_hydrophobic_area_density"]


def blocks(rows, feature):
    return grouped([r.get(feature, np.nan) for r in rows],
                   [cluster_of(r["name"]) for r in rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--coverage", required=True,
                        help="JSON written by excision.py, PDB entry to coverage")
    parser.add_argument("--threshold", type=float, default=0.90)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--resamples", type=int, default=4000)
    parser.add_argument("--features", nargs="+", default=FEATURES)
    parser.add_argument("--out", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args()

    rows = json.loads(Path(args.atlas).read_text(encoding="utf-8"))
    coverage = json.loads(Path(args.coverage).read_text(encoding="utf-8"))
    rng = np.random.default_rng(args.seed)

    if args.features == ["all"]:
        args.features = sorted({k for r in rows for k, v in r.items()
                                if k not in SKIP and isinstance(v, (int, float))})

    design = [r for r in rows if r["group"] == "design"]
    excised = [r for r in rows if r["group"] == "natural"]
    whole, fragment, unresolved = [], [], []
    for row in excised:
        value = coverage.get(cluster_of(row["name"]))
        if value is None:
            unresolved.append(row)
        elif value >= args.threshold:
            whole.append(row)
        else:
            fragment.append(row)

    entries = {"whole": {cluster_of(r["name"]) for r in whole},
               "fragment": {cluster_of(r["name"]) for r in fragment},
               "unresolved": {cluster_of(r["name"]) for r in unresolved}}
    print(f"excised arm: {len(excised)} records from "
          f"{len({cluster_of(r['name']) for r in excised})} entries")
    for label, subset in (("whole", whole), ("fragment", fragment),
                          ("unresolved", unresolved)):
        print(f"  {label:<12}{len(subset):>5} records "
              f"{len(entries[label]):>4} entries")
    # An entry can back several records, so the two medians differ. The entry
    # median describes the parents, the record median describes the arm.
    per_record = [coverage[cluster_of(r["name"])] for r in fragment]
    per_entry = sorted(coverage[name] for name in entries["fragment"])
    if per_record:
        print(f"  median fragment coverage {np.median(per_entry):.4f} per entry, "
              f"{np.median(per_record):.4f} per record")

    results, table = {}, []
    for feature in args.features:
        design_blocks = blocks(design, feature)
        entry = {}
        for label, subset in (("pooled", excised), ("fragment", fragment),
                              ("whole", whole)):
            control_blocks = blocks(subset, feature)
            if not control_blocks:
                continue
            delta = cliffs_delta(np.concatenate(design_blocks),
                                 np.concatenate(control_blocks))
            low, high, p = interval(design_blocks, control_blocks, rng,
                                    args.resamples)
            entry[label] = {"delta": float(delta), "ci": [low, high], "p": p,
                            "records": len(subset), "entries": len(control_blocks)}
            table.append({"feature": feature, "subset": label,
                          "records": len(subset), "entries": len(control_blocks),
                          "delta": round(float(delta), 6),
                          "ci_low": round(low, 6), "ci_high": round(high, 6),
                          "p": None if not np.isfinite(p) else round(float(p), 6)})
        results[feature] = entry
        line = "  ".join(
            f"{label}: {entry[label]['delta']:+.3f} "
            f"[{entry[label]['ci'][0]:+.2f}, {entry[label]['ci'][1]:+.2f}]"
            for label in entry)
        print(f"  {feature}\n    {line}")

    payload = {"threshold": args.threshold,
               "counts": {label: {"records": len(subset),
                                  "entries": len(entries[label])}
                          for label, subset in (("whole", whole),
                                                ("fragment", fragment),
                                                ("unresolved", unresolved))},
               "median_fragment_coverage_per_entry":
                   float(np.median(per_entry)) if per_entry else None,
               "median_fragment_coverage_per_record":
                   float(np.median(per_record)) if per_record else None,
               "features": results}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    path = Path(args.table)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["feature", "subset", "records",
                                                "entries", "delta", "ci_low",
                                                "ci_high", "p"])
        writer.writeheader()
        writer.writerows(table)
    print(f"\nwrote {out} and {path}")


if __name__ == "__main__":
    main()
