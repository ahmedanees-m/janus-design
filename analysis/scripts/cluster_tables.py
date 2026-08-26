"""Turn the cluster bootstrap into the shapes the tables and figures read.

Writes three files:

``cluster_stats.json``      one record per feature with a delta, interval and
                            Benjamini-Hochberg q per control arm, corrected
                            across the whole numeric panel separately for each
                            arm, which is the same correction the record
                            bootstrap uses.
``cluster_by_topology.json``one record per feature and topology family, with a
                            single correction across the pre-specified panel.
``Table_S21.csv``           both of the above flattened, plus the design-set
                            split.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def benjamini_hochberg(values):
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    finite = np.where(np.isfinite(values))[0]
    if not len(finite):
        return out
    order = finite[np.argsort(values[finite])]
    n = len(order)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        running = min(running, values[order[rank]] * n / (rank + 1))
        out[order[rank]] = running
    return out


def determined(interval):
    low, high = interval
    return bool(np.isfinite(low) and np.isfinite(high) and (low > 0) == (high > 0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", required=True, help="cluster bootstrap over every feature")
    parser.add_argument("--panel", required=True, help="cluster bootstrap over the pre-specified panel")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    full = json.loads(Path(args.full).read_text(encoding="utf-8"))
    panel = json.loads(Path(args.panel).read_text(encoding="utf-8"))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    arms = ["scrambled", "random", "natural", "whole_natural"]
    features = sorted(full["pooled"])

    stats = [{"feature": f} for f in features]
    for arm in arms:
        raw = []
        for record, feature in zip(stats, features, strict=True):
            entry = full["pooled"][feature].get(arm)
            if entry is None:
                record[f"delta_vs_{arm}"] = None
                record[f"ci_vs_{arm}"] = [None, None]
                raw.append(np.nan)
                continue
            record[f"delta_vs_{arm}"] = entry["delta"]
            record[f"ci_vs_{arm}"] = entry["ci"]
            raw.append(entry["p"] if entry["p"] is not None else np.nan)
        for record, q in zip(stats, benjamini_hochberg(raw), strict=True):
            record[f"q_vs_{arm}"] = None if not np.isfinite(q) else float(q)

    (out / "cluster_stats.json").write_text(json.dumps(stats, indent=1), encoding="utf-8")
    print(f"cluster_stats.json: {len(stats)} features, {len(arms)} arms")

    by_topology = {}
    for feature, split in panel["per_topology"].items():
        by_topology[feature] = [
            {"topology": r["topology"], "n_design": r["n_design"],
             "n_control": r["n_control_clusters"], "delta": r["delta"],
             "ci": r["ci"], "p": r["p"], "q": r.get("q"),
             "determined": determined(r["ci"])}
            for r in split
        ]
    (out / "cluster_by_topology.json").write_text(
        json.dumps(by_topology, indent=1), encoding="utf-8")
    tested = sum(1 for s in by_topology.values() for r in s
                 if r["p"] is not None and np.isfinite(r["p"]))
    print(f"cluster_by_topology.json: {len(by_topology)} features, {tested} tested cells")

    rows = []
    for feature in features:
        for arm in arms:
            entry = full["pooled"][feature].get(arm)
            if entry is None:
                continue
            record = next(s for s in stats if s["feature"] == feature)
            rows.append(["pooled", feature, arm, "", entry["delta"],
                         entry["ci"][0], entry["ci"][1], entry["p"],
                         record[f"q_vs_{arm}"]])
    for feature, split in by_topology.items():
        for r in split:
            rows.append(["per_topology", feature, "whole_natural", r["topology"],
                         r["delta"], r["ci"][0], r["ci"][1], r["p"], r["q"]])
    for feature, split in panel["per_design_set"].items():
        for r in split:
            rows.append(["per_design_set", feature, "whole_natural", r["design_set"],
                         r["delta"], r["ci"][0], r["ci"][1], r["p"], ""])

    target = out / "Table_S21_cluster_bootstrap.csv"
    with target.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["level", "feature", "control_arm", "stratum", "delta",
                         "ci_low", "ci_high", "p", "q"])
        for row in rows:
            writer.writerow(
                [c if not isinstance(c, float) else f"{c:.6f}" for c in row])
    print(f"{target.name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
