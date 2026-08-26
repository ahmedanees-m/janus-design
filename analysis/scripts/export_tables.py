"""Export the analysis outputs as supplementary tables in CSV.

Each JSON produced by the analysis becomes one flat CSV with a stable column
order. Nested records are flattened with a dotted key so nothing is dropped.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

EXPORTS = [
    ("S1_liability_atlas", "atlas.json", None),
    ("S2_effect_sizes_pooled", "atlas_stats.json", "features"),
    ("S3_effect_sizes_by_topology", "atlas_stats_by_topology.json", "flatten_by_key"),
    ("S4_effect_sizes_by_design_set", "atlas_stats_by_method.json", "flatten_by_key"),
    ("S5_entropy_budget", "entropy_budget.json", None),
    ("S6_liability_exchange", "aa_axis.json", None),
    ("S7_delta_sweep", "delta_sweep.json", None),
    ("S8_folding_frontier", "folding_sweep.json", None),
    ("S9_method_comparison", "baselines.json", None),
    ("S10_measured_stability", "ddg.json", None),
    ("S11_additivity", "epistasis.json", None),
    ("S12_predictive_features", "predictive.json", "results"),
    ("S13_cross_host_case_study", "case_study.json", None),
    ("S14_survivor_split", "survivors.json", None),
    ("S15_excision_check", "excision.json", None),
    ("S16_range_restriction", "restriction.json", "per_backbone"),
    ("S17_layer_decomposition", "decompose.json", None),
    ("S18_optimality_gap", "gap_coupling.json", None),
]


def flatten(record, prefix=""):
    out = {}
    for key, value in record.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten(value, f"{name}."))
        elif isinstance(value, (list, tuple)):
            if value and all(isinstance(v, (int, float, type(None))) for v in value):
                for index, item in enumerate(value):
                    out[f"{name}.{index}"] = item
            else:
                out[name] = json.dumps(value)
        else:
            out[name] = value
    return out


def rows_from(blob, mode):
    if mode == "flatten_by_key":
        rows = []
        for feature, records in blob.items():
            for record in records:
                rows.append({"feature": feature, **record})
        return rows
    if mode and isinstance(blob, dict):
        blob = blob[mode]
    if isinstance(blob, dict):
        return [{"key": key, **(value if isinstance(value, dict) else {"value": value})}
                for key, value in blob.items()]
    return blob


def write_table(rows, path):
    flat = [flatten(row) for row in rows]
    columns = []
    for row in flat:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat)
    return len(flat), len(columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    source = Path(args.processed)
    target = Path(args.out)
    target.mkdir(parents=True, exist_ok=True)

    manifest = []
    for name, filename, mode in EXPORTS:
        path = source / filename
        if not path.exists():
            print(f"skipped {name}: {filename} not present")
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        rows = rows_from(blob, mode)
        if not rows:
            print(f"skipped {name}: no rows")
            continue
        count, columns = write_table(rows, target / f"Table_{name}.csv")
        manifest.append((f"Table_{name}.csv", filename, count, columns))
        print(f"{name:<34}{count:>9} rows{columns:>5} columns")

    with (target / "Table_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "source", "rows", "columns"])
        writer.writerows(manifest)
    print(f"\nwrote {len(manifest)} tables and a manifest to {target}")


if __name__ == "__main__":
    main()
