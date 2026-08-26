"""Effect-size intervals that resample source proteins rather than records.

Both natural arms are matched to the design set with replacement, so 447 records
resolve to about a hundred distinct proteins and a record is not an independent
draw. Resampling records treats them as if they were and returns intervals that
are too narrow on the natural side.

This resamples the arms by cluster instead. A cluster is the source protein: the
PDB entry for the excised arm, the UniProt accession for the whole-protein arm,
the design itself for the design, scrambled and random arms, where every record
is already its own cluster and the two schemes coincide. Clusters are drawn with
replacement and every record belonging to a drawn cluster is carried in, so the
resample keeps the real dependence structure.

The point estimates are unchanged. Only the intervals move, and the p-values are
recomputed from the cluster distribution so that the interval and the test agree
about the same resampling unit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SKIP = {"name", "group", "topology", "n_end_class", "c_end_residue"}
ARMS = ("scrambled", "random", "natural", "whole_natural")

DESIGN_SETS = {
    "trRosetta hallucination": lambda name: "TrROS_Hall" in name,
    "EA set": lambda name: name.startswith("EA:"),
    "GG set": lambda name: name.startswith("GG:"),
}


def cluster_of(name: str) -> str:
    """The source protein a record was drawn from."""
    return name.split("_for_")[0]


def cliffs_delta(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    ordered = np.sort(b)
    greater = np.searchsorted(ordered, a, side="left").sum()
    less = (len(b) - np.searchsorted(ordered, a, side="right")).sum()
    return (greater - less) / (len(a) * len(b))


def grouped(values, clusters):
    """Values gathered by cluster, as a list of arrays."""
    order = {}
    for value, key in zip(values, clusters, strict=True):
        order.setdefault(key, []).append(value)
    return [np.asarray(v, dtype=float) for v in order.values()]


def draw(blocks, rng):
    """One cluster resample: draw clusters with replacement, keep all their records."""
    index = rng.integers(0, len(blocks), len(blocks))
    return np.concatenate([blocks[i] for i in index])


def interval(design_blocks, control_blocks, rng, resamples):
    """Percentile interval and a two-sided p-value from the cluster distribution."""
    total_design = sum(len(b) for b in design_blocks)
    total_control = sum(len(b) for b in control_blocks)
    if total_design < 3 or total_control < 3:
        return np.nan, np.nan, np.nan
    draws = np.empty(resamples)
    for i in range(resamples):
        draws[i] = cliffs_delta(draw(design_blocks, rng), draw(control_blocks, rng))
    draws = draws[np.isfinite(draws)]
    if len(draws) < resamples // 2:
        return np.nan, np.nan, np.nan
    if np.ptp(draws) == 0:
        # A feature with no spread on either side cannot be tested. Reporting a
        # p-value from a degenerate distribution would read as a determined
        # zero, which is not what a constant column establishes.
        value = float(draws[0])
        return value, value, np.nan
    low, high = np.percentile(draws, [2.5, 97.5])
    share = float(np.mean(draws > 0))
    p = 2 * min(share, 1 - share)
    return float(low), float(high), float(min(1.0, max(p, 1 / len(draws))))


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


def blocks_for(rows, feature, predicate):
    subset = [r for r in rows if predicate(r)]
    values = [r.get(feature, np.nan) for r in subset]
    clusters = [cluster_of(r["name"]) for r in subset]
    return grouped(values, clusters)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--resamples", type=int, default=4000)
    parser.add_argument("--features", nargs="+",
                        default=["all_degron_density_weighted", "low_complexity_fraction",
                                 "exposed_hydrophobic_area_density", "longest_repeat",
                                 "n_degron_density_weighted", "mean_rsa"])
    parser.add_argument("--control", default="whole_natural")
    args = parser.parse_args()

    rows = json.loads(Path(args.atlas).read_text(encoding="utf-8"))
    rng = np.random.default_rng(args.seed)

    if args.features == ["all"]:
        args.features = sorted({k for r in rows for k, v in r.items()
                                if k not in SKIP and isinstance(v, (int, float))})
        print(f"{len(args.features)} numeric features")

    sizes = {}
    for group in ("design",) + ARMS:
        names = [r["name"] for r in rows if r["group"] == group]
        if names:
            sizes[group] = (len(names), len({cluster_of(n) for n in names}))
    print("records and distinct source proteins per arm")
    for group, (records, clusters) in sizes.items():
        ratio = records / clusters
        print(f"  {group:<16}{records:>5} records{clusters:>6} proteins{ratio:>7.2f} per protein")

    results = {"arms": {g: {"records": r, "clusters": c} for g, (r, c) in sizes.items()},
               "pooled": {}, "per_topology": {}, "per_design_set": {}}

    print("\n=== pooled, cluster resampling ===")
    for feature in args.features:
        design_blocks = blocks_for(rows, feature, lambda r: r["group"] == "design")
        entry = {}
        for arm in ARMS:
            if arm not in sizes:
                continue
            control_blocks = blocks_for(rows, feature, lambda r, a=arm: r["group"] == a)
            delta = cliffs_delta(
                np.concatenate(design_blocks), np.concatenate(control_blocks))
            low, high, p = interval(design_blocks, control_blocks, rng, args.resamples)
            entry[arm] = {"delta": float(delta), "ci": [low, high], "p": p}
        results["pooled"][feature] = entry
        line = "  ".join(
            f"{arm}: {entry[arm]['delta']:+.3f} [{entry[arm]['ci'][0]:+.2f}, "
            f"{entry[arm]['ci'][1]:+.2f}]" for arm in entry)
        print(f"  {feature}\n    {line}")

    print(f"\n=== per topology, designs against {args.control} ===")
    families = sorted({r["topology"] for r in rows if r["group"] == "design"})
    flat = []
    for feature in args.features:
        split = []
        for family in families:
            design_blocks = blocks_for(
                rows, feature,
                lambda r, f=family: r["group"] == "design" and r["topology"] == f)
            control_blocks = blocks_for(
                rows, feature,
                lambda r, f=family: r["group"] == args.control and r["topology"] == f)
            if not design_blocks or not control_blocks:
                continue
            delta = cliffs_delta(
                np.concatenate(design_blocks), np.concatenate(control_blocks))
            low, high, p = interval(design_blocks, control_blocks, rng, args.resamples)
            record = {"topology": family,
                      "n_design": int(sum(len(b) for b in design_blocks)),
                      "n_control_clusters": len(control_blocks),
                      "delta": float(delta), "ci": [low, high], "p": p}
            split.append(record)
            flat.append(record)
        results["per_topology"][feature] = split

    adjusted = benjamini_hochberg([r["p"] for r in flat])
    for record, q in zip(flat, adjusted, strict=True):
        record["q"] = None if not np.isfinite(q) else float(q)

    for feature, split in results["per_topology"].items():
        print(f"\n  {feature}")
        print(f"    {'topology':<16}{'n':>5}{'clusters':>10}{'delta':>9}"
              f"{'95% CI':>22}{'q':>10}  verdict")
        for record in split:
            low, high = record["ci"]
            determined = np.isfinite(low) and (low > 0) == (high > 0)
            q = record.get("q")
            survives = determined and q is not None and q < 0.05
            verdict = ("determined" if survives
                       else "interval only" if determined
                       else "establishes nothing")
            print(f"    {record['topology']:<16}{record['n_design']:>5}"
                  f"{record['n_control_clusters']:>10}{record['delta']:>9.3f}"
                  f"{f'[{low:+.2f}, {high:+.2f}]':>22}"
                  f"{('%.2g' % q) if q is not None else '-':>10}  {verdict}")

    print(f"\n=== hallucination split by design set, against {args.control} ===")
    for feature in args.features:
        split = []
        for label, matches in DESIGN_SETS.items():
            design_blocks = blocks_for(
                rows, feature,
                lambda r, m=matches: r["group"] == "design" and m(r["name"]))
            control_blocks = blocks_for(
                rows, feature,
                lambda r, m=matches: (r["group"] == args.control
                                      and m(r["name"].split("_for_", 1)[-1])))
            if not design_blocks or not control_blocks:
                continue
            delta = cliffs_delta(
                np.concatenate(design_blocks), np.concatenate(control_blocks))
            low, high, p = interval(design_blocks, control_blocks, rng, args.resamples)
            split.append({"design_set": label,
                          "n_design": int(sum(len(b) for b in design_blocks)),
                          "n_control_clusters": len(control_blocks),
                          "delta": float(delta), "ci": [low, high], "p": p})
        results["per_design_set"][feature] = split
        if not split:
            continue
        print(f"\n  {feature}")
        print(f"    {'design set':<26}{'n':>5}{'clusters':>10}{'delta':>9}{'95% CI':>22}")
        for record in split:
            low, high = record["ci"]
            print(f"    {record['design_set']:<26}{record['n_design']:>5}"
                  f"{record['n_control_clusters']:>10}{record['delta']:>9.3f}"
                  f"{f'[{low:+.2f}, {high:+.2f}]':>22}")

    Path(args.out).write_text(json.dumps(results, indent=1), encoding="utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
