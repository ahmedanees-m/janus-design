"""Every aggregate the manuscript quotes from the run outputs, in one place.

Each block prints the statistic, how it was aggregated and the file it came
from, so a number in the text can be traced to a column without opening the
table.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import statistics
from pathlib import Path


def rows(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def summarise(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return {"n": len(values), "median": statistics.median(values),
            "mean": statistics.mean(values),
            "p10": statistics.quantiles(values, n=10)[0] if len(values) > 9 else None,
            "p90": statistics.quantiles(values, n=10)[8] if len(values) > 9 else None}


def show(label, stats, fmt="%.4f"):
    if stats is None:
        print(f"  {label:<34} no data")
        return
    body = f"n={stats['n']:<6} median={fmt % stats['median']:<12} mean={fmt % stats['mean']}"
    if stats["p10"] is not None:
        body += f"   p10={fmt % stats['p10']}  p90={fmt % stats['p90']}"
    print(f"  {label:<34}{body}")


def structure(path):
    print("\n== structural recovery ==", path)
    data = rows(path)
    by_arm = collections.defaultdict(list)
    for row in data:
        by_arm[row["arm"]].append(row)
    print(f"  {'arm':<22}{'n':>5}{'TM median':>12}{'TM mean':>10}"
          f"{'RMSD median':>14}{'pLDDT median':>14}{'TM>=0.5':>10}")
    for arm in sorted(by_arm):
        group = by_arm[arm]
        tm = [number(r["tm_score"]) for r in group]
        rmsd = [number(r["rmsd"]) for r in group]
        plddt = [number(r["mean_plddt"]) for r in group]
        share = sum(1 for v in tm if v is not None and v >= 0.5) / len(tm)
        print(f"  {arm:<22}{len(group):>5}{statistics.median(tm):>12.3f}"
              f"{statistics.mean(tm):>10.3f}{statistics.median(rmsd):>14.2f}"
              f"{statistics.median(plddt):>14.1f}{share:>10.3f}")


def budget(path):
    print("\n== budget-matched comparison ==", path)
    data = rows(path)
    by_key = collections.defaultdict(list)
    for row in data:
        by_key[(row["arm"], row["draws"])].append(row)
    print(f"  {'arm':<22}{'draws':>7}{'gain median':>14}{'burden':>10}"
          f"{'violations':>12}{'identity':>10}")
    for arm, draws in sorted(by_key, key=lambda k: (k[0], int(k[1]))):
        group = by_key[(arm, draws)]
        print(f"  {arm:<22}{draws:>7}"
              f"{statistics.median(number(r['gain_over_codon_dp']) for r in group):>14.4f}"
              f"{statistics.median(number(r['burden']) for r in group):>10.3f}"
              f"{statistics.mean(number(r['violations']) for r in group):>12.3f}"
              f"{statistics.mean(number(r['identity_to_design']) for r in group):>10.3f}")


def sensitivity(path):
    print("\n== weight sensitivity ==", path)
    data = rows(path)
    by_scheme = collections.defaultdict(list)
    for row in data:
        by_scheme[row["scheme"]].append(row)
    removed = [k for k in data[0] if k.startswith("removed_")]
    print(f"  {'scheme':<20}{'n':>5}{'gain median':>13}{'improved':>10}"
          f"{'burden drop':>13}{'viol drop':>11}{'identity':>10}")
    for scheme in sorted(by_scheme, key=lambda s: (s != "operating", s)):
        group = by_scheme[scheme]
        gains = [number(r["gain"]) for r in group]
        improved = sum(1 for g in gains if g > 1e-9) / len(gains)
        drop = [number(r["burden_codon_dp"]) - number(r["burden_joint"]) for r in group]
        viol = [number(r["violations_codon_dp"]) - number(r["violations_joint"])
                for r in group]
        print(f"  {scheme:<20}{len(group):>5}{statistics.median(gains):>13.4f}"
              f"{improved:>10.3f}{statistics.median(drop):>13.3f}"
              f"{statistics.mean(viol):>11.3f}"
              f"{statistics.mean(number(r['identity_to_design']) for r in group):>10.3f}")
    worst = min((number(r["gain"]) for r in data), default=None)
    print(f"  worst single gain across all schemes: {worst:.4f}")
    zero = all(abs(number(r[k])) < 1e-9 for r in data for k in removed
               if r.get(k) not in (None, ""))
    print(f"  codon-only removal is zero in every scheme: {zero}")


def scaling(paths):
    print("\n== length scaling ==", ", ".join(str(p) for p in paths))
    data = []
    for path in paths:
        data.extend(rows(path))
    by_set = collections.defaultdict(list)
    for row in data:
        by_set[(row["set"], row["delta"])].append(row)
    print(f"  {'set':<18}{'delta':>7}{'n':>5}{'len med':>9}{'len max':>9}"
          f"{'parse ms':>11}{'kbest ms':>11}{'ms/residue':>12}")
    for key in sorted(by_set):
        group = by_set[key]
        lengths = [number(r["length"]) for r in group]
        parse = [number(r["parse_ms"]) for r in group]
        kbest = [number(r["kbest_ms"]) for r in group]
        per = [number(r["parse_ms"]) / number(r["length"]) for r in group]
        print(f"  {key[0]:<18}{key[1]:>7}{len(group):>5}"
              f"{statistics.median(lengths):>9.0f}{max(lengths):>9.0f}"
              f"{statistics.median(parse):>11.3f}{statistics.median(kbest):>11.3f}"
              f"{statistics.median(per):>12.4f}")

    print("\n  by length band, delta 1.0, parse and gain")
    bands = [(26, 74), (75, 120), (121, 200), (201, 320), (321, 520)]
    for low, high in bands:
        group = [r for r in data if number(r["delta"]) == 1.0
                 and low <= number(r["length"]) <= high]
        if not group:
            continue
        parse = [number(r["parse_ms"]) for r in group]
        gains = [number(r.get("gain")) for r in group]
        gains = [g for g in gains if g is not None]
        moved = [number(r.get("residues_moved")) for r in group]
        moved = [m for m in moved if m is not None]
        line = (f"  {low:>4}-{high:<5}n={len(group):<5}"
                f"parse={statistics.median(parse):>8.3f} ms")
        if gains:
            line += (f"  gain={statistics.median(gains):>8.4f}"
                     f"  improved={sum(1 for g in gains if g > 1e-9) / len(gains):.3f}")
        if moved:
            line += f"  moved={statistics.median(moved):.1f}"
        print(line)

    index = {}
    deltas = sorted({number(r["delta"]) for r in data})
    ceiling = deltas[-1]
    for row in data:
        index[(row["set"], row["backbone"], number(row["delta"]))] = row
    for budget in deltas[:-1]:
        print(f"\n  share of the delta={ceiling:g} gain reached at "
              f"delta={budget:g}, by band")
        for low, high in bands:
            shares = []
            for (set_name, backbone, delta), row in index.items():
                if delta != budget or not (low <= number(row["length"]) <= high):
                    continue
                other = index.get((set_name, backbone, ceiling))
                if other is None:
                    continue
                part, full = number(row.get("gain")), number(other.get("gain"))
                if part is None or full is None or full <= 1e-9:
                    continue
                shares.append(part / full)
            if shares:
                print(f"  {low:>4}-{high:<5}n={len(shares):<5}"
                      f"median={statistics.median(shares):.4f}"
                      f"  mean={statistics.mean(shares):.4f}")


def coverage(path):
    print("\n== prefix coverage ==", path)
    data = rows(path)
    by_weight = collections.defaultdict(list)
    for row in data:
        by_weight[number(row["folding"])].append(row)
    depths = sorted({k for k in data[0] if k.startswith("recall_top1_at")},
                    key=lambda k: int(k.rsplit("at", 1)[1]))
    print(f"  {'folding':>8}{'n':>5}{'best from prefix':>18}{'median rank':>13}"
          f"{'shortfall':>12}   " + "".join(f"{d.rsplit('at', 1)[1]:>8}" for d in depths))
    for weight in sorted(by_weight):
        group = by_weight[weight]
        from_prefix = sum(1 for r in group if r["best_from_prefix"] in ("True", "true"))
        line = (f"  {weight:>8}{len(group):>5}{from_prefix / len(group):>18.3f}"
                f"{statistics.median(number(r['best_tier1_rank']) for r in group):>13.0f}"
                f"{statistics.median(number(r['shortfall_of_prefix']) for r in group):>12.4f}   ")
        line += "".join(f"{statistics.mean(number(r[d]) for r in group):>8.3f}"
                        for d in depths)
        print(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure")
    parser.add_argument("--budget")
    parser.add_argument("--sensitivity")
    parser.add_argument("--scaling", nargs="*")
    parser.add_argument("--coverage")
    parser.add_argument("--host")
    args = parser.parse_args()

    if args.structure:
        structure(args.structure)
    if args.budget:
        budget(args.budget)
    if args.sensitivity:
        sensitivity(args.sensitivity)
    if args.scaling:
        scaling([Path(p) for p in args.scaling])
    if args.coverage:
        coverage(args.coverage)
    if args.host:
        print("\n== host proteostasis ==", args.host)
        for entry in json.loads(Path(args.host).read_text(encoding="utf-8")):
            print(f"  {entry['source']:<20}{entry['arm']:<24}"
                  f"n={entry['n']:<6}primary={entry['primary_fraction']:.4f}"
                  f"  secondary={entry['secondary_fraction']:.4f}"
                  f"  cleaved={entry['cleaved_fraction']:.4f}")


if __name__ == "__main__":
    main()
