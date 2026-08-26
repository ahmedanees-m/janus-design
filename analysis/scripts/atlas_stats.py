"""Statistics for the atlas: effect sizes, normalised position, BH correction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

SKIP = {"name", "group", "topology", "n_end_class", "c_end_residue"}


def cliffs_delta(a, b):
    """Fraction of pairs where a exceeds b, minus the reverse."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    order = np.argsort(b)
    ordered = b[order]
    greater = np.searchsorted(ordered, a, side="left").sum()
    less = (len(b) - np.searchsorted(ordered, a, side="right")).sum()
    return (greater - less) / (len(a) * len(b))


def bootstrap_delta(a, b, rng, resamples=2000):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return (np.nan, np.nan)
    draws = np.empty(resamples)
    for i in range(resamples):
        draws[i] = cliffs_delta(rng.choice(a, len(a)), rng.choice(b, len(b)))
    return tuple(np.percentile(draws, [2.5, 97.5]))


def axis_width(random_set, natural):
    """Separation between the two endpoints, in units of their own scatter.

    The normalised position divides by the gap between the random and natural
    medians, so when that gap is small against the spread of the values the
    statistic is a ratio of two quantities that are both noise. Reporting the
    width lets a reader see which features the axis can carry rather than leaving
    the exclusion implicit.
    """
    pooled = np.concatenate([random_set, natural])
    pooled = pooled[~np.isnan(pooled)]
    if len(pooled) == 0:
        return np.nan
    spread = np.nanmedian(np.abs(pooled - np.nanmedian(pooled)))
    gap = abs(np.nanmedian(natural) - np.nanmedian(random_set))
    if not np.isfinite(gap):
        return np.nan
    return float(gap / spread) if spread > 1e-12 else np.inf


def normalised_position(design, random_set, natural, minimum=1.0):
    """Where designs sit on the random-to-natural axis, per the analysis plan.

    Returned only for features whose axis is at least ``minimum`` times the
    scatter of its endpoints. Elsewhere the effect size is the statistic to read.
    """
    if not np.isfinite(axis_width(random_set, natural)):
        return None
    if axis_width(random_set, natural) < minimum:
        return None
    d, r, n = (np.nanmedian(x) for x in (design, random_set, natural))
    return float((d - r) / (n - r))


def per_topology(rows, feature, control, rng, resamples, minimum=8):
    """Effect size and interval for one feature, split by topology family.

    A family that flipped one gate once, and a fold stratification with about a
    dozen designs per class, are both reasons to report this rather than a single
    pooled number. Families below ``minimum`` designs are reported as such rather
    than given an interval that would suggest more than the data holds.
    """
    families = sorted({r["topology"] for r in rows if r["group"] == "design"})
    out = []
    for name in families:
        a = np.array([r.get(feature, np.nan) for r in rows
                      if r["group"] == "design" and r["topology"] == name], dtype=float)
        b = np.array([r.get(feature, np.nan) for r in rows
                      if r["group"] == control and r["topology"] == name], dtype=float)
        record = {"topology": name, "n_design": int(np.sum(~np.isnan(a))),
                  "n_control": int(np.sum(~np.isnan(b)))}
        if record["n_design"] < minimum or record["n_control"] < minimum:
            record.update({"delta": None, "ci": None, "determined": False, "p": None})
        else:
            low, high = bootstrap_delta(a, b, rng, resamples)
            clean_a, clean_b = a[~np.isnan(a)], b[~np.isnan(b)]
            movable = np.ptp(clean_a) > 0 or np.ptp(clean_b) > 0
            record.update({"delta": float(cliffs_delta(a, b)), "ci": [low, high],
                           "determined": bool(np.isfinite(low) and np.isfinite(high)
                                              and (low > 0) == (high > 0)),
                           "p": float(mannwhitneyu(clean_a, clean_b).pvalue)
                           if movable else None})
        out.append(record)
    return out


# The hallucination family is not one method. Its members carry the design set's
# own identifiers, and the largest block is explicitly trRosetta hallucination.
DESIGN_METHODS = {
    "trRosetta hallucination": lambda name: "TrROS_Hall" in name,
    "EA set": lambda name: name.startswith("EA:"),
    "GG set": lambda name: name.startswith("GG:"),
}


def by_design_method(rows, feature, control, rng, resamples, minimum=8):
    """Split the hallucination family by the design set each member came from.

    A pooled family estimate is only worth quoting if the family is one thing.
    These are separate design campaigns that share a label here, so the split
    says whether the family effect belongs to all of them or to one.
    """
    def base(name):
        return name.split("_for_")[-1] if "_for_" in name else name

    out = []
    for label, predicate in DESIGN_METHODS.items():
        a = np.array([r.get(feature, np.nan) for r in rows
                      if r["group"] == "design" and predicate(r["name"])], dtype=float)
        b = np.array([r.get(feature, np.nan) for r in rows
                      if r["group"] == control and predicate(base(r["name"]))], dtype=float)
        record = {"method": label, "n_design": int(np.sum(~np.isnan(a))),
                  "n_control": int(np.sum(~np.isnan(b)))}
        if record["n_design"] < minimum or record["n_control"] < minimum:
            record.update({"delta": None, "ci": None, "determined": False})
        else:
            low, high = bootstrap_delta(a, b, rng, resamples)
            record.update({"delta": float(cliffs_delta(a, b)), "ci": [low, high],
                           "determined": bool(np.isfinite(low) and np.isfinite(high)
                                              and (low > 0) == (high > 0))})
        out.append(record)
    return out


def benjamini_hochberg(pvalues, rate=0.05):
    p = np.asarray(pvalues, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(len(p), np.nan)
    values = p[ok]
    order = np.argsort(values)
    ranked = values[order]
    m = len(ranked)
    adjusted = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    result = np.empty(m)
    result[order] = np.minimum(adjusted, 1.0)
    out[ok] = result
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--resamples", type=int, default=2000)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--axis-width", type=float, default=1.0,
                        help="minimum random-to-natural gap, in units of endpoint scatter")
    parser.add_argument("--topology-features", nargs="+",
                        default=["all_degron_density_weighted", "low_complexity_fraction",
                                 "exposed_hydrophobic_area_density", "longest_repeat",
                                 "gc", "cai"])
    args = parser.parse_args()

    from scipy.stats import mannwhitneyu

    rows = json.loads(Path(args.atlas).read_text(encoding="utf-8"))
    rng = np.random.default_rng(args.seed)
    groups = {g: [r for r in rows if r["group"] == g] for g in
              ("design", "scrambled", "random", "natural", "whole_natural")}
    groups = {g: v for g, v in groups.items() if v}
    print({g: len(v) for g, v in groups.items()})

    features = sorted({k for r in rows for k in r if k not in SKIP
                       and isinstance(r.get(k), (int, float))})
    print(f"{len(features)} numeric features\n")

    results = []
    for feature in features:
        values = {g: np.array([r.get(feature, np.nan) for r in v], dtype=float)
                  for g, v in groups.items()}
        if np.all(np.isnan(values["design"])):
            continue
        record = {"feature": feature}
        for control in [g for g in groups if g != "design"]:
            delta = cliffs_delta(values["design"], values[control])
            low, high = bootstrap_delta(values["design"], values[control], rng, args.resamples)
            record[f"delta_vs_{control}"] = delta
            record[f"ci_vs_{control}"] = [low, high]
            a = values["design"][~np.isnan(values["design"])]
            b = values[control][~np.isnan(values[control])]
            if len(a) > 2 and len(b) > 2 and (np.ptp(a) > 0 or np.ptp(b) > 0):
                record[f"p_vs_{control}"] = float(mannwhitneyu(a, b).pvalue)
            else:
                record[f"p_vs_{control}"] = np.nan
        record["axis_width"] = axis_width(values["random"], values["natural"])
        record["normalised_position"] = normalised_position(
            values["design"], values["random"], values["natural"], args.axis_width)
        if "whole_natural" in values:
            record["axis_width_whole"] = axis_width(values["random"],
                                                    values["whole_natural"])
            record["normalised_position_whole"] = normalised_position(
                values["design"], values["random"], values["whole_natural"],
                args.axis_width)
        record["median_design"] = float(np.nanmedian(values["design"]))
        record["median_random"] = float(np.nanmedian(values["random"]))
        record["median_natural"] = float(np.nanmedian(values["natural"]))
        if "whole_natural" in values:
            record["median_whole_natural"] = float(np.nanmedian(values["whole_natural"]))
        results.append(record)

    for control in [g for g in groups if g != "design"]:
        adjusted = benjamini_hochberg([r[f"p_vs_{control}"] for r in results])
        for record, q in zip(results, adjusted, strict=True):
            record[f"q_vs_{control}"] = None if np.isnan(q) else float(q)

    Path(args.out).write_text(json.dumps(results), encoding="utf-8")

    print("=== largest effects, designs against natural ===")
    ranked = sorted(results, key=lambda r: -abs(r["delta_vs_natural"] or 0))
    header = f"{'feature':<34}{'vs nat':>9}{'vs rand':>9}{'vs scram':>10}{'pos':>8}{'q(nat)':>10}"
    print(header)
    for record in ranked[: args.top]:
        pos = record["normalised_position"]
        q = record["q_vs_natural"]
        print(f"{record['feature']:<34}"
              f"{record['delta_vs_natural']:>9.3f}"
              f"{record['delta_vs_random']:>9.3f}"
              f"{record['delta_vs_scrambled']:>10.3f}"
              f"{('%.2f' % pos) if pos is not None else '-':>8}"
              f"{('%.2g' % q) if q is not None else '-':>10}")

    print("\n=== degron features, the Gate 4 question ===")
    degron = [r for r in results if "degron" in r["feature"]]
    print(f"{'feature':<34}{'vs nat':>9}{'CI':>20}{'vs rand':>9}{'q(nat)':>10}")
    for record in sorted(degron, key=lambda r: r["feature"]):
        low, high = record["ci_vs_natural"]
        q = record["q_vs_natural"]
        print(f"{record['feature']:<34}"
              f"{record['delta_vs_natural']:>9.3f}"
              f"{f'[{low:+.2f}, {high:+.2f}]':>20}"
              f"{record['delta_vs_random']:>9.3f}"
              f"{('%.2g' % q) if q is not None else '-':>10}")

    if "whole_natural" in groups:
        print("\n=== degron features against whole natural proteins ===")
        print("the excised arm cannot answer this: a motif outside the cut is "
              "unseen, not absent")
        print(f"{'feature':<34}{'vs whole':>10}{'CI':>20}{'vs excised':>12}{'q':>10}")
        for record in sorted(degron, key=lambda r: r["feature"]):
            low, high = record["ci_vs_whole_natural"]
            q = record["q_vs_whole_natural"]
            print(f"{record['feature']:<34}"
                  f"{record['delta_vs_whole_natural']:>10.3f}"
                  f"{f'[{low:+.2f}, {high:+.2f}]':>20}"
                  f"{record['delta_vs_natural']:>12.3f}"
                  f"{('%.2g' % q) if q is not None else '-':>10}")

        print("\n=== largest effects against whole natural proteins ===")
        ranked = sorted(results, key=lambda r: -abs(r["delta_vs_whole_natural"] or 0))
        print(f"{'feature':<34}{'vs whole':>10}{'vs excised':>12}{'design':>11}"
              f"{'whole':>11}{'q':>10}")
        for record in ranked[: args.top]:
            q = record["q_vs_whole_natural"]
            print(f"{record['feature']:<34}"
                  f"{record['delta_vs_whole_natural']:>10.3f}"
                  f"{record['delta_vs_natural']:>12.3f}"
                  f"{record['median_design']:>11.3f}"
                  f"{record.get('median_whole_natural', float('nan')):>11.3f}"
                  f"{('%.2g' % q) if q is not None else '-':>10}")

    print("\n=== normalised position, every feature whose axis is wide enough ===")
    print(f"axis width is the random-to-natural gap over the scatter of those two "
          f"groups; the threshold is {args.axis_width:.1f}")
    carried = [r for r in results if r["normalised_position"] is not None]
    excluded = len(results) - len(carried)
    print(f"{len(carried)} of {len(results)} features clear it, {excluded} do not "
          f"and are read by effect size instead")
    print(f"{'feature':<34}{'position':>10}{'width':>9}{'vs nat':>9}{'q(nat)':>10}")
    for record in sorted(carried, key=lambda r: -r["normalised_position"]):
        q = record["q_vs_natural"]
        print(f"{record['feature']:<34}{record['normalised_position']:>10.2f}"
              f"{record['axis_width']:>9.2f}{record['delta_vs_natural']:>9.3f}"
              f"{('%.2g' % q) if q is not None else '-':>10}")

    print("\n=== per-topology effect sizes ===")
    print("a pooled effect can be carried by one family, and the fold "
          "stratification is thin, so the split is the primary reporting")
    control = "whole_natural" if "whole_natural" in groups else "natural"
    topology_results = {}
    for feature in args.topology_features:
        if feature not in {r["feature"] for r in results}:
            continue
        topology_results[feature] = per_topology(rows, feature, control, rng,
                                                 args.resamples)

    # One correction across the whole per-family panel rather than inside each
    # feature. Splitting every feature by family multiplies the tests, and the
    # panel as a whole is what gets read.
    flat = [record for split in topology_results.values() for record in split]
    adjusted = benjamini_hochberg([r["p"] if r["p"] is not None else np.nan
                                   for r in flat])
    for record, q in zip(flat, adjusted, strict=True):
        record["q"] = None if not np.isfinite(q) else float(q)
    tested = sum(1 for r in flat if r["p"] is not None)
    families = len({r["topology"] for r in flat})
    print(f"\n  Benjamini-Hochberg across the whole per-family panel: {tested} "
          f"tests over {len(topology_results)} features and {families} families, "
          f"false discovery rate 0.05")

    for feature, split in topology_results.items():
        print(f"\n  {feature}, designs against {control}")
        print(f"    {'topology':<16}{'n':>6}{'delta':>9}{'95% CI':>22}{'q':>10}  verdict")
        for record in split:
            if record["delta"] is None:
                print(f"    {record['topology']:<16}{record['n_design']:>6}"
                      f"{'-':>9}{'-':>22}{'-':>10}  too few to determine")
                continue
            low, high = record["ci"]
            q = record.get("q")
            survives = record["determined"] and q is not None and q < 0.05
            verdict = ("determined" if survives
                       else "interval only" if record["determined"]
                       else "establishes nothing")
            print(f"    {record['topology']:<16}{record['n_design']:>6}"
                  f"{record['delta']:>9.3f}{f'[{low:+.2f}, {high:+.2f}]':>22}"
                  f"{('%.2g' % q) if q is not None else '-':>10}  {verdict}")

    # The hallucination family is the largest and the only one that comes out
    # enriched, so whether it is one effect or one campaign's effect matters.
    method_results = {}
    for feature in args.topology_features:
        if feature not in {r["feature"] for r in results}:
            continue
        split = by_design_method(rows, feature, control, rng, args.resamples)
        method_results[feature] = split
        print(f"\n  {feature}, hallucination split by design set")
        print(f"    {'design set':<26}{'n':>6}{'delta':>9}{'95% CI':>22}  verdict")
        for record in split:
            if record["delta"] is None:
                print(f"    {record['method']:<26}{record['n_design']:>6}"
                      f"{'-':>9}{'-':>22}  too few to determine")
                continue
            low, high = record["ci"]
            verdict = "determined" if record["determined"] else "establishes nothing"
            print(f"    {record['method']:<26}{record['n_design']:>6}"
                  f"{record['delta']:>9.3f}{f'[{low:+.2f}, {high:+.2f}]':>22}  {verdict}")

    weighted = [r for r in degron if r["feature"].endswith("_weighted")]
    raw = [r for r in degron if r["feature"].endswith("_raw")]
    print(f"\n  max |delta| against natural, weighted: "
          f"{max(abs(r['delta_vs_natural']) for r in weighted):.3f}")
    print(f"  max |delta| against natural, unweighted: "
          f"{max(abs(r['delta_vs_natural']) for r in raw):.3f}")
    topology_path = Path(args.out).with_name(Path(args.out).stem + "_by_topology.json")
    topology_path.write_text(json.dumps(topology_results), encoding="utf-8")
    method_path = Path(args.out).with_name(Path(args.out).stem + "_by_method.json")
    method_path.write_text(json.dumps(method_results), encoding="utf-8")
    print(f"\nwrote {args.out} ({len(results)} features) and {topology_path}")


if __name__ == "__main__":
    main()
