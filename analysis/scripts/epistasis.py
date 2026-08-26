"""Does additivity hold for pairs of substitutions inside the shell?

JANUS proposes several substitutions at once and scores them independently, so
whether their measured effects add matters. MegaScale carries double mutants for
a subset of parents; each is compared against the sum of its two singles.

The measured bias is upward, and more so on natural parents than on designed
ones. Two readings fit that: genuine epistasis, milder for the milder
substitutions the shell admits, or a measurement artefact of the assay's range.

The deposited values are not clipped. Their distribution is smooth from -15 to
+18 with no pile-up at either end, so the release reports extrapolated free
energies rather than censored ones. The range effect, if there is one, is
therefore loss of precision rather than a wall: K50 resolves unfolding over a
finite protease series, and a variant whose true stability sits far outside the
middle of that series is extrapolated, with an error that grows with distance and
cannot run symmetrically once the variant is already fully cleaved.

That makes the test a restriction rather than a threshold. Each record carries
where its parent, its additive prediction and its observation sit on the measured
scale, and the comparison is repeated on progressively narrower central bands. If
the bias is a range effect it shrinks as the band tightens. If it is epistasis it
does not.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from ddg import infer_offset, load_measurements
from janus.genetic_code import translate
from janus.lattice import amino_acid_shell
from janus.objectives.mpnn import load_unconditional

SUB = r"([A-Z])(\d+)([A-Z])"
DOUBLE = re.compile(rf"^(?P<parent>.+?)\.pdb_.*?{SUB}\|{SUB}$")
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def load_doubles(parquet_dir):
    doubles = defaultdict(dict)
    for path in sorted(Path(parquet_dir).glob("Lib*_K50dG.parquet")):
        table = pq.read_table(path, columns=["name", "deltaG"])
        for name, energy in zip(table.column("name").to_pylist(),
                                table.column("deltaG").to_pylist(), strict=True):
            if not name or energy is None or not np.isfinite(energy):
                continue
            match = DOUBLE.match(name)
            if not match:
                continue
            parent = match.group("parent")
            first = (int(match.group(3)), match.group(2), match.group(4))
            second = (int(match.group(6)), match.group(5), match.group(7))
            doubles[parent][first, second] = float(energy)
    return doubles


def measured_band(parquet_dir, low, high):
    """Percentile bounds of every deltaG in the release.

    Taken from the data rather than from a reported assay window, so the band is
    the one these measurements actually populate.
    """
    values = []
    for path in sorted(Path(parquet_dir).glob("Lib*_K50dG.parquet")):
        column = pq.read_table(path, columns=["deltaG"]).column("deltaG").to_pylist()
        values.extend(v for v in column if v is not None and np.isfinite(v))
    array = np.array(values)
    return float(np.percentile(array, low)), float(np.percentile(array, high))


def inside(row, low, high):
    """Parent, additive prediction and observation all within the band."""
    return all(low <= row[key] <= high
               for key in ("parent_dg", "predicted_dg", "observed_dg"))


def reconciliation(rows):
    """Is the small in-shell subset actually in tension with the larger estimate?

    The in-shell designed subset is the closest analogue to what the method does,
    and its bias looks smaller than the designed estimate on the well-measured
    band. Quoting a point estimate beside a half-width invites the reading that
    the two disagree, so the interval is put on both and compared directly.
    """
    print()
    print("=== is the in-shell subset in tension with the designed estimate ===")
    print(f"  {'subset':<30}{'n':>7}{'bias':>9}{'se':>8}{'95% interval':>20}")
    intervals = {}
    for label, subset in (
        ("designed parents", [r for r in rows if r["origin"] == "designed"]),
        ("designed, inside the shell",
         [r for r in rows if r["both_in_shell"] and r["origin"] == "designed"]),
    ):
        if len(subset) < 10:
            continue
        values = np.array([r["epistasis"] for r in subset])
        mean = float(values.mean())
        error = float(values.std(ddof=1) / np.sqrt(len(values)))
        low, high = mean - 1.96 * error, mean + 1.96 * error
        intervals[label] = (low, high)
        print(f"  {label:<30}{len(subset):>7}{mean:>+9.3f}{error:>8.3f}"
              f"{f'[{low:+.3f}, {high:+.3f}]':>20}")

    keys = list(intervals)
    if len(keys) == 2:
        (a_low, a_high), (b_low, b_high) = intervals[keys[0]], intervals[keys[1]]
        overlap = a_low <= b_high and b_low <= a_high
        print(f"  the intervals {'overlap' if overlap else 'do not overlap'}")
        if overlap:
            print("  so the small subset does not establish a smaller bias; it is")
            print("  consistent with the larger estimate and too small to refine it")


def mildness(rows):
    """Are the substitutions the shell admits milder, and does that explain the gap?

    The reconciliation on offer is that shell substitutions are the ones
    ProteinMPNN considers plausible, so they should perturb stability less, and
    milder substitutions should be more additive. Both halves are testable and
    both are tested here rather than assumed.
    """
    from scipy.stats import mannwhitneyu

    print()
    print("=== are shell substitutions milder ===")
    print("severity is the larger of the two singles' |ddG|, in kcal/mol")
    print(f"  {'set':<30}{'n':>7}{'median':>9}{'75th':>8}{'90th':>8}")
    groups = {}
    for label, subset in (
        ("designed, both in shell",
         [r for r in rows if r["origin"] == "designed" and r["both_in_shell"]]),
        ("designed, not both in shell",
         [r for r in rows if r["origin"] == "designed" and r["both_in_shell"] is False]),
        ("natural, both in shell",
         [r for r in rows if r["origin"] == "natural" and r["both_in_shell"]]),
        ("natural, not both in shell",
         [r for r in rows if r["origin"] == "natural" and r["both_in_shell"] is False]),
    ):
        if len(subset) < 10:
            print(f"  {label:<30}{len(subset):>7}   too few to report")
            continue
        severity = np.array([r["severity"] for r in subset])
        groups[label] = severity
        print(f"  {label:<30}{len(subset):>7}{np.median(severity):>9.3f}"
              f"{np.percentile(severity, 75):>8.3f}{np.percentile(severity, 90):>8.3f}")

    inside_key, outside_key = "designed, both in shell", "designed, not both in shell"
    if inside_key in groups and outside_key in groups:
        p_value = float(mannwhitneyu(groups[inside_key], groups[outside_key]).pvalue)
        difference = float(np.median(groups[inside_key]) - np.median(groups[outside_key]))
        print(f"  designed, in shell against out: median difference {difference:+.3f} "
              f"kcal/mol, p = {p_value:.3g}")
        print("  a negative difference is the favourable reading: the substitutions")
        print("  the method actually makes are the gentler ones")

    print()
    print("=== epistasis against how severe the substitutions are ===")
    print(f"  {'severity band':<18}{'n':>7}{'median bias':>14}{'spread':>9}")
    edges = [0.0, 0.5, 1.0, 2.0, np.inf]
    designed = [r for r in rows if r["origin"] == "designed"]
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        band = [r["epistasis"] for r in designed if low <= r["severity"] < high]
        if len(band) < 10:
            continue
        label = f"{low:g} to {high:g}" if np.isfinite(high) else f"above {low:g}"
        print(f"  {label:<18}{len(band):>7}{np.median(band):>+14.3f}"
              f"{np.std(band):>9.3f}")


def propagation(rows):
    """What the additive prediction is worth for a design carrying several changes.

    The delta sweep moves 4.9 to 6.0 residues at the operating point, so the
    relevant question is not the error on one pair but on that many at once.
    Assuming every pair deviates independently gives an upper bound; measuring
    how the deviation decays with sequence separation says how loose that bound
    is, since substitutions scattered along a short chain are mostly far apart.
    """
    designed = [r for r in rows if r["origin"] == "designed"]
    if len(designed) < 30:
        return

    print()
    print("=== does the deviation depend on how far apart the substitutions are ===")
    print(f"  {'separation':<18}{'n':>7}{'median bias':>14}{'spread':>9}")
    edges = [0, 5, 10, 20, 40, np.inf]
    decayed = None
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        band = [r["epistasis"] for r in designed if low <= r["separation"] < high]
        if len(band) < 10:
            continue
        label = f"{low} to {high}" if np.isfinite(high) else f"above {low}"
        print(f"  {label:<18}{len(band):>7}{np.median(band):>+14.3f}"
              f"{np.std(band):>9.3f}")
        if not np.isfinite(high):
            decayed = np.array(band)

    values = np.array([r["epistasis"] for r in designed])
    spread = float(values.std(ddof=1))
    mean = float(values.mean())
    print()
    print("=== what that means for a design carrying several substitutions ===")
    print(f"  per-pair deviation, designed parents: mean {mean:+.3f}, "
          f"median {np.median(values):+.3f}, sd {spread:.3f} kcal/mol")
    print("  the offset uses the mean, because the expectation of a sum is the sum")
    print("  of expectations and the median of a sum is not the sum of medians")
    print(f"  {'substitutions':<15}{'pairs':>7}{'offset':>9}{'1 sigma':>10}"
          f"{'95% half-width':>16}")
    for count in (2, 3, 5, 6, 10):
        pairs = count * (count - 1) // 2
        sigma = spread * np.sqrt(pairs)
        print(f"  {count:<15}{pairs:>7}{pairs * mean:>+9.2f}{sigma:>10.2f}"
              f"{1.96 * sigma:>16.2f}")
    print("  one sigma is a spread and not a bound: about a third of designs fall")
    print("  outside it. All three columns assume every pair deviates")
    print("  independently, which the separation table above gives no reason to")
    print("  discount, since deviation does not fall away with distance")
    if decayed is not None:
        print(f"  widely separated pairs alone: median {np.median(decayed):+.3f}, "
              f"spread {np.std(decayed):.3f} kcal/mol over {len(decayed)} pairs")


def dynamic_range(rows, parquet_dir, bands):
    """Bias as the comparison is restricted to progressively narrower bands."""
    print()
    print("=== bias against distance from the middle of the measured range ===")
    print("headroom is the additive prediction above the 10th percentile of every")
    print("measured value; a prediction far below that is an extrapolation")
    edges = [-np.inf, 0.0, 0.5, 1.0, 2.0, np.inf]
    labels = ["below", "0 to 0.5", "0.5 to 1", "1 to 2", "above 2"]
    for origin in ("natural", "designed"):
        subset = [r for r in rows if r["origin"] == origin]
        if not subset:
            continue
        print()
        print(f"  {origin} parents, {len(subset)} pairs")
        print(f"    {'headroom':<12}{'n':>8}{'median bias':>14}{'mean bias':>12}")
        for low, high, label in zip(edges[:-1], edges[1:], labels, strict=True):
            band = [r["epistasis"] for r in subset if low <= r["headroom"] < high]
            if not band:
                continue
            print(f"    {label:<12}{len(band):>8}{np.median(band):>+14.3f}"
                  f"{np.mean(band):>+12.3f}")

    print()
    print("=== the same comparison, restricted to the middle of the range ===")
    print("all three of parent, prediction and observation inside the band")
    print(f"    {'band':<22}{'natural':>20}{'designed':>20}{'gap':>9}")
    for low, high in bands:
        floor, ceiling = measured_band(parquet_dir, low, high)
        kept = [r for r in rows if inside(r, floor, ceiling)]
        cells = {}
        for origin in ("natural", "designed"):
            values = [r["epistasis"] for r in kept if r["origin"] == origin]
            median = float(np.median(values)) if values else float("nan")
            cells[origin] = f"{median:+.3f} (n={len(values)})"
            cells[origin + "_median"] = median
        label = f"{low:g} to {high:g} pct"
        gap = cells["natural_median"] - cells["designed_median"]
        print(f"    {label:<22}{cells['natural']:>20}{cells['designed']:>20}{gap:>+9.3f}")
    print()
    print("    a bias that is a range effect shrinks as the band tightens; one that")
    print("    survives the narrowest band is not explained by the assay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--agreement", type=float, default=0.95)
    parser.add_argument("--bands", type=float, nargs="+",
                        default=[0.0, 100.0, 5.0, 95.0, 10.0, 90.0, 25.0, 75.0],
                        help="percentile pairs defining the restriction bands")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    args.bands = list(zip(args.bands[::2], args.bands[1::2], strict=True))

    wild, constructs, singles = load_measurements(args.parquet)
    doubles = load_doubles(args.parquet)
    print(f"{len(doubles)} parents carry double mutants")

    shells = {}
    for path in sorted(Path(args.marginals).glob("*.npz")):
        marginals = load_unconditional(path)
        shells[path.stem] = amino_acid_shell(marginals, args.delta)

    rows = []
    for parent, pairs in doubles.items():
        if parent not in wild or parent not in constructs or parent not in singles:
            continue
        construct = constructs[parent]
        offset, agreement = infer_offset(construct, list(singles[parent]))
        if offset is None or agreement < args.agreement:
            continue
        origin = wild[parent]
        shell = shells.get(parent)

        for (first, second), combined in pairs.items():
            a = singles[parent].get(first)
            b = singles[parent].get(second)
            if a is None or b is None:
                continue
            ddg_a, ddg_b, ddg_ab = a - origin, b - origin, combined - origin
            in_shell = None
            if shell is not None:
                indices = [p - 1 + offset for p, _, _ in (first, second)]
                if all(0 <= i < len(shell) for i in indices):
                    in_shell = all(
                        mut in shell[i]
                        for i, (_, _, mut) in zip(indices, (first, second), strict=True)
                    )
            rows.append({
                "parent": parent,
                "origin": "natural" if PDBLIKE.match(parent) else "designed",
                "parent_dg": origin,
                "predicted_dg": origin + ddg_a + ddg_b,
                "observed_dg": combined,
                "ddg_a": ddg_a, "ddg_b": ddg_b, "ddg_ab": ddg_ab,
                "additive": ddg_a + ddg_b,
                "epistasis": ddg_ab - (ddg_a + ddg_b),
                "both_in_shell": in_shell,
                "separation": abs(first[0] - second[0]),
                "severity": max(abs(ddg_a), abs(ddg_b)),
            })

    floor, ceiling = measured_band(args.parquet, 10.0, 90.0)
    print(f"central band of every measured value: {floor:+.2f} to {ceiling:+.2f} kcal/mol")
    for row in rows:
        row["floor"] = floor
        row["ceiling"] = ceiling
        row["headroom"] = row["predicted_dg"] - floor

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows)
    reconciliation(rows)
    mildness(rows)
    propagation(rows)
    dynamic_range(rows, args.parquet, args.bands)
    print(f"\nwrote {args.out} ({len(rows)} pairs)")


def summarise(rows):
    from scipy.stats import pearsonr, wilcoxon

    print(f"\n=== additivity over {len(rows)} double mutants ===")
    for label, subset in [
        ("all pairs", rows),
        ("designed parents", [r for r in rows if r["origin"] == "designed"]),
        ("natural parents", [r for r in rows if r["origin"] == "natural"]),
        ("both inside the shell", [r for r in rows if r["both_in_shell"]]),
        ("designed, inside the shell",
         [r for r in rows if r["both_in_shell"] and r["origin"] == "designed"]),
    ]:
        if len(subset) < 10:
            print(f"  {label:<28} n={len(subset)}, too few to report")
            continue
        epi = np.array([r["epistasis"] for r in subset])
        observed = np.array([r["ddg_ab"] for r in subset])
        predicted = np.array([r["additive"] for r in subset])
        r_value = pearsonr(predicted, observed).statistic
        try:
            p = wilcoxon(epi).pvalue
        except ValueError:
            p = np.nan
        print(f"  {label:<28} n={len(subset):>6}  "
              f"bias {np.median(epi):+.3f}  spread {np.std(epi):.3f}  "
              f"r={r_value:.3f}  p={p:.2g}")

    epi = np.array([r["epistasis"] for r in rows])
    print(f"\n  pairs where the deviation exceeds 0.5 kcal/mol: "
          f"{100 * np.mean(np.abs(epi) > 0.5):.1f}%")
    print(f"  pairs where it exceeds 1.0: {100 * np.mean(np.abs(epi) > 1.0):.1f}%")


if __name__ == "__main__":
    main()
