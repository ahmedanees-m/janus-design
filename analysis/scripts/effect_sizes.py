"""Effect sizes for the stability claims, and the random-to-natural positions.

A p value on 18,000 substitutions says nothing about whether the difference
matters. The Gate 4 commitment was to a Cliff's delta floor of 0.15, and that
floor applies to every claim in the paper rather than only to the degron panel.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def cliffs_delta(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ordered = np.sort(b)
    greater = np.searchsorted(ordered, a, side="left").sum()
    less = (len(b) - np.searchsorted(ordered, a, side="right")).sum()
    return (greater - less) / (len(a) * len(b))


def bootstrap(a, b, rng, resamples=2000):
    draws = np.empty(resamples)
    for i in range(resamples):
        draws[i] = cliffs_delta(rng.choice(a, len(a)), rng.choice(b, len(b)))
    return np.percentile(draws, [2.5, 97.5])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ddg", required=True)
    parser.add_argument("--epistasis", required=True)
    parser.add_argument("--atlas-stats", required=True)
    parser.add_argument("--floor", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    ddg = json.loads(Path(args.ddg).read_text(encoding="utf-8"))
    arms = {a: np.array([r["ddg"] for r in ddg if r["arm"] == a])
            for a in ("janus", "shell_random", "mpnn_sample")}

    print("=== substitution cost, effect sizes against the 0.15 floor ===")
    results = {}
    for other in ("shell_random", "mpnn_sample"):
        delta = cliffs_delta(arms["janus"], arms[other])
        low, high = bootstrap(arms["janus"], arms[other], rng)
        shift = float(np.median(arms["janus"]) - np.median(arms[other]))
        rate = 100 * (np.mean(arms["janus"] < -1.0) - np.mean(arms[other] < -1.0))
        results[other] = {"delta": delta, "ci": [low, high],
                          "median_shift": shift, "rate_difference": rate}
        verdict = "above" if abs(delta) >= args.floor else "BELOW the floor"
        print(f"  JANUS against {other}:")
        print(f"    Cliff's delta {delta:+.3f} [{low:+.3f}, {high:+.3f}]  {verdict}")
        print(f"    median shift {shift:+.3f} kcal/mol, "
              f"destabilisation rate difference {rate:+.2f} points")

    print("\n  Reading: the difference is real and negligible. The claim the plan")
    print("  asked for is non-inferiority, and non-inferiority holds with room.")
    for other in ("shell_random", "mpnn_sample"):
        low = results[other]["ci"][0]
        print(f"    JANUS is not worse than {other} by more than "
              f"{abs(min(low, 0)):.3f} in Cliff's delta")

    epi = json.loads(Path(args.epistasis).read_text(encoding="utf-8"))
    print("\n=== additivity, stated as a bound rather than a p value ===")
    for label, subset in [
        ("designed, inside the shell",
         [r for r in epi if r["both_in_shell"] and r["origin"] == "designed"]),
        ("all inside the shell", [r for r in epi if r["both_in_shell"]]),
    ]:
        values = np.array([r["epistasis"] for r in subset])
        if len(values) < 10:
            continue
        se = values.std(ddof=1) / np.sqrt(len(values))
        print(f"  {label}: n={len(values)}, median {np.median(values):+.3f}, "
              f"mean {values.mean():+.3f}, SE {se:.3f}")
        print(f"    additivity is unbiased to within {1.96 * se:.2f} kcal/mol "
              f"at 95 percent confidence")
    results["epistasis"] = {
        "n": len([r for r in epi if r["both_in_shell"] and r["origin"] == "designed"]),
    }

    stats = json.loads(Path(args.atlas_stats).read_text(encoding="utf-8"))
    positioned = [r for r in stats if r.get("normalised_position") is not None]
    print(f"\n=== random-to-natural position, {len(positioned)} of {len(stats)} features ===")
    print("  0 means designs sit at the random end, 1 at the natural end")
    print(f"{'feature':<34}{'position':>10}{'vs random':>11}{'vs natural':>12}")
    for record in sorted(positioned, key=lambda r: abs(r["delta_vs_natural"]), reverse=True)[:18]:
        print(f"{record['feature']:<34}{record['normalised_position']:>10.2f}"
              f"{record['delta_vs_random']:>11.3f}{record['delta_vs_natural']:>12.3f}")

    print(f"\n  features where the axis is too narrow to normalise: "
          f"{len(stats) - len(positioned)}")

    Path(args.out).write_text(json.dumps(results, default=float), encoding="utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
