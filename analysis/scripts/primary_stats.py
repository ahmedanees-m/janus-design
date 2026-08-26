"""The two statistics the plan committed to that were still missing.

Section 8.3 named the normalised position on the random-to-natural axis, with
bootstrap intervals, as H1's primary statistic. Every table so far has reported
design-against-natural deltas instead.

The H5 fold stratification is the most quotable result there and the most likely
to be underpowered: 614 designs across 53 classes is about a dozen each, and an
AUC of 0.5 on a dozen designs is indistinguishable from noise. Per-class
intervals decide whether the dissociation can be stated at all.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from features import protein_features
from janus.objectives.proteostasis import load_classes

HEADLINE = [
    "n_degron_density_weighted", "n_degron_density_raw",
    "all_degron_density_weighted", "all_degron_density_raw",
    "internal_degron_density_weighted", "c_degron_density_weighted",
    "exposed_hydrophobic_area_density", "longest_exposed_hydrophobic_run",
    "low_complexity_fraction", "longest_repeat", "longest_at_homopolymer",
    "initiation_dg", "transcript_mfe", "gc", "max_gc_20", "codon_pair_score",
]
DROP = {"n_end_class", "c_end_residue", "length"}


def normalised_position(design, random_set, natural):
    d, r, n = (np.nanmedian(x) for x in (design, random_set, natural))
    if not np.isfinite(n - r) or abs(n - r) < 1e-12:
        return np.nan
    return (d - r) / (n - r)


def position_interval(atlas, feature, rng, resamples=2000):
    groups = {g: np.array([r.get(feature, np.nan) for r in atlas if r["group"] == g],
                          dtype=float)
              for g in ("design", "random", "natural")}
    if any(np.all(np.isnan(v)) for v in groups.values()):
        return np.nan, (np.nan, np.nan), np.nan
    point = normalised_position(groups["design"], groups["random"], groups["natural"])
    draws = np.empty(resamples)
    for i in range(resamples):
        draws[i] = normalised_position(
            rng.choice(groups["design"], len(groups["design"])),
            rng.choice(groups["random"], len(groups["random"])),
            rng.choice(groups["natural"], len(groups["natural"])),
        )
    draws = draws[np.isfinite(draws)]
    if len(draws) < 100:
        return point, (np.nan, np.nan), np.nan
    axis = abs(np.nanmedian(groups["natural"]) - np.nanmedian(groups["random"]))
    return point, tuple(np.percentile(draws, [2.5, 97.5])), axis


def auc_interval(labels, scores, rng, resamples=2000):
    labels, scores = np.asarray(labels), np.asarray(scores)
    if len(set(labels)) < 2:
        return np.nan, (np.nan, np.nan)
    point = roc_auc_score(labels, scores)
    draws = []
    index = np.arange(len(labels))
    for _ in range(resamples):
        pick = rng.choice(index, len(index))
        if len(set(labels[pick])) < 2:
            continue
        draws.append(roc_auc_score(labels[pick], scores[pick]))
    if len(draws) < 100:
        return point, (np.nan, np.nan)
    return point, tuple(np.percentile(draws, [2.5, 97.5]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    atlas = json.loads(Path(args.atlas).read_text(encoding="utf-8"))

    print("=== H1 primary statistic: position on the random-to-natural axis ===")
    print("  0 sits at random, 1 at natural; beyond 1 is past natural")
    print(f"{'feature':<36}{'position':>10}{'95% interval':>22}{'axis width':>12}")
    positions = {}
    for feature in HEADLINE:
        point, (low, high), axis = position_interval(atlas, feature, rng)
        positions[feature] = {"position": point, "ci": [low, high], "axis": axis}
        if np.isnan(point):
            continue
        span = f"[{low:+.2f}, {high:+.2f}]" if np.isfinite(low) else "unstable"
        print(f"{feature:<36}{point:>10.2f}{span:>22}{axis:>12.4f}")

    print("\n  A wide or unstable interval means the random and natural medians are")
    print("  close, so the axis is short and the ratio is poorly determined. Those")
    print("  features are reported as effect sizes only.")

    print("\n=== H5 per-class discrimination, with intervals ===")
    classes = load_classes(args.elm)
    with Path(args.benchmark).open(newline="", encoding="utf-8", errors="replace") as fh:
        raw = list(csv.DictReader(fh))

    records = []
    for row in raw:
        sequence = (row.get("sequence") or "").strip().upper()
        label = (row.get("Exp_success") or "").strip().upper()
        plddt = row.get("ESMFold_pLDDT")
        if not sequence or label not in ("TRUE", "FALSE") or not plddt:
            continue
        if set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        computed = protein_features(sequence, classes, None, None)
        computed = {k: v for k, v in computed.items()
                    if k not in DROP and isinstance(v, (int, float))}
        for key in list(computed):
            if key.endswith("_raw"):
                computed[f"{key}_density"] = 100.0 * computed[key] / len(sequence)
        computed["seq_len"] = len(sequence)
        records.append({"fold": row.get("Fold", ""), "success": label == "TRUE",
                        "plddt": float(plddt), **computed})

    names = sorted({k for r in records for k in r
                    if isinstance(r.get(k), (int, float)) and k not in ("success", "plddt")})
    matrix = np.nan_to_num(np.array([[float(r.get(k, 0.0)) for k in names] for r in records]))
    labels = np.array([r["success"] for r in records])
    plddt = np.array([r["plddt"] for r in records])

    predictions = np.zeros(len(records))
    for train, test in StratifiedKFold(5, shuffle=True, random_state=0).split(matrix, labels):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1))
        model.fit(matrix[train], labels[train])
        predictions[test] = model.predict_proba(matrix[test])[:, 1]

    strata = {}
    for i, record in enumerate(records):
        strata.setdefault(record["fold"], []).append(i)

    print(f"{'fold class':<20}{'n':>5}{'success':>9}{'pLDDT AUC':>24}{'liability AUC':>24}")
    per_class = {}
    for fold, index in sorted(strata.items(), key=lambda kv: -len(kv[1]))[:10]:
        if len(index) < 10:
            continue
        y = labels[index]
        if len(set(y)) < 2:
            print(f"{fold:<20}{len(index):>5}{'single class':>9}")
            continue
        a, (al, ah) = auc_interval(y, plddt[index], rng)
        b, (bl, bh) = auc_interval(y, predictions[index], rng)
        per_class[fold] = {"n": len(index), "plddt": [a, al, ah], "liability": [b, bl, bh]}
        print(f"{fold:<20}{len(index):>5}{100 * y.mean():>8.0f}%"
              f"{a:>10.3f} [{al:.2f},{ah:.2f}]{b:>10.3f} [{bl:.2f},{bh:.2f}]")

    print("\n  Intervals this wide mean the per-class ordering is not established.")
    print("  The pooled comparison is the claimable result.")

    Path(args.out).write_text(json.dumps(
        {"positions": positions, "per_class": per_class}, default=float), encoding="utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
