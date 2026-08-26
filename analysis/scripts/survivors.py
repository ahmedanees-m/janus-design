"""Does the atlas describe designs, or only the ones that worked?

The 447 backbones in the atlas come from published, experimentally characterised
sets, so every one of them survived whatever filtering its authors applied. If
liability load is part of what selection removes, the atlas understates what a
generator actually emits, and the gap between designs and natural proteins is
larger than measured.

The Garcia benchmark answers this directly, since it carries both the designs
that worked and the designs that did not. Splitting its 614 on the experimental
outcome and comparing liability features between the two arms measures how much
of the liability load selection is removing. Structures come from the ESMFold run
so the accessibility-weighted features are available, with a pLDDT floor.

A caveat this cannot escape: a design that failed may have failed for reasons
that also produce a bad model, so the structure-derived features are entangled
with the outcome. The sequence-only features are not, and both are reported.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from atlas_stats import bootstrap_delta, benjamini_hochberg, cliffs_delta
from features import protein_features
from janus.objectives import liability
from janus.objectives.proteostasis import load_classes
from predictive import read_model

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
DROP = {"n_end_class", "c_end_residue"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--structures")
    parser.add_argument("--plddt-floor", type=float, default=70.0)
    parser.add_argument("--resamples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    classes = load_classes(args.elm)
    rng = np.random.default_rng(args.seed)

    with Path(args.benchmark).open(newline="", encoding="utf-8", errors="replace") as fh:
        raw = list(csv.DictReader(fh))

    rows = []
    for record in raw:
        sequence = (record.get("sequence") or "").strip().upper()
        label = (record.get("Exp_success") or "").strip().upper()
        if not sequence or label not in ("TRUE", "FALSE") or set(sequence) - STANDARD:
            continue
        rsa = sse = None
        if args.structures:
            path = Path(args.structures) / f"{record['Name']}.pdb"
            if path.exists():
                model = read_model(path, sequence)
                if model is not None and model[2] >= args.plddt_floor:
                    rsa, sse = model[0], model[1]
        computed = protein_features(sequence, classes, rsa, sse)
        computed = {k: v for k, v in computed.items()
                    if k not in DROP and isinstance(v, (int, float))}
        for key in list(computed):
            if key.endswith("_raw"):
                computed[f"{key}_density"] = 100.0 * computed[key] / len(sequence)
        computed["protein_repeat"] = liability.score("protein_repeat", sequence)
        if rsa is not None:
            for name in ("exposed_hydrophobic", "exposed_hydrophobic_run"):
                computed[name] = liability.score(name, sequence, rsa)
        rows.append({"name": record["Name"], "success": label == "TRUE",
                     "modelled": rsa is not None, **computed})

    successes = [r for r in rows if r["success"]]
    failures = [r for r in rows if not r["success"]]
    print(f"{len(rows)} designs, {len(successes)} succeeded and {len(failures)} failed; "
          f"{sum(r['modelled'] for r in rows)} have a model above the pLDDT floor")

    features = sorted({k for r in rows for k in r
                       if isinstance(r.get(k), (int, float)) and k not in ("success",)})
    results = []
    for feature in features:
        a = np.array([r.get(feature, np.nan) for r in failures], dtype=float)
        b = np.array([r.get(feature, np.nan) for r in successes], dtype=float)
        if np.all(np.isnan(a)) or np.all(np.isnan(b)):
            continue
        low, high = bootstrap_delta(a, b, rng, args.resamples)
        results.append({
            "feature": feature,
            "delta_failed_vs_succeeded": float(cliffs_delta(a, b)),
            "ci": [low, high],
            "median_failed": float(np.nanmedian(a)),
            "median_succeeded": float(np.nanmedian(b)),
            "determined": bool(np.isfinite(low) and np.isfinite(high) and (low > 0) == (high > 0)),
        })

    from scipy.stats import mannwhitneyu
    pvalues = []
    for record in results:
        a = np.array([r.get(record["feature"], np.nan) for r in failures], dtype=float)
        b = np.array([r.get(record["feature"], np.nan) for r in successes], dtype=float)
        a, b = a[~np.isnan(a)], b[~np.isnan(b)]
        pvalues.append(float(mannwhitneyu(a, b).pvalue)
                       if len(a) > 2 and len(b) > 2 and (np.ptp(a) > 0 or np.ptp(b) > 0)
                       else np.nan)
    for record, q in zip(results, benjamini_hochberg(pvalues), strict=True):
        record["q"] = None if np.isnan(q) else float(q)

    ranked = sorted(results, key=lambda r: -abs(r["delta_failed_vs_succeeded"]))
    print()
    print("=== liability in failed designs against successful ones ===")
    print("positive means the failures carry more of it, so selection removes it")
    print(f"{'feature':<36}{'delta':>9}{'95% CI':>20}{'failed':>10}{'passed':>10}{'q':>10}")
    for record in ranked[: args.top]:
        low, high = record["ci"]
        q = record["q"]
        print(f"{record['feature']:<36}{record['delta_failed_vs_succeeded']:>9.3f}"
              f"{f'[{low:+.2f}, {high:+.2f}]':>20}"
              f"{record['median_failed']:>10.3f}{record['median_succeeded']:>10.3f}"
              f"{('%.2g' % q) if q is not None else '-':>10}")

    determined = [r for r in results if r["determined"]]
    positive = [r for r in determined if r["delta_failed_vs_succeeded"] > 0]
    print()
    print(f"{len(determined)} of {len(results)} features are determined; "
          f"{len(positive)} of those are higher in the failures")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
