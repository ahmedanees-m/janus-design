"""Is the local collapse in rank correlation just range restriction?

The surrogate correlates with the conditional model at 0.74 across the shell and
at 0.06 inside the parser's top 500. Those are the same relationship seen
through two window widths unless the data say otherwise, so the comparison to
make is against what restriction alone predicts.

Thorndike's case II correction, applied forward, gives the correlation expected
inside a narrower window from the wide-window correlation and the two spreads.
The quantity that does not depend on window width is the residual scatter of the
conditional score about its regression on the surrogate, and that is what should
be reported instead of a correlation.

Writes the numbers only. Plotting lives outside this repository.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr


def load_set(index_path, score_dir):
    index = json.loads(Path(index_path).read_text(encoding="utf-8"))
    score_dir = Path(score_dir)
    out = {}
    for name, record in index.items():
        conditional = []
        for i in range(1, record["n_sequences"] + 1):
            path = score_dir / f"{name}_fasta_{i}.npz"
            if not path.exists():
                conditional = None
                break
            # Per-residue mean; multiply by length for a total in nats.
            conditional.append(-float(np.mean(np.load(path)["score"])) * record["length"])
        if conditional:
            out[name] = (np.array(record["unconditional"]), np.array(conditional),
                         record["length"])
    return out


def restricted_rho(rho, wide_sd, narrow_sd):
    """Thorndike case II, applied forward from a wide window to a narrow one."""
    ratio = narrow_sd / wide_sd
    return rho * ratio / np.sqrt(1 - rho**2 * (1 - ratio**2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--narrow-index", required=True)
    parser.add_argument("--narrow-scores", required=True)
    parser.add_argument("--wide-index", required=True)
    parser.add_argument("--wide-scores", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    narrow = load_set(args.narrow_index, args.narrow_scores)
    wide = load_set(args.wide_index, args.wide_scores)
    shared = sorted(set(narrow) & set(wide))
    print(f"{len(shared)} backbones in both sets")

    rows = []
    for name in shared:
        u_wide, c_wide, length = wide[name]
        u_narrow, c_narrow, _ = narrow[name]

        rho_wide = spearmanr(u_wide, c_wide).statistic
        r_wide = pearsonr(u_wide, c_wide).statistic
        rho_narrow = spearmanr(u_narrow, c_narrow).statistic

        slope, intercept = np.polyfit(u_wide, c_wide, 1)
        residual = c_wide - (slope * u_wide + intercept)

        predicted = restricted_rho(r_wide, np.std(u_wide), np.std(u_narrow))
        rows.append({
            "backbone": name,
            "length": length,
            "rho_wide": float(rho_wide),
            "pearson_wide": float(r_wide),
            "rho_narrow": float(rho_narrow),
            "sd_unconditional_wide": float(np.std(u_wide)),
            "sd_unconditional_narrow": float(np.std(u_narrow)),
            "sd_conditional_wide": float(np.std(c_wide)),
            "residual_sd": float(np.std(residual)),
            "predicted_narrow": float(predicted),
            "band_width": float(u_narrow.max() - u_narrow.min()),
        })

    observed = np.array([r["rho_narrow"] for r in rows])
    predicted = np.array([r["predicted_narrow"] for r in rows])
    residual = np.array([r["residual_sd"] for r in rows])
    band = np.array([r["band_width"] for r in rows])

    print("\n=== restriction accounting ===")
    print(f"  wide-window Spearman:     median {np.median([r['rho_wide'] for r in rows]):.3f}")
    print(f"  wide-window Pearson:      median {np.median([r['pearson_wide'] for r in rows]):.3f}")
    print(f"  narrow-window observed:   median {np.median(observed):.3f}")
    print(f"  narrow-window predicted:  median {np.median(predicted):.3f}")
    print(f"  observed minus predicted: median {np.median(observed - predicted):+.3f}")
    print(f"  backbones where observed exceeds predicted: "
          f"{100 * np.mean(observed > predicted):.0f}%")

    print("\n=== the restriction-immune quantity ===")
    print(f"  residual scatter of conditional about its regression on the surrogate:")
    print(f"    median {np.median(residual):.2f} nats, "
          f"10th {np.percentile(residual, 10):.2f}, 90th {np.percentile(residual, 90):.2f}")
    print(f"  width of the parser's top-500 band: median {np.median(band):.2f} nats")
    print(f"  ratio, scatter to band width: median {np.median(residual / band):.1f}x")

    # rho against k, from prefixes of the same top-500 list.
    ks = [10, 25, 50, 100, 200, 300, 400, 500]
    curve = {"k": [], "observed": [], "predicted": [], "sd": []}
    for k in ks:
        obs, pred, sds = [], [], []
        for name in shared:
            u_narrow, c_narrow, _ = narrow[name]
            if len(u_narrow) < k:
                continue
            u, c = u_narrow[:k], c_narrow[:k]
            if np.std(u) == 0:
                continue
            obs.append(spearmanr(u, c).statistic)
            u_wide, c_wide, _ = wide[name]
            pred.append(restricted_rho(pearsonr(u_wide, c_wide).statistic,
                                       np.std(u_wide), np.std(u)))
            sds.append(np.std(u))
        curve["k"].append(k)
        curve["observed"].append(float(np.nanmedian(obs)))
        curve["predicted"].append(float(np.nanmedian(pred)))
        curve["sd"].append(float(np.median(sds)))

    print("\n=== rank correlation against candidates kept ===")
    print(f"{'k':>6}{'band SD':>10}{'observed':>11}{'predicted':>11}")
    for i, k in enumerate(curve["k"]):
        print(f"{k:>6}{curve['sd'][i]:>10.3f}{curve['observed'][i]:>11.3f}"
              f"{curve['predicted'][i]:>11.3f}")

    Path(args.out).write_text(json.dumps({"per_backbone": rows, "curve": curve}),
                              encoding="utf-8")


if __name__ == "__main__":
    main()
