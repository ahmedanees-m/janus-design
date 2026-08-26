"""Gate 8: how far the surrogate ranking is from the conditional one.

The parser optimises the sum of unconditional single-pass marginals. The model's
own likelihood is autoregressive and order dependent. This measures the distance
between the two rankings on the same sequences, and how deep the Tier-1 list has
to go before it contains the conditional best.

ProteinMPNN's score_only path reports mean negative log-likelihood per residue,
once per random decoding order. Sequences for one backbone share a length, so
negating it ranks identically to total conditional log-likelihood. The spread
across orders is reported as well: it bounds how well any conditional ranking is
defined in the first place.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


def read_scores(score_dir: Path, name: str, count: int):
    """Conditional scores per sequence, as (mean over orders, spread over orders)."""
    means, spreads = [], []
    for index in range(1, count + 1):
        path = score_dir / f"{name}_fasta_{index}.npz"
        if not path.exists():
            return None, None
        per_order = np.load(path)["score"]
        means.append(-float(np.mean(per_order)))
        spreads.append(float(np.std(per_order)))
    return np.array(means), np.array(spreads)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    index = json.loads(Path(args.index).read_text(encoding="utf-8"))
    score_dir = Path(args.scores)

    rows = []
    for name, record in sorted(index.items()):
        count = record["n_sequences"]
        conditional, spread = read_scores(score_dir, name, count)
        if conditional is None:
            print(f"  {name}: scores missing, skipped")
            continue
        unconditional = np.array(record["unconditional"])

        rho = spearmanr(unconditional, conditional).statistic
        best = int(np.argmax(conditional))
        rows.append({
            "backbone": name,
            "length": record["length"],
            "n_sequences": count,
            "spearman": float(rho),
            "conditional_best_rank": best + 1,
            "top1_is_conditional_best": best == 0,
            "order_spread_median": float(np.median(spread)),
            "order_spread_max": float(np.max(spread)),
            "conditional_range": float(conditional.max() - conditional.min()),
            "conditional_at_rank1": float(conditional[0]),
            "conditional_best": float(conditional.max()),
        })

    if not rows:
        raise SystemExit("no backbones scored")

    rho = np.array([r["spearman"] for r in rows])
    ranks = np.array([r["conditional_best_rank"] for r in rows])
    spread = np.array([r["order_spread_median"] for r in rows])
    gap = np.array([r["conditional_best"] - r["conditional_at_rank1"] for r in rows])

    print(f"=== surrogate gap over {len(rows)} backbones ===")
    print(f"  Spearman rho: median {np.median(rho):.3f}, "
          f"10th {np.percentile(rho, 10):.3f}, 90th {np.percentile(rho, 90):.3f}, "
          f"min {rho.min():.3f}")
    print(f"  conditional best at Tier-1 rank: median {int(np.median(ranks))}, "
          f"90th {int(np.percentile(ranks, 90))}, max {int(ranks.max())}")
    print(f"  Tier-1 top path is the conditional best on "
          f"{100 * np.mean(ranks == 1):.1f}% of backbones")
    print(f"  conditional likelihood gained from rank 1 to the best: "
          f"median {np.median(gap):.4f} per-residue nats")

    print("\n=== how well the conditional ranking is defined at all ===")
    print(f"  spread across 10 decoding orders, per sequence: "
          f"median {np.median(spread):.4f}, max {np.max([r['order_spread_max'] for r in rows]):.4f}")
    print(f"  median conditional range across the candidate set: "
          f"{np.median([r['conditional_range'] for r in rows]):.4f}")
    ratio = np.median(spread) / np.median([r["conditional_range"] for r in rows])
    print(f"  order noise as a fraction of the range being ranked: {ratio:.3f}")

    print("\n=== recall: conditional best inside the Tier-1 top k ===")
    for k in [1, 5, 10, 25, 50, 100, 200, 300, 400, 500]:
        covered = np.mean(ranks <= k)
        print(f"  k={k:>4}: {100 * covered:5.1f}%")

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    print(f"\nwrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
