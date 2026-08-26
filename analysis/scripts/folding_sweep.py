"""The folding-weight frontier, in normalised units.

The Tier-1 score is in nats and the initiation term in kcal/mol, so an
unnormalised weight is an accidental unit conversion and will not transfer
between backbones of different length. Each term is divided by its spread across
the candidate pool for that backbone first, which makes the weight dimensionless
and reads as standard deviations of one term traded for the other.

The candidate pool is the k-best prefix plus uniform draws from the same shell,
because the prefix alone stops being the better pool once the folding weight is
high.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from janus import Weights, design, hosts
from janus.objectives import mrna
from janus.objectives.mpnn import load_unconditional
from janus.rescore import FoldingWeights, Scales, pool_scales, rescore
from janus.sample import shell_samples

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
GRID = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--k", type=int, default=200)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--limit", type=int, default=120)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not mrna.available():
        raise SystemExit("ViennaRNA is required")

    host = hosts.load(args.host)
    rng = np.random.default_rng(args.seed)
    # The exchange rate is a statement about designs, so the natural domains
    # sharing this directory are excluded rather than silently averaged in.
    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)][: args.limit]

    rows = []
    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        pool = design(marginals, host, weights=TIER1, delta=args.delta, k=args.k)
        pool += shell_samples(marginals, host, weights=TIER1, delta=args.delta,
                              count=args.samples, rng=rng)

        scales = pool_scales(pool, host)
        reference = max(pool, key=lambda d: d.score)

        for lam in GRID:
            ranked = rescore(pool, host, FoldingWeights(initiation=lam), scales)
            best = ranked[0]
            clean = [r for r in ranked if r.design.synthesisable]
            rows.append({
                "backbone": path.stem,
                "length": len(marginals),
                "lam": lam,
                "tier1": best.tier1,
                "initiation": best.initiation_energy,
                "tier1_reference": reference.score,
                "initiation_reference": mrna.initiation_energy(reference.cds, host),
                "scale_tier1": scales.tier1,
                "scale_initiation": scales.initiation,
                "synthesisable": bool(best.design.synthesisable),
                "clean_fraction": len(clean) / len(ranked),
                "identity_to_reference": sum(
                    a == b for a, b in zip(best.protein, reference.protein, strict=True)
                ) / len(reference.protein),
            })
        if number % 20 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows)
    print(f"\nwrote {args.out} ({len(rows)} rows)")


def summarise(rows):
    backbones = sorted({r["backbone"] for r in rows})
    print(f"\n=== pool scales, {len(backbones)} backbones ===")
    at_zero = [r for r in rows if r["lam"] == 0.0]
    print(f"  Tier-1 spread across the pool:      median {np.median([r['scale_tier1'] for r in at_zero]):.3f} nats")
    print(f"  initiation spread across the pool:  median {np.median([r['scale_initiation'] for r in at_zero]):.3f} kcal/mol")
    ratio = np.median([r["scale_tier1"] / r["scale_initiation"] for r in at_zero])
    print(f"  ratio: one kcal/mol of initiation is worth {ratio:.2f} nats of Tier-1")
    print("  an unnormalised weight of 1.0 was therefore implicitly weighting")
    print(f"  initiation at about {1 / ratio:.2f} of its normalised value")

    print("\n=== frontier ===")
    print("  medians over backbones, since what the shell can open is right-skewed")
    print(f"{'weight':>8}{'Tier-1 given up':>18}{'initiation gained':>19}"
          f"{'exchange':>11}{'identity':>10}{'selected':>10}")
    print(f"{'':>8}{'nats':>18}{'kcal/mol':>19}{'kcal/nat':>11}{'':>10}{'clean %':>10}")
    previous = None
    curve = []
    for lam in sorted({r["lam"] for r in rows}):
        subset = [r for r in rows if r["lam"] == lam]
        # Medians, not means. What the shell can open is right-skewed, so a mean
        # sits above the typical backbone and drifts from the same quantity
        # reported anywhere else.
        cost = float(np.median([r["tier1_reference"] - r["tier1"] for r in subset]))
        gain = float(np.median([r["initiation"] - r["initiation_reference"] for r in subset]))
        identity = float(np.median([r["identity_to_reference"] for r in subset]))
        # Fraction of the pool that is clean does not depend on the weight;
        # what matters is whether the candidate actually selected is clean.
        clean = 100 * float(np.mean([r["synthesisable"] for r in subset]))
        rate = ""
        if previous and cost - previous[0] > 1e-9:
            rate = f"{(gain - previous[1]) / (cost - previous[0]):.2f}"
        curve.append({"lam": lam, "cost": cost, "gain": gain,
                      "identity": identity, "clean": clean})
        print(f"{lam:>8.3f}{cost:>18.4f}{gain:>19.4f}{rate:>11}{identity:>10.3f}{clean:>10.1f}")
        previous = (cost, gain)

    positive = [c for c in curve if c["gain"] > 0]
    if positive:
        best = max(positive, key=lambda c: c["gain"] / max(c["cost"], 1e-9))
        print(f"\n  best exchange rate at weight {best['lam']}: "
              f"{best['gain']:.3f} kcal/mol for {best['cost']:.4f} nats "
              f"= {best['gain'] / max(best['cost'], 1e-9):.1f} kcal/mol per nat")
    knee = None
    for i in range(1, len(curve)):
        if curve[i]["cost"] > 0 and curve[i]["gain"] > 0:
            marginal = ((curve[i]["gain"] - curve[i - 1]["gain"])
                        / max(curve[i]["cost"] - curve[i - 1]["cost"], 1e-9))
            if knee is None or marginal < knee[1]:
                knee = (curve[i]["lam"], marginal)
    print(f"  diminishing returns set in by weight {knee[0]}" if knee else "")


if __name__ == "__main__":
    main()
