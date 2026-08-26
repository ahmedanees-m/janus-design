"""Which layer buys the initiation term: residues or codons?

At the defensible operating weight the selected candidate is 98.7 percent
identical to the Tier-1 optimum in amino acids, which on a 43-residue design is
about half a substitution. If the gain is coming almost entirely from synonymous
codon changes then the amino-acid axis, which is what distinguishes this work
from fixed-protein codon optimisation, is contributing little at the operating
point. That has to be measured rather than assumed.

Three arms, all ranked under the same normalised objective:

- ``synonymous``: the protein is pinned to the Tier-1 optimum and only codons
  move, which is what a fixed-protein codon optimiser can reach
- ``residue``: residues move and each candidate's codons are the exact Tier-1
  optimum for that protein, so any change in the initiation term is attributable
  to the residue choice
- ``joint``: both move, which is the full method
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from janus import Weights, design, hosts
from janus.genetic_code import CODON_TO_AA
from janus.objectives import mrna
from janus.objectives.mpnn import load_unconditional
from janus.rescore import FoldingWeights, pool_scales, rescore
from janus.sample import shell_samples

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
WEIGHTS = [0.125, 0.5]

PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def codon_changes(a, b):
    return sum(
        a[i : i + 3] != b[i : i + 3] for i in range(0, min(len(a), len(b)), 3)
    )


def residue_changes(a, b):
    return sum(x != y for x, y in zip(a, b, strict=False))


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
    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)][: args.limit]

    rows = []
    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)

        prefix = design(marginals, host, weights=TIER1, delta=args.delta, k=args.k)
        reference = prefix[0]
        pinned = {i: r for i, r in enumerate(reference.protein)}

        # Codons only: the protein is held at the Tier-1 optimum.
        synonymous = design(marginals, host, weights=TIER1, delta=args.delta,
                            k=args.k, fixed=pinned)
        # Residues only: codons are the exact Tier-1 optimum for each protein.
        residue = shell_samples(marginals, host, weights=TIER1, delta=args.delta,
                                count=args.samples, rng=rng)
        joint = prefix + residue

        scales = pool_scales(joint, host)
        base_initiation = mrna.initiation_energy(reference.cds, host)

        for lam in WEIGHTS:
            for arm, pool in (("synonymous", synonymous),
                              ("residue", residue),
                              ("joint", joint)):
                best = rescore(pool, host, FoldingWeights(initiation=lam), scales)[0]
                rows.append({
                    "backbone": path.stem,
                    "length": len(marginals),
                    "lam": lam,
                    "arm": arm,
                    "initiation_gain": best.initiation_energy - base_initiation,
                    "tier1_cost": reference.score - best.tier1,
                    "residue_changes": residue_changes(best.protein, reference.protein),
                    "codon_changes": codon_changes(best.cds, reference.cds),
                    "synthesisable": bool(best.design.synthesisable),
                })
        if number % 20 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows)
    print(f"\nwrote {args.out} ({len(rows)} rows)")


def summarise(rows):
    lengths = {r["backbone"]: r["length"] for r in rows}
    print(f"\n=== decomposition over {len(lengths)} backbones, "
          f"median length {int(np.median(list(lengths.values())))} residues ===")

    for lam in sorted({r["lam"] for r in rows}):
        print(f"\n--- folding weight {lam} ---")
        print(f"{'arm':<14}{'initiation gained':>19}{'Tier-1 cost':>14}"
              f"{'residue changes':>17}{'codon changes':>15}{'clean':>8}")
        joint_gain = None
        for arm in ("synonymous", "residue", "joint"):
            subset = [r for r in rows if r["lam"] == lam and r["arm"] == arm]
            if not subset:
                continue
            gain = float(np.mean([r["initiation_gain"] for r in subset]))
            cost = float(np.mean([r["tier1_cost"] for r in subset]))
            res = float(np.mean([r["residue_changes"] for r in subset]))
            cod = float(np.mean([r["codon_changes"] for r in subset]))
            clean = 100 * float(np.mean([r["synthesisable"] for r in subset]))
            if arm == "joint":
                joint_gain = gain
            print(f"{arm:<14}{gain:>19.4f}{cost:>14.4f}{res:>17.2f}{cod:>15.2f}{clean:>7.1f}%")

        if joint_gain:
            syn = float(np.mean([r["initiation_gain"] for r in rows
                                 if r["lam"] == lam and r["arm"] == "synonymous"]))
            res = float(np.mean([r["initiation_gain"] for r in rows
                                 if r["lam"] == lam and r["arm"] == "residue"]))
            print(f"\n  synonymous alone recovers {100 * syn / joint_gain:.1f}% "
                  f"of the joint gain")
            print(f"  residue alone recovers     {100 * res / joint_gain:.1f}%")

    print("\n=== how often the joint arm actually changes a residue ===")
    for lam in sorted({r["lam"] for r in rows}):
        subset = [r for r in rows if r["lam"] == lam and r["arm"] == "joint"]
        changed = [r["residue_changes"] for r in subset]
        print(f"  weight {lam}: {100 * np.mean([c > 0 for c in changed]):.1f}% of backbones, "
              f"median {int(np.median(changed))} substitutions, max {int(np.max(changed))}")

    print("\n=== enumeration check ===")
    median_length = int(np.median(list(lengths.values())))
    print(f"  the initiation window covers codons 0 to 12, so exhaustive search over")
    print(f"  one substitution there is 13 positions x 19 alternatives x their codons,")
    print(f"  on the order of a thousand candidates. The parser is not needed to find")
    print(f"  that; it is needed to know the exact optimum the substitution is measured")
    print(f"  against, over a space of roughly 3.5^{median_length} residue assignments.")


if __name__ == "__main__":
    main()
