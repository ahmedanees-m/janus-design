"""Coverage of the Tier-2 optimum by the Tier-1 k-best prefix.

Rescoring can only reorder what the parse returned, so what the two-tier design
rests on is whether the prefix already contains the candidates the second tier
would choose. This scores the Tier-1 k-best prefix and a large uniform draw from
the same shell under the full objective, and reports where the best candidate in
the union sits in the Tier-1 ordering, how much of the Tier-2 top is inside the
prefix at each depth, and how both move with the folding weight.
"""

from __future__ import annotations

import argparse
import csv
import zlib
from pathlib import Path

import numpy as np
from baselines import PDBLIKE, TIER1, make_objective, read_backbone

from janus import design, hosts
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights
from janus.sample import shell_samples

DEPTHS = (10, 50, 100, 500)
TOPS = (1, 10, 50)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--prefix", type=int, default=500)
    parser.add_argument("--draws", type=int, default=5000)
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--foldings", nargs="+", type=float,
                        default=[0.125, 0.25, 0.5, 1.0])
    parser.add_argument("--penalty", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int, default=120)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    if args.shards > 1:
        files = files[args.shard::args.shards]
    print(f"{len(files)} backbones, prefix {args.prefix}, {args.draws} draws",
          flush=True)

    rows = []
    for number, path in enumerate(files, start=1):
        name = path.stem
        marginals = load_unconditional(path)
        try:
            native, rsa = read_backbone(Path(args.backbones) / f"{name}.pdb")
        except Exception:
            continue
        if len(rsa) != len(marginals) or "X" in native:
            continue

        seed = args.seed + zlib.crc32(name.encode())
        prefix = design(marginals, host, weights=TIER1, delta=args.delta,
                        k=args.prefix, anchor=native)
        drawn = shell_samples(marginals, host, weights=TIER1, delta=args.delta,
                              count=args.draws,
                              rng=np.random.default_rng(seed + 1))
        union = prefix + drawn
        order = np.argsort([-c.score for c in union])

        for folding in args.foldings:
            weights = FoldingWeights(initiation=folding, liability=args.penalty)
            total, _ = make_objective(marginals, host, rsa, classes, args.delta,
                                      native, args.pool, weights,
                                      np.random.default_rng(seed))
            values = np.array([total(c) for c in union])
            ranked = np.argsort(-values)
            tier1_rank = {int(index): position
                          for position, index in enumerate(order)}
            best = int(ranked[0])
            row = {"backbone": name, "length": len(native), "folding": folding,
                   "prefix": len(prefix), "draws": len(drawn),
                   "best_in_prefix": float(values[:len(prefix)].max()),
                   "best_in_union": float(values[best]),
                   "best_tier1_rank": tier1_rank[best] + 1,
                   "best_from_prefix": bool(best < len(prefix))}
            row["shortfall_of_prefix"] = row["best_in_union"] - row["best_in_prefix"]
            for top in TOPS:
                winners = set(int(i) for i in ranked[:top])
                for depth in DEPTHS:
                    if depth > len(prefix):
                        continue
                    kept = set(int(i) for i in order[:depth])
                    row[f"recall_top{top}_at{depth}"] = len(winners & kept) / top
            rows.append(row)

        if number % 5 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
