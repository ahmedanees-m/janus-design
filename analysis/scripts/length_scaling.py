"""Parse cost and attainable gain as a function of chain length.

The lattice is linear in length and the measurements in this work are made on
miniproteins, so the two have to be separated: what scales is a property of the
algorithm, and what was measured is a property of the set. This times the parse
and the k-best extraction over a length ladder, records the shell the budget
admits at each length, and runs the shipped search so the gain over the
fixed-protein codon dynamic program can be read against length as well.

Timing is the median of repeated solves on a warm interpreter, so it reports the
solve rather than the first call into it.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import time
import zlib
from pathlib import Path

import numpy as np
from baselines import PDBLIKE, TIER1, coordinate_descent, make_objective, read_backbone

from janus import design, hosts
from janus.lattice import amino_acid_shell
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights
from janus.sample import ShellSearch


def timed(call, repeats):
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        call()
        samples.append(time.perf_counter() - started)
    return statistics.median(samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--set", required=True, help="label for the backbone set")
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--deltas", nargs="+", type=float, default=[0.5, 1.0])
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--kbest", type=int, default=50)
    parser.add_argument("--folding", type=float, default=0.25)
    parser.add_argument("--penalty", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-length", type=int,
                        help="skip chains longer than this")
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--search", action="store_true",
                        help="also run the residue search and record the gain")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm) if args.search else None
    weights = FoldingWeights(initiation=args.folding, liability=args.penalty)

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if args.set != "design" or not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    if args.shards > 1:
        files = files[args.shard::args.shards]
    print(f"{len(files)} backbones in set {args.set}", flush=True)

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
        if args.max_length and len(native) > args.max_length:
            continue

        entropy = float(sum(-(np.exp(row) * row).sum() for row in marginals))
        # One objective for every budget on this backbone. The scales are
        # measured on the widest pool, so gains at different budgets are on one
        # scale and can be compared.
        common = None
        if args.search:
            rng = np.random.default_rng(args.seed + zlib.crc32(name.encode()))
            common = make_objective(marginals, host, rsa, classes,
                                    max(args.deltas), native, args.pool,
                                    weights, rng)
        for delta in args.deltas:
            shell = amino_acid_shell(marginals, delta)
            admitted = statistics.mean(len(s) for s in shell)
            best = timed(lambda m=marginals, d=delta: design(
                m, host, weights=TIER1, delta=d, k=1), args.repeats)
            kbest = timed(lambda m=marginals, d=delta: design(
                m, host, weights=TIER1, delta=d, k=args.kbest), args.repeats)
            row = {"set": args.set, "backbone": name, "length": len(native),
                   "delta": delta, "total_entropy": entropy,
                   "mean_shell": admitted,
                   "parse_ms": best * 1000, "kbest_ms": kbest * 1000,
                   "k": args.kbest}
            if args.search:
                total, burden = common
                search = ShellSearch(marginals, host, TIER1, delta, anchor=native)
                pinned = design(marginals, host, weights=TIER1, delta=delta, k=1,
                                fixed=dict(enumerate(native)), anchor=native)[0]
                started = time.perf_counter()
                moved, _, evaluations = coordinate_descent(search, pinned, total,
                                                           marginals, args.rounds)
                row.update({
                    "objective_codon_dp": total(pinned),
                    "objective_joint": total(moved),
                    "gain": total(moved) - total(pinned),
                    "burden_codon_dp": burden(pinned.protein),
                    "burden_joint": burden(moved.protein),
                    "violations_codon_dp": len(pinned.violations),
                    "violations_joint": len(moved.violations),
                    "residues_moved": sum(a != b for a, b in
                                          zip(moved.protein, pinned.protein, strict=True)),
                    "identity_to_native": sum(
                        a == b for a, b in zip(moved.protein, native, strict=True)
                    ) / len(native),
                    "search_seconds": time.perf_counter() - started,
                    "evaluations": evaluations,
                })
            rows.append(row)

        if number % 10 == 0:
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
