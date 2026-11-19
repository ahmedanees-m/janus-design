"""How the result moves when the objective weights move.

The weights express a preference between terms that measure different things, so
a result that only holds at one setting is a result about the setting. This runs
the shipped search under alternative weightings, one term at a time and then at
random perturbations of all of them, and records what changes and what does not.

Three statements are checked at every scheme rather than argued: the codon-only
arm removes no protein-level liability, the joint arm attains more of the
objective than the fixed-protein codon dynamic program, and synthesis violations
fall. The rest is reported as a spread.
"""

from __future__ import annotations

import argparse
import csv
import zlib
from pathlib import Path

import numpy as np

from baselines import (LIABILITIES, PDBLIKE, coordinate_descent,
                       make_objective, read_backbone)
from janus import Weights, design, hosts
from janus.objectives import liability, mrna
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights
from janus.sample import ShellSearch

OPERATING = {"mpnn": 1.0, "cai": 0.5, "cpb": 0.3,
             "initiation": 0.25, "liability": 0.25}


def schemes(count, seed):
    out = [("operating", dict(OPERATING))]
    for name in OPERATING:
        for factor in (0.5, 2.0):
            scheme = dict(OPERATING)
            scheme[name] = OPERATING[name] * factor
            out.append((f"{name} x{factor:g}", scheme))
    out.append(("equal", {k: 1.0 for k in OPERATING}))
    out.append(("initiation off", {**OPERATING, "initiation": 0.0}))
    out.append(("liability off", {**OPERATING, "liability": 0.0}))
    rng = np.random.default_rng(seed)
    for i in range(count):
        factors = np.exp(rng.uniform(np.log(0.5), np.log(2.0), len(OPERATING)))
        out.append((f"random {i + 1}",
                    {k: OPERATING[k] * float(f)
                     for k, f in zip(OPERATING, factors)}))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--random-schemes", type=int, default=12)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int, default=120)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    grid = schemes(args.random_schemes, args.seed)

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    if args.shards > 1:
        files = files[args.shard::args.shards]
    print(f"{len(files)} backbones, {len(grid)} weight schemes", flush=True)

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

        for label, scheme in grid:
            tier1 = Weights(mpnn=scheme["mpnn"], cai=scheme["cai"],
                            cpb=scheme["cpb"])
            weights = FoldingWeights(initiation=scheme["initiation"],
                                     liability=scheme["liability"])
            rng = np.random.default_rng(args.seed + zlib.crc32(name.encode()))
            total, burden = make_objective(marginals, host, rsa, classes,
                                           args.delta, native, args.pool,
                                           weights, rng, tier1=tier1)
            search = ShellSearch(marginals, host, tier1, args.delta,
                                 anchor=native)
            pinned = design(marginals, host, weights=tier1, delta=args.delta,
                            k=1, fixed=dict(enumerate(native)), anchor=native)[0]
            moved, _, evaluations = coordinate_descent(search, pinned, total,
                                                       marginals, args.rounds)
            removed = {
                term: liability.score(term, pinned.protein, rsa, classes)
                - liability.score(term, moved.protein, rsa, classes)
                for term in LIABILITIES
            }
            rows.append({
                "backbone": name, "scheme": label, "length": len(native),
                **{f"w_{k}": v for k, v in scheme.items()},
                "objective_codon_dp": total(pinned),
                "objective_joint": total(moved),
                "gain": total(moved) - total(pinned),
                "burden_codon_dp": burden(pinned.protein),
                "burden_joint": burden(moved.protein),
                "violations_codon_dp": len(pinned.violations),
                "violations_joint": len(moved.violations),
                "initiation_codon_dp": mrna.initiation_energy(pinned.cds, host),
                "initiation_joint": mrna.initiation_energy(moved.cds, host),
                "identity_to_design": sum(a == b for a, b in
                                          zip(moved.protein, native)) / len(native),
                "residues_moved": sum(a != b for a, b in
                                      zip(moved.protein, pinned.protein)),
                "evaluations": evaluations,
                **{f"removed_{k}": v for k, v in removed.items()},
            })

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
