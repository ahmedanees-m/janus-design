"""Every arm at a matched budget of full-objective evaluations.

A sampler given one draw and a search given thousands are not matched in
expenditure, so a side-by-side table reflects the budget as much as the method.
This scores a batch of draws from each stochastic model under the objective the
search arms optimise and reports the best of the first n for a ladder of n,
which places every arm on one axis: objective attained against evaluations of
the full objective spent.

The codon layer is solved exactly for every draw, so a sampled residue sequence
is scored at its own best coding sequence rather than at an arbitrary one.
"""

from __future__ import annotations

import argparse
import csv
import json
import zlib
from pathlib import Path

import numpy as np
from baselines import PDBLIKE, TIER1, make_objective, read_backbone

from janus import design, hosts
from janus.objectives import mrna
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights
from janus.sample import ShellSearch

LADDER = (1, 2, 4, 8, 16, 32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--samples", nargs="+", required=True,
                        help="arm=path pairs of sampled residue sequences")
    parser.add_argument("--search", required=True,
                        help="CSV of search-arm sequences from redesign_sequences.py")
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--folding", type=float, default=0.25)
    parser.add_argument("--penalty", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    weights = FoldingWeights(initiation=args.folding, liability=args.penalty)

    drawn = {}
    for item in args.samples:
        arm, _, path = item.partition("=")
        drawn[arm] = json.loads(Path(path).read_text(encoding="utf-8"))

    search_rows = {}
    with Path(args.search).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            search_rows.setdefault(row["backbone"], {})[row["arm"]] = row["protein"]

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    if args.shards > 1:
        files = files[args.shard::args.shards]
    print(f"{len(files)} designed backbones, {len(drawn)} sampled arms", flush=True)

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

        rng = np.random.default_rng(args.seed + zlib.crc32(name.encode()))
        total, burden = make_objective(marginals, host, rsa, classes, args.delta,
                                       native, args.pool, weights, rng)
        search = ShellSearch(marginals, host, TIER1, args.delta, anchor=native)
        pinned = design(marginals, host, weights=TIER1, delta=args.delta, k=1,
                        fixed=dict(enumerate(native)), anchor=native)[0]
        reference = total(pinned)

        def solve(protein, marginals=marginals, native=native, search=search):
            if len(protein) != len(native) or "X" in protein:
                return None
            if any(r not in search.slots[i] for i, r in enumerate(protein)):
                return design(marginals, host, weights=TIER1, delta=args.delta,
                              k=1, fixed=dict(enumerate(protein)),
                              anchor=protein)[0]
            return search.best(protein)

        def emit(arm, temperature, draws, candidate, evaluations, name=name,
                 native=native, total=total, reference=reference,
                 burden=burden):
            rows.append({
                "backbone": name, "arm": arm, "temperature": temperature,
                "draws": draws, "length": len(native),
                "objective": total(candidate),
                "gain_over_codon_dp": total(candidate) - reference,
                "tier1": candidate.score,
                "burden": burden(candidate.protein),
                "violations": len(candidate.violations),
                "initiation": mrna.initiation_energy(candidate.cds, host),
                "identity_to_design": sum(
                    a == b for a, b in zip(candidate.protein, native, strict=True)
                ) / len(native),
                "evaluations": evaluations,
            })

        for arm, table in drawn.items():
            entry = table.get(name)
            if entry is None:
                continue
            temperature = entry.get("temperature")
            scored = []
            for protein in entry["proteins"]:
                candidate = solve(protein)
                if candidate is not None:
                    scored.append((total(candidate), candidate))
            if not scored:
                continue
            for n in LADDER:
                if n > len(scored):
                    continue
                best = max(scored[:n], key=lambda item: item[0])[1]
                emit(arm, temperature, n, best, n)

        for arm, evaluations in (("design", 1), ("coordinate_descent", 437),
                                 ("joint", 3121)):
            protein = search_rows.get(name, {}).get(arm)
            if protein is None:
                continue
            candidate = pinned if arm == "design" else solve(protein)
            if candidate is not None:
                emit(arm, None, 1, candidate, evaluations)

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
