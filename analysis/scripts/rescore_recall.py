"""Does the parser's ordering help find good folding candidates, or would
sampling the same shell do as well?

mRNA folding cannot enter the parser: base pairs form at arbitrary range, so a
path's folding energy is not a sum over its positions. The two-tier design has
the parser propose candidates and folding rank them. That is only worth doing if
a Tier-1 prefix of size N contains better folding candidates than N draws from
the same shell, so both are given the same number of folding evaluations.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np

from janus import hosts
from janus.design import design
from janus.genetic_code import AA_ALPHABET, CODON_TO_AA
from janus.lattice import build
from janus.objectives import Weights, assemble, mrna
from janus.objectives.mpnn import load_unconditional
from janus.parse import kbest

AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}
PREFIXES = [1, 5, 10, 25, 50, 100, 250, 500, 1000]
TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def combined(tier1, energy, lam, scale_tier1=1.0, scale_initiation=1.0):
    """Tier-1 score plus a weighted initiation term, optionally on the spread scale.

    With unit scales this is the raw combination the first run used, where lambda
    silently carries the nats-per-kcal conversion and so does not transfer between
    backbones. Dividing each term by its spread over the ranked pool makes lambda
    dimensionless and puts it on the same axis as the folding sweep.
    """
    return tier1 / scale_tier1 + lam * energy / scale_initiation


def sample_candidates(lattice, marginals, node, edge, host, count, temperature, rng):
    """Draw residue sequences from the shell and optimise their codons exactly."""
    table = []
    for codons in lattice.codons:
        by_residue = {}
        for index, codon in enumerate(codons):
            by_residue.setdefault(CODON_TO_AA[codon], []).append(index)
        table.append(by_residue)

    weights = []
    for position, admitted in enumerate(lattice.amino_acids):
        logits = np.array([marginals[position][AA_INDEX[a]] for a in admitted]) / temperature
        logits -= logits.max()
        probability = np.exp(logits)
        weights.append(probability / probability.sum())

    out = []
    for _ in range(count):
        keep = []
        for position, admitted in enumerate(lattice.amino_acids):
            choice = admitted[rng.choice(len(admitted), p=weights[position])]
            keep.append(table[position][choice])
        sub_node = [node[i][keep[i]] for i in range(len(node))]
        sub_edge = [edge[i][np.ix_(keep[i], keep[i + 1])] for i in range(len(edge))]
        path = kbest(sub_node, sub_edge, k=1)[0]
        cds = "".join(lattice.codons[i][keep[i][s]] for i, s in enumerate(path.states))
        out.append((path.score, cds))
    return out


def run(marginal_dir, host, delta, limit, budget, lambdas, seed, normalised=False):
    rng = np.random.default_rng(seed)
    files = [p for p in sorted(Path(marginal_dir).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if limit:
        files = files[:limit]
    rows = []

    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        lattice = build(marginals, delta=delta)
        node, edge = assemble(lattice, host, marginals, TIER1)

        start = time.perf_counter()
        ranked = design(marginals, host, weights=TIER1, delta=delta, k=budget)
        parse_seconds = time.perf_counter() - start
        drawn = sample_candidates(lattice, marginals, node, edge, host, budget, 0.1, rng)

        # Folding is the cost here, so evaluate each candidate once and reuse the
        # energies across the whole lambda grid.
        parser_tier1 = np.array([d.score for d in ranked])
        sample_tier1 = np.array([s for s, _ in drawn])
        parser_energy = np.array([mrna.initiation_energy(d.cds, host) for d in ranked])
        sample_energy = np.array([mrna.initiation_energy(c, host) for _, c in drawn])

        if normalised:
            pooled_tier1 = np.concatenate([parser_tier1, sample_tier1])
            pooled_energy = np.concatenate([parser_energy, sample_energy])
            scale_tier1 = max(float(pooled_tier1.std()), 1e-9)
            scale_initiation = max(float(pooled_energy.std()), 1e-9)
        else:
            scale_tier1 = scale_initiation = 1.0

        for lam in lambdas:
            parser_scores = combined(parser_tier1, parser_energy, lam,
                                     scale_tier1, scale_initiation)
            sample_scores = combined(sample_tier1, sample_energy, lam,
                                     scale_tier1, scale_initiation)

            best_index = int(np.argmax(parser_scores))
            row = {
                "backbone": path.stem,
                "length": len(marginals),
                "lambda_initiation": lam,
                "normalised": normalised,
                "scale_tier1": scale_tier1,
                "scale_initiation": scale_initiation,
                "budget": budget,
                "parse_seconds": parse_seconds,
                "tier1_top": float(parser_tier1[0]),
                "tier1_of_winner": float(parser_tier1[best_index]),
                "initiation_top": float(parser_energy[0]),
                "initiation_winner": float(parser_energy[best_index]),
                "winner_rank": best_index + 1,
                "parser_best": float(parser_scores.max()),
                "sample_best": float(sample_scores.max()),
            }
            for prefix in PREFIXES:
                if prefix <= len(parser_scores):
                    row[f"parser_at_{prefix}"] = float(parser_scores[:prefix].max())
                if prefix <= len(sample_scores):
                    row[f"sample_at_{prefix}"] = float(sample_scores[:prefix].max())
            rows.append(row)

        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)
    return rows


def summarise(rows, lambdas):
    if rows and rows[0]["normalised"]:
        ratio = np.array([r["scale_tier1"] / r["scale_initiation"] for r in rows])
        print(f"spread ratio tier1/initiation: median {np.median(ratio):.4f}, "
              f"10th {np.percentile(ratio, 10):.4f}, 90th {np.percentile(ratio, 90):.4f}")
        print("  a raw lambda of L sits at L / that ratio on this axis, so the "
              f"earlier 0.1 / 0.3 / 1.0 grid lands near "
              f"{0.1 / np.median(ratio):.2f} / {0.3 / np.median(ratio):.2f} / "
              f"{1.0 / np.median(ratio):.2f}")

    for lam in lambdas:
        subset = [r for r in rows if r["lambda_initiation"] == lam]
        if not subset:
            continue
        print(f"\n=== lambda_initiation = {lam} ({len(subset)} backbones) ===")
        print(f"{'budget':>8} {'parser':>14} {'sampling':>14} {'parser wins':>13}")
        for prefix in PREFIXES:
            key_p, key_s = f"parser_at_{prefix}", f"sample_at_{prefix}"
            have = [r for r in subset if key_p in r and key_s in r]
            if not have:
                continue
            parser = np.array([r[key_p] for r in have])
            sampled = np.array([r[key_s] for r in have])
            print(f"{prefix:>8} {np.median(parser):>14.4f} {np.median(sampled):>14.4f} "
                  f"{100 * np.mean(parser > sampled):>12.1f}%")

        ranks = np.array([r["winner_rank"] for r in subset])
        print(f"\n  rank of the folding winner in the Tier-1 list: "
              f"median {int(np.median(ranks))}, 90th pct {int(np.percentile(ranks, 90))}, "
              f"max {int(ranks.max())}")
        print(f"  winner is the Tier-1 top path on {100 * np.mean(ranks == 1):.1f}% of backbones")

        gained = np.array([r["initiation_winner"] - r["initiation_top"] for r in subset])
        paid = np.array([r["tier1_top"] - r["tier1_of_winner"] for r in subset])
        print(f"  initiation energy opened by: median {np.median(gained):+.3f} kcal/mol, "
              f"90th pct {np.percentile(gained, 90):+.3f}")
        print(f"  Tier-1 score paid for it   : median {np.median(paid):.4f} nats, "
              f"90th pct {np.percentile(paid, 90):.4f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--budget", type=int, default=1000)
    parser.add_argument("--lambdas", type=float, nargs="+", default=[0.1, 0.3, 1.0])
    parser.add_argument("--normalised", action="store_true",
                        help="read lambda on the pool-spread scale, as rescore does")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not mrna.available():
        raise SystemExit("ViennaRNA is required for this experiment")

    host = hosts.load(args.host)
    rows = run(args.marginals, host, args.delta, args.limit,
               args.budget, args.lambdas, args.seed, args.normalised)
    summarise(rows, args.lambdas)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows), encoding="utf-8")
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
