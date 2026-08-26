"""How much of the objective does cheap search leave on the table?

The parser gives the exact optimum; the same DP over negated weights gives the
worst attainable path. The two fix the range the objective spans on a backbone,
so heuristics are reported as a fraction of that range and backbones of
different lengths stay comparable.

Greedy, posterior sampling and annealing all optimise the identical function.
Shared setup is excluded from the budget; search time is what is matched.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from janus import hosts
from janus.genetic_code import AA_ALPHABET, CODON_TO_AA
from janus.lattice import build
from janus.objectives import Weights, assemble
from janus.objectives.mpnn import load_unconditional
from janus.parse import kbest

AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}

ACCUMULATING = {
    "mpnn": Weights(mpnn=1.0),
    "mpnn+cai": Weights(mpnn=1.0, cai=0.5),
    "mpnn+cai+gc": Weights(mpnn=1.0, cai=0.5, gc=0.5),
    "mpnn+cai+gc+cpb": Weights(mpnn=1.0, cai=0.5, gc=0.5, cpb=0.3),
}

# Only the codon-pair term couples adjacent positions; the rest are node-local.
# Sweeping its weight separates what the parser buys from how much the objective
# is coupled in the first place.
COUPLING = {
    f"cpb={value}": Weights(mpnn=1.0, cai=0.5, gc=0.5, cpb=value)
    for value in (0.0, 0.1, 0.3, 1.0, 3.0)
}

METHODS = ["greedy", "sample_x1", "sample_x10", "sample_x100",
           "anneal_x1", "anneal_x10", "anneal_x100"]


def path_value(node, edge, states):
    total = sum(node[i][s] for i, s in enumerate(states))
    total += sum(edge[i][states[i], states[i + 1]] for i in range(len(states) - 1))
    return float(total)


def best_path(node, edge):
    return kbest(node, edge, k=1)[0].score


def worst_path(node, edge):
    return -kbest([-x for x in node], [-x for x in edge], k=1)[0].score


def greedy(node, edge):
    """Left to right, taking the best next state given the one already chosen."""
    states = [int(np.argmax(node[0]))]
    for i in range(1, len(node)):
        states.append(int(np.argmax(node[i] + edge[i - 1][states[-1]])))
    return states


def residue_states(lattice):
    """Per position, the state indices belonging to each admitted residue."""
    table = []
    for codons in lattice.codons:
        by_residue = {}
        for index, codon in enumerate(codons):
            by_residue.setdefault(CODON_TO_AA[codon], []).append(index)
        table.append(by_residue)
    return table


def sample_best(node, edge, lattice, marginals, table, temperature, seconds, rng):
    """Posterior sampling with exact codons, the strongest form of Baseline 0.

    Residues are drawn from the tempered marginals within the shell, then codons
    for that fixed sequence are optimised exactly. Sampling codons too would be a
    weaker comparator.
    """
    weights = []
    for position, admitted in enumerate(lattice.amino_acids):
        logits = np.array([marginals[position][AA_INDEX[a]] for a in admitted]) / temperature
        logits -= logits.max()
        probability = np.exp(logits)
        weights.append(probability / probability.sum())

    best = -math.inf
    draws = 0
    deadline = time.perf_counter() + seconds
    while time.perf_counter() < deadline:
        keep = []
        for position, admitted in enumerate(lattice.amino_acids):
            choice = admitted[rng.choice(len(admitted), p=weights[position])]
            keep.append(table[position][choice])
        sub_node = [node[i][keep[i]] for i in range(len(node))]
        sub_edge = [edge[i][np.ix_(keep[i], keep[i + 1])] for i in range(len(edge))]
        best = max(best, best_path(sub_node, sub_edge))
        draws += 1
    return best, draws


def anneal(node, edge, seconds, rng):
    """Single-position moves scored incrementally, with geometric cooling."""
    states = [int(np.argmax(x)) for x in node]
    current = path_value(node, edge, states)
    best = current
    length = len(node)
    hot, cold = 1.0, 0.01
    steps = 0
    deadline = time.perf_counter() + seconds
    while time.perf_counter() < deadline:
        for _ in range(256):
            i = int(rng.integers(length))
            width = len(node[i])
            if width < 2:
                continue
            candidate = int(rng.integers(width))
            if candidate == states[i]:
                continue
            delta = node[i][candidate] - node[i][states[i]]
            if i > 0:
                delta += (edge[i - 1][states[i - 1], candidate]
                          - edge[i - 1][states[i - 1], states[i]])
            if i + 1 < length:
                delta += (edge[i][candidate, states[i + 1]]
                          - edge[i][states[i], states[i + 1]])
            temperature = max(cold, hot * (cold / hot) ** min(1.0, steps / 20000))
            if delta >= 0 or rng.random() < math.exp(delta / temperature):
                states[i] = candidate
                current += delta
                best = max(best, current)
            steps += 1
    return best, steps


def run(marginal_dir, host, delta, limit, seed, objectives):
    rng = np.random.default_rng(seed)
    files = sorted(Path(marginal_dir).glob("*.npz"))
    if limit:
        files = files[:limit]
    rows = []

    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        lattice = build(marginals, delta=delta)
        table = residue_states(lattice)

        for label, weights in objectives.items():
            node, edge = assemble(lattice, host, marginals, weights)

            start = time.perf_counter()
            optimum = best_path(node, edge)
            budget = time.perf_counter() - start
            floor = worst_path(node, edge)
            span = optimum - floor

            row = {
                "backbone": path.stem,
                "length": len(marginals),
                "objective": label,
                "terms": sum(1 for v in vars(weights).values() if v),
                "optimum": optimum,
                "worst": floor,
                "span": span,
                "janus_seconds": budget,
                "branching": lattice.branching,
                "greedy": path_value(node, edge, greedy(node, edge)),
            }
            for factor in (1, 10, 100):
                value, draws = sample_best(node, edge, lattice, marginals, table,
                                           0.1, budget * factor, rng)
                row[f"sample_x{factor}"] = value
                row[f"sample_x{factor}_draws"] = draws
                value, steps = anneal(node, edge, budget * factor, rng)
                row[f"anneal_x{factor}"] = value
                row[f"anneal_x{factor}_steps"] = steps

            for method in METHODS:
                row[f"shortfall_{method}"] = (optimum - row[method]) / span if span > 0 else 0.0
            rows.append(row)

        if number % 25 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)
    return rows


def summarise(rows, objectives):
    print("\n=== shortfall as a fraction of the attainable objective range ===")
    header = f"{'objective':<18} {'terms':>5}"
    for method in ["greedy", "sample x100", "anneal x1", "anneal x100"]:
        header += f" {method:>24}"
    print(header)
    for label in objectives:
        subset = [r for r in rows if r["objective"] == label]
        if not subset:
            continue
        line = f"{label:<18} {subset[0]['terms']:>5}"
        for method in ["greedy", "sample_x100", "anneal_x1", "anneal_x100"]:
            values = np.array([r[f"shortfall_{method}"] for r in subset])
            line += f" {np.median(values):11.6f} [{np.percentile(values, 90):9.6f}]"
        print(line)
    print("\n  median [90th percentile]; zero means the heuristic found the optimum")

    print("\n=== proportion of backbones on which each method is exactly optimal ===")
    header = f"{'objective':<18}"
    for method in METHODS:
        header += f" {method:>13}"
    print(header)
    for label in objectives:
        subset = [r for r in rows if r["objective"] == label]
        if not subset:
            continue
        line = f"{label:<18}"
        for method in METHODS:
            hit = np.mean([r[f"shortfall_{method}"] < 1e-9 for r in subset])
            line += f" {100 * hit:12.1f}%"
        print(line)

    print("\n=== budget actually used, on the last objective ===")
    subset = [r for r in rows if r["objective"] == list(objectives)[-1]]
    if subset:
        median_ms = 1000 * np.median([r["janus_seconds"] for r in subset])
        print(f"  parser solve time      : median {median_ms:.2f} ms")
        for factor in (1, 10, 100):
            draws = [r[f"sample_x{factor}_draws"] for r in subset]
            steps = [r[f"anneal_x{factor}_steps"] for r in subset]
            print(f"  at x{factor:<4}: median {int(np.median(draws)):>6} posterior draws, "
                  f"{int(np.median(steps)):>8} annealing steps")
        print(f"  mean branching factor  : {np.mean([r['branching'] for r in subset]):.2f}")
        print(f"  backbones              : {len(subset)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--mode", choices=["accumulating", "coupling"], default="accumulating",
                        help="add objectives one at a time, or sweep the coupling weight")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    objectives = ACCUMULATING if args.mode == "accumulating" else COUPLING
    host = hosts.load(args.host)
    rows = run(args.marginals, host, args.delta, args.limit, args.seed, objectives)
    summarise(rows, objectives)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows), encoding="utf-8")
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
