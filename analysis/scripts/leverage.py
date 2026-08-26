"""H2b: does design freedom sit where the gene-level objective can use it?

Leverage at position i is the range of a gene-level term reachable by varying
that position alone, with every other position held at the Tier-1 optimum. H2b
predicts leverage and entropy coincide over the N-terminal region. Entropy is
measured in the same run so the two profiles are directly comparable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from janus import Weights, design, hosts
from janus.genetic_code import CODON_TO_AA, SYNONYMOUS, gc_fraction
from janus.lattice import amino_acid_shell, build
from janus.objectives import mrna
from janus.objectives.mpnn import load_unconditional

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)


def substitute(cds, position, codon):
    return cds[: 3 * position] + codon + cds[3 * position + 3 :]


def position_leverage(cds, shells, host):
    """Reachable range of the initiation term and of GC, one position at a time."""
    initiation, gc = [], []
    for position, admitted in enumerate(shells):
        codons = [c for residue in admitted for c in SYNONYMOUS[residue]]
        values = [mrna.initiation_energy(substitute(cds, position, c), host) for c in codons]
        fractions = [gc_fraction(substitute(cds, position, c)) for c in codons]
        initiation.append(float(max(values) - min(values)))
        gc.append(float(max(fractions) - min(fractions)))
    return np.array(initiation), np.array(gc)


def permutation_null(entropy, leverage, window, rng, draws=2000):
    """Shuffle positions within the design, keeping the window fixed."""
    observed = spearmanr(entropy[:window], leverage[:window]).statistic
    null = np.empty(draws)
    order = np.arange(len(entropy))
    for i in range(draws):
        rng.shuffle(order)
        null[i] = spearmanr(entropy[order][:window], leverage[order][:window]).statistic
    return observed, null


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--limit", type=int, default=150)
    parser.add_argument("--window", type=int, default=15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not mrna.available():
        raise SystemExit("ViennaRNA is required")

    host = hosts.load(args.host)
    rng = np.random.default_rng(args.seed)
    files = sorted(Path(args.marginals).glob("*.npz"))[: args.limit]

    rows = []
    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        if len(marginals) <= args.window + 5:
            continue
        shells = amino_acid_shell(marginals, args.delta)
        entropy = -(np.exp(marginals) * marginals).sum(axis=1)
        top = design(marginals, host, weights=TIER1, delta=args.delta, k=1)[0]

        initiation, gc = position_leverage(top.cds, shells, host)
        observed, null = permutation_null(entropy, initiation, args.window, rng)

        rows.append({
            "backbone": path.stem,
            "length": len(marginals),
            "entropy": entropy.tolist(),
            "initiation_leverage": initiation.tolist(),
            "gc_leverage": gc.tolist(),
            "rho_window": None if np.isnan(observed) else float(observed),
            "null_mean": float(np.nanmean(null)),
            "null_sd": float(np.nanstd(null)),
            "rho_remainder": float(
                spearmanr(entropy[args.window :], initiation[args.window :]).statistic
            ) if not np.all(initiation[args.window :] == initiation[args.window]) else None,
        })
        if number % 25 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    Path(args.out).write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows, args.window)
    print(f"\nwrote {args.out} ({len(rows)} rows)")


def summarise(rows, window):
    length = min(r["length"] for r in rows)
    entropy = np.array([r["entropy"][:length] for r in rows])
    initiation = np.array([r["initiation_leverage"][:length] for r in rows])
    gc = np.array([r["gc_leverage"][:length] for r in rows])

    print(f"\n=== profiles along the chain, {len(rows)} backbones, first {length} positions ===")
    print(f"{'pos':>4} {'entropy':>9} {'5prime lev':>11} {'GC lev':>9}")
    for i in list(range(min(16, length))) + [length - 1]:
        print(f"{i:>4} {entropy[:, i].mean():>9.3f} {initiation[:, i].mean():>11.4f} "
              f"{gc[:, i].mean():>9.4f}")

    print("\n=== where the initiation term can be moved at all ===")
    movable = initiation.mean(axis=0) > 1e-6
    last = int(np.max(np.flatnonzero(movable))) if movable.any() else -1
    print(f"  last position with any initiation leverage: {last}")
    print(f"  mean leverage over positions 0 to {last}: "
          f"{initiation[:, : last + 1].mean():.4f} kcal/mol")
    print(f"  mean leverage beyond it: {initiation[:, last + 1 :].mean():.6f} kcal/mol")

    print("\n=== entropy, N-terminal window against the rest ===")
    head = entropy[:, :window].mean()
    tail = entropy[:, window:].mean()
    print(f"  positions 0 to {window - 1}: {head:.3f} nats")
    print(f"  beyond:                  {tail:.3f} nats")
    print(f"  difference:              {head - tail:+.3f} nats")

    print("\n=== H2b: entropy against leverage inside the window ===")
    rho = np.array([r["rho_window"] for r in rows if r["rho_window"] is not None])
    null_mean = np.array([r["null_mean"] for r in rows if r["rho_window"] is not None])
    null_sd = np.array([r["null_sd"] for r in rows if r["rho_window"] is not None])
    print(f"  observed rho: median {np.median(rho):+.3f}, "
          f"10th {np.percentile(rho, 10):+.3f}, 90th {np.percentile(rho, 90):+.3f}")
    print(f"  permutation null: mean {np.mean(null_mean):+.3f}, sd {np.mean(null_sd):.3f}")
    z = (rho - null_mean) / np.where(null_sd > 0, null_sd, np.nan)
    print(f"  z against the null: median {np.nanmedian(z):+.3f}")
    print(f"  backbones with rho above the null 95th percentile: "
          f"{100 * np.nanmean(z > 1.645):.1f}%")


if __name__ == "__main__":
    main()
