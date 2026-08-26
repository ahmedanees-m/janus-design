"""Candidate sets spanning a wide range of unconditional score, as a control.

The k-best set answers what the parser's own ranking is worth inside the band it
searches. It cannot distinguish a surrogate that carries no information from one
whose information is exhausted within 0.3 nats. This draws sequences from across
the whole shell so the unconditional axis spans nats rather than tenths, which
is the comparison that separates those two readings.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from janus.genetic_code import AA_ALPHABET
from janus.lattice import amino_acid_shell
from janus.objectives.mpnn import load_unconditional
from janus.parse import kbest

AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}


def node_scores(marginals, shells):
    return [
        np.array([marginals[i][AA_INDEX[a]] for a in shell], dtype=np.float64)
        for i, shell in enumerate(shells)
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--top", type=int, default=20, help="sequences taken from the k-best head")
    parser.add_argument("--sampled", type=int, default=180, help="uniform draws from the shell")
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    (out / "fasta").mkdir(parents=True, exist_ok=True)

    files = sorted(Path(args.marginals).glob("*.npz"))[: args.limit]
    index = {}
    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        shells = amino_acid_shell(marginals, args.delta)
        node = node_scores(marginals, shells)
        edge = [np.zeros((len(node[i]), len(node[i + 1]))) for i in range(len(node) - 1)]

        collected = {}
        for item in kbest(node, edge, k=args.top):
            protein = "".join(shells[i][s] for i, s in enumerate(item.states))
            collected[protein] = item.score

        tries = 0
        while len(collected) < args.top + args.sampled and tries < 40 * args.sampled:
            tries += 1
            picked = [int(rng.integers(len(shell))) for shell in shells]
            protein = "".join(shells[i][s] for i, s in enumerate(picked))
            if protein in collected:
                continue
            collected[protein] = float(sum(node[i][s] for i, s in enumerate(picked)))

        proteins = list(collected)
        scores = [collected[p] for p in proteins]

        name = path.stem
        lines = []
        for i, protein in enumerate(proteins, start=1):
            lines.append(f">{name}_{i}")
            lines.append(protein)
        (out / "fasta" / f"{name}.fa").write_text("\n".join(lines) + "\n", encoding="utf-8")

        index[name] = {
            "length": len(marginals),
            "n_sequences": len(proteins),
            "unconditional": scores,
        }
        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    (out / "index.json").write_text(json.dumps(index), encoding="utf-8")
    spans = [max(v["unconditional"]) - min(v["unconditional"]) for v in index.values()]
    total = sum(v["n_sequences"] for v in index.values())
    print(f"{len(index)} backbones, {total} sequences")
    print(f"unconditional span: median {np.median(spans):.2f} nats, min {min(spans):.2f}")


if __name__ == "__main__":
    main()
