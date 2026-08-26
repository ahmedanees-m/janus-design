"""Write the top-K residue sequences by summed unconditional marginal.

The parser optimises this sum, so this is its ranking by construction. Whatever
reorders it under conditional scoring is the surrogate gap.

The lattice used here carries one state per admitted residue rather than one per
codon. Over the full codon lattice the k-best paths are synonymous variants of a
single protein, which is the right answer for gene design and the wrong
candidate set for a question about residue choice.
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


def top_proteins(marginals, delta, k):
    shells = amino_acid_shell(marginals, delta)
    node = [
        np.array([marginals[i][AA_INDEX[a]] for a in shell], dtype=np.float64)
        for i, shell in enumerate(shells)
    ]
    edge = [np.zeros((len(node[i]), len(node[i + 1]))) for i in range(len(node) - 1)]
    out = []
    for path in kbest(node, edge, k=k):
        protein = "".join(shells[i][s] for i, s in enumerate(path.states))
        out.append((protein, path.score))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--k", type=int, default=500)
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = Path(args.out)
    (out / "fasta").mkdir(parents=True, exist_ok=True)

    files = sorted(Path(args.marginals).glob("*.npz"))[: args.limit]
    index = {}
    for number, path in enumerate(files, start=1):
        marginals = load_unconditional(path)
        ranked = top_proteins(marginals, args.delta, args.k)

        name = path.stem
        lines = []
        for i, (protein, _) in enumerate(ranked, start=1):
            lines.append(f">{name}_{i}")
            lines.append(protein)
        (out / "fasta" / f"{name}.fa").write_text("\n".join(lines) + "\n", encoding="utf-8")

        index[name] = {
            "length": len(marginals),
            "n_sequences": len(ranked),
            "unconditional": [score for _, score in ranked],
            "identity_to_top": [
                sum(a == b for a, b in zip(protein, ranked[0][0], strict=True)) / len(protein)
                for protein, _ in ranked
            ],
        }
        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    (out / "index.json").write_text(json.dumps(index), encoding="utf-8")
    total = sum(v["n_sequences"] for v in index.values())
    spans = [max(v["unconditional"]) - min(v["unconditional"]) for v in index.values()]
    print(f"{len(index)} backbones, {total} sequences to score")
    print(f"unconditional score span across the top-K: median {np.median(spans):.3f} nats")
    ident = [np.mean(v["identity_to_top"]) for v in index.values()]
    print(f"mean identity to the top path: {np.mean(ident):.3f}")


if __name__ == "__main__":
    main()
