"""Objective terms and their assembly into lattice weights.

Only terms additive over nodes or over adjacent-codon edges can enter the DP.
Fold compatibility, codon adaptation and GC are position-local; codon-pair bias
is first-order. Wider-context terms (mRNA folding, sites spanning three codons,
repeats, windowed GC) do not decompose here and are handled elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..hosts import Host
from ..lattice import Lattice
from .codon import codon_pair_matrix, gc_content_term, relative_adaptiveness_term
from .mpnn import marginal_term

__all__ = ["Weights", "assemble", "score_terms"]


@dataclass(frozen=True)
class Weights:
    """Multipliers on the decomposable objective terms.

    ``gc`` is a linear bias: negative pushes AT-rich, positive GC-rich. The
    host's hard GC range is a constraint, not a term, and is checked afterwards.
    """

    mpnn: float = 1.0
    cai: float = 0.0
    cpb: float = 0.0
    gc: float = 0.0


def assemble(
    lattice: Lattice,
    host: Host,
    log_marginals: np.ndarray,
    weights: Weights,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Return per-position node scores and adjacent-position transition matrices."""
    node_scores = []
    for position, codons in enumerate(lattice.codons):
        scores = np.zeros(len(codons))
        if weights.mpnn:
            scores += weights.mpnn * marginal_term(log_marginals[position], codons)
        if weights.cai:
            scores += weights.cai * relative_adaptiveness_term(host, codons)
        if weights.gc:
            scores += weights.gc * gc_content_term(codons)
        node_scores.append(scores)

    edge_scores = []
    for position in range(len(lattice) - 1):
        left, right = lattice.codons[position], lattice.codons[position + 1]
        if weights.cpb:
            edge_scores.append(weights.cpb * codon_pair_matrix(host, left, right))
        else:
            edge_scores.append(np.zeros((len(left), len(right))))

    return node_scores, edge_scores


def score_terms(cds: str, host: Host, log_marginals: np.ndarray) -> dict[str, float]:
    """Report each decomposable term for a finished coding sequence, unweighted."""
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    adaptiveness = host.relative_adaptiveness
    pair_scores = host.codon_pair_scores

    log_w = [np.log(adaptiveness[c]) for c in codons if adaptiveness.get(c, 0) > 0]
    junctions = [pair_scores.get(pair, 0.0) for pair in zip(codons, codons[1:], strict=False)]

    return {
        "mpnn": float(sum(marginal_term(log_marginals[i], (c,))[0] for i, c in enumerate(codons))),
        "cai": float(np.exp(np.mean(log_w))) if log_w else 0.0,
        "cpb": float(np.mean(junctions)) if junctions else 0.0,
        "gc": float(sum(b in "GC" for b in cds) / len(cds)) if cds else 0.0,
    }
