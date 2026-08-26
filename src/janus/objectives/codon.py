"""Codon adaptation, codon-pair bias and GC content as lattice weights."""

from __future__ import annotations

import numpy as np

from ..hosts import Host

_LOG_FLOOR = -10.0


def relative_adaptiveness_term(host: Host, codons: tuple[str, ...]) -> np.ndarray:
    """Log relative adaptiveness per codon.

    The sum over a sequence divided by its length is log CAI, so maximising the
    sum at fixed length maximises CAI. Unobserved codons are floored rather than
    excluded, keeping them reachable but heavily penalised.
    """
    adaptiveness = host.relative_adaptiveness
    return np.array(
        [
            np.log(adaptiveness[codon]) if adaptiveness.get(codon, 0.0) > 0 else _LOG_FLOOR
            for codon in codons
        ],
        dtype=np.float64,
    )


def codon_pair_matrix(
    host: Host,
    left: tuple[str, ...],
    right: tuple[str, ...],
) -> np.ndarray:
    """Codon-pair score for every admitted junction between two positions.

    Pairs absent from the reference set score zero. Absence in one genome is weak
    evidence, and a floor would act as a constraint the data cannot support.
    """
    scores = host.codon_pair_scores
    matrix = np.zeros((len(left), len(right)))
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            matrix[i, j] = scores.get((first, second), 0.0)
    return matrix


def gc_content_term(codons: tuple[str, ...]) -> np.ndarray:
    """Fraction of GC bases in each codon."""
    return np.array([sum(b in "GC" for b in codon) / 3.0 for codon in codons], dtype=np.float64)
