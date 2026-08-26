"""Open shell, marginal term only: the solver must return the marginal argmax.

Isolates the amino-acid layer. The LinearDesign test holds the protein fixed, so
a fault in amino-acid branching would pass it silently. Only jointly sufficient.
"""

from __future__ import annotations

import numpy as np
import pytest

from janus import Weights, design
from janus.genetic_code import AA_ALPHABET, CODON_TO_AA, translate
from janus.lattice import amino_acid_shell


def argmax_protein(marginals: np.ndarray) -> str:
    return "".join(AA_ALPHABET[i] for i in marginals.argmax(axis=1))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_open_shell_marginal_only_recovers_argmax(ecoli, marginals, seed):
    log_marginals = marginals(length=60, seed=seed)
    result = design(
        log_marginals,
        ecoli,
        weights=Weights(mpnn=1.0),
        delta=np.inf,
        k=1,
    )[0]

    assert result.protein == argmax_protein(log_marginals)
    assert translate(result.cds) == result.protein


def test_recovery_holds_on_a_near_flat_posterior(ecoli, marginals):
    """Flat posterior: small argmax margins expose off-by-one node weighting."""
    log_marginals = marginals(length=40, seed=7, concentration=0.05)
    result = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=np.inf, k=1)[0]
    assert result.protein == argmax_protein(log_marginals)


def test_score_equals_the_summed_marginals(ecoli, marginals):
    log_marginals = marginals(length=30, seed=11)
    result = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=np.inf, k=1)[0]

    expected = sum(
        log_marginals[i, AA_ALPHABET.index(CODON_TO_AA[result.cds[3 * i : 3 * i + 3]])]
        for i in range(len(result.protein))
    )
    assert result.score == pytest.approx(expected)


def test_open_shell_admits_every_residue(marginals):
    shells = amino_acid_shell(marginals(length=25, seed=3), delta=np.inf)
    assert all(len(shell) == len(AA_ALPHABET) for shell in shells)


def test_shell_always_contains_the_argmax(marginals):
    log_marginals = marginals(length=25, seed=5)
    for delta in (0.0, 0.5, 1.0, 2.0, 3.0):
        shells = amino_acid_shell(log_marginals, delta=delta)
        for position, shell in enumerate(shells):
            assert AA_ALPHABET[log_marginals[position].argmax()] in shell


def test_shells_grow_monotonically_with_delta(marginals):
    log_marginals = marginals(length=25, seed=5)
    previous = amino_acid_shell(log_marginals, delta=0.0)
    for delta in (0.5, 1.0, 2.0, 4.0):
        current = amino_acid_shell(log_marginals, delta=delta)
        for earlier, later in zip(previous, current, strict=True):
            assert set(earlier) <= set(later)
        previous = current
