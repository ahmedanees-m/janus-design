"""Invariants of k-best extraction.

A k-best list that is merely good, or that repeats paths, still yields plausible
sequences and plausible Pareto fronts, so these compare against exhaustive
enumeration on lattices small enough to enumerate.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from janus.parse import Path, kbest


def random_lattice(widths, seed=0):
    rng = np.random.default_rng(seed)
    node_scores = [rng.normal(size=width) for width in widths]
    edge_scores = [
        rng.normal(size=(widths[i], widths[i + 1])) for i in range(len(widths) - 1)
    ]
    return node_scores, edge_scores


def enumerate_paths(node_scores, edge_scores):
    scored = []
    for states in itertools.product(*(range(len(s)) for s in node_scores)):
        total = sum(node_scores[i][s] for i, s in enumerate(states))
        total += sum(edge_scores[i][states[i], states[i + 1]] for i in range(len(states) - 1))
        scored.append(Path(score=total, states=states))
    scored.sort(key=lambda p: -p.score)
    return scored


@pytest.mark.parametrize("widths", [(3, 4, 2, 5), (2, 2, 2, 2, 2, 2), (6, 1, 6), (4,)])
@pytest.mark.parametrize("k", [1, 3, 10])
def test_kbest_matches_exhaustive_enumeration(widths, k):
    node_scores, edge_scores = random_lattice(widths, seed=sum(widths) + k)
    expected = enumerate_paths(node_scores, edge_scores)[:k]
    found = kbest(node_scores, edge_scores, k=k)

    assert len(found) == len(expected)
    for got, want in zip(found, expected, strict=True):
        assert got.score == pytest.approx(want.score)
    assert [p.states for p in found] == [p.states for p in expected]


def test_kbest_paths_are_distinct():
    node_scores, edge_scores = random_lattice((5, 5, 5, 5), seed=42)
    found = kbest(node_scores, edge_scores, k=50)
    assert len({p.states for p in found}) == len(found)


def test_kbest_scores_are_non_increasing():
    node_scores, edge_scores = random_lattice((4, 6, 3, 7, 2), seed=17)
    scores = [p.score for p in kbest(node_scores, edge_scores, k=40)]
    assert all(a >= b for a, b in zip(scores, scores[1:], strict=False))


def test_kbest_returns_every_path_when_k_exceeds_the_lattice():
    widths = (2, 3, 2)
    node_scores, edge_scores = random_lattice(widths, seed=1)
    found = kbest(node_scores, edge_scores, k=1000)
    assert len(found) == 2 * 3 * 2


def test_reported_score_matches_the_returned_path():
    node_scores, edge_scores = random_lattice((4, 5, 4, 5), seed=23)
    for path in kbest(node_scores, edge_scores, k=15):
        total = sum(node_scores[i][s] for i, s in enumerate(path.states))
        total += sum(
            edge_scores[i][path.states[i], path.states[i + 1]] for i in range(len(path.states) - 1)
        )
        assert path.score == pytest.approx(total)


def test_ties_do_not_collapse_distinct_paths():
    """Ties are common once weights are zeroed; a merge keyed on score rather
    than provenance drops paths silently."""
    node_scores = [np.zeros(3), np.zeros(3), np.zeros(3)]
    edge_scores = [np.zeros((3, 3)), np.zeros((3, 3))]
    found = kbest(node_scores, edge_scores, k=27)
    assert len({p.states for p in found}) == 27


def test_rejects_mismatched_transition_matrices():
    node_scores = [np.zeros(3), np.zeros(4)]
    with pytest.raises(ValueError, match="shape"):
        kbest(node_scores, [np.zeros((3, 3))], k=1)

    with pytest.raises(ValueError, match="transition matrices"):
        kbest(node_scores, [], k=1)
