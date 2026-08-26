"""Weighted lattice parsing and k-best extraction."""

from __future__ import annotations

import heapq
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Path:
    score: float
    states: tuple[int, ...]


def kbest(
    node_scores: list[np.ndarray],
    edge_scores: list[np.ndarray],
    k: int = 1,
) -> list[Path]:
    """Return up to ``k`` highest-scoring paths, best first.

    Exact for the scores it is given: the list is the true top ``k`` under the
    sum of node and edge scores, with no beam and no pruning. That says nothing
    about whether those scores are the right objective, which is the caller's
    question and is qualified where the objective is assembled.

    ``node_scores[i]`` holds one score per state at position ``i``.
    ``edge_scores[i]`` is the transition matrix from position ``i`` to ``i + 1``,
    shaped ``(len(node_scores[i]), len(node_scores[i + 1]))``.
    """
    if k < 1:
        raise ValueError("k must be at least one")
    length = len(node_scores)
    if length == 0:
        return []
    if len(edge_scores) != length - 1:
        raise ValueError(
            f"expected {length - 1} transition matrices for {length} positions, "
            f"got {len(edge_scores)}"
        )
    for i, matrix in enumerate(edge_scores):
        expected = (len(node_scores[i]), len(node_scores[i + 1]))
        if matrix.shape != expected:
            raise ValueError(f"transition matrix {i} has shape {matrix.shape}, expected {expected}")

    # prefixes[i][s] is a descending list of (score, predecessor state, rank).
    prefixes: list[list[list[tuple[float, int, int]]]] = [
        [[] for _ in range(len(scores))] for scores in node_scores
    ]
    for state, score in enumerate(node_scores[0]):
        prefixes[0][state] = [(float(score), -1, -1)]

    for i in range(1, length):
        transitions = edge_scores[i - 1]
        previous = prefixes[i - 1]
        for state, node in enumerate(node_scores[i]):
            column = transitions[:, state]
            heap = [
                (-(previous[p][0][0] + column[p]), p, 0)
                for p in range(len(previous))
                if previous[p]
            ]
            heapq.heapify(heap)

            merged: list[tuple[float, int, int]] = []
            while heap and len(merged) < k:
                negated, p, rank = heapq.heappop(heap)
                merged.append((-negated + float(node), p, rank))
                if rank + 1 < len(previous[p]):
                    heapq.heappush(
                        heap,
                        (-(previous[p][rank + 1][0] + column[p]), p, rank + 1),
                    )
            prefixes[i][state] = merged

    final = prefixes[-1]
    heap = [(-final[s][0][0], s, 0) for s in range(len(final)) if final[s]]
    heapq.heapify(heap)

    paths: list[Path] = []
    while heap and len(paths) < k:
        negated, state, rank = heapq.heappop(heap)
        paths.append(Path(score=-negated, states=_backtrack(prefixes, state, rank)))
        if rank + 1 < len(final[state]):
            heapq.heappush(heap, (-final[state][rank + 1][0], state, rank + 1))

    return paths


def _backtrack(
    prefixes: list[list[list[tuple[float, int, int]]]],
    state: int,
    rank: int,
) -> tuple[int, ...]:
    states = [state]
    for i in range(len(prefixes) - 1, 0, -1):
        _, state, rank = prefixes[i][state][rank]
        states.append(state)
    states.reverse()
    return tuple(states)
