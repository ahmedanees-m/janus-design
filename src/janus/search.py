"""Searching the shell under terms the parser cannot carry.

The parser solves the codon layer exactly for any fixed protein, and it does so
in milliseconds. What it cannot do is optimise a term that fails to decompose
over lattice nodes: mRNA folding, because base pairs form at arbitrary range, and
the protein-level liabilities, because they are functions of the whole residue
sequence. Those terms have to be searched.

Two searches are here and they differ only in acceptance order. ``optimise``
visits positions in a fixed order and takes the best admitted substitution at
each before moving on. ``best_improvement`` finds the best move over the whole
chain before taking any. Both re-solve the codon layer exactly after every move,
so the codon layer is never the approximate part.

``optimise`` is the default, chosen on measurement rather than taste. Over 447
backbones the two reach objectives within 1.2 percent of each other, and
``optimise`` gets there on 437 evaluations against 3121. The landscape is benign
enough that the cheaper ordering suffices. What the parser contributes is not the
ordering but the exact codon layer inside either one, and the exact optimum the
result is measured against.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from .design import Design, design
from .hosts import Host
from .objectives import Weights
from .sample import ShellSearch


def _entropy_order(log_marginals) -> list[int]:
    """Positions by descending marginal entropy, so the freest move first."""
    entropy = [-(np.exp(row) * row).sum() for row in log_marginals]
    return [int(i) for i in np.argsort(entropy)[::-1]]


def optimise(
    log_marginals,
    host: Host,
    score: Callable[[Design], float],
    weights: Weights | None = None,
    delta: float = 1.0,
    start: Design | None = None,
    anchor: str | None = None,
    order: Sequence[int] | None = None,
    rounds: int = 6,
) -> tuple[Design, float, int]:
    """Coordinate descent over the shell, codons re-solved exactly after each move.

    ``score`` takes a ``Design`` and returns a value to maximise; it is where the
    non-decomposable terms live. ``start`` defaults to the shell's Tier-1 optimum,
    and ``anchor`` admits a named residue sequence alongside the shell so a search
    can begin from a sequence the shell would not otherwise contain.

    Returns the design, its score, and the number of times ``score`` was called,
    since that count is the honest budget when the terms inside it are expensive.
    """
    search = ShellSearch(log_marginals, host, weights, delta, anchor=anchor)
    current = start if start is not None else design(
        log_marginals, host, weights=weights, delta=delta, k=1, anchor=anchor)[0]
    value = score(current)
    evaluations = 1
    positions = list(order) if order is not None else _entropy_order(log_marginals)

    for _ in range(rounds):
        improved = False
        for position in positions:
            best, best_value = None, value
            for residue in search.admitted[position]:
                if residue == current.protein[position]:
                    continue
                candidate = search.best(
                    current.protein[:position] + residue + current.protein[position + 1:]
                )
                evaluations += 1
                candidate_value = score(candidate)
                if candidate_value > best_value + 1e-12:
                    best, best_value = candidate, candidate_value
            if best is not None:
                current, value, improved = best, best_value, True
        if not improved:
            break
    return current, value, evaluations


def best_improvement(
    log_marginals,
    host: Host,
    score: Callable[[Design], float],
    weights: Weights | None = None,
    delta: float = 1.0,
    start: Design | None = None,
    anchor: str | None = None,
    max_steps: int = 40,
) -> tuple[Design, float, int]:
    """Take the best move over the whole chain at each step, rather than per position.

    Thorough and about seven times more expensive for the same answer on the
    backbones measured here. Kept because the comparison between the two is what
    established that the search strategy does not matter much, and a result that
    rests on a comparison should leave both arms runnable.
    """
    search = ShellSearch(log_marginals, host, weights, delta, anchor=anchor)
    current = start if start is not None else design(
        log_marginals, host, weights=weights, delta=delta, k=1, anchor=anchor)[0]
    value = score(current)
    evaluations, steps = 1, 0

    while steps < max_steps:
        best, best_value = None, value
        for position, residues in enumerate(search.admitted):
            for residue in residues:
                if residue == current.protein[position]:
                    continue
                candidate = search.best(
                    current.protein[:position] + residue + current.protein[position + 1:]
                )
                evaluations += 1
                candidate_value = score(candidate)
                if candidate_value > best_value + 1e-12:
                    best, best_value = candidate, candidate_value
        if best is None:
            break
        current, value, steps = best, best_value, steps + 1
    return current, value, evaluations
