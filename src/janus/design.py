"""Top-level entry point: backbone marginals and a host in, coding sequences out."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import lattice as lattice_module
from . import parse
from .genetic_code import translate
from .hosts import Host
from .objectives import Weights, assemble, score_terms
from .objectives.codon import codon_pair_matrix, gc_content_term, relative_adaptiveness_term
from .objectives.mpnn import marginal_term
from .objectives.synthesis import Violation, violations


@dataclass(frozen=True)
class Design:
    protein: str
    cds: str
    score: float
    terms: dict[str, float]
    violations: list[Violation]

    @property
    def synthesisable(self) -> bool:
        return not self.violations


def evaluate(cds: str, host: Host, log_marginals: np.ndarray,
             weights: Weights | None = None) -> Design:
    """Score a coding sequence the lattice did not produce.

    Baselines return finished genes, and comparing them needs the same objective
    the parser maximises rather than the unweighted report ``score_terms`` gives.
    A trailing stop codon is dropped, since the marginals cover the chain only.
    """
    weights = weights or Weights()
    body = cds[:-3] if len(cds) == 3 * (len(log_marginals) + 1) else cds
    codons = [body[i : i + 3] for i in range(0, len(body), 3)]
    if len(codons) != len(log_marginals):
        raise ValueError(
            f"expected {len(log_marginals)} codons, got {len(codons)}"
        )

    score = 0.0
    for position, codon in enumerate(codons):
        if weights.mpnn:
            score += weights.mpnn * float(marginal_term(log_marginals[position], (codon,))[0])
        if weights.cai:
            score += weights.cai * float(relative_adaptiveness_term(host, (codon,))[0])
        if weights.gc:
            score += weights.gc * float(gc_content_term((codon,))[0])
    if weights.cpb:
        for left, right in zip(codons, codons[1:], strict=False):
            score += weights.cpb * float(codon_pair_matrix(host, (left,), (right,))[0][0])

    return Design(
        protein=translate(body),
        cds=body,
        score=score,
        terms=score_terms(body, host, log_marginals),
        violations=violations(body, host),
    )


def design(
    log_marginals: np.ndarray,
    host: Host,
    weights: Weights | None = None,
    delta: float = 1.0,
    k: int = 1,
    fixed: dict[int, str] | None = None,
    anchor: str | None = None,
) -> list[Design]:
    """Return the ``k`` highest-scoring (protein, gene) pairs, best first.

    The score is exact with respect to the decomposable objective built from
    ``weights``, defined on unconditional marginals rather than the conditional
    posterior. It is not a fold likelihood under the full autoregressive model
    and carries no folding, cis-element or proteostasis term. Synthesis
    constraints are reported per design, not optimised.
    """
    weights = weights or Weights()
    built = lattice_module.build(log_marginals, delta=delta, fixed=fixed, anchor=anchor)
    node_scores, edge_scores = assemble(built, host, log_marginals, weights)

    designs = []
    for path in parse.kbest(node_scores, edge_scores, k=k):
        cds = "".join(built.codons[i][state] for i, state in enumerate(path.states))
        designs.append(
            Design(
                protein=translate(cds),
                cds=cds,
                score=path.score,
                terms=score_terms(cds, host, log_marginals),
                violations=violations(cds, host),
            )
        )
    return designs
