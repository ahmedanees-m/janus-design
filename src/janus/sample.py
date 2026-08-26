"""Candidate generation over the shell with exact codons.

The k-best prefix enumerates small perturbations around the Tier-1 optimum, so
it searches a narrow neighbourhood deeply. Once a rescoring term is strong enough
to move the optimum out of that neighbourhood the prefix stops being the better
candidate pool and breadth pays instead. Drawing residue sequences from the shell
and optimising their codons exactly gives that breadth while keeping the codon
layer optimal.

Building the lattice and the term tables costs far more than optimising the
codons for one residue assignment, and it does not depend on the residues. A
caller that visits many assignments, such as a local search, should build once
through ``ShellSearch`` and re-optimise per assignment.
"""

from __future__ import annotations

import numpy as np

from .design import Design
from .genetic_code import AA_ALPHABET, CODON_TO_AA, translate
from .hosts import Host
from .lattice import build
from .objectives import Weights, assemble, score_terms
from .objectives.synthesis import violations
from .parse import kbest

AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}


def _residue_slots(lattice):
    table = []
    for codons in lattice.codons:
        by_residue: dict[str, list[int]] = {}
        for index, codon in enumerate(codons):
            by_residue.setdefault(CODON_TO_AA[codon], []).append(index)
        table.append(by_residue)
    return table


class ShellSearch:
    """A built shell lattice that optimises codons for any admitted residue string."""

    def __init__(self, log_marginals, host: Host, weights: Weights | None = None,
                 delta: float = 1.0, anchor: str | None = None):
        self.log_marginals = log_marginals
        self.host = host
        self.weights = weights or Weights()
        self.lattice = build(log_marginals, delta=delta, anchor=anchor)
        self.node, self.edge = assemble(self.lattice, host, log_marginals, self.weights)
        self.slots = _residue_slots(self.lattice)

    @property
    def admitted(self):
        """Residues the shell allows at each position, in lattice order."""
        return self.lattice.amino_acids

    def best(self, residues) -> Design:
        """Highest-scoring gene encoding ``residues``, which must lie in the shell.

        Exact with respect to the same decomposable objective the parser uses,
        with the residues held fixed. It is the codon layer's optimum for that
        protein, not the lattice's optimum.
        """
        keep = [self.slots[i][r] for i, r in enumerate(residues)]
        sub_node = [self.node[i][keep[i]] for i in range(len(self.node))]
        sub_edge = [self.edge[i][np.ix_(keep[i], keep[i + 1])] for i in range(len(self.edge))]
        path = kbest(sub_node, sub_edge, k=1)[0]
        cds = "".join(self.lattice.codons[i][keep[i][s]] for i, s in enumerate(path.states))
        return Design(
            protein=translate(cds),
            cds=cds,
            score=path.score,
            terms=score_terms(cds, self.host, self.log_marginals),
            violations=violations(cds, self.host),
        )


def shell_samples(
    log_marginals,
    host: Host,
    weights: Weights | None = None,
    delta: float = 1.0,
    count: int = 200,
    temperature: float | None = None,
    rng=None,
) -> list[Design]:
    """Draw ``count`` residue sequences from the shell, codons optimised exactly.

    ``temperature`` samples residues from the tempered marginal; ``None`` draws
    uniformly over the shell, which explores it more evenly and is what the
    two-tier comparison used.
    """
    rng = rng or np.random.default_rng()
    search = ShellSearch(log_marginals, host, weights, delta)

    picker = []
    for position, admitted in enumerate(search.admitted):
        if temperature is None:
            picker.append(None)
            continue
        logits = np.array([log_marginals[position][AA_INDEX[a]] for a in admitted]) / temperature
        logits -= logits.max()
        probability = np.exp(logits)
        picker.append(probability / probability.sum())

    out = []
    for _ in range(count):
        residues = []
        for position, admitted in enumerate(search.admitted):
            if picker[position] is None:
                residues.append(admitted[int(rng.integers(len(admitted)))])
            else:
                residues.append(admitted[int(rng.choice(len(admitted), p=picker[position]))])
        out.append(search.best(residues))
    return out
