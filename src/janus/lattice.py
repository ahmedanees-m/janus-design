"""Amino-acid-degenerate codon lattice.

LinearDesign encodes the synonymous coding sequences of a fixed protein as a
codon DFA. Here each position additionally admits the codons of every amino acid
within a log-probability shell of the inverse-folding marginal, so one path
picks a protein and a gene together.

Shells are defined on ProteinMPNN unconditional marginals. The conditional
posterior depends on residues chosen elsewhere and does not decompose over
nodes, so it cannot be optimised here; it belongs in the rescoring tier.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .genetic_code import AA_ALPHABET, SYNONYMOUS


@dataclass(frozen=True)
class Lattice:
    """Admitted (amino acid, codon) branches at each position of the chain."""

    codons: tuple[tuple[str, ...], ...]
    amino_acids: tuple[tuple[str, ...], ...]

    def __len__(self) -> int:
        return len(self.codons)

    @property
    def branching(self) -> float:
        """Mean number of amino acids admitted per position."""
        return float(np.mean([len(a) for a in self.amino_acids]))

    @property
    def log_size(self) -> float:
        """Natural log of the number of paths, for reporting search-space size."""
        return float(sum(np.log(len(c)) for c in self.codons))


def amino_acid_shell(log_marginals: np.ndarray, delta: float) -> list[list[str]]:
    """Amino acids within ``delta`` nats of the per-position marginal argmax.

    Zero collapses each position to its argmax and recovers the fixed-protein
    codon lattice. Above roughly three nats most of the alphabet is admitted at
    unconstrained positions.
    """
    if delta < 0:
        raise ValueError("delta must be non-negative")
    if log_marginals.ndim != 2 or log_marginals.shape[1] != len(AA_ALPHABET):
        raise ValueError(
            f"expected marginals of shape (L, {len(AA_ALPHABET)}), got {log_marginals.shape}"
        )

    cutoffs = log_marginals.max(axis=1) - delta
    shells = []
    for row, cutoff in zip(log_marginals, cutoffs, strict=True):
        admitted = [AA_ALPHABET[j] for j in np.flatnonzero(row >= cutoff)]
        shells.append(admitted)
    return shells


def build(
    log_marginals: np.ndarray,
    delta: float,
    fixed: dict[int, str] | None = None,
    anchor: str | None = None,
) -> Lattice:
    """Build the lattice for one backbone.

    ``fixed`` pins positions to a given amino acid, for catalytic residues,
    interfaces and tags.

    ``anchor`` additionally admits one residue sequence wherever the shell does
    not already contain it. The shell is centred on the marginal argmax, so a
    given protein, such as a design's own sequence, need not lie inside it; a
    comparison that has to start from that protein and then use the freedom
    needs a lattice that holds both. ``fixed`` still wins where the two meet.
    """
    shells = amino_acid_shell(log_marginals, delta)

    for position, residue in enumerate(anchor or ""):
        if residue not in SYNONYMOUS:
            raise ValueError(f"anchor position {position} has unknown residue {residue!r}")
        if residue not in shells[position]:
            shells[position] = shells[position] + [residue]

    for position, residue in (fixed or {}).items():
        if residue not in SYNONYMOUS:
            raise ValueError(f"position {position} pinned to unknown residue {residue!r}")
        shells[position] = [residue]

    codons = tuple(
        tuple(codon for aa in shell for codon in SYNONYMOUS[aa]) for shell in shells
    )
    amino_acids = tuple(tuple(shell) for shell in shells)
    return Lattice(codons=codons, amino_acids=amino_acids)
