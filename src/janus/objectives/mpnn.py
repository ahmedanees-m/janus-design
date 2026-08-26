"""Inverse-folding marginals as lattice node weights.

ProteinMPNN run with ``--unconditional_probs_only`` returns one probability
vector per position from a single forward pass. These are additive over
positions, which is what makes the DP well posed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..genetic_code import AA_TO_INDEX, CODON_TO_AA, MPNN_ALPHABET

_LOG_FLOOR = -30.0


def load_unconditional(path: str | Path) -> np.ndarray:
    """Read a ProteinMPNN unconditional-probability archive as log marginals.

    ``log_p`` has shape ``(batch, length, 21)``. Batches hold repeated decodes of
    one backbone, so they are averaged in probability space; the unknown-token
    column is dropped and the rest renormalised.
    """
    archive = np.load(Path(path))
    key = "log_p" if "log_p" in archive else archive.files[0]
    log_p = np.asarray(archive[key], dtype=np.float64)

    if log_p.ndim == 3:
        log_p = np.log(np.exp(log_p).mean(axis=0))
    if log_p.ndim != 2:
        raise ValueError(f"expected a 2- or 3-dimensional array in {path}, got {log_p.shape}")
    if log_p.shape[1] != len(MPNN_ALPHABET):
        raise ValueError(
            f"expected {len(MPNN_ALPHABET)} alphabet columns in {path}, got {log_p.shape[1]}"
        )

    probabilities = np.exp(log_p[:, : len(AA_TO_INDEX)])
    totals = probabilities.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError(f"{path} contains a position with no probability mass on the 20 residues")
    return np.log(probabilities / totals)


def marginal_term(position_marginals: np.ndarray, codons: tuple[str, ...]) -> np.ndarray:
    """Node weight for each admitted codon: the marginal of the residue it encodes."""
    return np.array(
        [_lookup(position_marginals, CODON_TO_AA[codon]) for codon in codons],
        dtype=np.float64,
    )


def _lookup(position_marginals: np.ndarray, residue: str) -> float:
    index = AA_TO_INDEX.get(residue)
    if index is None:
        return _LOG_FLOOR
    return float(position_marginals[index])
