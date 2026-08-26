"""Protein-level liabilities: the terms a codon optimiser cannot reach.

Every score here is a function of the residue sequence and, where structure
matters, of a fixed per-residue accessibility. None of them depends on the
codons. That is the point: for these terms a fixed-protein codon designer has an
achievable improvement of exactly zero, and the joint lattice is necessary rather
than marginal.

Higher is worse throughout, so the rescoring tier subtracts them.
"""

from __future__ import annotations

import re
from collections import Counter

import numpy as np

from .proteostasis import load_by_family, scan

# Kyte and Doolittle 1982.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
EXPOSED = 0.25
HYDROPHOBIC = 1.8
LC_WINDOW = 12
LC_ENTROPY = 1.5


def low_complexity(protein: str) -> float:
    """Fraction of windows whose residue-identity entropy falls below threshold."""
    if len(protein) < LC_WINDOW:
        return 0.0
    flagged = 0
    for start in range(len(protein) - LC_WINDOW + 1):
        counts = np.array(list(Counter(protein[start : start + LC_WINDOW]).values()), dtype=float)
        p = counts / counts.sum()
        if -(p * np.log(p)).sum() < LC_ENTROPY:
            flagged += 1
    return flagged / (len(protein) - LC_WINDOW + 1)


def longest_repeat(protein: str, minimum: int = 3) -> int:
    """Length of the longest residue substring occurring more than once."""
    for length in range(len(protein) // 2, minimum - 1, -1):
        seen = set()
        for start in range(len(protein) - length + 1):
            window = protein[start : start + length]
            if window in seen:
                return length
            seen.add(window)
    return 0


def exposed_hydrophobic(protein: str, accessibility) -> float:
    """Summed relative accessibility over exposed hydrophobic residues.

    Node-decomposable given a fixed backbone, unlike the other terms here, but
    still untouchable by codon choice.
    """
    if accessibility is None:
        return float(sum(1 for r in protein if HYDROPATHY.get(r, 0.0) > HYDROPHOBIC))
    total = 0.0
    for residue, area in zip(protein, accessibility, strict=True):
        if area is not None and area > EXPOSED and HYDROPATHY.get(residue, 0.0) > HYDROPHOBIC:
            total += area
    return total


def longest_exposed_hydrophobic_run(protein: str, accessibility) -> int:
    if accessibility is None:
        return 0
    best = current = 0
    for residue, area in zip(protein, accessibility, strict=True):
        if area is not None and area > EXPOSED and HYDROPATHY.get(residue, 0.0) > HYDROPHOBIC:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def degron_load(protein: str, accessibility, classes) -> float:
    """Accessibility-weighted degron motif load.

    Returns zero when no motif classes are supplied, so a caller that has not
    loaded ELM gets a term that is absent rather than one that is wrong.
    """
    if not classes:
        return 0.0
    hits = scan(protein, classes, accessibility)
    return float(load_by_family(hits, weighted=True)["all_degron"])


TERMS = {
    "low_complexity": lambda p, a, c: low_complexity(p),
    "protein_repeat": lambda p, a, c: float(longest_repeat(p)),
    "exposed_hydrophobic": lambda p, a, c: exposed_hydrophobic(p, a),
    "exposed_hydrophobic_run": lambda p, a, c: float(longest_exposed_hydrophobic_run(p, a)),
    "degron": lambda p, a, c: degron_load(p, a, c),
}


def score(name: str, protein: str, accessibility=None, classes=None) -> float:
    """One named liability for one protein. Higher is worse."""
    if name not in TERMS:
        raise KeyError(f"unknown liability {name!r}; known: {sorted(TERMS)}")
    return float(TERMS[name](protein, accessibility, classes))


def profile(protein: str, accessibility=None, classes=None) -> dict[str, float]:
    """Every liability at once, for reporting."""
    return {name: score(name, protein, accessibility, classes) for name in TERMS}
