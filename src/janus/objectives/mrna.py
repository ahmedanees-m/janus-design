"""mRNA folding terms.

These do not decompose over lattice nodes: base pairs form at arbitrary range,
so a path's folding energy is not a sum over its positions. They are evaluated
on finished transcripts in the rescoring tier.

The 5' UTR is held fixed at a host-appropriate sequence and reported as a
controlled variable. Joint UTR and CDS design is a different axis, already
occupied by LinearDesign2, and is out of scope here.
"""

from __future__ import annotations

from functools import lru_cache

from ..hosts import Host

_RNA = None


def _vienna():
    global _RNA
    if _RNA is None:
        import RNA

        _RNA = RNA
    return _RNA


def available() -> bool:
    try:
        _vienna()
    except ImportError:
        return False
    return True


def transcript(cds: str, host: Host) -> str:
    """Fixed 5' UTR, a start codon if the design lacks one, then the design."""
    body = cds.upper()
    if not body.startswith("ATG"):
        body = "ATG" + body
    return (host.initiation.utr + body).replace("T", "U")


def initiation_energy(cds: str, host: Host) -> float:
    """Folding free energy of the initiation window, kcal/mol.

    Kudla 2009 and Goodman 2013 both find the folding stability of this region
    to be the dominant coding-sequence determinant of expression in *E. coli*,
    ahead of codon usage. More negative means more structure over the region the
    30S subunit must load onto, so the objective wants this near zero.
    """
    window = _initiation_window(cds, host)
    return _fold(window)[1] if window else 0.0


def global_mfe(cds: str, host: Host) -> float:
    """Whole-transcript minimum free energy, kcal/mol."""
    return _fold(transcript(cds, host))[1]


def initiation_structure(cds: str, host: Host) -> str:
    return _fold(_initiation_window(cds, host))[0]


def _initiation_window(cds: str, host: Host) -> str:
    sequence = transcript(cds, host)
    origin = len(host.initiation.utr)
    low, high = host.initiation.window
    return sequence[max(0, origin + low) : min(len(sequence), origin + high)]


@lru_cache(maxsize=200_000)
def _fold(sequence: str) -> tuple[str, float]:
    if not sequence:
        return "", 0.0
    structure, energy = _vienna().fold(sequence)
    return structure, float(energy)
