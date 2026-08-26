"""Standard genetic code and the alphabets used throughout the package.

Codons are enumerated in NCBI order (bases TCAG, first position slowest), which
is the order the standard-code amino-acid string below is written against.
"""

from __future__ import annotations

BASES = "TCAG"

# NCBI translation table 1, in NCBI codon order.
_TABLE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODONS: tuple[str, ...] = tuple(
    b1 + b2 + b3 for b1 in BASES for b2 in BASES for b3 in BASES
)

CODON_TO_AA: dict[str, str] = dict(zip(CODONS, _TABLE1, strict=True))

STOP_CODONS: tuple[str, ...] = tuple(c for c in CODONS if CODON_TO_AA[c] == "*")

SENSE_CODONS: tuple[str, ...] = tuple(c for c in CODONS if CODON_TO_AA[c] != "*")

# The 20 proteinogenic amino acids, in the order ProteinMPNN emits them. Its
# probability tensors carry a 21st column for the unknown token 'X', which the
# loader in objectives.mpnn drops.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
MPNN_ALPHABET = AA_ALPHABET + "X"

AA_TO_INDEX: dict[str, int] = {a: i for i, a in enumerate(AA_ALPHABET)}

SYNONYMOUS: dict[str, tuple[str, ...]] = {
    aa: tuple(c for c in SENSE_CODONS if CODON_TO_AA[c] == aa) for aa in AA_ALPHABET
}


def translate(cds: str) -> str:
    """Translate a coding sequence, stopping at the first in-frame stop codon."""
    if len(cds) % 3:
        raise ValueError(f"coding sequence length {len(cds)} is not a multiple of three")
    protein = []
    for i in range(0, len(cds), 3):
        aa = CODON_TO_AA[cds[i : i + 3].upper()]
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(b in "GC" for b in seq.upper()) / len(seq)
