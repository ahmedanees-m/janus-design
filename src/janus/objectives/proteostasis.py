"""Degron and linear-motif load, weighted by structural accessibility.

Motif definitions come from the Eukaryotic Linear Motif resource rather than
being transcribed from the primary papers, so they are versioned and match what
other groups scan with.

Accessibility weighting matters here. The peptide screens that defined the
C-degron classes fuse a short peptide to a reporter, leaving a free solvent
exposed C terminus, which is the precondition for recognition. In a folded
miniprotein the same residues may be packed against the core. Counting raw motif
hits transfers a peptide result to a folded domain without that correction and
overstates the risk; both counts are therefore reported.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

FAMILIES = {
    "DEG_Cend": "c_degron",
    "DEG_Nend": "n_degron",
    "DEG": "internal_degron",
    "CLV": "protease_site",
    "TRG": "targeting",
    "MOD": "modification",
}


@dataclass(frozen=True)
class MotifClass:
    identifier: str
    family: str
    description: str
    pattern: re.Pattern


@dataclass(frozen=True)
class MotifHit:
    identifier: str
    family: str
    start: int
    end: int
    accessibility: float


def family_of(identifier: str) -> str | None:
    for prefix in ("DEG_Cend", "DEG_Nend", "DEG", "CLV", "TRG", "MOD"):
        if identifier.startswith(prefix):
            return FAMILIES[prefix]
    return None


def load_classes(path: str | Path, families: set[str] | None = None) -> list[MotifClass]:
    """Read the ELM classes table, keeping the families we score."""
    classes = []
    with Path(path).open(encoding="utf-8", newline="") as fh:
        rows = [line for line in fh if not line.startswith("#")]
    for row in csv.DictReader(rows, delimiter="\t"):
        identifier = row["ELMIdentifier"]
        family = family_of(identifier)
        if family is None or (families and family not in families):
            continue
        try:
            pattern = re.compile(f"(?=({row['Regex']}))")
        except re.error:
            continue
        classes.append(MotifClass(identifier, family, row["Description"], pattern))
    return classes


def scan(protein: str, classes: list[MotifClass], accessibility=None) -> list[MotifHit]:
    """Every match of every class, including overlaps.

    ``accessibility`` is one relative solvent accessibility per residue. When it
    is absent every hit is weighted one, which reproduces the unweighted count.
    """
    hits = []
    for motif in classes:
        for match in motif.pattern.finditer(protein):
            text = match.group(1)
            start = match.start()
            end = start + len(text)
            if accessibility is None:
                weight = 1.0
            else:
                window = [a for a in accessibility[start:end] if a is not None]
                weight = float(sum(window) / len(window)) if window else 0.0
            hits.append(MotifHit(motif.identifier, motif.family, start, end, weight))
    return hits


def load_by_family(hits: list[MotifHit], weighted: bool) -> dict[str, float]:
    """Total motif load per family, weighted by accessibility or not."""
    totals = {family: 0.0 for family in set(FAMILIES.values())}
    for hit in hits:
        totals[hit.family] += hit.accessibility if weighted else 1.0
    totals["all_degron"] = sum(
        v for k, v in totals.items() if k.endswith("degron")
    )
    return totals
