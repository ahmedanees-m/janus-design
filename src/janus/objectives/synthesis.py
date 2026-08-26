"""Synthesisability constraints on a finished coding sequence.

Vendors reject outright on repeats, skewed GC, long homopolymers and internal
restriction sites, so these are constraints rather than weighted penalties.

None decomposes over adjacent codon pairs: a six-base site spans up to three
codons, a homopolymer limit of nine spans four, and repeat and windowed-GC
checks are global. They are checked against the finished sequence, with the
k-best list supplying fallbacks. Folding them into the search needs the
nucleotide-level lattice.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..genetic_code import gc_fraction
from ..hosts import GCWindow, Host


@dataclass(frozen=True)
class Violation:
    kind: str
    detail: str
    start: int
    end: int


def violations(cds: str, host: Host) -> list[Violation]:
    """Every hard constraint the sequence breaks, in order of first position."""
    sequence = cds.upper()
    limits = host.constraints
    found: list[Violation] = []

    overall = gc_fraction(sequence)
    low, high = limits.gc_range
    if not low <= overall <= high:
        found.append(
            Violation("gc_content", f"{overall:.3f} outside [{low}, {high}]", 0, len(sequence))
        )

    for window in limits.gc_windows:
        found.extend(_gc_window(sequence, window))
    found.extend(_homopolymers(sequence, limits.max_homopolymer))
    found.extend(_forbidden_sites(sequence, limits.forbidden_sites))
    found.extend(_repeats(sequence, limits.max_repeat_length))

    return sorted(found, key=lambda v: (v.start, v.kind))


def _gc_window(sequence: str, window: GCWindow) -> list[Violation]:
    width = window.width
    if width <= 0 or len(sequence) < width:
        return []
    found = []
    for start in range(len(sequence) - width + 1):
        value = gc_fraction(sequence[start : start + width])
        if not window.minimum <= value <= window.maximum:
            found.append(
                Violation(
                    "gc_window",
                    f"{width} bp window at {start} is {value:.3f}, outside "
                    f"[{window.minimum}, {window.maximum}]",
                    start,
                    start + width,
                )
            )
    return _merge_overlapping(found)


def _homopolymers(sequence: str, limits: dict[str, int]) -> list[Violation]:
    found = []
    start = 0
    while start < len(sequence):
        end = start
        while end < len(sequence) and sequence[end] == sequence[start]:
            end += 1
        run = end - start
        allowed = limits.get(sequence[start])
        if allowed is not None and run > allowed:
            found.append(
                Violation("homopolymer", f"{run} consecutive {sequence[start]} (max {allowed})", start, end)
            )
        start = end
    return found


def _forbidden_sites(sequence: str, sites: dict[str, str]) -> list[Violation]:
    found = []
    for name, motif in sites.items():
        start = sequence.find(motif)
        while start != -1:
            found.append(Violation("forbidden_site", f"{name} ({motif})", start, start + len(motif)))
            start = sequence.find(motif, start + 1)
    return found


def _repeats(sequence: str, minimum: int) -> list[Violation]:
    if minimum <= 0 or len(sequence) < 2 * minimum:
        return []
    seen: dict[str, int] = {}
    found = []
    for start in range(len(sequence) - minimum + 1):
        kmer = sequence[start : start + minimum]
        first = seen.setdefault(kmer, start)
        if first != start:
            found.append(
                Violation("repeat", f"{minimum} bp repeat of {kmer} first seen at {first}", start, start + minimum)
            )
    return _merge_overlapping(found)


def _merge_overlapping(found: list[Violation]) -> list[Violation]:
    """Collapse runs of overlapping violations of the same kind into one report."""
    if not found:
        return []
    merged = [found[0]]
    for violation in found[1:]:
        last = merged[-1]
        if violation.kind == last.kind and violation.start <= last.end:
            merged[-1] = Violation(last.kind, last.detail, last.start, max(last.end, violation.end))
        else:
            merged.append(violation)
    return merged
