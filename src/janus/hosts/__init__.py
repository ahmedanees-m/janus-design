"""Host definitions.

A host is a data file rather than a code path: ``hosts/<name>.yaml`` plus
sidecar count tables. Adding an organism must not require touching the lattice
or the parser.

The YAML records observed counts from a named reference CDS set. Relative
adaptiveness and codon-pair scores are derived at load time.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from ..genetic_code import CODON_TO_AA, SENSE_CODONS, SYNONYMOUS

HOST_DIR = Path(__file__).parent

__all__ = [
    "Constraints",
    "GCWindow",
    "Host",
    "Initiation",
    "available",
    "codon_pair_scores",
    "count_reference_cds",
    "load",
    "read_codon_counts",
    "read_pair_counts",
    "relative_adaptiveness",
    "write_counts",
]


@dataclass(frozen=True)
class GCWindow:
    width: int
    minimum: float
    maximum: float


@dataclass(frozen=True)
class Constraints:
    """Hard sequence requirements, checked against a finished coding sequence."""

    gc_range: tuple[float, float]
    gc_windows: tuple[GCWindow, ...]
    max_homopolymer: dict[str, int]
    max_repeat_length: int = 10
    forbidden_sites: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class Initiation:
    """Host-class-specific 5' model.

    ``shine_dalgarno`` penalises local structure over ``window``, in nucleotides
    relative to the A of the start codon. ``cap_scanning`` uses Kozak compliance
    plus start-proximal openness and flips the sign of the global folding term.
    Separate models rather than one sign flip: the mechanisms differ and only the
    prokaryotic one has large library support.
    """

    model: str
    window: tuple[int, int]
    global_mfe_sign: int
    utr: str = ""
    kozak: str | None = None


@dataclass(frozen=True)
class Host:
    name: str
    taxid: int
    source: dict
    constraints: Constraints
    initiation: Initiation
    codon_counts: dict[str, int]
    codon_pair_counts: dict[tuple[str, str], int]
    relative_adaptiveness: dict[str, float]
    codon_pair_scores: dict[tuple[str, str], float]


def available() -> list[str]:
    return sorted(path.stem for path in HOST_DIR.glob("*.yaml"))


def load(name_or_path: str | Path) -> Host:
    """Load a host by shipped name (``ecoli_bl21``) or by path to a YAML file."""
    path = Path(name_or_path)
    if not path.suffix:
        path = HOST_DIR / f"{path}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"no host file at {path}; shipped hosts: {', '.join(available())}")

    with path.open(encoding="utf-8") as fh:
        spec = yaml.safe_load(fh)

    codon_counts = read_codon_counts(path.parent / spec["codon_counts"])
    unknown = set(codon_counts) - set(SENSE_CODONS)
    if unknown:
        raise ValueError(f"{spec['codon_counts']} lists non-sense codons: {sorted(unknown)}")

    pair_counts = read_pair_counts(path.parent / spec["codon_pair_counts"])

    return Host(
        name=spec["name"],
        taxid=int(spec["taxid"]),
        source=spec.get("source", {}),
        constraints=_constraints(spec["constraints"]),
        initiation=_initiation(spec["initiation"]),
        codon_counts=codon_counts,
        codon_pair_counts=pair_counts,
        relative_adaptiveness=relative_adaptiveness(codon_counts),
        codon_pair_scores=codon_pair_scores(codon_counts, pair_counts),
    )


def relative_adaptiveness(codon_counts: dict[str, int]) -> dict[str, float]:
    """Sharp and Li w values: each codon's count over the best count in its family."""
    values = {}
    for family in SYNONYMOUS.values():
        best = max(codon_counts.get(codon, 0) for codon in family)
        for codon in family:
            values[codon] = (codon_counts.get(codon, 0) / best) if best else 0.0
    return values


def codon_pair_scores(
    codon_counts: dict[str, int],
    pair_counts: dict[tuple[str, str], int],
) -> dict[tuple[str, str], float]:
    """Coleman codon-pair score: log of observed over expected pair frequency.

    The expectation holds amino-acid pair usage and codon bias fixed. Unobserved
    pairs are omitted rather than floored.
    """
    total_pairs = sum(pair_counts.values())
    total_codons = sum(codon_counts.values())
    if not total_pairs or not total_codons:
        return {}

    aa_counts: Counter[str] = Counter()
    for codon, n in codon_counts.items():
        aa_counts[CODON_TO_AA[codon]] += n

    aa_pair_counts: Counter[tuple[str, str]] = Counter()
    for (first, second), n in pair_counts.items():
        aa_pair_counts[CODON_TO_AA[first], CODON_TO_AA[second]] += n

    scores = {}
    for (first, second), observed in pair_counts.items():
        x, y = CODON_TO_AA[first], CODON_TO_AA[second]
        f_first = codon_counts.get(first, 0) / total_codons
        f_second = codon_counts.get(second, 0) / total_codons
        f_x = aa_counts[x] / total_codons
        f_y = aa_counts[y] / total_codons
        f_xy = aa_pair_counts[x, y] / total_pairs
        if not (f_first and f_second and f_x and f_y and f_xy):
            continue
        expected = (f_first * f_second) / (f_x * f_y) * f_xy
        scores[first, second] = math.log((observed / total_pairs) / expected)
    return scores


def count_reference_cds(records) -> tuple[Counter[str], Counter[tuple[str, str]]]:
    """Count codons and adjacent codon pairs over an iterable of coding sequences.

    Skips sequences that are not a multiple of three or carry an internal stop;
    genome-wide CDS dumps contain partial entries and pseudogenes.
    """
    codons: Counter[str] = Counter()
    pairs: Counter[tuple[str, str]] = Counter()

    for record in records:
        seq = record.upper().replace("U", "T")
        if len(seq) < 6 or len(seq) % 3:
            continue
        triplets = [seq[i : i + 3] for i in range(0, len(seq), 3)]
        if triplets[-1] in ("TAA", "TAG", "TGA"):
            triplets = triplets[:-1]
        if any(triplet not in SENSE_CODONS for triplet in triplets):
            continue
        codons.update(triplets)
        pairs.update(zip(triplets, triplets[1:], strict=False))

    return codons, pairs


def _constraints(spec: dict) -> Constraints:
    windows = tuple(
        GCWindow(
            width=int(window["width"]),
            minimum=float(window.get("minimum", 0.0)),
            maximum=float(window.get("maximum", 1.0)),
        )
        for window in spec.get("gc_windows", [])
    )
    return Constraints(
        gc_range=tuple(spec["gc_range"]),
        gc_windows=windows,
        max_homopolymer={base.upper(): int(n) for base, n in spec["max_homopolymer"].items()},
        max_repeat_length=int(spec.get("max_repeat_length", 10)),
        forbidden_sites={name: site.upper() for name, site in spec.get("forbidden_sites", {}).items()},
    )


def _initiation(spec: dict) -> Initiation:
    return Initiation(
        model=spec["model"],
        window=tuple(spec["window"]),
        global_mfe_sign=int(spec["global_mfe_sign"]),
        utr=spec.get("utr", "").upper().replace(" ", ""),
        kozak=spec.get("kozak"),
    )


def read_codon_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    for fields in _read_table(path, columns=2):
        counts[fields[0].upper()] = int(fields[1])
    return counts


def read_pair_counts(path: Path) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    for fields in _read_table(path, columns=3):
        counts[fields[0].upper(), fields[1].upper()] = int(fields[2])
    return counts


def write_counts(path: Path, rows) -> None:
    """Write a count table, most frequent first, with a provenance header."""
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("# generated by janus build-host; see the accompanying YAML for the source\n")
        for row in rows:
            fh.write("\t".join(str(value) for value in row) + "\n")


def _read_table(path: Path, columns: int):
    if not path.exists():
        raise FileNotFoundError(f"count table {path} is missing")
    with path.open(encoding="utf-8") as fh:
        for number, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) != columns:
                raise ValueError(f"{path}:{number} has {len(fields)} fields, expected {columns}")
            yield fields
