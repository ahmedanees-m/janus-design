"""The feature panel of analysis_plan.md, computed for one record at a time.

A record is a residue sequence, optionally a backbone for accessibility and
secondary structure, and optionally a coding sequence. Controls have no
backbone, so structure-dependent features are absent rather than zero and the
accessibility weighting falls back to unweighted.
"""

from __future__ import annotations

import math
import re
from collections import Counter

import numpy as np

from janus.genetic_code import gc_fraction
from janus.objectives import mrna
from janus.objectives.proteostasis import load_by_family, scan
from janus.objectives.synthesis import violations

# Kyte and Doolittle 1982.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
# Bjellqvist side-chain pKa values, as used by ProtParam.
PKA_POS = {"K": 10.0, "R": 12.0, "H": 5.98}
PKA_NEG = {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
N_TERM_PKA, C_TERM_PKA = 7.5, 3.55

EXPOSED = 0.25
HYDROPHOBIC = 1.8
LC_WINDOW = 12
LC_ENTROPY = 1.5

# N-end rule classes, Varshavsky. Glycine is separated out after Timms 2019.
N_END = {
    **{r: "type1_primary" for r in "RKH"},
    **{r: "type2_primary" for r in "FWYLI"},
    **{r: "secondary" for r in "NQDE"},
    "C": "tertiary",
    "G": "glycine",
}


def charge_at(protein, ph=7.4):
    total = 1.0 / (1.0 + 10 ** (ph - N_TERM_PKA))
    total -= 1.0 / (1.0 + 10 ** (C_TERM_PKA - ph))
    counts = Counter(protein)
    for residue, pka in PKA_POS.items():
        total += counts[residue] / (1.0 + 10 ** (ph - pka))
    for residue, pka in PKA_NEG.items():
        total -= counts[residue] / (1.0 + 10 ** (pka - ph))
    return total


def isoelectric_point(protein):
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        if charge_at(protein, mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def low_complexity_fraction(protein):
    if len(protein) < LC_WINDOW:
        return 0.0
    flagged = 0
    for start in range(len(protein) - LC_WINDOW + 1):
        window = protein[start : start + LC_WINDOW]
        counts = np.array(list(Counter(window).values()), dtype=float)
        p = counts / counts.sum()
        if -(p * np.log(p)).sum() < LC_ENTROPY:
            flagged += 1
    return flagged / (len(protein) - LC_WINDOW + 1)


def longest_run(flags):
    best = current = 0
    for flag in flags:
        current = current + 1 if flag else 0
        best = max(best, current)
    return best


def longest_homopolymer(cds, bases):
    best = 0
    for match in re.finditer(r"(.)\1*", cds):
        if match.group(1) in bases:
            best = max(best, len(match.group(0)))
    return best


def longest_repeat(cds, cap=40):
    for length in range(min(cap, len(cds) // 2), 3, -1):
        seen = set()
        for start in range(len(cds) - length + 1):
            kmer = cds[start : start + length]
            if kmer in seen:
                return length
            seen.add(kmer)
    return 0


def repeat_count(cds, length=20):
    seen, hits = set(), 0
    for start in range(len(cds) - length + 1):
        kmer = cds[start : start + length]
        if kmer in seen:
            hits += 1
        seen.add(kmer)
    return hits


def max_gc_window(cds, width):
    if len(cds) < width:
        return gc_fraction(cds)
    return max(gc_fraction(cds[i : i + width]) for i in range(len(cds) - width + 1))


def internal_sd(cds, skip=40):
    """AGGAGG allowing one mismatch, past the initiation window."""
    target = "AGGAGG"
    hits = 0
    for start in range(skip, len(cds) - len(target) + 1):
        window = cds[start : start + len(target)]
        if sum(a != b for a, b in zip(window, target, strict=True)) <= 1:
            hits += 1
    return hits


def protein_features(protein, classes, accessibility=None, sse=None):
    out = {}
    hits = scan(protein, classes, accessibility)
    for weighted, suffix in [(False, "_raw"), (True, "_weighted")]:
        for family, value in load_by_family(hits, weighted).items():
            out[f"{family}{suffix}"] = value

    out["n_end_class"] = N_END.get(protein[0], "stabilising") if protein else "none"
    out["c_end_residue"] = protein[-1] if protein else "none"

    if accessibility is not None:
        usable = [a for a in accessibility if a is not None]
        out["mean_rsa"] = float(np.mean(usable)) if usable else math.nan
        out["n_term_rsa"] = accessibility[0] if accessibility[0] is not None else math.nan
        out["c_term_rsa"] = accessibility[-1] if accessibility[-1] is not None else math.nan
        exposed_hydrophobic = [
            a is not None and a > EXPOSED and HYDROPATHY.get(r, 0) > HYDROPHOBIC
            for r, a in zip(protein, accessibility, strict=True)
        ]
        out["longest_exposed_hydrophobic_run"] = longest_run(exposed_hydrophobic)
        out["exposed_hydrophobic_area"] = float(sum(
            a for r, a in zip(protein, accessibility, strict=True)
            if a is not None and a > EXPOSED and HYDROPATHY.get(r, 0) > HYDROPHOBIC
        ))
        out["exposed_hydrophobic_patches"] = sum(
            1 for match in re.finditer(r"1{2,}", "".join("1" if f else "0" for f in exposed_hydrophobic))
        )

    if sse is not None:
        out["helix_fraction"] = sse.count("a") / len(sse)
        out["strand_fraction"] = sse.count("b") / len(sse)
        out["coil_fraction"] = sse.count("c") / len(sse)

    out["gravy"] = float(np.mean([HYDROPATHY.get(r, 0.0) for r in protein]))
    out["net_charge"] = charge_at(protein)
    out["isoelectric_point"] = isoelectric_point(protein)
    out["free_cysteines"] = protein.count("C")
    out["low_complexity_fraction"] = low_complexity_fraction(protein)
    out["length"] = len(protein)
    return out


def gene_features(cds, host):
    out = {
        "cai": 0.0,
        "codon_pair_score": 0.0,
        "gc": gc_fraction(cds),
        "max_gc_20": max_gc_window(cds, 20),
        "max_gc_100": max_gc_window(cds, 100),
        "gc_first_60": gc_fraction(cds[:60]),
        "longest_at_homopolymer": longest_homopolymer(cds, "AT"),
        "longest_gc_homopolymer": longest_homopolymer(cds, "GC"),
        "longest_repeat": longest_repeat(cds),
        "repeats_over_20": repeat_count(cds, 20),
        "internal_sd_hits": internal_sd(cds),
        "restriction_hits": sum(
            1 for v in violations(cds, host) if v.kind == "forbidden_site"
        ),
    }
    adaptiveness = host.relative_adaptiveness
    pair_scores = host.codon_pair_scores
    codons = [cds[i : i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    weights = [math.log(adaptiveness[c]) for c in codons if adaptiveness.get(c, 0) > 0]
    if weights:
        out["cai"] = float(math.exp(sum(weights) / len(weights)))
    junctions = [pair_scores.get(p, 0.0) for p in zip(codons, codons[1:], strict=False)]
    if junctions:
        out["codon_pair_score"] = float(np.mean(junctions))
    if mrna.available():
        out["initiation_dg"] = mrna.initiation_energy(cds, host)
        out["transcript_mfe"] = mrna.global_mfe(cds, host)
    return out
