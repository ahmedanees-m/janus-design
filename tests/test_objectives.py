"""Objective terms and the synthesis constraint checker."""

from __future__ import annotations

import math

import numpy as np
import pytest

from janus.genetic_code import (
    CODON_TO_AA,
    CODONS,
    SENSE_CODONS,
    STOP_CODONS,
    SYNONYMOUS,
    translate,
)
from janus.hosts import codon_pair_scores, relative_adaptiveness
from janus.objectives import Weights, score_terms
from janus.objectives.codon import gc_content_term, relative_adaptiveness_term
from janus.objectives.synthesis import violations


def test_genetic_code_is_complete_and_consistent():
    assert len(CODONS) == 64
    assert len(SENSE_CODONS) == 61
    assert set(STOP_CODONS) == {"TAA", "TAG", "TGA"}
    assert sum(len(family) for family in SYNONYMOUS.values()) == 61
    assert CODON_TO_AA["ATG"] == "M"
    assert CODON_TO_AA["TGG"] == "W"
    for aa, family in SYNONYMOUS.items():
        assert all(CODON_TO_AA[codon] == aa for codon in family)


def test_translate_stops_at_the_first_stop_codon():
    assert translate("ATGGCTTAAATGGCT") == "MA"
    with pytest.raises(ValueError, match="multiple of three"):
        translate("ATGG")


def test_relative_adaptiveness_peaks_at_one_per_family():
    counts = {codon: 1 for codon in SENSE_CODONS}
    counts["CTG"] = 500
    values = relative_adaptiveness(counts)
    assert values["CTG"] == 1.0
    assert all(values[c] < 1.0 for c in SYNONYMOUS["L"] if c != "CTG")
    assert max(values[c] for c in SYNONYMOUS["P"]) == 1.0


def test_codon_pair_score_is_zero_when_pairs_are_used_as_expected():
    """Under independence every score should be zero. A term that is not is
    measuring codon bias twice rather than pair-specific preference."""
    codons = SYNONYMOUS["L"] + SYNONYMOUS["P"]
    counts = {codon: 100 for codon in codons}
    total = sum(counts.values())
    pairs = {
        (first, second): counts[first] * counts[second] // total
        for first in codons
        for second in codons
    }
    scores = codon_pair_scores(counts, pairs)
    assert scores
    assert all(abs(value) < 1e-9 for value in scores.values())


def test_codon_pair_score_is_positive_for_an_over_used_pair(ecoli):
    scores = ecoli.codon_pair_scores
    assert len(scores) > 3000
    assert any(value > 0 for value in scores.values())
    assert any(value < 0 for value in scores.values())
    assert all(math.isfinite(value) for value in scores.values())


def test_gc_term_counts_bases():
    assert gc_content_term(("AAA", "GCG", "GCC")).tolist() == [0.0, 1.0, 1.0]
    assert gc_content_term(("ATG",))[0] == pytest.approx(1 / 3)


def test_unobserved_codons_are_floored_not_dropped(ecoli):
    values = relative_adaptiveness_term(ecoli, SENSE_CODONS)
    assert len(values) == len(SENSE_CODONS)
    assert np.all(np.isfinite(values))
    assert np.all(values <= 0)


def test_score_terms_reports_cai_as_a_geometric_mean(ecoli):
    best = {aa: max(family, key=lambda c: ecoli.relative_adaptiveness[c]) for aa, family in SYNONYMOUS.items()}
    cds = "".join(best[aa] for aa in "MAKLVGT")
    marginals = np.log(np.full((len(cds) // 3, 20), 1 / 20))
    assert score_terms(cds, ecoli, marginals)["cai"] == pytest.approx(1.0)


def test_weights_default_to_fold_compatibility_only():
    weights = Weights()
    assert weights.mpnn == 1.0
    assert (weights.cai, weights.cpb, weights.gc) == (0.0, 0.0, 0.0)


def test_synthesis_flags_a_forbidden_site(ecoli):
    cds = "ATG" + "GAATTC" + "AAAGCTGTT"
    kinds = {v.kind for v in violations(cds, ecoli)}
    assert "forbidden_site" in kinds


def test_synthesis_flags_a_long_homopolymer(ecoli):
    cds = "ATG" + "A" * 12 + "GCTGTTAAC"
    found = [v for v in violations(cds, ecoli) if v.kind == "homopolymer"]
    assert found and "12 consecutive A" in found[0].detail


def test_synthesis_flags_a_direct_repeat(ecoli):
    unit = "ATGGCTGTTAACCTG"
    found = [v for v in violations(unit * 2, ecoli) if v.kind == "repeat"]
    assert found


def test_synthesis_flags_extreme_gc(ecoli):
    assert any(v.kind == "gc_content" for v in violations("GCGGCGGCGGCGGCG", ecoli))
    assert any(v.kind == "gc_content" for v in violations("ATAATAATAATAATA", ecoli))


def test_a_balanced_sequence_passes(ecoli):
    cds = "ATGGCTAAACTGGTTGACATCAACCTGTATGAAGGCACCCGTTGGAGCCAGTTCAAAGATATTCCGGTGAACCTGCATTCTGAA"
    assert violations(cds, ecoli) == []


def test_anchor_admits_a_protein_the_shell_would_exclude(ecoli, marginals):
    """A sequence outside the shell has to be reachable for a matched comparison."""
    from janus import design
    from janus.genetic_code import AA_ALPHABET
    from janus.lattice import amino_acid_shell, build

    log_marginals = marginals(length=25, seed=41)
    shells = amino_acid_shell(log_marginals, 0.5)
    outside = "".join(
        next(a for a in AA_ALPHABET if a not in shell and a != "*") for shell in shells
    )
    anchored = build(log_marginals, delta=0.5, anchor=outside)
    for position, residue in enumerate(outside):
        assert residue in anchored.amino_acids[position]
        for admitted in shells[position]:
            assert admitted in anchored.amino_acids[position]

    pinned = design(log_marginals, ecoli, delta=0.5, k=1,
                    fixed=dict(enumerate(outside)), anchor=outside)[0]
    assert pinned.protein == outside
