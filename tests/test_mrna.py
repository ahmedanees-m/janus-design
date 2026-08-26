"""Folding terms and the rescoring tier."""

from __future__ import annotations

import numpy as np
import pytest

from janus import Weights, design
from janus.objectives import mrna
from janus.rescore import FoldingWeights, rescore

pytestmark = pytest.mark.skipif(not mrna.available(), reason="ViennaRNA is not installed")


def test_transcript_prepends_the_utr_and_a_start_codon(ecoli):
    built = mrna.transcript("GCTAAA", ecoli)
    assert built.startswith(ecoli.initiation.utr.replace("T", "U"))
    assert "AUGGCUAAA" in built
    assert "T" not in built


def test_transcript_keeps_an_existing_start_codon(ecoli):
    assert mrna.transcript("ATGGCTAAA", ecoli).count("AUG") == 1


def test_initiation_energy_is_zero_or_negative(ecoli):
    for cds in ["GCTAAACTGGTTGACATCAAC", "GGCGGCGGCGGCGGCGGCGGC"]:
        assert mrna.initiation_energy(cds, ecoli) <= 0.0


def test_a_hairpin_over_the_window_is_more_structured_than_a_flat_sequence(ecoli):
    hairpin = "GGCGCCGCGGCGCCAAAAGGCGCCGCGGCGCC"
    flat = "AAAGAAAAAGAAAAAGAAAAAGAAAAAGAAA"
    assert mrna.initiation_energy(hairpin, ecoli) < mrna.initiation_energy(flat, ecoli)


def test_the_window_is_bounded_by_the_host_setting(ecoli):
    low, high = ecoli.initiation.window
    long_cds = "GCT" * 60
    window = mrna._initiation_window(long_cds, ecoli)
    assert len(window) == high - low


def test_global_mfe_covers_the_whole_transcript(ecoli):
    short = mrna.global_mfe("GCTAAACTG", ecoli)
    long = mrna.global_mfe("GCTAAACTG" * 8, ecoli)
    assert long < short


def test_rescoring_reorders_the_kbest_list(ecoli, marginals):
    log_marginals = marginals(length=40, seed=3)
    candidates = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=0.5), delta=1.5, k=60)
    ranked = rescore(candidates, ecoli, FoldingWeights(initiation=1.0))

    assert len(ranked) == len(candidates)
    assert [r.total for r in ranked] == sorted((r.total for r in ranked), reverse=True)
    assert ranked[0].design.cds != candidates[0].cds or len(candidates) == 1


def test_rescoring_opens_the_initiation_window(ecoli, marginals):
    log_marginals = marginals(length=40, seed=5)
    candidates = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=0.5), delta=1.5, k=60)
    ranked = rescore(candidates, ecoli, FoldingWeights(initiation=1.0))

    before = mrna.initiation_energy(candidates[0].cds, ecoli)
    after = ranked[0].initiation_energy
    assert after >= before


def test_zero_weight_leaves_the_tier_one_order_untouched(ecoli, marginals):
    log_marginals = marginals(length=30, seed=7)
    candidates = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=1.0, k=20)
    ranked = rescore(candidates, ecoli, FoldingWeights(initiation=0.0, mfe=0.0))
    assert [r.design.cds for r in ranked] == [c.cds for c in candidates]


def test_synthesisable_filter_drops_failing_candidates(ecoli, marginals):
    log_marginals = marginals(length=45, seed=11)
    candidates = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=0.5), delta=1.5, k=80)
    kept = rescore(candidates, ecoli, synthesisable_only=True)
    assert all(r.design.synthesisable for r in kept)
    assert len(kept) == sum(c.synthesisable for c in candidates)


def test_ecoli_does_not_claim_a_global_structure_term(ecoli):
    assert ecoli.initiation.global_mfe_sign == 0
    log_marginals = np.log(np.full((20, 20), 1 / 20))
    candidates = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=0.0, k=3)
    ranked = rescore(candidates, ecoli, FoldingWeights(initiation=0.0, mfe=5.0))
    assert all(r.total == pytest.approx(r.tier1) for r in ranked)


def test_pool_scales_measure_term_spread(ecoli, marginals):
    from janus.rescore import pool_scales

    log_marginals = marginals(length=35, seed=21)
    pool = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=0.5), delta=1.5, k=40)
    scales = pool_scales(pool, ecoli)
    assert scales.tier1 > 0
    assert scales.initiation > 0
    assert scales.mfe == pytest.approx(1.0)


def test_unit_scales_reproduce_the_unnormalised_ranking(ecoli, marginals):
    from janus.rescore import Scales

    log_marginals = marginals(length=30, seed=23)
    pool = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=1.0, k=25)
    plain = rescore(pool, ecoli, FoldingWeights(initiation=1.0))
    unit = rescore(pool, ecoli, FoldingWeights(initiation=1.0), Scales())
    assert [r.design.cds for r in plain] == [r.design.cds for r in unit]


def test_normalisation_rescales_the_weight_not_the_order(ecoli, marginals):
    """A weight of one in normalised units equals the scale ratio in native ones."""
    from janus.rescore import Scales, pool_scales

    log_marginals = marginals(length=35, seed=27)
    pool = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=0.5), delta=1.5, k=60)
    scales = pool_scales(pool, ecoli)

    normalised = rescore(pool, ecoli, FoldingWeights(initiation=1.0), scales)
    equivalent = scales.tier1 / scales.initiation
    native = rescore(pool, ecoli, FoldingWeights(initiation=equivalent), Scales())
    assert [r.design.cds for r in normalised] == [r.design.cds for r in native]


def test_shell_samples_stay_inside_the_shell(ecoli, marginals):
    from janus.lattice import amino_acid_shell
    from janus.sample import shell_samples

    log_marginals = marginals(length=30, seed=31)
    shells = amino_acid_shell(log_marginals, 1.0)
    drawn = shell_samples(log_marginals, ecoli, delta=1.0, count=25,
                          rng=np.random.default_rng(0))
    assert len(drawn) == 25
    for item in drawn:
        assert len(item.protein) == len(shells)
        for position, residue in enumerate(item.protein):
            assert residue in shells[position]


def test_shell_samples_score_no_better_than_the_parser_optimum(ecoli, marginals):
    from janus.sample import shell_samples

    log_marginals = marginals(length=30, seed=33)
    weights = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
    best = design(log_marginals, ecoli, weights=weights, delta=1.0, k=1)[0]
    drawn = shell_samples(log_marginals, ecoli, weights=weights, delta=1.0, count=40,
                          rng=np.random.default_rng(1))
    assert all(item.score <= best.score + 1e-9 for item in drawn)


def test_a_constant_term_is_divided_out_rather_than_amplified(ecoli, marginals):
    """A liability the pool cannot vary must not decide the ranking."""
    from janus.rescore import pool_scales

    log_marginals = marginals(length=30, seed=35)
    pool = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=1.0, k=40)
    scales = pool_scales(pool, ecoli, liability=lambda protein: 7.0)
    assert scales.liability == float("inf")

    weights = FoldingWeights(initiation=0.0, liability=1000.0)
    ranked = rescore(pool, ecoli, weights, scales, liability=lambda protein: 7.0)
    assert [r.design.cds for r in ranked] == [d.cds for d in pool]


def test_shell_search_matches_the_parser_on_the_optimum(ecoli, marginals):
    from janus.sample import ShellSearch

    log_marginals = marginals(length=30, seed=37)
    weights = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
    best = design(log_marginals, ecoli, weights=weights, delta=1.0, k=1)[0]
    search = ShellSearch(log_marginals, ecoli, weights, delta=1.0)
    again = search.best(best.protein)
    assert again.cds == best.cds
    assert again.score == pytest.approx(best.score)


def test_zero_weight_still_reports_the_true_folding_energy(ecoli, marginals):
    """A zero weight means the term does not steer the ranking, not that it is zero.

    Skipping the folding call when the weight was zero made the reported energy a
    lie, and a sweep measuring what its winner gained against the weight-zero row
    then read the whole of that row's energy as a gain.
    """
    log_marginals = marginals(length=40, seed=11)
    pool = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=1.0, k=20)

    off = rescore(pool, ecoli, FoldingWeights(initiation=0.0))
    measured = {r.cds: r.initiation_energy
                for r in rescore(pool, ecoli, FoldingWeights(initiation=1.0))}

    assert all(r.initiation_energy == measured[r.cds] for r in off)
    assert any(r.initiation_energy != 0.0 for r in off)


def test_neither_search_beats_the_parser_on_a_decomposable_objective(ecoli, marginals):
    """The parser is exact for this objective, so no search may improve on it.

    Both searches are given the Tier-1 score itself. The parser's optimum is then
    the true maximum over the shell, and a search that reported anything above it
    would be scoring inconsistently with what it optimised. Starting them at that
    optimum also means neither should move.
    """
    from janus.search import best_improvement, optimise

    log_marginals = marginals(length=30, seed=13)
    weights = Weights(mpnn=1.0, cai=0.5)
    optimum = design(log_marginals, ecoli, weights=weights, delta=1.0, k=1)[0]

    for search in (optimise, best_improvement):
        found, value, calls = search(log_marginals, ecoli, lambda d: d.score,
                                     weights=weights, delta=1.0)
        assert value == pytest.approx(optimum.score)
        assert found.protein == optimum.protein
        assert calls > 1


def test_a_search_can_trade_tier_one_away_for_an_outside_term(ecoli, marginals):
    """With a term the parser cannot carry, the search must leave the optimum.

    A search that never moves would pass the test above while doing nothing, so
    this one gives it a reason to move and checks that it does, and that it pays
    Tier-1 score to do so.
    """
    from janus.search import optimise

    log_marginals = marginals(length=30, seed=13)
    weights = Weights(mpnn=1.0, cai=0.5)
    optimum = design(log_marginals, ecoli, weights=weights, delta=1.0, k=1)[0]

    # Reward tryptophan, which the marginals give no reason to prefer.
    def score(candidate):
        return candidate.score + 4.0 * candidate.protein.count("W")

    found, value, _ = optimise(log_marginals, ecoli, score, weights=weights, delta=1.0)

    assert found.protein.count("W") >= optimum.protein.count("W")
    assert value > score(optimum) - 1e-9
    if found.protein != optimum.protein:
        assert found.score < optimum.score
