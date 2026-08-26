# Does the atlas describe designs, or only the ones that worked

614 benchmark designs, 269 successful and 345 failed, ESMFold structures with a
pLDDT floor of 70, 20 August 2026.

The 447 backbones in the atlas come from published, experimentally characterised
sets, so every one survived whatever filtering its authors applied. If liability
load is part of what selection removes, the atlas understates what a generator
emits and the gap between designs and natural proteins is larger than measured.

The Garcia benchmark answers this directly, because it kept the failures. Splitting
it on the experimental outcome measures how much liability selection is removing.
Positive means the failed designs carry more of it.

## Where selection is doing work

| feature | delta | 95% CI | failed | passed | q |
|---|---|---|---|---|---|
| modification motifs, weighted | **+0.375** | [+0.29, +0.45] | 5.000 | 3.242 | 3.1e-14 |
| coil fraction | +0.199 | [+0.09, +0.31] | 0.250 | 0.226 | 0.0028 |
| net charge | +0.169 | [+0.08, +0.26] | -1.337 | -2.329 | 0.0028 |
| isoelectric point | +0.165 | [+0.08, +0.25] | 6.354 | 5.671 | 0.0028 |
| exposed hydrophobic area | +0.140 | [+0.03, +0.25] | 0.894 | 0.629 | 0.038 |
| all degrons, weighted | +0.120 | [+0.03, +0.21] | 0.260 | 0.179 | 0.024 |
| protease sites, raw | -0.189 | [-0.27, -0.10] | 6.000 | 7.000 | 0.00071 |
| targeting motifs, raw | -0.165 | [-0.25, -0.08] | 3.000 | 4.000 | 0.0028 |
| helix fraction | -0.159 | [-0.27, -0.05] | 0.486 | 0.505 | 0.022 |
| length | -0.155 | [-0.25, -0.06] | 94 | 101 | 0.0044 |

22 of 40 features are determined and 15 of those are higher in the failures. So
selection does remove liability load, and the atlas is biased in the direction
the concern predicts. The effects are small: the largest is 0.375 and most sit
between 0.10 and 0.20.

Three of these are not liabilities at all. Coil and helix fraction, and length,
are properties of what was attempted rather than of what was designed into it,
and shorter designs succeeding more often is a fact about the benchmark's
composition. Protease sites and targeting motifs run the wrong way, being higher
in the *successes*, which is a reminder that a motif count is not a liability
until something acts on it.

## Where selection is not doing work, which is the part that matters

| feature | delta | 95% CI | q | verdict |
|---|---|---|---|---|
| low-complexity fraction | -0.033 | [-0.11, +0.05] | 0.56 | establishes nothing |
| protein repeat | +0.077 | [-0.01, +0.16] | 0.15 | establishes nothing |
| longest exposed hydrophobic run | +0.003 | [-0.08, +0.08] | 0.98 | establishes nothing |
| all degrons, raw | +0.031 | [-0.05, +0.11] | 0.59 | establishes nothing |
| GRAVY | +0.040 | [-0.05, +0.13] | 0.52 | establishes nothing |
| free cysteines | -0.000 | [-0.04, +0.04] | 0.98 | establishes nothing |

**Low-complexity content and repeat structure are not selected on.** Their
intervals straddle zero and are narrow enough to exclude anything approaching the
effects above. Designs that failed carry them at the same rate as designs that
succeeded, so nothing in the filtering that produced the atlas's 447 backbones
removed them.

That is the specific answer the atlas needs. Its two robust findings, that
designs carry more low-complexity content and more repeat structure than whole
natural proteins, are the two features this test shows selection does not touch.
Those results are not survivor artefacts. They are also, separately, the two
liabilities the amino-acid axis clears outright for about a fifth of a nat.

The features that *are* selected on include exposed hydrophobic area, at +0.140,
and accessibility-weighted degron load, at +0.120. Both are atlas features, so
those two comparisons are survivor-biased and the true design-to-natural gap on
them is larger than the atlas measures. For exposed hydrophobic area that makes
the reported result conservative, since designs already look cleaner than natural
proteins there. For degrons it means the pooled null is measured on a filtered
population, which is one more reason the degron claim is made per topology or not
at all.

## What this cannot settle

A design that failed may have failed for reasons that also produce a poor
structural model, so the structure-derived features are entangled with the
outcome. That shows up directly: whether a design has a model above the pLDDT
floor is itself the second-largest effect in the panel, at -0.357. The
sequence-only features, which include low-complexity fraction and protein repeat,
carry no such entanglement, and those are the ones the conclusion rests on.

This is also one benchmark's notion of success, from one group's pipeline. It
measures the filtering that produced these 614, not filtering in general.

## Files

- `analysis/scripts/survivors.py`
- `janus-data/processed/survivors.json`
- `janus-data/interim/garcia_esmfold`, 614 models
