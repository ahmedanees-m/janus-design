# Additivity, and what the assay's range does to it

230992 double mutants across 541 MegaScale parents, delta 1.0 nats,
20 August 2026.

JANUS proposes several substitutions at once and scores them independently, so
whether their measured effects add is a real question for the method. Each double
mutant is compared against the sum of its two singles on the same parent.

## The bias, as first measured

| set | n | median bias | spread | r |
|---|---|---|---|---|
| all pairs | 230992 | +0.787 | 1.785 | 0.775 |
| natural parents | 208654 | +0.850 | 1.767 | 0.770 |
| designed parents | 22338 | +0.323 | 1.737 | 0.628 |
| both inside the shell | 3120 | +0.597 | 1.884 | 0.774 |
| designed, inside the shell | 220 | +0.186 | 1.943 | 0.734 |

Deviations exceed 0.5 kcal/mol on 71.7% of pairs and 1.0 on 50.4%.

The pattern that mattered was the gradient: +0.850 on natural parents down to
+0.186 on designed parents inside the shell. That was read as epistasis being
milder for the milder substitutions JANUS makes, which is what an additivity
bound in H4 would need. It was equally consistent with a measurement artefact,
and this note settles which.

## The mechanism is imprecision, not clipping

The deposited free energies are not censored. Across all 1868872 measured values
the distribution runs smoothly from -15.28 to +17.76 with no pile-up at either
end, so the release reports extrapolated values rather than truncated ones. There
is no wall to hit.

What there is instead is loss of precision. K50 resolves unfolding over a finite
protease series, and a variant whose stability sits far outside the middle of
that series is extrapolated, with an error that grows with distance and cannot
run symmetrically once the variant is already fully cleaved. So the test is a
restriction, not a threshold: repeat the comparison on progressively narrower
central bands of the measured distribution and see whether the effect survives.

Central band of every measured value, 10th to 90th percentile: -1.34 to +4.28
kcal/mol.

## Bias against distance from the measured range

Headroom is the additive prediction above the bottom of that band.

| headroom | natural n | natural bias | designed n | designed bias |
|---|---|---|---|---|
| below | 63224 | +2.333 | 1512 | +1.650 |
| 0 to 0.5 | 14612 | +1.389 | 855 | +0.782 |
| 0.5 to 1 | 17238 | +1.100 | 1204 | +0.696 |
| 1 to 2 | 40491 | +0.719 | 3591 | +0.622 |
| above 2 | 73089 | +0.200 | 15176 | +0.091 |

Monotone in both arms, and an order of magnitude across the range. The bias is
largest exactly where the prediction is furthest into extrapolated territory.

## The test

Restricting so that the parent, the additive prediction and the observation all
fall inside a central band of the measured distribution.

| band | natural | designed | gap |
|---|---|---|---|
| 0 to 100 pct | +0.850 (n=208648) | +0.323 (n=22338) | +0.527 |
| 5 to 95 pct | +0.726 (n=159469) | +0.375 (n=19553) | +0.351 |
| 10 to 90 pct | +0.606 (n=113085) | +0.409 (n=15441) | +0.196 |
| 25 to 75 pct | +0.333 (n=26021) | +0.369 (n=3907) | -0.036 |

The gap closes monotonically and is gone in the narrowest band, where the
designed arm is if anything marginally higher than the natural one.

## What this changes

**The natural-against-designed gradient does not survive.** The claim that
epistasis is milder on designed parents, and milder still inside the shell, is a
range effect. In the well-measured middle the two arms are indistinguishable, and
the n=220 shell subset that carried the strongest version of the claim was never
significant on its own (p = 0.48). The pattern argument is withdrawn.

**A real non-additivity remains, and it is larger than the earlier reading
suggested.** In the narrowest band both arms sit at about +0.35 kcal/mol of
median positive deviation. That is not an artefact; it is present where the assay
measures best. So the substitutions JANUS makes are not additive to within
measurement error, and any bound in H4 has to carry roughly +0.35 kcal/mol of
systematic optimism per pair rather than the +0.186 the earlier subset suggested.

The direction matters for how it is used. Positive deviation means a double
mutant is more stable than the sum of its singles predicts, so a design scored
additively is being underestimated rather than overestimated. That makes the
independence assumption conservative for stability, which is the safer direction
to be wrong in, but it is an assumption with a measured size now rather than a
hoped-for one.

## Files

- `analysis/scripts/epistasis.py`
- `janus-data/processed/epistasis.json`, 230992 pairs

## Is the in-shell subset in tension with the designed estimate

It looked as though it might be. The in-shell designed subset gives a bias near
zero on n=220, while designed parents on the well-measured band give +0.37. Two
numbers for the same quantity, one apparently null.

They are not in tension. Putting an interval on both, on means so the two are the
same statistic:

| subset | n | bias | standard error | 95% interval |
|---|---|---|---|---|
| designed parents | 22338 | +0.166 | 0.012 | [+0.143, +0.189] |
| designed, inside the shell | 220 | -0.014 | 0.131 | [-0.272, +0.243] |

The intervals overlap. The in-shell interval is eleven times wider and contains
the designed estimate outright, so **the subset does not establish a smaller bias;
it is consistent with the larger one and too small to refine it**.

The earlier phrasing, that the in-shell subset was unbiased to plus or minus 0.26
kcal/mol, invited exactly the wrong reading. A wide interval around a small point
estimate is not evidence of absence, and stating it that way made an underpowered
subset look like a favourable finding. It was not one.

## Are shell substitutions milder

The favourable reconciliation on offer was that shell substitutions are the ones
ProteinMPNN considers plausible, so they should perturb stability less, and milder
substitutions should be more additive. Both halves are testable.

Severity is the larger of the two singles' absolute ddG, in kcal/mol.

| set | n | median | 75th | 90th |
|---|---|---|---|---|
| designed, both in shell | 220 | 1.336 | 2.312 | 3.265 |
| designed, not both in shell | 8663 | 1.470 | 2.076 | 3.298 |
| natural, both in shell | 2900 | 2.128 | 3.522 | 5.301 |
| natural, not both in shell | 149942 | 2.125 | 3.232 | 4.878 |

**Shell substitutions are not measurably milder.** On designed parents the median
difference is -0.134 kcal/mol in the favourable direction, p = 0.168. On natural
parents there is no difference at all. The first half of the reconciliation fails.

The second half holds. Deviation from additivity grows with how severe the
substitutions are:

| severity | n | median bias | spread |
|---|---|---|---|
| 0 to 0.5 | 3474 | +0.181 | 1.370 |
| 0.5 to 1 | 5745 | +0.107 | 1.255 |
| 1 to 2 | 8725 | +0.403 | 1.119 |
| above 2 | 4394 | +0.661 | 3.012 |

So milder substitutions are more additive, but the shell does not select milder
substitutions, and the mechanism that would have made this a point in the method's
favour is not there. It is stated as measured rather than left as the plausible
story it was.

## What the additive prediction is worth for a real design

The delta sweep moves 4.9 to 6.0 residues at the operating point, so the question
is the error on that many substitutions at once, not on one pair.

Deviation does not decay with how far apart the substitutions sit:

| separation, residues | n | median bias | spread |
|---|---|---|---|
| 0 to 5 | 10747 | +0.326 | 1.797 |
| 5 to 10 | 3642 | +0.398 | 1.582 |
| 10 to 20 | 4294 | +0.231 | 1.848 |
| 20 to 40 | 3653 | +0.347 | 1.498 |

That closes the escape route. If deviation fell away with distance, substitutions
scattered along a 43-residue chain would mostly not interact and the error would
grow far more slowly than the pair count. It does not fall away.

On designed parents the per-pair deviation has mean +0.166, median +0.323 and
standard deviation 1.737 kcal/mol. Two quantities follow and they are different
things: a systematic offset, and a spread around it.

The offset uses the **mean**, because the expectation of a sum is the sum of
expectations while the median of a sum is not the sum of medians. Multiplying the
median would overstate it by about a factor of two.

| substitutions | pairs | offset | 1 sigma | 95% half-width |
|---|---|---|---|---|
| 2 | 1 | +0.17 | 1.74 | 3.41 |
| 3 | 3 | +0.50 | 3.01 | 5.90 |
| 5 | 10 | +1.66 | 5.49 | 10.77 |
| 6 | 15 | +2.49 | 6.73 | 13.19 |

**One sigma is a spread, not a bound.** About a third of designs fall outside it,
so a five-substitution design's additive prediction sits within roughly plus or
minus 11 kcal/mol at 95%, around an offset of +1.7. An earlier version of this
note called the 5.49 figure a bound; it is not, and that wording is withdrawn.

All three columns assume every pair deviates independently. The separation table
above gives no reason to discount that, since deviation does not fall away with
distance.

## What H4 can and cannot say

The two halves have to be stated separately, because one is measured directly and
the other is propagated.

**H4 holds where it is measured.** Single substitutions inside the shell are
non-inferior to random draws from the same shell on measured stability, on
n = 58456 single mutants. That is a direct comparison on the assay's own numbers
and it stands at full strength.

**Per-design aggregate stability prediction is imprecise.** For the five or six
substitutions this method makes at the operating point, the additive prediction
carries an offset near +1.7 to +2.5 kcal/mol and a 95% half-width near 11 to 13.
On a domain whose whole folding free energy is a few kcal/mol, that is not a
footnote.

The direction is the safer one. Positive deviation means a design scored
additively is underestimated rather than overestimated, so the error runs toward
more stable than predicted. **JANUS can say its substitutions are individually
non-inferior. It cannot promise a stability-neutral design**, and the paper says
both.
