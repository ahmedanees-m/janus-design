# Audit of the normalisation defect

A term with zero spread over its candidate pool was divided by a floor of 1e-9
rather than dropped, so a constant term could dominate every other by nine orders
of magnitude. The delta sweep was the first experiment to trip it, which means it
was live for an unknown period before that. This is the audit of what it reached.

## What was on which side of the fix

Every experiment that calls `pool_scales` was placed on one side of the fix by
the date of the run that produced its output file.

| experiment | side of the fix | re-run |
|---|---|---|
| folding-weight frontier | before | yes |
| layer decomposition | before | yes |
| amino-acid exchange rate | at the fix | no |
| case study, three arms | after | no |
| delta sweep | after | no |
| baseline comparison | after | no |

The two-tier bound never used `pool_scales` at all; it combined terms directly,
which is a separate defect handled below.

## What the re-runs found

**Nothing.** Both pre-fix experiments were re-run against the fixed code and
reproduce their published values to four decimal places. The frontier at weight
0.125 gives 0.0536 nats and 2.1900 kcal/mol before and after; the decomposition
gives 87.6% synonymous recovery before and after.

The reason is mechanical. The defect only bites when a term's spread over the
pool is numerically zero, and initiation-window folding energy over a pool of 400
distinct coding sequences never is. The delta sweep tripped it because it was the
first experiment to score a liability term that genuinely takes one value on every
candidate: a design carrying no low-complexity window has none anywhere in its
shell either.

The defect was real, was fixed, and reached nothing that was reported. The two
re-runs are what establish the last of those.

## Two other defects the audit turned up

**A weight of zero was reported as an energy of zero.** `rescore` skipped the
folding call when the initiation weight was zero, as an optimisation, and returned
0.0 for the reported energy. `Rescored.initiation_energy` is a reported field as
well as a scoring term, so a sweep measuring what its winner gained against the
weight-zero row read that row's entire energy as a gain: the frontier briefly
showed 6.69 kcal/mol opened for nothing at weight zero. Introduced alongside the fix above, and caught here. The energy is now always
measured, and a test asserts that a zero weight changes the ranking rather than
the reported energy.

**Three scripts averaged natural domains into a statement about designs.** The
frontier, the decomposition and the two-tier bound globbed a marginals directory
holding both the 447 designs and 415 natural domains, and took the first N in
sorted order. All three now filter to designs. In practice the runs had been
pointed at design-only inputs and the numbers do not move, but the scripts no
longer depend on being invoked carefully.

## The exchange rate had three values, and now has one

Three numbers were in circulation for the same quantity:

| source | weight | nats paid | kcal opened |
|---|---|---|---|
| frontier table | 0.125 | 0.0536 | 2.1900 |
| expression table | 0.125 | 0.0067 | 0.60 |
| expression table | 0.250 | 0.0344 | 2.10 |

None of them was wrong and none of them was comparable, for two reasons that had
nothing to do with the defect above.

**They rank different pools.** The frontier ranks the 200-path k-best prefix plus
200 shell draws. The expression table was derived from the two-tier experiment,
which ranks a 1000-path prefix alone. Shell draws reach further into protein
space, so at the same weight they open the window more and pay more for it.

**They took different averages.** The frontier reported means and the expression
table medians. What the shell can open is right-skewed, so the mean sits above the
typical backbone.

Both are now medians, both are reported with the pool named, and the pool is
recorded in the output file. The reconciled figures, 120 designed backbones:

| pool | weight | nats paid | kcal opened | fold change | budget |
|---|---|---|---|---|---|
| k-best prefix | 0.125 | 0.0078 | 1.450 | 1.92x | 0.008% |
| k-best prefix | **0.250** | **0.0362** | **2.000** | **2.46x** | **0.039%** |
| k-best prefix | 0.500 | 0.0735 | 2.500 | 3.08x | 0.079% |
| prefix and draws | 0.125 | 0.0492 | 2.000 | 2.46x | 0.053% |
| prefix and draws | 0.250 | 0.0566 | 2.000 | 2.46x | 0.061% |
| prefix and draws | 0.500 | 0.0606 | 2.000 | 2.46x | 0.065% |

**The headline is the k-best prefix at weight 0.25: 2.0 kcal/mol of
initiation-window opening for 0.036 nats, 0.039% of the entropy budget, a
predicted 2.5-fold change in initiation rate.** The prefix is the pool the
two-tier bound endorses at that weight, and the bound is what says so. Adding
shell draws buys the same 2.0 kcal/mol for 0.057 nats, which is the same exchange
reached less efficiently, and both are reported.

Above weight 1.0 the hybrid pool's cost jumps to 10.1 nats, because there the
shell draws start winning. That is the two-tier bound again, in exchange-rate
form.

## Files

- `src/janus/rescore.py`, `_spread` and the initiation reporting
- `tests/test_mrna.py`, the two regression tests
- `analysis/scripts/folding_sweep.py`, `decompose.py`, `rescore_recall.py`,
  `expression.py`
