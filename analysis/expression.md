
# The exchange rate in expression units

120 designed backbones, delta 1.0 nats, host *E. coli* BL21(DE3), 21 August 2026.

The exchange rate is measured in kcal/mol of initiation-window folding energy
bought per nat of unconditional fold log-likelihood given up. Neither unit is
directly interpretable as expression, so the folding energy is converted into a
predicted fold-change in initiation rate.

## The conversion

The Ribosome Binding Site Calculator relates translation initiation rate to the
total free energy of initiation as rate proportional to exp(-beta dG_total). In
an ideal system beta would be 1/RT, about 1.62 mol/kcal at 37 C. Fitted against
measured expression inside *E. coli* it is **0.45 plus or minus 0.05 mol/kcal**
(Salis, Mirsky and Voigt, *Nature Biotechnology* 2009), the difference reflecting
the non-ideal crowded interior of the cell.

## What the shell buys

The exchange rate depends on the candidate pool being ranked as well as on the
weight, so each table names its pool. Weight is on the spread-normalised scale,
where the operating region is 0.125 to 0.5. Budget is the Tier-1 score given up
as a share of the backbone's total marginal entropy, 92.58 nats at the median.
All values are medians over the 120 backbones.

The k-best prefix, which is the pool the two-tier bound endorses across the
operating region and the pool the reported rate is quoted on:

| weight | window opened | Tier-1 paid | predicted fold change | range | budget spent |
|---|---|---|---|---|---|
| 0.125 | +1.45 kcal/mol | 0.0078 nats | 1.92x | 1.79 to 2.06 | 0.008% |
| 0.250 | +2.00 kcal/mol | 0.0362 nats | **2.46x** | 2.23 to 2.72 | **0.039%** |
| 0.500 | +2.50 kcal/mol | 0.0735 nats | 3.08x | 2.72 to 3.49 | 0.079% |
| 1.000 | +2.50 kcal/mol | 0.0985 nats | 3.08x | 2.72 to 3.49 | 0.106% |
| 2.000 | +2.60 kcal/mol | 0.1153 nats | 3.22x | 2.83 to 3.67 | 0.125% |
| 4.000 | +2.60 kcal/mol | 0.1160 nats | 3.22x | 2.83 to 3.67 | 0.125% |

That prefix plus uniform shell draws, which is the pool the folding-weight sweep
ranks:

| weight | window opened | Tier-1 paid | predicted fold change | range | budget spent |
|---|---|---|---|---|---|
| 0.125 | +2.00 kcal/mol | 0.0492 nats | 2.46x | 2.23 to 2.72 | 0.053% |
| 0.250 | +2.00 kcal/mol | 0.0566 nats | 2.46x | 2.23 to 2.72 | 0.061% |
| 0.500 | +2.00 kcal/mol | 0.0606 nats | 2.46x | 2.23 to 2.72 | 0.065% |
| 1.000 | +3.35 kcal/mol | 0.1016 nats | 4.52x | 3.82 to 5.34 | 0.110% |
| 2.000 | +5.10 kcal/mol | 10.1112 nats | 9.92x | 7.69 to 12.81 | 10.922% |
| 4.000 | +5.30 kcal/mol | 10.9403 nats | 10.86x | 8.33 to 14.15 | 11.818% |
| 8.000 | +5.35 kcal/mol | 11.6674 nats | 11.11x | 8.50 to 14.51 | 12.603% |

At the operating weight of 0.25 on the prefix pool, **0.039% of the backbone's
entropy budget buys a predicted 2.46-fold increase in initiation rate**. The
window saturates beyond a weight of 0.5, holding at 2.50 kcal/mol through 1.0 and
reaching 2.60 at the top of the sweep while the Tier-1 cost keeps rising.

The two pools reach the same 2.00 kcal/mol at the operating weight, the prefix
for 0.0362 nats and the hybrid for 0.0566. Above a weight of 1.0 the hybrid
pool's cost rises to 10.1 nats, where the shell draws begin to win, which is the
two-tier boundary expressed as an exchange rate.

## Three ways this is an extrapolation, not a measurement

**The energy term is not the same term.** The RBS Calculator's dG_total is a
composite of ribosome-mRNA hybridisation, spacing, standby site and the cost of
unfolding mRNA structure. What is measured here is the folding energy of the
initiation window alone. Applying beta to it assumes that window's energy enters
dG_total roughly one for one. That is what the model intends, but it is not an
equivalence that has been measured for this term, and if the window's
contribution is partial the fold-change is an overestimate.

**The coefficient is fitted in *E. coli*.** Nothing here transfers to HEK293,
whose initiation is cap-scanning rather than Shine-Dalgarno, and no equivalent
fitted coefficient is used for that host.

**Initiation rate is not yield.** Initiation is rate-limiting for expression over
a wide range but not everywhere. Where transcription, folding, degradation or
toxicity becomes limiting, a predicted 2.46-fold gain in initiation delivers less
than 2.46-fold of protein, and possibly none.

The figure is reported as a prediction of initiation rate and labelled as such
wherever it appears.

## What it does and does not support

It supports the claim that the amino-acid freedom inside a one-nat shell is
cheap relative to what it buys at the gene layer: a fraction of a percent of the
available sequence entropy for a predicted several-fold change in a quantity the
literature agrees is the dominant coding-sequence determinant of bacterial
expression.

It does not support a claim about protein yield, and it does not by itself
justify the joint lattice. The decomposition attributes 87.6% of this particular
gain to synonymous change alone, so most of the 2.46-fold is an argument for a
good codon optimiser rather than for joint design. The case for the amino-acid
axis rests on the protein-level liabilities, where the codon layer's achievable
gain is exactly zero.

## Files

- `analysis/scripts/expression.py`
- `janus-data/processed/expression_prefix.json`, the k-best prefix pool
- `janus-data/processed/expression_hybrid.json`, prefix plus shell draws
- `janus-data/processed/folding_sweep.json`, the per-backbone rows behind the
  hybrid pool
