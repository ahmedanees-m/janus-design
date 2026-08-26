
# The exchange rate in expression units

120 designed backbones, delta 1.0 nats, host *E. coli* BL21(DE3), 20 August 2026.

The exchange rate is measured in kcal/mol of initiation-window folding energy
bought per nat of unconditional fold log-likelihood given up. Neither unit is
useful to someone deciding whether to order a gene. This converts the first.

## The conversion

The Ribosome Binding Site Calculator relates translation initiation rate to the
total free energy of initiation as rate proportional to exp(-beta dG_total). In
an ideal system beta would be 1/RT, about 1.62 mol/kcal at 37 C. Fitted against
measured expression inside *E. coli* it is **0.45 plus or minus 0.05 mol/kcal**
(Salis, Mirsky and Voigt, *Nature Biotechnology* 2009), the difference reflecting
the non-ideal crowded interior of the cell.

## What the shell buys

Median over 120 backbones. Weight is on the spread-normalised scale, where the
defensible operating region is 0.125 to 0.5. Budget is the Tier-1 score given up
as a share of the backbone's total marginal entropy.

| weight | window opened | Tier-1 paid | predicted fold change | range | budget spent |
|---|---|---|---|---|---|
| 0.125 | +0.60 kcal/mol | 0.0067 nats | 1.31x | 1.27 to 1.35 | 0.007% |
| 0.250 | +2.10 kcal/mol | 0.0344 nats | **2.57x** | 2.32 to 2.86 | **0.037%** |
| 0.500 | +2.65 kcal/mol | 0.0778 nats | 3.30x | 2.89 to 3.76 | 0.084% |
| 1.000 | +2.65 kcal/mol | 0.1003 nats | 3.30x | 2.89 to 3.76 | 0.108% |
| 4.000 | +2.65 kcal/mol | 0.1036 nats | 3.30x | 2.89 to 3.76 | 0.112% |

At the 90th percentile of what the shell can open:

| weight | window opened | Tier-1 paid | predicted fold change |
|---|---|---|---|
| 0.125 | +4.31 kcal/mol | 0.0915 nats | 6.96x |
| 0.250 | +5.01 kcal/mol | 0.1490 nats | 9.53x |
| 0.500 | +5.60 kcal/mol | 0.1875 nats | 12.43x |

At the operating point, **0.037% of the backbone's entropy budget buys a
predicted 2.6-fold increase in initiation rate**, and on the most tractable tenth
of backbones 0.15% buys a predicted 9.5-fold.

The window saturates at 2.65 kcal/mol. Past a normalised weight of 0.5 the extra
weight buys no further energy and only costs Tier-1 score, which is a better
argument for the operating region than the sweep alone gives.

## Three ways this is an extrapolation, not a measurement

**The energy term is not the same term.** The RBS Calculator's dG_total is a
composite of ribosome-mRNA hybridisation, spacing, standby site and the cost of
unfolding mRNA structure. What is measured here is the folding energy of the
initiation window alone. Applying beta to it assumes that window's energy enters
dG_total roughly one for one. That is what the model intends, but it is not an
equivalence anyone has measured for this particular term, and if the window's
contribution is partial the fold-change is an overestimate.

**It is *E. coli* only.** beta was fitted there. Nothing here transfers to
HEK293, whose initiation is cap-scanning rather than Shine-Dalgarno, and no
equivalent fitted coefficient is used for that host.

**Initiation rate is not yield.** Initiation is rate-limiting for expression over
a wide range but not everywhere. Where transcription, folding, degradation or
toxicity becomes limiting, a predicted 2.6-fold gain in initiation delivers less
than 2.6-fold of protein, and possibly none.

The number belongs in the discussion with these caveats attached, not on a figure
axis as though it had been measured.

## What it does and does not support

It supports the claim that the amino-acid freedom inside a one-nat shell is
cheap relative to what it buys at the gene layer: a fraction of a percent of the
available sequence entropy for a predicted several-fold change in a quantity the
literature agrees is the dominant coding-sequence determinant of bacterial
expression.

It does not support a claim about protein yield, and it does not by itself
justify the joint lattice. The decomposition showed 87.6% of this particular gain
is reachable by synonymous change alone, so most of the 2.6-fold is an argument
for a good codon optimiser rather than for joint design. The case for the
amino-acid axis rests on the protein-level liabilities, where the codon layer's
achievable gain is exactly zero.

## Files

- `analysis/scripts/expression.py`
- `janus-data/processed/expression.json`
- `janus-data/processed/rescore_recall_normalised.json`
