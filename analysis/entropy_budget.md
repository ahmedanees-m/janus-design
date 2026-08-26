# Entropy budget

All 862 MegaScale AlphaFold backbones, 447 designed and 415 natural,
26 to 74 residues. Leverage measured on 150 designed backbones. 19 August 2026.

Secondary structure is biotite's P-SEA rather than DSSP, which needs a binary
unavailable in the container. Relative accessibility is Shrake and Rupley with
the theoretical maxima of Tien et al. 2013. The substitution is recorded here
rather than left implicit.

## How much freedom there is

| | designed | natural |
|---|---|---|
| backbones | 447 | 415 |
| median length | 43 | 57 |
| mean per-position entropy | 2.072 nats | 1.931 nats |
| positions effectively fixed (H below 0.1) | 0.0% | 0.0% |
| positions above 2 nats | 67.2% | 58.3% |
| mean shell at delta 0.5 | 2.04 | 1.92 |
| mean shell at delta 1.0 | 3.54 | 3.23 |
| mean shell at delta 2.0 | 7.04 | 6.31 |
| mean shell at delta 3.0 | 10.91 | 9.72 |
| median total freedom per chain | 93 nats | 110 nats |

H2 holds at corpus scale. A shell of 3.54 residues per position at delta 1.0
lands inside the 2 to 5 range the plan predicted, and the preliminary figure of
3.21 from 60 backbones of one topology was representative. No position on any
backbone is effectively fixed.

De novo backbones carry more sequence freedom per position than natural domains
of the same size class, on every measure: higher mean entropy, a larger shell at
every budget, and nine percentage points more positions above two nats. Natural
chains carry more freedom in total only because they are longer. As far as we
can establish this has not been measured before.

## Where the freedom sits

| | mean entropy |
|---|---|
| helix | 2.142 |
| strand | 1.960 |
| coil | 1.880 |
| buried, relative SASA below 0.15 | 1.300 |
| exposed, relative SASA above 0.40 | 2.293 |

Freedom tracks exposure strongly, Spearman 0.539 across 44,000 positions. Buried
positions carry a little over half the entropy of exposed ones, which is what an
inverse-folding posterior should do and is a useful check that the marginals
behave.

Positional dependence is weak. Spearman between entropy and distance from the
nearest terminus is -0.120; positions within three residues of a terminus
average 2.229 nats against 1.968 for the interior.

## H2b: freedom and leverage do not coincide

Leverage at a position is the range of the initiation term reachable by varying
that position alone with the rest held at the Tier-1 optimum.

Leverage is hard bounded. The initiation window runs from four bases before the
start codon to 37 after it, which covers the start codon and 34 bases of coding
sequence, so no position past codon 12 can move the term at all. Measured mean
leverage is 2.83 kcal/mol over positions 0 to 12 and exactly zero beyond.

Entropy over the same chains is flat: 2.095 nats over positions 0 to 14 against
2.051 beyond, a difference of 0.045 nats.

Inside the window the correlation between entropy and leverage has a median of
+0.221, against a permutation null with mean +0.083 and standard deviation
0.221. The median z is +0.69 and 6.7% of backbones exceed the null's 95th
percentile, against 5% expected by chance.

**H2b as stated fails.** It predicted that freedom and gene-level leverage would
coincide at the N terminus. Leverage is concentrated there because the window is
fixed at the 5' end; freedom is not concentrated anywhere.

The conclusion H2b was meant to support survives, and by a cleaner route. The
exchange rate does not need freedom to be concentrated where leverage is: it
needs enough freedom there, and uniform freedom at roughly two nats per position
supplies it. The separately measured exchange rate bears this out, at about
2.5 kcal/mol of initiation opening for a tenth of a nat of fold log-likelihood.
The prediction, the null, and the reason the conclusion holds anyway are all
reportable.
