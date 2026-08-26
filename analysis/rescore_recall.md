# Rescoring recall

120 designed backbones, delta 1.0 nats, host *E. coli* BL21(DE3), budget 1000
folding evaluations per arm. First run 19 August 2026, re-run on the normalised
weight scale 20 August 2026.

mRNA folding cannot enter the parser: base pairs form at arbitrary range, so a
path's folding energy is not a sum over its positions. The two-tier design has
the parser propose and folding rank. That is worth doing only if a Tier-1 prefix
of size N holds better folding candidates than N draws from the same shell.

Both arms use the parser. The sampling arm draws a residue sequence from the
tempered marginals within the shell and then optimises its codons exactly with
the same DP. The comparison is therefore between two ways of exploring the
amino-acid layer, not between using the DP and not using it.

## The weight scale, corrected

The first run combined the Tier-1 score in nats with the initiation energy in
kcal/mol directly, so `lambda` silently carried the unit conversion and its
numeric value meant nothing outside this one experiment. The rescoring tier was
normalised by pool spread afterwards, and the folding sweep and the operating
region of 0.125 to 0.5 are quoted on that normalised scale. The two tables were
therefore on different axes.

Measured over the ranked pool, the spread ratio Tier-1 to initiation has median
0.2752 (10th percentile 0.1577, 90th 0.4738). The first run's grid of 0.1 / 0.3 /
1.0 lands at roughly **0.36 / 1.09 / 3.63** on the normalised scale. Its
strongest setting sat about seven times above the top of the defensible operating
region, and only its weakest setting was inside it. Everything below is re-run
with `lambda` read on the normalised scale.

## Parser prefix against sampling, at matched folding evaluations

Median combined score, and the share of backbones on which the parser prefix
wins. Scores are unitless because both terms are divided by their pool spread,
so compare within a column and not across the table.

| budget | 0.125 | | 0.25 | | 0.5 | | 1.0 | | 2.0 | | 4.0 | |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| | median | wins | median | wins | median | wins | median | wins | median | wins | median | wins |
| 1 | -111.79 | 100.0% | -112.51 | 100.0% | -113.60 | 95.8% | -115.63 | 85.0% | -119.80 | 71.7% | -129.99 | 64.2% |
| 10 | -111.69 | 95.8% | -112.40 | 94.2% | -113.60 | 84.2% | -115.44 | 56.7% | -118.88 | 35.0% | -128.24 | 29.2% |
| 100 | -111.65 | 87.5% | -112.30 | 86.7% | -113.60 | 62.5% | -115.10 | 45.8% | -118.75 | 25.0% | -126.19 | 17.5% |
| 1000 | -111.65 | 67.5% | -112.30 | 70.0% | -113.25 | 60.0% | -114.93 | 42.5% | -117.46 | 27.5% | -122.84 | 20.0% |

## Where the two-tier design holds

Across the whole defensible operating region, 0.125 to 0.5, the parser prefix is
the better candidate pool at every budget from one evaluation to a thousand,
winning on 60% to 100% of backbones. Its margin narrows as the budget grows, as
it must, since both arms converge on the same shell.

The crossover is at a normalised weight of about 1. At 1.0 the parser is ahead up
to roughly 25 evaluations and sampling takes over beyond that. At 2.0 and 4.0
sampling is clearly better from ten evaluations on.

The mechanism is unchanged: k-best enumerates small perturbations around the
Tier-1 optimum, so it searches a narrow neighbourhood deeply, while sampling
searches the shell broadly and shallowly. Once the non-decomposable term is
strong enough to move the optimum away from the Tier-1 neighbourhood, breadth is
what pays.

This bounds the two-tier architecture rather than refuting it, and the corrected
scale moves the bound outward rather than inward. Generate and rank is sound
while the terms outside the parser are secondary to those inside it, which is the
entire operating region. It is the wrong shape when they dominate, which begins
at roughly twice the top of that region, and there the options are broad sampling
with exact codon optimisation, or annealing on the combined objective.

The earlier text described the architecture as failing at `lambda` 1.0 without
saying that 1.0 was 3.63 in the units the rest of the analysis uses. The failure
is real and stays in the limitations, at the weight where it actually happens.

## How deep the list has to be

| lambda | median rank of the winner | 90th pct | winner is the Tier-1 top path |
|---|---|---|---|
| 0.125 | 2 | 289 | 48.3% |
| 0.25 | 11 | 612 | 29.2% |
| 0.5 | 124 | 718 | 19.2% |
| 1.0 | 196 | 792 | 13.3% |
| 2.0 | 215 | 831 | 11.7% |
| 4.0 | 222 | 831 | 11.7% |

Taking the Tier-1 best and folding it is not enough anywhere. Even at the bottom
of the operating region the top path is the folding winner on under half of
backbones, and by the top of it on a fifth.

## Exchange rate

What opening the initiation window costs in fold compatibility.

| lambda | initiation opened | Tier-1 paid |
|---|---|---|
| 0.125 | +0.600 kcal/mol | 0.0067 nats |
| 0.25 | +2.100 kcal/mol | 0.0344 nats |
| 0.5 | +2.650 kcal/mol | 0.0778 nats |
| 1.0 | +2.650 kcal/mol | 0.1003 nats |
| 4.0 | +2.650 kcal/mol | 0.1036 nats |

Median values. The window saturates at about 2.65 kcal/mol of opening: past a
normalised weight of 0.5 the extra weight buys almost no further energy and only
costs Tier-1 score. At the operating point 2.1 kcal/mol costs 0.034 nats, and the
90th percentile opens 5.01 kcal/mol for 0.149. The shell is large enough to buy
the gene-level term almost free, which is what H2 asserts, and the saturation
point is a better argument for the operating region than the sweep alone.

## Files

- `analysis/scripts/rescore_recall.py`, `--normalised`
- `janus-data/processed/rescore_recall.json`, first run, raw units
- `janus-data/processed/rescore_recall_normalised.json`, 720 rows
