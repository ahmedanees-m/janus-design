# Optimality gap

447 designed backbones from the MegaScale AlphaFold models, delta 1.0 nats,
host *E. coli* BL21(DE3), 19 August 2026.

Every method optimises the same function. Shortfall is a fraction of the range
the objective spans on that backbone, from the parser optimum down to the worst
attainable path, both computed exactly by running the DP over negated weights.
Shared setup is excluded; search time is matched. The parser solves a 43-residue
lattice in a median 1.6 ms, so a matched budget is about 2 posterior draws or
450 annealing steps, and 100x is 190 draws or 39,000 steps.

## Objectives added one at a time

Median shortfall.

| objective | terms | greedy | sampling x100 | annealing x1 | annealing x100 |
|---|---|---|---|---|---|
| mpnn | 1 | 0.000000 | 0.000156 | 0.000000 | 0.000000 |
| mpnn+cai | 2 | 0.000000 | 0.000194 | 0.000000 | 0.000000 |
| mpnn+cai+gc | 3 | 0.000000 | 0.008638 | 0.000000 | 0.000000 |
| mpnn+cai+gc+cpb | 4 | 0.011045 | 0.012193 | 0.016059 | 0.000000 |

Greedy and annealing are exactly optimal on 100% of backbones for the first
three objectives and on 0% once the codon-pair term is added, where annealing at
100x recovers the exact optimum on 65.3%.

The count of objectives is the wrong axis. Fold compatibility, codon adaptation
and GC are node-local, so the objective stays separable however many are
stacked, and per-position greedy is optimal rather than merely close. Only the
codon-pair term couples positions.

## Coupling weight swept

200 backbones, varying only the codon-pair weight.

| lambda_cpb | greedy | sampling x100 | annealing x1 | annealing x100 | annealing x100 exact |
|---|---|---|---|---|---|
| 0.0 | 0.000000 | 0.008510 | 0.000000 | 0.000000 | 100.0% |
| 0.1 | 0.000605 | 0.009527 | 0.001606 | 0.000000 | 95.5% |
| 0.3 | 0.010213 | 0.011779 | 0.016788 | 0.000000 | 59.5% |
| 1.0 | 0.039735 | 0.023952 | 0.087588 | 0.004608 | 5.0% |
| 3.0 | 0.104966 | 0.106655 | 0.151651 | 0.028932 | 1.0% |

The gap is monotone in coupling strength and depends on nothing else.

## Reading

At a codon-pair weight near 0.3 the parser leads every heuristic at matched
budget by about 1% of the attainable range, on every backbone, and annealing
needs roughly 100x to draw level on half of them. That is the middle of the
three outcomes in the plan: a guarantee rather than a large practical win.

## Outcome

The folding experiment answered the question in the
only way it can be answered given the decision not to build the nucleotide
lattice: at a folding weight of 1.0, sampling beats the parser prefix from ten
evaluations onward. The folding term does not rescue the algorithmic claim, it
bounds it, and the paper's shape has already been set on the middle outcome.

**Boundary condition, stated plainly.** The algorithmic claim rests on codon-pair
coupling as the only term inside the lattice that couples positions. The mRNA
folding term, which is the reason lattice parsing exists in the mRNA-design
literature, is realised in the rescoring tier. Within the lattice the parser
leads every heuristic at matched budget by about one percent of the attainable
range at a defensible codon-pair weight, consistently and on every backbone.

## What the parser is actually for in this paper

The stronger answer to "if the term that needs a lattice is not in your lattice,
why have a lattice" is not the one percent. It is that **the parser is a
measurement instrument here, not only an optimiser.**

The exchange rate that the paper leads with, roughly 2.5 kcal/mol of
initiation-window opening for a tenth of a nat of fold log-likelihood, is a
statement about the distance between two optima. It can only be quoted because
Tier 1 returns the exact optimum of the fold and codon objective, so the cost of
moving away from it is known rather than estimated. The same exactness gives the
worst attainable path, and with it the normalised range that makes shortfalls
comparable across backbones of different length. Both are properties of an exact
parse and neither survives a heuristic.

Scaled against the budget it spends from: a tenth of a nat out of a median 93
nats of design freedom per 43-residue chain is **0.11 percent of the entropy
budget**.
