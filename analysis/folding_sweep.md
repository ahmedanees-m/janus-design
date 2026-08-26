# Folding-weight frontier

120 designed backbones, host *E. coli* BL21(DE3), delta 1.0 nats. Candidate pool
is the 200-path k-best prefix plus 200 uniform draws from the same shell with
codons optimised exactly.

## The weight was doing unit conversion

The Tier-1 score is in nats and the initiation term in kcal/mol, so adding them
directly makes the weight an implicit exchange rate rather than a preference.
Measured across the candidate pools:

| | median spread across the pool |
|---|---|
| Tier-1 score | 6.317 nats |
| initiation energy | 2.376 kcal/mol |

One kcal/mol of initiation is worth 2.67 nats of Tier-1 in spread terms. An
unnormalised weight of 1.0 was therefore weighting initiation at about **0.37 of
parity**, which is why earlier runs at that setting did not buy the term. Each
term is now divided by its pool spread before weighting, so a weight of one
trades one standard deviation for one standard deviation and transfers between
backbones of different length.

## Frontier

Cost and gain are measured against the Tier-1 optimum, which is the weight-zero
selection.

| weight | Tier-1 given up (nats) | initiation opened (kcal/mol) | rate (kcal/nat) | identity to Tier-1 optimum | selected candidate clean |
|---|---|---|---|---|---|
| 0 | 0.0000 | 0.0000 | | 1.000 | 50.0% |
| 0.125 | 0.0536 | 2.1900 | 40.9 | 0.987 | 48.3% |
| 0.25 | 0.0615 | 2.2108 | 2.63 | 0.986 | 47.5% |
| 0.5 | 0.1128 | 2.2583 | 0.93 | 0.983 | 47.5% |
| 1 | 3.8458 | 4.1933 | 0.52 | 0.828 | 53.3% |
| 2 | 8.0583 | 5.3183 | 0.27 | 0.658 | 65.8% |
| 4 | 10.3326 | 5.6750 | 0.16 | 0.567 | 75.0% |
| 8 | 11.3321 | 5.7525 | 0.08 | 0.532 | 78.3% |

**The rate column is not uniform and should not be read as one series.** The
entry at weight 0.125 is the average rate from the origin, because there is no
preceding grid point; every entry below it is the marginal rate between adjacent
rows. The apparent forty-fold drop between the first two rows is that change of
definition plus a genuinely steep knee, not an error. Read the first row as
average-from-origin and the rest as marginal.

The knee is sharp and early. A weight of 0.125 buys **2.19 kcal/mol of
initiation-window opening for 0.054 nats**, averaged from the origin at
40.9 kcal/mol per nat, while leaving 98.7 percent of the Tier-1 optimum's
residues unchanged. Going from 0.5 to 1.0 costs a further 3.7 nats to gain
1.9 kcal/mol.

The defensible operating region is a weight between 0.125 and 0.5, chosen from
the curve rather than after the fact. That brackets the region a reviewer
predicted from the unnormalised table.

## The exchange rate, restated

Against the entropy budget it spends from, 0.054 nats out of a median 93 nats of
design freedom per 43-residue chain is **0.058 percent of the budget** for
2.19 kcal/mol of window opening. An earlier draft quoted roughly 2.5 kcal/mol for
0.1 nats from an unnormalised run; the normalised frontier gives a better rate at
a principled operating point.

Converting that into an expression fold-change would make the number quotable,
and Kudla 2009 and Goodman 2013 both supply the necessary relationship between
5' folding energy and protein output. That extraction has not been done, and a
figure is not asserted here without it.

## Two things the frontier does not say

**Synthesisability is a filter, by design, not a term.** The selected candidate
is clean about half the time throughout the defensible region and the rise at
high weight is incidental, so the folding weight neither helps nor hurts here.
That is the intended architecture rather than a shortcoming: a six-base
restriction site spans three codons and a repeat check is global, so neither
decomposes over adjacent codon pairs and neither can enter the parse without the
nucleotide-level lattice. Applying them as a filter over the ranked pool is the
correct treatment of a hard constraint, and the cross-host case study reaches 100
percent that way.

The cost of that filter is list depth. About half the pool passes every
constraint at the operating weight, so the first clean candidate sits within the
first few ranks and a k-best list of a few hundred is ample. It is not free: a
method that returned a single path would fail the constraint on roughly half of
backbones, which is why k-best rather than best-path is the interface.

**This is measured against the Tier-1 optimum, not against a vendor baseline.**
The cross-host case study compares against the native design sequence with
CAI-maximal codons, and there JANUS still lost on initiation. Both are true and
they measure different things: moving from the native protein to the
marginal-optimal protein costs initiation energy, and the folding term recovers
part but not all of it. A matched-protein arm, holding the residue sequence at
the native design and optimising codons only, is needed before the vendor
comparison can be read as a statement about the folding term.
