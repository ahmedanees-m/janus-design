# What the amino-acid freedom is worth, against how much of it there is

447 designed backbones, host *E. coli* BL21(DE3), Tier-1 weights mpnn 1.0 / cai
0.5 / cpb 0.3, folding weight 0.25 and liability weight 0.25 per term on the
spread-normalised scale, 20 August 2026.

Delta zero collapses the shell to the per-position marginal argmax, so the
lattice degenerates to a fixed-protein codon DP and the search has no residue
move available. Every delta above zero is the joint lattice. Holding the whole
objective fixed and sweeping delta prices the amino-acid axis directly.

The objective is the complete one: Tier 1 in the parser, exact; initiation-window
folding and the five protein-level liabilities optimised by hill climbing from
the Tier-1 optimum of the same shell, one residue at a time. Synthesis and
cis-element constraints are checked rather than weighted and reported as
violation counts. Every term is divided by its spread over one candidate pool
drawn from the widest shell in the sweep, so the arms share one scale.

## Attained objective

Gain is the spread-normalised objective above the fixed-protein arm, so it is
unitless.

| delta | mean shell | median gain | 10th | 90th | residues moved | wins | violations |
|---|---|---|---|---|---|---|---|
| 0.0 | 1.00 | 0.0000 | 0.0000 | 0.0000 | 0.00 | - | 0.83 |
| 0.5 | 2.05 | 1.5239 | 0.2869 | 14.04 | 4.89 | **100%** | 0.40 |
| 1.0 | 3.56 | 1.6915 | 0.4143 | 14.23 | 5.73 | **100%** | 0.37 |
| 2.0 | 7.08 | 1.7043 | 0.4679 | 14.29 | 5.94 | 100% | 0.38 |
| 3.0 | 10.96 | 1.7043 | 0.4751 | 14.29 | 6.00 | 100% | 0.38 |
| unbounded | 20.00 | 1.7043 | 0.4818 | 14.29 | 6.02 | 100% | 0.37 |

**The joint lattice beats the fixed-protein codon DP on every one of the 447
backbones, at every delta above zero.** That is the direct measurement the whole
project turns on, and it is not close.

**The freedom saturates early.** Half a nat captures 89% of the achievable gain
and one nat captures 99.2%. Going from delta 1 to an unbounded shell enlarges the
mean shell from 3.56 residues to all 20 and adds 0.013 to the median gain, under
one percent. The operating point chosen for the atlas is the right one, and it
was not chosen with this measurement in hand.

**The spread is enormous.** At delta 1 the tenth percentile gains 0.41 and the
ninetieth gains 14.23, a factor of 34. On some backbones the residue freedom is
worth very little and on others a great deal. A median is a poor summary of this
distribution and the paper should not lead with one alone.

## Where the freedom is spent

Median change from the fixed-protein arm, in native units.

| delta | Tier-1, nats | initiation, kcal/mol | low complexity | repeat | exposed hydrophobic | run | degron |
|---|---|---|---|---|---|---|---|
| 0.5 | -0.692 | +2.800 | 0.000 | -3.000 | -0.546 | 0.000 | 0.000 |
| 1.0 | -1.355 | +3.800 | 0.000 | -3.000 | -0.813 | 0.000 | 0.000 |
| 2.0 | -1.772 | +4.000 | 0.000 | -3.000 | -0.889 | 0.000 | 0.000 |
| 3.0 | -1.919 | +4.000 | 0.000 | -3.000 | -0.900 | -1.000 | 0.000 |
| unbounded | -1.969 | +4.000 | 0.000 | -3.000 | -0.900 | -1.000 | 0.000 |

At the operating point the freedom pays 1.36 nats of Tier-1 score and buys 3.8
kcal/mol of initiation window, a cleared repeat, and 0.81 of exposed hydrophobic
accessibility.

The zeros in the low-complexity and degron columns are medians, not absences. 57%
of designs carry a low-complexity window and 44% carry a degron, so the median
design carries neither and its median change is zero. What those terms are worth
on the designs that carry them is in `aa_axis.md`, where both clear outright for
about a fifth of a nat.

**Synthesis violations more than halve**, from 0.83 per design under the
fixed-protein DP to 0.37 at delta 0.5 and above. That is not a weighted term at
all; the constraints are checked, not optimised. The residue freedom removes them
as a side effect of moving off a protein whose codon options were too constrained
to avoid them. It is the most practical of the benefits here and the least
expected.

## Reading

The initiation column dominates the decomposition in raw magnitude, and 87.6% of
that gain is reachable by synonymous change alone. So the fixed-protein DP is not
being beaten mainly on initiation; it is being beaten on it because at delta 0
the protein is pinned to the marginal argmax rather than chosen, and a different
protein has a different set of codons available to open the window with.

The terms where the fixed-protein arm cannot compete at all, by construction, are
the protein-level liabilities. Those are worth less in normalised units but they
are the ones no codon optimiser reaches.

## Files

- `analysis/scripts/delta_sweep.py`
- `janus-data/processed/delta_sweep.json`, 2682 rows
