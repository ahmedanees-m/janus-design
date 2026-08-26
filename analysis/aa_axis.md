# What the amino-acid axis buys

447 designed backbones, delta 1.0 nats, host *E. coli* BL21(DE3), Tier-1 weights
mpnn 1.0 / cai 0.5 / cpb 0.3, ELM classes v1.4, 20 August 2026.

The 5' initiation exchange is real but mostly synonymous: 87.6% of the cheap gain
there is reachable without changing a residue, so it prices the codon layer more
than the joint lattice. This experiment asks the complementary question on terms
where the codon layer has no reach at all.

Five protein-level liabilities, all functions of the residue sequence and, where
structure enters, of a fixed per-residue accessibility from the backbone:

| term | what it counts |
|---|---|
| `low_complexity` | fraction of 12-residue windows below 1.5 nats of residue entropy |
| `protein_repeat` | length of the longest residue substring occurring twice |
| `exposed_hydrophobic` | summed relative accessibility over exposed hydrophobic residues |
| `exposed_hydrophobic_run` | longest consecutive run of them |
| `degron` | accessibility-weighted ELM degron motif load |

None depends on the codons. A fixed-protein codon optimiser's achievable
improvement on any of them is zero, and the design of the experiment is to
measure that zero rather than assert it.

## Two arms, one search

Both arms hill-climb the same shell from the same Tier-1 optimum, taking the
admitted single-residue substitution with the largest burden reduction per nat of
Tier-1 score given up. The joint arm may use every residue the shell admits. The
codon-only arm pins each position to the residue the Tier-1 optimum chose, so
only the codons can move. It runs the same code.

The first attempt ranked a fixed pool of uniform shell draws under a rising
weight instead. That produced a step function: the winner sat at the Tier-1
optimum until some distant random draw overtook it, then jumped 27 residues and
17 nats at once. The apparent exchange rate was a property of what happened to be
sampled. Searching from the optimum traces a path along the exchange, and each
weight reads off a point on that trace. The trace is greedy, so it is an
attainable exchange and a lower bound on what the shell holds, not the true
Pareto frontier; nothing below depends on it being optimal, since the comparison
is against an arm that reaches zero.

## What a liability costs to remove

Weight is on the spread-normalised scale. `removed` and `codon arm` are mean
burden removed; `nats` is mean Tier-1 score given up; `budget` is the mean share
of the shell's total marginal log-probability that the changed positions spend.

### low_complexity, carried by 257 of 447 designs, removable on 256 of those

| weight | removed | codon arm | nats | residues | budget | moved |
|---|---|---|---|---|---|---|
| 0.125 | 0.0646 | 0.0000 | 0.127 | 1.05 | 0.3% | 55% |
| 0.250 | 0.0710 | 0.0000 | 0.163 | 1.19 | 0.5% | 56% |
| 0.500 | 0.0745 | 0.0000 | 0.195 | 1.28 | 0.6% | 57% |
| 1.000 | 0.0764 | 0.0000 | 0.213 | 1.32 | 0.6% | 57% |
| 4.000 | 0.0766 | 0.0000 | 0.217 | 1.33 | 0.6% | 57% |

Exhausting the frontier clears **100% of the burden on the median carrier**, 10th
and 90th percentiles both 100%, for a median 0.21 nats and 90th percentile 0.93.

### degron, carried by 197 of 447, removable on 164

| weight | removed | codon arm | nats | residues | budget | moved |
|---|---|---|---|---|---|---|
| 0.125 | 0.1423 | 0.0000 | 0.046 | 0.30 | 0.2% | 27% |
| 0.250 | 0.1832 | 0.0000 | 0.095 | 0.39 | 0.3% | 34% |
| 0.500 | 0.1946 | 0.0000 | 0.113 | 0.42 | 0.4% | 36% |
| 1.000 | 0.1964 | 0.0000 | 0.119 | 0.43 | 0.4% | 37% |
| 4.000 | 0.1968 | 0.0000 | 0.121 | 0.43 | 0.4% | 37% |

Median 100% of the burden cleared, 10th and 90th percentiles both 100%, for a
median 0.22 nats.

### protein_repeat, carried by 335 of 447, removable on 228

| weight | removed | codon arm | nats | residues | budget | moved |
|---|---|---|---|---|---|---|
| 0.125 | 1.3400 | 0.0000 | 0.062 | 0.53 | 0.2% | 46% |
| 0.250 | 1.4944 | 0.0000 | 0.093 | 0.60 | 0.3% | 50% |
| 0.500 | 1.5280 | 0.0000 | 0.106 | 0.61 | 0.3% | 51% |
| 4.000 | 1.5280 | 0.0000 | 0.106 | 0.61 | 0.3% | 51% |

Median 100% cleared but a long lower tail, 10th percentile 25%, for 0.14 nats.
A third of carriers cannot be improved at all: the repeat is load-bearing under
the marginals, and no admitted substitution breaks it.

### exposed_hydrophobic, carried by 443 of 447, removable on 431

| weight | removed | codon arm | nats | residues | budget | moved |
|---|---|---|---|---|---|---|
| 0.125 | 0.3822 | 0.0000 | 0.067 | 0.74 | 0.0% | 53% |
| 0.250 | 0.6897 | 0.0000 | 0.219 | 1.37 | 0.4% | 76% |
| 0.500 | 1.1270 | 0.0000 | 0.645 | 2.32 | 1.7% | 90% |
| 1.000 | 1.4804 | 0.0000 | 1.284 | 3.19 | 3.7% | 96% |
| 4.000 | 1.5668 | 0.0000 | 1.540 | 3.46 | 4.5% | 96% |

The expensive one, and the only term whose exchange is still moving at weight 1.
Median 73% of the burden cleared, 10th percentile 43%, for a median 1.40 nats and
90th percentile 3.03. Unlike the others it is nearly universal and nearly always
partially removable, but rarely removable in full: exposed hydrophobic surface is
distributed across the chain rather than concentrated in one defect.

### exposed_hydrophobic_run, carried by 443 of 447, removable on only 142

| weight | removed | codon arm | nats | residues | budget | moved |
|---|---|---|---|---|---|---|
| 0.125 | 0.2729 | 0.0000 | 0.045 | 0.23 | 0.1% | 22% |
| 0.500 | 0.4094 | 0.0000 | 0.123 | 0.36 | 0.3% | 32% |
| 4.000 | 0.4161 | 0.0000 | 0.130 | 0.37 | 0.3% | 32% |

The weakest of the five. Only 32% of carriers can be improved and the median
improvement is half the run, for 0.31 nats.

## The codon-only arm

Across 3129 runs per term, seven weights on 447 backbones, the largest burden the
codon-only arm removed on any term was **0.000e+00**. That is 15645 measured
zeros. The arm is not a different objective or a different search; it is the same
hill climb with the residue moves withheld.

## Reading

Two terms are cleared outright at the operating point for about a fifth of a nat:
low-complexity on 100% of its 256 removable carriers, degrons on 100% of its 164.
That is 0.6% and 0.4% of the shell's entropy budget respectively, on roughly one
residue. Against an arm that removes exactly nothing.

This is the regime the joint lattice exists for. Where the initiation result made
the amino-acid axis look like a small addition on top of a codon optimiser, these
terms are not on the codon axis at all, and the exchange rate for clearing them
is cheap in fold compatibility.

## How often the condition is met

Term by term the coverage looks partial: 57% of designs carry a low-complexity
window, 44% a degron. Taken one at a time that reads as a conditional result with
a narrow condition. The union is the right denominator, because a design needs
only one removable liability for the codon layer's reach to be the thing that
matters.

| | designs | share |
|---|---|---|
| carry at least one of the five | 446 / 447 | 99.8% |
| at least one is removable | **443 / 447** | **99.1%** |
| at least two removable | 402 / 447 | 89.9% |
| at least three removable | 257 / 447 | 57.5% |
| at least one of the four besides exposed hydrophobic | 409 / 447 | 91.5% |

**99.1% of designs carry a liability that is present, removable, and outside the
codon layer's reach.** Excluding exposed hydrophobic surface, which is nearly
universal and might be argued to be a property of small proteins rather than a
defect, it is still 91.5%.

The condition is real and it is almost always met. The claim is not that the
amino-acid layer always pays; it is that when a protein-level liability is present
and removable the codon layer's achievable gain is zero and the joint lattice's is
not, at a price of well under a nat, and that this describes 443 of 447 designs
rather than a favourable subset.

## Files

- `analysis/scripts/aa_axis.py`
- `janus-data/processed/aa_axis.json`, 31290 rows
- `src/janus/objectives/liability.py`
- `src/janus/sample.py`, `ShellSearch`
