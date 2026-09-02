# Cross-host redesign, decomposed by layer

12 designs, three per topology family, two hosts, delta 1.0 nats, folding weight
0.25 and liability weight 0.25 per term on the spread-normalised scale,
20 August 2026.

Three arms on the same backbones:

1. **vendor**: the design's own residue sequence with the highest
   relative-adaptiveness codon at every position
2. **codon**: the same residue sequence with the JANUS codon layer, every
   position pinned
3. **joint**: a hill climb from arm 2 over the shell, one residue at a time,
   on the same combined objective

Arm 2 minus arm 1 is the codon layer with the protein held fixed. Arm 3 minus
arm 2 is what the residue freedom adds on top of that.

## Why the arms had to be rebuilt

The first version compared arm 1 against a joint search whose shell is centred on
the ProteinMPNN marginal argmax. That shell does not generally contain the
design's own sequence, so the two arms differed in two ways at once: the residue
freedom, and a different starting protein. The initiation term went the wrong way
under that comparison and was read as a cost of the amino-acid move.

Arm 3 now starts from arm 2 and searches a lattice anchored on the native
sequence, so it contains arm 2 by construction and the difference is attributable.
Under the matched arms initiation improves rather than degrades, and every
protein-level liability moves in the direction the objective pushes it.

## *E. coli* BL21(DE3)

| feature | vendor | codon | joint | codon buys | residues buy |
|---|---|---|---|---|---|
| initiation window dG | -2.183 | -1.958 | -0.783 | +0.225 | +1.175 |
| transcript MFE | -31.718 | -30.867 | -29.175 | +0.852 | +1.692 |
| CAI | 1.000 | 0.971 | 0.967 | -0.029 | -0.004 |
| codon-pair score | -0.014 | 0.096 | 0.130 | +0.111 | +0.033 |
| max GC, 20 bp | 0.742 | 0.708 | 0.675 | -0.033 | -0.033 |
| longest repeat, nt | 9.000 | 8.667 | 7.667 | -0.333 | -1.000 |
| synthesis violations | 0.500 | 0.417 | 0.083 | -0.083 | -0.333 |
| low complexity | 0.058 | 0.058 | 0.000 | **+0.000** | -0.058 |
| protein repeat | 1.667 | 1.667 | 0.000 | **+0.000** | -1.667 |
| exposed hydrophobic | 0.645 | 0.645 | 0.445 | **+0.000** | -0.200 |
| exposed hydrophobic run | 0.750 | 0.750 | 0.583 | **+0.000** | -0.167 |
| degron load | 0.191 | 0.191 | 0.000 | **+0.000** | -0.191 |

Passing every synthesis constraint: vendor 58%, codon 67%, joint 92%.

## HEK293

| feature | vendor | codon | joint | codon buys | residues buy |
|---|---|---|---|---|---|
| initiation window dG | -4.633 | -5.592 | -1.725 | -0.958 | +3.867 |
| transcript MFE | -44.558 | -48.633 | -47.125 | -4.075 | +1.508 |
| CAI | 1.000 | 0.963 | 0.962 | -0.037 | -0.002 |
| codon-pair score | 0.005 | 0.213 | 0.251 | +0.208 | +0.038 |
| longest repeat, nt | 9.083 | 9.000 | 8.167 | -0.083 | -0.833 |
| synthesis violations | 0.500 | 0.500 | 0.167 | +0.000 | -0.333 |
| low complexity | 0.058 | 0.058 | 0.000 | **+0.000** | -0.058 |
| protein repeat | 1.667 | 1.667 | 0.000 | **+0.000** | -1.667 |
| exposed hydrophobic | 0.645 | 0.645 | 0.421 | **+0.000** | -0.224 |
| exposed hydrophobic run | 0.750 | 0.750 | 0.500 | **+0.000** | -0.250 |
| degron load | 0.191 | 0.191 | 0.000 | **+0.000** | -0.191 |

Passing every synthesis constraint: vendor 58%, codon 58%, joint 83%.

## Reading the split

The codon layer's column is exactly zero on all five protein-level liabilities in
both hosts, for the same reason as in the exchange-rate experiment: they are not
functions of the codons. What it does buy is gene-level, and it buys it by giving
up CAI. In *E. coli* it trades 0.029 of CAI for +0.111 of codon-pair score, 0.225
kcal/mol of initiation window and a third of a synthesis violation. That is the
whole case for a multi-objective codon DP over a max-CAI vendor table, and it is
a modest case.

The residue layer clears the protein-level terms outright: degron load,
low-complexity fraction and repeat all go to zero in both hosts, exposed
hydrophobic falls by about a third. It also improves initiation further, by more
than the codon layer managed, which is the result the mismatched arms had hidden.

## What the residues cost, and the circularity in that number

| residues moved | objective | Tier-1 | initiation | liability burden |
|---|---|---|---|---|
| 0 | +0.0000 | +0.000 | +0.000 | +0.000 |
| 1 | +1.2306 | +2.035 | -0.150 | -3.658 |
| 2 | +1.7515 | +3.767 | +0.083 | -4.557 |
| 3 | +2.1620 | +5.900 | -0.217 | -4.989 |
| 5 | +2.7963 | +8.977 | +0.408 | -5.361 |
| 8 | +3.4426 | +12.445 | +0.517 | -5.656 |

*E. coli*, mean over the 12 designs, as change from the codon arm.

The Tier-1 column rises rather than falls, so on this baseline the residue moves
are not paying for anything: they gain fold compatibility and clear liabilities
at the same time. That is because the starting protein is the design's own
sequence, which sits well below the ProteinMPNN optimum inside its own shell,
not because the amino-acid axis is free. Where the exchange-rate experiment
starts from the Tier-1 optimum and measures a genuine price of 0.2 nats, this one
starts from a point with slack in every direction.

The Tier-1 gain is also measured on the same marginals the objective maximises,
so on its own it is a statement about the model rather than about the protein.
The independent check on that is the measured-stability arm, not this table.

The practical claim lives in the first rows. The first residue moved is worth
+1.23 on the combined objective and removes most of the liability burden; three
residues reach +2.16. Running to a local optimum takes 21 residues in *E. coli*
and 22 in HEK293, which leaves identity to the design at 0.52 and 0.49. At that
point it is a redesign rather than a codon optimisation of an existing design,
and the small-budget rows are the defensible operating point.

## Files

- `analysis/scripts/case_study.py`
- `janus-data/processed/case_study.json`, 72 records
- `src/janus/lattice.py`, the `anchor` argument

## The trajectory, decomposed by term

The combined-objective number is contaminated for one of its components. Tier-1
score rises along this trajectory because the starting protein is the design's own
sequence, which sits below the ProteinMPNN optimum inside its own shell, and it is
measured on the same marginals the objective maximises. A gain there is the model
agreeing with itself.

Initiation energy and liability burden are not contaminated that way. Folding
energy comes from ViennaRNA and the liabilities from motif and composition
definitions, neither of which the search consults through ProteinMPNN. The
`independent` column below is the objective with the model-internal term removed
in matching units, and it is the column to read.

*E. coli* BL21(DE3), 12 designs, means over backbones:

| residues moved | independent | initiation dG | burden | Tier-1, model-internal | objective |
|---|---|---|---|---|---|
| 0 | +0.0000 | +0.000 | +0.000 | +0.000 | +0.0000 |
| **1** | **+0.9006** | -0.150 | **-3.658** | +2.035 | +1.2306 |
| 2 | +1.1470 | +0.083 | -4.557 | +3.767 | +1.7515 |
| **3** | **+1.2201** | -0.217 | **-4.989** | +5.900 | +2.1620 |
| 5 | +1.3741 | +0.408 | -5.361 | +8.977 | +2.7963 |
| 8 | +1.4750 | +0.517 | -5.656 | +12.445 | +3.4426 |

HEK293, same designs:

| residues moved | independent | initiation dG | burden | Tier-1, model-internal | objective |
|---|---|---|---|---|---|
| **1** | **+0.9859** | +0.417 | **-3.803** | +1.851 | +1.2986 |
| 2 | +1.2308 | +1.033 | -4.589 | +3.787 | +1.8563 |
| **3** | **+1.3177** | +1.033 | -4.937 | +5.869 | +2.2703 |
| 5 | +1.6054 | +2.542 | -5.490 | +8.504 | +2.9812 |

**One residue substitution removes 3.7 units of liability burden**, and three
remove 5.0. Those are the independently measured quantities and they carry most of
the trajectory's value: the independent column reaches 73% of its eventual level
by one residue and 90% by three.

**The gain is almost entirely liability, not initiation.** In *E. coli* the
initiation window barely moves over the first three residues and moves slightly
the wrong way, -0.15 and -0.22 kcal/mol. The residue layer is not being used to
open the 5' window here; the codon layer already did that, +0.225 kcal/mol on its
own. In HEK293 initiation does improve, +0.42 and +1.03, so the split between the
two terms is host-dependent.

The first three rows are the operating regime and are what the paper quotes.
Running to a local optimum takes 21.33 residues and leaves identity at 0.519,
which is a redesign rather than an optimisation, and the +4.6755 combined-objective
figure for the full run should not be read as what the method offers a user with a
design in hand.

Two further things the full-run row hides. Constraint passing improves from 58% to
92%, which is a real and independently checked gain. And the degron, low-complexity
and repeat burdens go to exactly zero, which the codon arm cannot do at all: every
one of those columns reads +0.000 for what the codon layer buys.
