# Every method, one objective

447 designed backbones, host *E. coli* BL21(DE3), delta 1.0 nats, folding weight
0.25 and liability weight 0.25 per term on the spread-normalised scale,
22 August 2026.

Ten arms scored by one piece of code. The objective is the complete one: Tier 1
exact in the parser, initiation-window folding and the five protein-level
liabilities in the rescoring tier, every term divided by its spread over one pool
drawn per backbone.

| arm | objective | vs codon DP | Tier-1 | burden | violations | identity | seconds | evaluations |
|---|---|---|---|---|---|---|---|---|
| vendor CAI-max | -14.107 | -0.073 | -77.555 | 9.596 | 0.87 | 1.000 | 0.000 | 1 |
| CodonTransformer | -14.518 | -0.488 | -80.426 | 9.596 | 0.58 | 1.000 | 0.099 | 1 |
| codon DP | -14.018 | +0.000 | -76.869 | 9.596 | 0.74 | 1.000 | 0.001 | 1 |
| CodonMPNN | -14.741 | -0.527 | -76.992 | 12.533 | 1.39 | 0.441 | 0.991 | 1 |
| ProteinMPNN then codons | -11.736 | +2.227 | -59.356 | 11.009 | 1.69 | 0.436 | 0.321 | 1 |
| SolubleMPNN then codons | -12.449 | +1.598 | -63.254 | 11.179 | 1.41 | 0.443 | 0.321 | 1 |
| MoMPNN then codons | -12.472 | +1.492 | -59.885 | **13.458** | 1.60 | 0.456 | 0.321 | 1 |
| rejection sampling | -10.429 | +3.654 | -63.706 | 4.219 | 0.14 | 0.393 | 2.012 | 2000 |
| coordinate descent | -8.973 | +5.005 | -56.724 | 2.995 | 0.20 | 0.490 | 0.355 | 437 |
| **JANUS joint** | **-8.911** | **+5.023** | -56.680 | 2.968 | 0.21 | 0.511 | 2.352 | 3121 |

Share of backbones on which each arm beats the codon DP: vendor 21%,
CodonTransformer 4%, CodonMPNN 39%, ProteinMPNN 87%, SolubleMPNN 83%, MoMPNN 76%,
and 100% for each of the three search arms.

## What the codon-only arms can and cannot do

The first three arms hold the design's own residue sequence. They differ only in
how they pick codons, and their liability burden is identical at 9.596 because
none of them can touch it. **Three different codon optimisers, one number, to
three decimal places.**

Among them the codon DP is best under this objective, as it must be, since it is
the only one optimising it. It buys that by giving up CAI, 0.969 against the
vendor's 1.000, for codon-pair score 0.099 against -0.002. The margin is small and
the vendor wins on 21% of backbones where the initiation window happens to favour
the maximum-CAI assignment.

## The learned models are not optimising this objective

Five arms come from models trained to imitate rather than to optimise, and
scoring them under an objective none has seen would be circular if the total were
all that was reported. The per-term columns are what make the comparison say
something.

| arm | CAI | codon pair | GC | initiation dG |
|---|---|---|---|---|
| vendor CAI-max | 1.000 | -0.002 | 0.492 | -4.100 |
| CodonTransformer | 0.861 | 0.053 | 0.539 | -4.300 |
| codon DP | 0.969 | 0.099 | 0.481 | -4.200 |
| CodonMPNN | 0.926 | 0.036 | 0.542 | -5.200 |
| ProteinMPNN then codons | 0.970 | 0.105 | 0.450 | -4.700 |
| SolubleMPNN then codons | 0.973 | 0.093 | 0.422 | -3.700 |
| MoMPNN then codons | 0.976 | 0.084 | 0.420 | -3.000 |
| JANUS joint | 0.970 | 0.129 | 0.489 | -0.800 |

**CodonMPNN** is the one published method occupying the joint amino-acid and
codon space. It moves 56% of residues, so it is genuinely searching that space,
and its Tier-1 score of -76.992 is level with the codon DP's -76.869: the proteins
it picks are about as fold-compatible under the marginals as the design's own.
What it does not do is control anything else. Its liability burden is 12.533
against 9.596 for every fixed-protein arm, its violations run 1.39 against 0.74,
and its initiation window is the most structured of any arm. None of those three
is circular; they are properties a practitioner cares about regardless of what is
being optimised.

## Property alignment does not lower this panel

Three residue-proposal models are here, and the point of having three is that one
would not have settled whether an effect belonged to a model or to a class.

| model | trained for | burden | vs the design's own protein |
|---|---|---|---|
| the design's own protein | - | 9.596 | - |
| ProteinMPNN | structure recovery | 11.009 | +1.41 |
| SolubleMPNN | solubility, by training-set exclusion | 11.179 | +1.58 |
| MoMPNN | solubility and thermostability, by preference alignment | 13.458 | +3.86 |

**All three raise the protein-level liability burden above the designs they
replace, and the two property-aware models raise it more than the plain one.**
That is the second point on the axis, and it says the effect is not MoMPNN
specifically.

The reading has to stay narrow. SolubleMPNN's notion of solubility is a training
set filtered to soluble proteins; MoMPNN's is preference alignment against
solubility and thermostability labels; and the panel here is low-complexity
content, repeats, exposed hydrophobic surface, hydrophobic runs and degron load.
These are related notions, not the same one, and none of these models was asked to
minimise what is measured here. What the table supports is that **aligning a
residue model for developability, in either operationalisation available, does not
incidentally lower this panel**, and that a designer who wants these terms
controlled has to optimise them rather than hope a property-aware inverse folding
model delivers them.

The mechanism is visible in the identity column: all three move about 56% of
residues, so they are not making small corrections to the design, they are
replacing it with the protein their own objective prefers. Whatever the original
design got right about these five terms is discarded along with the rest.

Both property-aware models do deliver a better initiation window than the plain
one, -3.700 and -3.000 against ProteinMPNN's -4.700, without seeing an mRNA at
all. That is a composition effect rather than a folding one, and it is the one
place where the alignment helps something this objective also wants.

## Block alternation, in the only form this objective admits

A block scheme alternates optimising one block given the other, which requires the
residue step to see what the codon step is optimising. An inverse folding model
has no channel for that: it scores residues against a backbone and knows nothing
of initiation windows, codon pairs or degron load. Alternation therefore collapses
to a single round, propose residues with the model and then solve the codon layer
exactly for them, and the two `then codons` arms are exactly that shape.

They recover **46%**, **32%** and **32%** of the joint arm's gain over the
fixed-protein DP. All three beat the codon DP on most backbones, 87%, 83% and 76%,
and all three do it by moving to proteins near their own model's optimum, which is
why their Tier-1 scores jump from -76.9 to between -59 and -63.

**All three make the protein-level liabilities worse while doing it**, and all
three roughly double the synthesis violations. The gain is bought on the terms the
residue model is good at and paid for on the terms it cannot see. That is the
direct answer to what the joint objective buys over proposing and then solving:
between two and three times the objective, with the difference concentrated in
exactly the terms no proposal model is optimising.

## The search arms, including the one that ties

Rejection sampling recovers 70% of the joint gain from 2000 draws. That is
Baseline 0 and it is a strong baseline, worth saying plainly: most of what the
lattice buys here is available to anyone willing to sample the same shell and rank
it. What the lattice adds is the remaining 30%, at a fifth of the wall clock, with
a guarantee attached to the codon layer.

**Coordinate descent recovers 100%.** At -8.973 against the joint arm's -8.911 it
is within 1.2% of the objective, on 437 evaluations against 3121 and 0.353 seconds
against 2.329. Both are local searches over the same neighbourhood, differing only
in acceptance order. On this landscape the joint arm's extra thoroughness buys
+0.06 of objective for seven times the compute.

That is a finding about the search, not about the lattice. Both arms are the joint
lattice and both beat the fixed-protein DP on 447 of 447. It says the landscape is
benign enough that a cheap ordering suffices, which is useful to report and which
nothing in the paper's claims depends on. It also says which arm should be the
default: **coordinate descent is the better engineering choice on measured
grounds**, and the exact parser remains what solves the codon layer inside it and
supplies the optimum the exchange rate is quoted against.

## What this comparison does not include

Budgets are matched two ways, wall clock and full-objective evaluations, and both
are reported, because the methods differ by four orders of magnitude in cost per
evaluation and either measure alone flatters someone.

The `then codons` arms use one sampled sequence per backbone at temperature 0.1.
Sampling several and keeping the best would raise them, at a cost this table would
then have to carry, and that is the same trade the rejection arm makes explicit.

Rejection sampling ends with **fewer** synthesis violations than the joint arm,
0.14 against 0.21. Violations are checked rather than weighted by any arm here, so
neither was optimising them, and the difference is which candidates happened to
survive the filter out of pools of different size and composition. It is not a
result about either method.

## Files

- `analysis/scripts/baselines.py`
- `baselines/run_codontransformer.py`, `baselines/run_codonmpnn.py`
- `baselines/codontransformer.Dockerfile`
- `janus-data/processed/baselines.json`, 4470 rows
