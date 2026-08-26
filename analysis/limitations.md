# Limitations

What was measured, what it bounds, and where the measurement lives.

## The objective is a surrogate, and the gap is large

The parser is exact with respect to ProteinMPNN's unconditional single-pass
marginals. It is not exact with respect to the conditional autoregressive
posterior, which is what ProteinMPNN actually models. The conditional
log-probability at one position depends on residues chosen elsewhere, so it does
not decompose over lattice nodes and no dynamic program can optimise it.

The gap was measured, not assumed. Across the candidate pool the conditional
score scatters by **2.42 nats** against a **0.30-nat** band of surrogate values.
So the surrogate orders candidates the conditional model would order differently,
and a design that is optimal under the lattice is not the conditional optimum.
Reporting the collapse of the correlation as an r value would misdescribe it,
because the surrogate's range is restricted by construction. The two spreads are
the statement that survives that restriction. See `analysis/surrogate_gap.md`.

## Generate and rank has a measured breaking point

Terms that do not decompose over lattice nodes are handled by ranking what the
parser proposes. That is sound only while those terms are secondary to the ones
inside the parser. Once a Tier-2 term is strong enough to move the optimum out of
the Tier-1 neighbourhood, drawing broadly from the shell finds better candidates
than a deep k-best prefix does.

The crossover is at a spread-normalised weight of about 1.0, which is roughly
twice the top of the operating region these terms are used at, 0.125 to 0.5.
Inside that region the parser prefix wins on 60 to 100 percent of backbones at
every budget from one evaluation to a thousand. Outside it, the architecture is
the wrong shape and the alternatives are broad sampling with exact codon
optimisation, or annealing on the combined objective. See
`analysis/rescore_recall.md`.

An earlier version of this table reported the crossover at a weight of 1.0 on an
unnormalised scale, where the two terms were combined in nats and kcal/mol
directly. That value was 3.63 in normalised units. The bound is real; it sits
further out than the earlier table implied.

## The folding term is in Tier 2, not in the lattice

mRNA folding energy is not a sum over positions, because base pairs form at
arbitrary range. It is therefore evaluated on finished transcripts and used to
rank, not to search. The claim of exactness never covers it. A method that put
folding into the lattice would be solving a different and harder problem, which
is what LinearDesign does with the protein held fixed.

## There is no nucleotide-level lattice

Synthesis constraints, restriction sites, homopolymer runs, local GC windows and
cryptic regulatory elements span more than adjacent codons. They are applied as a
filter over the k-best list rather than as terms in the parse, so a design that
violates one is reported and skipped rather than avoided by construction. In
practice this is not costly, and the amino-acid freedom more than halves the
violation rate as a side effect, from 0.83 to 0.37 per design. But the guarantee
does not extend to them.

## Substitutions are scored independently and they are not additive

The method proposes several substitutions at once and scores them independently.
Measured on 230992 MegaScale double mutants, deviation from additivity is
positive and, in the central band of the assay where measurement is best, about
**+0.35 kcal/mol** on both natural and designed parents, with a per-pair spread
near 1.9 kcal/mol.

An earlier reading held that epistasis was milder on designed parents inside the
shell, +0.186 against +0.850. That gradient does not survive restriction to the
well-measured band: the gap closes monotonically from +0.527 to -0.036 and is
gone. **The pattern argument is withdrawn**, and the n=220 subset that carried
its strongest form was never significant on its own (p = 0.48).

Nor is that subset evidence of anything else. Its 95% interval is [-0.272,
+0.243], eleven times wider than the designed estimate's and containing it
outright, so it is consistent with the larger bias rather than contradicting it.
An earlier note described it as unbiased to plus or minus 0.26 kcal/mol, which
presented a wide interval around a small point estimate as though it established
absence. That wording is withdrawn.

The reconciliation that would have helped is not available either. Shell
substitutions would be more additive if they were milder, and deviation does grow
with severity, from +0.181 below 0.5 kcal/mol to +0.661 above 2. But the shell
does not select milder substitutions: median severity is 1.336 inside against
1.470 outside, a difference of -0.134 at p = 0.168.

**And deviation does not decay with sequence separation**, +0.326 for pairs five
residues apart against +0.347 for pairs forty apart, so substitutions scattered
along a short chain cannot be assumed independent.

Per-pair deviation on designed parents has mean +0.166 and standard deviation
1.737 kcal/mol. For the five or six substitutions the method makes at the
operating point, treating pairs as independent gives a systematic offset of +1.66
to +2.49 kcal/mol and a 95% half-width of 10.77 to 13.19. Those are two different
quantities and both are quoted: the offset uses the mean, since the expectation of
a sum is the sum of expectations, and the interval is a spread rather than a
bound.

An earlier version of this section called the one-sigma figure of 5.49 kcal/mol a
bound. It is not; about a third of designs fall outside it. That wording is
withdrawn.

The claim therefore splits in two. **H4 holds where it is measured**: single
substitutions inside the shell are non-inferior to random draws from the same
shell on measured stability, n = 58456, and that stands at full strength.
**Per-design aggregate stability prediction is imprecise**, by the offset and
spread above, on domains whose whole folding free energy is a few kcal/mol. The
direction is the safer one, since positive deviation means a design scored
additively is underestimated. JANUS can say its substitutions are individually
non-inferior; it cannot promise a stability-neutral design. See
`analysis/epistasis.md`.

## The atlas describes designs that survived

The 447 backbones come from published, experimentally characterised sets, so each
survived whatever filtering its authors applied. Splitting a separate benchmark
of 614 designs on experimental outcome shows selection does remove liability: 15
of 22 determined features are higher in the failures, including exposed
hydrophobic area at +0.140 and weighted degron load at +0.120. Both are atlas
features, so those comparisons are survivor-biased.

The two findings the paper rests on are not. Low-complexity fraction, at -0.033
[-0.11, +0.05], and protein repeat, at +0.077 [-0.01, +0.16], establish nothing:
designs that failed carry them at the same rate as designs that succeeded.
Selection does not filter on the liabilities the amino-acid axis clears. See
`analysis/survivors.md`.

This is one benchmark's notion of success from one group's pipeline. It measures
the filtering that produced those 614, not filtering in general, and it is not a
substitute for running the atlas on unselected generator output.

## The natural controls are two different things

The like-for-like natural arm is drawn from the same MegaScale release and is
mostly excised domains: of 104 distinct entries only 15 are whole UniProt
proteins, and the median covers 0.11 of its parent. A motif outside the cut is
unseen, not absent, and that is enough to flip the sign of a comparison.

A separate arm of whole natural proteins was assembled for the comparisons that
need it, matched from a pool of 2235 reviewed UniProt entries. The two are
reported separately and never pooled. The whole-protein arm carries its own
limitation: designs and short natural proteins differ in fold as well as in
origin, with mean accessibility differing by -0.812 and strand fraction by
+0.460, and matching on length and composition does not fix that. See
`analysis/atlas.md`.

Both natural arms are matched with replacement, so 447 records resolve to 105
distinct proteins in the whole-protein arm and 104 in the excised arm. The
bootstrap intervals resample records and so do not price that repetition; they
are narrower than an interval over distinct proteins would be.

## The degron panel is null and one family carries what is left

Against whole natural proteins the pooled accessibility-weighted degron effect is
-0.107 [-0.30, +0.11] and does not exclude zero. Per topology only EEHEE clears
correction, at -0.603 [-0.92, -0.12] and q = 0.045, itself marginal against the
0.05 threshold. No family is enriched. Record resampling gave a narrower pooled
interval that excluded zero and put the hallucination family and its trRosetta
subset above threshold; both fall back inside zero once the natural side is
resampled by protein, and no degron enrichment is claimed.

## The H5 benchmark is contaminated asymmetrically

ESMFold's training data includes PDB structures that overlap the benchmark's
designs, and the liability score's motif definitions do not. Contamination
therefore favours the baseline and not the challenger, which makes the reported
improvement conservative. It cannot be stratified away, because the overlap is
not recorded per design.

Separately, accessibility weighting does not help here. It appeared worth +0.106
in AUC and is worth **+0.009** once the comparison is restricted to designs that
all have a usable ESMFold model. The apparent gain was a restatement of pLDDT,
which is the comparator.

The same restriction reverses which score wins, and that cuts the other way. On
all 614 pLDDT beats the sequence-only liability score, 0.731 against 0.672,
p = 0.049; on the 411 with a usable model the liability score wins, 0.688 against
0.591, p = 0.014. The 614 comparison lets "has no confident model" act as a covert
predictor on the pLDDT side, so it cannot settle which score is better, and the
411 cohort is primary for that question. The 614 figure remains the right one for
deployment, where designs without confident models do arrive and do have to be
filtered. See `analysis/predictive.md`.

## The eukaryotic 5' model is weaker than the prokaryotic one

The *E. coli* initiation term rests on Kudla 2009 and Goodman 2013, which both
find folding stability over the initiation window to be the dominant
coding-sequence determinant of expression. The HEK293 term rests on cap-scanning
arguments with no comparable quantitative anchor. The expression-unit conversion
is applied only to *E. coli*, because the fitted coefficient it uses was fitted
there.

## The expression conversion is an extrapolation

The predicted fold-change uses the apparent Boltzmann factor from the Ribosome
Binding Site Calculator, 0.45 plus or minus 0.05 mol/kcal. That coefficient
relates initiation rate to a composite free energy of which the initiation-window
folding term is one part, so applying it to that term alone assumes a one-for-one
contribution that has not been measured. And a fold-change in initiation rate is
not a fold-change in yield. The number belongs in the discussion, not on a figure
axis. See `analysis/expression.md`.

## MegaScale's readout is proteolysis-derived

Stability comes from K50 protease resistance on small domains, converted to a
free energy, not from calorimetry. It is a large and internally consistent
dataset and it is not a thermodynamic gold standard. The deposited values are
extrapolated rather than clipped, so their tails are imprecise rather than
censored, which is why the additivity analysis restricts to a central band.

## The Tier-1 gain in the case study is measured on its own model

Along the case study's trajectory the Tier-1 score rises rather than falls,
because the starting protein is the design's own sequence, which sits below the
ProteinMPNN optimum inside its own shell. That gain is measured on the same
marginals the objective maximises, so by itself it is a statement about the model
rather than about the protein. The independent check is the measured-stability
arm.

## The searches outside the parser are heuristic

Everything in Tier 2 is hill climbing over single-residue substitutions from the
Tier-1 optimum. The traces are greedy, so they are attainable exchanges and lower
bounds on what the shell holds, not Pareto frontiers. Nothing in the central
comparison depends on their optimality, since the arm they are measured against
achieves exactly zero, but the reported costs are upper bounds on what a better
search would pay.

The coordinate-descent baseline is a different acceptance order over the same
neighbourhood, not the conditional-ProteinMPNN alternation that LinearDesign2
used. That comparison is not implemented.

## No wet lab

Nothing here was expressed, purified or measured. Every claim about expression is
a prediction from a folding model and a literature coefficient, and every claim
about stability is a claim about someone else's measurements on someone else's
proteins.
