# Liability atlas

447 de novo designed backbones against four control sets. Feature panel,
thresholds and statistics fixed beforehand in `analysis_plan.md`. The first three
are below; the whole-protein arm was added later and has its own section.

Controls: scrambled, a within-design permutation; random, draws from the
aggregate residue composition of the design corpus; natural, MegaScale natural
domains matched to each design on length within 10 percent and then, inside that
band, by nearest neighbour in amino-acid composition space. Both strata of D8
are used. Scrambled and random keep the design's own accessibility profile, so
weighted features compare the same structure carrying a different sequence.

Gene-layer features are computed on coding sequences generated identically for
every group by taking the highest relative-adaptiveness codon per residue, so
differences reflect the residue sequence and not the codon policy.

## The degron reading moved twice before settling

Worth recording in order, because two of the three readings were wrong and both
would have been reportable.

| controls | max weighted degron delta vs natural | reading |
|---|---|---|
| 200 designs, length-matched | -0.35 | large effect |
| 447 designs, length-matched | +0.02 | below the floor |
| 447 designs, length and composition matched | **-0.36** | large effect |

The first failed because the 200 were taken in name order and dominated by one
topology family. The second failed because matching on length alone leaves
amino-acid composition free, and composition drives both motif content and, under
a fixed codon policy, GC content. The third is what the analysis plan specified.

The largest accessibility-weighted degron effect is 0.356, well above the 0.15
floor, so the degron panel is retained.

## Degron burden

| feature | unweighted | weighted | 95% interval, weighted | q |
|---|---|---|---|---|
| N-degron | -0.305 | **-0.356** | [-0.42, -0.29] | 1.7e-20 |
| all degrons | -0.238 | -0.293 | [-0.36, -0.22] | 5.5e-14 |
| internal degron | -0.071 | -0.094 | [-0.15, -0.04] | 0.0028 |
| C-degron | -0.015 | -0.015 | [-0.05, +0.02] | 0.56 |

On the pooled natural arm, designs carry fewer degron motifs than matched
natural domains, and accessibility weighting enlarges the effect rather than
shrinking it, from -0.305 to -0.356 for N-degrons.

**That effect does not survive the excision check and is withdrawn.** See below.

C-degrons show nothing at all, in any comparison. At 43 residues fewer than one
chain in ten carries a C-end match under any ELM class.

## What separates designs, with composition controlled

| feature | vs natural | vs random | position on the random-to-natural axis |
|---|---|---|---|
| exposed hydrophobic area, density | -0.323 | -0.902 | 1.50 |
| N-degron, weighted density | -0.356 | -0.212 | - |
| low-complexity fraction | +0.275 | +0.235 | - |
| longest repeat | +0.232 | +0.285 | 1.00 |
| initiation window dG | +0.157 | +0.100 | - |
| codon-pair score | +0.152 | +0.049 | - |

The GC and transcript-structure differences reported before composition matching
have gone. They were composition artefacts, which is exactly the failure D8's
second stratum exists to prevent.

Exposed hydrophobic area is the largest effect in the panel against random,
-0.902, and its normalised position is 1.50, past natural. On surface
hydrophobicity these designs are cleaner than the natural domains they were
matched to. They come from stability-selected sets, so this reports what survived
selection rather than what a generator emits.

Low-complexity content and repeat length are higher in designs than in every
control including composition-matched random, so those are not composition
effects. They are the compositional liabilities selection removes and a design
pipeline does not.

## The natural arm is mostly excised domains, and that explains the degron result

MegaScale's natural set is drawn from the PDB at 40 to 80 residues, which is the
size range where excised domains live. A domain cut out of a larger chain has
termini that were never solvent exposed in the parent, so the excision protocol
can manufacture an accessibility difference on its own.

Classifying each natural entry by how much of its UniProt reference the deposited
entity covers: of 104 distinct entries, 88 resolved, and only **15 are whole
proteins** at 90 percent coverage or better. **73 are fragments, with a median
coverage of 0.11**, so the typical natural control covers about a ninth of its
parent protein.

Splitting the comparison, with cluster bootstrap intervals over source entries:

| feature | all natural | whole proteins only | excised only |
|---|---|---|---|
| N-terminal accessibility | -0.397 | **-0.487** [-0.75, -0.23] | -0.415 [-0.69, -0.10] |
| C-terminal accessibility | -0.373 | **-0.409** [-0.63, -0.07] | -0.428 [-0.61, -0.21] |
| mean accessibility | -0.429 | **-0.403** [-0.71, +0.27] | -0.481 [-0.62, -0.29] |
| N-degron, weighted density | -0.356 | **+0.040** [-0.46, +0.44] | -0.331 [-0.52, -0.08] |
| all degrons, weighted density | -0.293 | **+0.014** [-0.41, +0.41] | -0.299 [-0.52, +0.02] |

Terminal burial survives on the whole-protein subset, with intervals excluding
zero and a point estimate if anything larger: designs do bury their N and C
termini more than whole natural proteins of matched length and composition. Mean
accessibility no longer separates once the interval is taken over 15 entries
instead of 27 records. The degron depletion does not survive: its interval
straddles zero on whole proteins and stays negative on the fragments.

The whole-protein arm is 27 records drawn from 15 distinct entries, so the
interval is wide and cannot exclude a moderate effect in either direction. What it
does exclude is the pooled effect.

## The primary statistic, with intervals

Section 8.3 named the normalised position on the random-to-natural axis as H1's
primary statistic. It is only well determined where that axis is wide, and for
most of the panel it is not.

| feature | axis width | position | 95% interval |
|---|---|---|---|
| exposed hydrophobic area, density | 4.39 | **1.50** | [+1.30, +1.56] |
| longest exposed hydrophobic run | 1.00 | 1.00 | [+1.00, +1.00] |
| longest repeat | 1.00 | 1.00 | [+0.00, +1.00] |
| transcript MFE | 4.60 | -0.07 | [-0.68, +0.31] |
| all degrons, weighted density | 0.29 | -1.25 | [-2.73, -0.59] |
| GC | 0.005 | -0.62 | [-13.8, +13.4] |
| initiation window dG | 0.20 | -4.00 | [-11.0, +11.0] |

Where the random and natural medians sit close together the denominator is near
zero and the ratio is meaningless, which the intervals show directly: GC has an
axis 0.005 wide and an interval spanning plus and minus fourteen. Those features
are reported as effect sizes only, and that decision is now justified by the
intervals rather than by a suppression rule.

Exposed hydrophobic area is the one headline feature where the statistic is
sharp: designs sit at **1.50 [1.30, 1.56]**, past the natural end of the axis.
The degron positions in this table are computed against the pooled natural arm
and inherit the excision problem described above; they are not claimable.


## A whole-protein natural arm, built from scratch

The 15 whole entries above were too thin to settle H1, so a second natural arm
was assembled: reviewed Swiss-Prot entries of 26 to 74 residues with
protein-level evidence, not fragments, carrying no signal peptide, propeptide,
transit peptide or peptide feature, and no transmembrane or intramembrane
annotation. The cleavage filters are what make them whole, since an entry with a
cleaved region has a precursor sequence and a different mature protein. The
membrane filter keeps the band from filling with recently annotated microproteins
whose helices are buried in a bilayer the monomeric model does not contain.

**2235 candidates**, against 415 in the excised arm, structures from AlphaFold DB
so accessibility is computed the same way as for the designed backbones. Matched
to each design on length within 10 percent, then nearest neighbour in composition
space, exactly as the excised arm is matched.

Matching is with replacement, so the arm is 447 records over 105 distinct
proteins, against 104 distinct entries in the excised arm. Every comparison below
is 447 against 447 with repeats on the natural side; the effective sample size
there is the distinct count, not the record count.

### Degrons

| feature | vs whole | 95% CI | vs excised | q |
|---|---|---|---|---|
| all degrons, weighted | **-0.103** | [-0.18, -0.03] | -0.275 | 0.0094 |
| all degrons, weighted density | -0.107 | [-0.18, -0.03] | -0.293 | 0.0069 |
| all degrons, raw | +0.050 | [-0.01, +0.12] | -0.176 | 0.18 |
| N-degron, raw | +0.112 | [+0.04, +0.17] | -0.208 | 0.0013 |
| N-degron, weighted | +0.012 | [-0.06, +0.08] | -0.321 | 0.74 |
| C-degron, weighted | -0.044 | [-0.08, -0.00] | -0.016 | 0.042 |

Designs are not degron-enriched against whole natural proteins. On raw counts
they carry slightly more N-degrons, and that difference vanishes once
accessibility is applied: the N-degrons designs carry are buried. On
accessibility-weighted total load they are marginally *lower*, -0.103 with an
interval that excludes zero but only just.

The pooled degron effect against the excised arm, -0.275, is roughly two and a
half times the whole-protein estimate. Part of that is excision and part is the
matching: with 2235 candidates instead of 415 the nearest composition neighbour
is a much closer match, so composition-driven differences shrink.

### What separates designs from whole natural proteins

| feature | delta | design | whole | q |
|---|---|---|---|---|
| exposed hydrophobic area, density | **-0.877** | 3.25 | 10.94 | 1.5e-112 |
| mean accessibility | -0.812 | 0.386 | 0.516 | 7.5e-97 |
| N-terminal accessibility | -0.628 | 0.669 | 0.977 | 2.6e-58 |
| C-terminal accessibility | -0.523 | 0.761 | 0.987 | 9.8e-41 |
| strand fraction | +0.460 | 0.234 | 0.000 | 3.2e-34 |
| protease site density, raw | +0.346 | 6.06 | 2.44 | 8.8e-19 |
| longest repeat | +0.340 | 9.00 | 8.00 | 1.3e-18 |
| low-complexity fraction | +0.339 | 0.000 | 0.000 | 4.2e-30 |

Designs are worse on low-complexity content, repeat structure and protease-site
density. They are cleaner on exposed hydrophobic surface, by a very large margin,
and they bury their termini more.

**The two populations are not structurally matched, and the accessibility
comparisons inherit that.** Mean accessibility differs by -0.812 and strand
fraction by +0.460: natural proteins in this size band are almost entirely
helical or extended, while the design set includes beta topologies, and the
designs are more compact. Matching on length and composition does not fix that.
The scrambled and random arms do control for structure, since they keep the
design's own backbone, and against those the exposed-hydrophobic result holds.
Against the natural arms it should be read as a difference between two
populations that differ in fold as well as in origin.

## Per-topology effect sizes

One topology family carries the pooled number, so the pooled number is not the
primary reporting. Designs against the whole-protein arm, bootstrap intervals, families
below 8 designs not reported.

### Degron load, weighted density: the pooled number is a mix of signs

| topology | n | delta | 95% CI | verdict |
|---|---|---|---|---|
| EEHEE | 65 | -0.603 | [-0.76, -0.44] | determined |
| HEEH | 69 | -0.208 | [-0.39, -0.02] | determined |
| hallucination | 148 | **+0.148** | [+0.01, +0.28] | determined |
| HHH | 64 | -0.176 | [-0.37, +0.02] | establishes nothing |
| other | 47 | -0.183 | [-0.41, +0.04] | establishes nothing |
| EHEE | 54 | +0.003 | [-0.21, +0.22] | establishes nothing |

The pooled -0.103 averages over families whose intervals do not overlap.
Hallucinated designs are degron **enriched** relative to whole natural proteins;
EEHEE designs are strongly depleted. Reporting the pooled estimate alone would
describe neither. This is the clearest case in the analysis for per-topology
reporting, and the degron claim is now made per family or not at all.

### Low-complexity fraction: robust across families

| topology | n | delta | 95% CI | verdict |
|---|---|---|---|---|
| hallucination | 148 | +0.472 | [+0.37, +0.56] | determined |
| HEEH | 69 | +0.447 | [+0.32, +0.57] | determined |
| EHEE | 54 | +0.347 | [+0.22, +0.48] | determined |
| other | 47 | +0.324 | [+0.17, +0.49] | determined |
| EEHEE | 65 | +0.170 | [+0.03, +0.31] | determined |
| HHH | 64 | +0.102 | [-0.01, +0.21] | establishes nothing |

Positive in all six, determined in five. The effect is largest on hallucinated
designs and smallest on three-helix bundles.

### Exposed hydrophobic area density: robust and very large

Determined in all six families, from -0.691 on hallucinated designs to -1.000 on
EEHEE and HHH, where every design is below every matched natural control.

### Longest repeat: positive everywhere, determined in four

+0.546 hallucination, +0.528 HEEH, +0.368 HHH, +0.217 EHEE, all determined;
EEHEE +0.008 and other +0.181 establish nothing.

### GC: nothing, in five of six families

Only EEHEE is determined, at -0.331. This is consistent with the earlier finding
that the pooled GC differences were composition artefacts.

CAI comes out at exactly 0.000 with a zero-width interval in every family. That
is not a result: both arms are given coding sequences by the same maximum
relative-adaptiveness rule, so their CAI is identical by construction. It is
reported here only because it confirms the gene-layer control is doing what it
should.

## The primary statistic as a systematic panel

Section 8.3 named the normalised position on the random-to-natural axis as H1's
primary statistic. It is applied here to every feature, with a stated threshold:
the gap between the random and natural medians must be at least one median
absolute deviation of those two groups pooled. **14 of 61 features clear it; the
other 47 are read by effect size instead.**

| feature | position | axis width | delta vs natural | q |
|---|---|---|---|---|
| exposed hydrophobic area, density | **1.50** | 1.63 | -0.323 | 4.6e-16 |
| exposed hydrophobic area | **1.48** | 1.57 | -0.305 | 1.8e-14 |
| longest exposed hydrophobic run, density | 1.02 | 2.07 | +0.065 | 0.16 |
| exposed hydrophobic patches | 1.00 | 1.00 | +0.066 | 0.051 |
| longest repeat | 1.00 | 1.00 | +0.232 | 4.2e-09 |
| targeting motifs, raw | 1.00 | 1.00 | +0.012 | 0.87 |
| strand fraction | -0.00 | 2.01 | +0.453 | 4e-33 |
| mean accessibility | 0.00 | 1.19 | -0.429 | 3.6e-27 |
| helix fraction | 0.00 | 1.24 | -0.272 | 8.4e-12 |
| C-terminal accessibility | 0.00 | 1.26 | -0.373 | 7.2e-21 |
| N-terminal accessibility | 0.00 | 1.09 | -0.397 | 2.1e-23 |
| longest A/T homopolymer | -0.00 | 1.00 | +0.111 | 0.0081 |

A second limitation of the statistic is visible in the column: it piles up at
exactly 0.00 and 1.00. Those are integer-valued or near-integer features whose
design median coincides with one endpoint's median, so the ratio collapses to a
boundary and carries no information about position. Only the two
exposed-hydrophobic-area features, which are continuous, give an interior value.
The axis-width threshold does not catch this, and the statistic should be read as
informative on continuous features only.

Benjamini-Hochberg is applied across all 61 numeric features, separately for each
control arm, at a false discovery rate of 0.05. Adjusted values are in the q
columns throughout this note and in `atlas_stats.json`.

## Consequence for H1

H1 predicted that designs carry systematically higher gene-level and proteostatic
liability than matched natural proteins.

**On the proteostatic half the answer is null.** An earlier draft reported
designs as carrying measurably less degron burden; that was an artefact of
comparing against excised domains, and it is withdrawn. Against the whole-protein
arm the pooled weighted degron effect is -0.107 [-0.30, +0.11] and does not
exclude zero. Per family, only EEHEE clears correction, at -0.603 [-0.92, -0.12]
and q = 0.045. No family is enriched: the hallucination family sits at +0.148
[-0.11, +0.45] and its trRosetta subset at +0.232 [-0.09, +0.53]. H5's feature
attribution reaches the same place from the other direction, finding that degron
features carry none of the discrimination between successful and failed
designs.

**What survives is the sequence-composition half.** Designs carry more
low-complexity content and more repeat structure than whole natural proteins,
positive in all six topology families; low-complexity content clears correction
in four of them and repeat structure in one, with two more at q = 0.053 and
0.083. Selection removes these from natural sequences and no step of a
design pipeline addresses them. They are also the two liabilities the amino-acid
axis can clear outright, for about a fifth of a nat, against a codon-only arm
that removes exactly nothing.

**On surface hydrophobicity designs are cleaner than natural proteins**, by the
largest effect in the panel, determined in every family. Part of that is real and
part is a fold difference the matching does not control, since natural proteins
of this length are largely helical while the design set is not.

## The hallucination family is one design campaign, not a family

The per-topology split leaves hallucinated designs as the only family enriched in
accessibility-weighted degron density against matched whole natural proteins,
+0.148 [+0.018, +0.274] on n=148. That family is not one method. Its members carry
their own design sets' identifiers, and it is three campaigns pooled under one
label: 97 named `_TrROS_Hall`, 27 prefixed `GG:`, 24 prefixed `EA:`.

Splitting on the design set, against the same matched whole-protein controls:

| design set | n | delta | 95% CI | verdict |
|---|---|---|---|---|
| trRosetta hallucination | 97 | **+0.232** | [+0.08, +0.39] | determined |
| GG set | 27 | +0.167 | [-0.17, +0.47] | establishes nothing |
| EA set | 24 | +0.142 | [-0.20, +0.48] | establishes nothing |

**The enrichment belongs to the trRosetta hallucination designs**, where it is
larger than the pooled family estimate and determined on its own. The other two
sets point the same way and are too small to say. So the sharper statement is
available: degron load tracks the design campaign, and the campaign carrying it is
network hallucination.

Set beside EEHEE designs at -0.422 [-0.50, -0.34] against the same controls, the
range across design methods is about 0.65 of Cliff's delta, with intervals nowhere
near overlapping. The pooled estimate over all designs, -0.103, describes neither
end and should not be quoted alone.

Two limits on how hard to lean on this. The three subsets are compared against
their own matched controls, so the comparison is internally consistent, but n=97
against n=24 and n=27 means only one of the three could have been determined
whatever the truth. And these are the hallucinated designs that reached an
experimental characterisation, so the survivor caveat in `survivors.md` applies
here too, and degron load is one of the features that caveat bites on.

Exposed hydrophobic area density runs the other way in every set and is
determined in all three, from -0.560 to -0.909, so the campaigns differ on
degrons and agree on hydrophobic exposure.

## Resampling unit

Both natural arms are matched with replacement. The excised arm's 447 records
come from 104 distinct PDB entries and the whole-protein arm's from 105 UniProt
accessions, about 4.3 records per protein in each. A record is therefore not an
independent draw, and resampling records gives intervals that are too narrow on
the natural side.

Every interval in this note against a natural arm is a cluster bootstrap:
clusters are the source proteins, they are drawn with replacement, and all
records belonging to a drawn cluster are carried in. The p-value is read off the
same cluster distribution as the interval so that both refer to one resampling
unit. Features with no spread on either side cannot be tested and are dropped
from the correction, which is what takes the per-family panel to 30 tests.
Scrambled and random carry one record per design, so their intervals are
unchanged.

Two conclusions moved. The pooled degron effect against whole proteins stops
excluding zero, and the apparent degron enrichment in the hallucination family
and its trRosetta subset stops excluding zero. Both had cleared correction under
record resampling. The composition results are unaffected.

`cluster_bootstrap.py` regenerates all of it.
