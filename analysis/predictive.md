# Predictive filtering

H5. Garcia, Dixit and Rocklin benchmark, bioRxiv 2025.07.29.667290 v2
supplementary CSV, CC-BY. 614 designs, 269 successful (43.8 percent).

An earlier note recorded this as blocked behind a publisher paywall. That was
wrong: the journal version is paywalled, the preprint supplement is not, and it
carries sequences, labels, ESMFold pLDDT, ProteinMPNN scores and fold class.

## Scope, narrowed by what the benchmark contains

The benchmark supplies sequences but no backbones. All 614 were folded with
ESMFold so the accessibility weighting could be tested; the section below reports
what it is worth, which is close to nothing once the obvious circularity is
controlled. The headline score remains the sequence-only one.

28 sequence-level features, five-fold cross-validated, logistic regression with
standardisation. Protein-level repeat content was added to the set, since it
needs no structure; the accessibility-dependent liabilities still cannot be
computed here.

## Result

| score | ROC AUC | PR AUC |
|---|---|---|
| ESMFold pLDDT alone | 0.731 | 0.634 |
| liability score alone | 0.672 | 0.582 |
| pLDDT and liability | **0.788** | **0.704** |

DeLong against pLDDT alone: the liability score by itself is *worse*
(p = 0.0492), and the combination is better (p = 0.000164).

The improvement is +0.057 in AUC over pLDDT alone. Modest, consistent with the
source paper's finding that combining metrics gives only moderate improvement,
and on the right side of zero.

## Which liabilities carry it

A single AUC over a bag of features says nothing about which liability is doing
the work, and there was a tension that needed it. The atlas finds designs are not
degron-enriched against natural proteins, while this score adds AUC over pLDDT.
Each family reported alone, and with itself removed from the full set:

| family | features | alone | without | drop |
|---|---|---|---|---|
| degron motifs | 12 | 0.521 | 0.673 | **-0.001** |
| other ELM motifs | 9 | 0.606 | 0.651 | +0.022 |
| complexity and repeats | 2 | 0.512 | 0.634 | +0.038 |
| bulk composition | 4 | 0.579 | 0.628 | +0.045 |
| length | 1 | 0.574 | 0.671 | +0.001 |

**Degrons carry none of it.** Alone they are at chance, 0.521, and removing all
twelve of them changes the score by -0.001. The tension dissolves: H1 finds no
degron enrichment between groups and H5 finds no degron discrimination within
designs, which is one consistent finding rather than two awkward ones. The paper
should not describe this score as degron-driven.

What does carry it is bulk composition, the largest single contribution at
+0.045, then complexity and repeat content at +0.038, then the non-degron ELM
motif classes at +0.022. None of these is strong alone: the best family on its
own reaches 0.606 against a full score of 0.672. The discrimination is a
composite, and no single liability in it would be worth reporting by itself.

That is the through-line the atlas and the exchange rate also point at. The
liabilities that separate designs from natural proteins, and the ones the
amino-acid axis can remove cheaply, are complexity and repeat content, not
degrons.

## The fold stratification is mostly underpowered

614 designs across 53 fold classes is about a dozen each. With bootstrap
intervals:

| fold class | n | success | pLDDT AUC | liability AUC |
|---|---|---|---|---|
| Foldit | 172 | 33% | **0.882 [0.82, 0.93]** | **0.617 [0.53, 0.70]** |
| HBI_b | 43 | 21% | 0.595 [0.33, 0.83] | 0.608 [0.43, 0.78] |
| R2x3_BP1_A | 18 | 39% | 0.610 [0.31, 0.88] | 0.662 [0.32, 0.93] |
| R3x3_BP2 | 17 | 65% | 0.500 [0.19, 0.83] | 0.758 [0.42, 1.00] |
| NF1 | 16 | 38% | 0.517 [0.19, 0.83] | 0.583 [0.27, 0.87] |
| Di-III | 12 | 33% | 0.656 [0.27, 1.00] | 0.781 [0.49, 1.00] |

Only Foldit is determined. There the two intervals do not overlap and pLDDT is
genuinely the stronger discriminator.

**The dissociation claimed in an earlier draft is withdrawn.** That draft read
pLDDT at 0.500 on R3x3_BP2 against a liability score of 0.758 as the two being
informative on different topologies. On 17 designs those intervals are
[0.19, 0.83] and [0.42, 1.00]; they overlap across almost their whole range and
establish nothing. The pooled comparison is the claimable result.

## Contamination cuts in our favour, and cannot be stratified away

R13 asked for a PDB release-date stratification separating designs inside
ProteinMPNN and ESMFold training data from those outside it. That split cannot be
made here: the benchmark is drawn from 11 studies published between 2012 and
2021, so essentially every design predates both models' training cutoffs. There
is no clean stratum.

That is worth stating positively rather than leaving as a pending analysis. The
baseline is contaminated and the challenger is not. ESMFold pLDDT is a learned
score evaluated on designs that were almost certainly in its training data, which
inflates it. The liability score is hand-built motif and composition features
with no learned component and no exposure to these sequences. An uncontaminated
challenger adding 0.041 AUC to an inflated baseline is a conservative estimate of
what it would add to an uncontaminated one.

What remains genuinely missing is accessibility weighting. The benchmark supplies
no backbones, and weighting enlarged the degron effect substantially in the
atlas, so the liability-alone figure of 0.634 is a lower bound. Folding the 614
sequences would supply the structures; it is specified and not run.

## Accessibility weighting is not worth what it looks like

All 614 sequences were folded with ESMFold, mean pLDDT 72.9, and the motif
features recomputed with accessibility weighting. 411 models clear a pLDDT floor
of 70; the other 203 do not supply accessibility.

| cohort | n | pLDDT | sequence only | weighted | gain | combined |
|---|---|---|---|---|---|---|
| all designs | 614 | 0.731 | 0.672 | 0.779 | **+0.106** | 0.796 |
| usable models | 411 | 0.591 | 0.688 | 0.697 | **+0.009** | 0.721 |

The first row looks like a large gain and is not one. The structures come from
ESMFold and the comparator is ESMFold pLDDT, so a feature read off the model
carries the confidence signal: a design the model folds badly gets an extended,
highly accessible structure, and its weighted liabilities restate a low pLDDT.
Worse, the 203 designs without a usable model receive a systematically different
feature vector, and having no confident model is itself a predictor of failure.

Restricting to designs that all have usable models removes both effects, and the
weighting is then worth **+0.009**. That is the number to report. The expectation
that a structure-aware score would be strictly better than the sequence-only one,
carried since the atlas showed weighting enlarged the degron effect, is not borne
out on this benchmark.

The restricted cohort is a harder setting for everything, since filtering on
pLDDT removes most of pLDDT's own discriminative range: it falls from 0.731 to
0.591 there. That is a reason the restricted row cannot be compared to the
unrestricted one across columns, and it does not affect the within-row comparison
that matters.

## The ordering reverses between the cohorts, and that is the result

Read down the columns rather than across, and the two rows disagree about which
score is better. DeLong, two-sided, computed inside each cohort rather than on
the pooled set:

| cohort | n | pLDDT | sequence only | p | weighted | p | combined | p |
|---|---|---|---|---|---|---|---|---|
| all designs | 614 | 0.731 | 0.672 | 0.049 | 0.779 | 0.013 | 0.796 | 6.7e-05 |
| usable models | 411 | 0.591 | 0.688 | **0.014** | 0.697 | 0.0055 | 0.721 | 7.2e-05 |

On all 614, pLDDT beats the sequence-only liability score, 0.731 against 0.672.
**On the 411 with a usable model, the liability score beats pLDDT, 0.688 against
0.591, p = 0.014.** The ordering is not just narrowed, it is reversed, and both
directions clear significance.

The mechanism is the same one that inflates the weighting gain. Restricting to
confident models strips out the easy cases, the designs a folding model fails on
so badly that low pLDDT flags them without any liability being involved. What is
left is the subset where the structure looks fine and the design fails anyway, and
there a composition-and-complexity score is the better discriminator.

**The 411 cohort is primary for the comparison between scores.** The 614 figure
lets "has no confident model" act as a covert predictor, and it favours the
pLDDT-derived side of the comparison, so it cannot settle which score is better.
The 614 figure is the right one for the deployment question, because in practice
designs without confident models do arrive and do have to be filtered.

Combined beats both components in both cohorts, 0.796 and 0.721, at p below 1e-4.
That is the defensible headline: the liability score carries information pLDDT
does not, most visibly where pLDDT has least to say.

Attribution on the restricted cohort tells the same story as the sequence-only
version, with the model-derived features reported as their own family:

| family | features | alone | without | drop |
|---|---|---|---|---|
| bulk composition | 4 | 0.593 | 0.651 | +0.046 |
| complexity and repeats | 2 | 0.573 | 0.667 | +0.030 |
| model geometry | 11 | 0.587 | 0.685 | +0.013 |
| other ELM motifs | 9 | 0.628 | 0.691 | +0.006 |
| length | 1 | 0.496 | 0.698 | -0.001 |
| degron motifs | 12 | 0.448 | 0.701 | -0.003 |

Degrons are below chance alone and negative on removal. Everything read off the
model contributes +0.013 between eleven features. Composition and complexity
remain the whole of it.
