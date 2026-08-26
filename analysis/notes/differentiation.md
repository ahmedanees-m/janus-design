# Differentiation from the nearest prior art

Sources read 19 August 2026.

## R1: CodonMPNN

*Stark, Padia, Balla, Diao and Church, arXiv:2409.17265.*

CodonMPNN generates a codon sequence conditioned on a backbone structure and an
organism label, so the amino acid is implied by the codon and the method already
occupies the joint amino-acid and codon space. Any claim to have been first to
formulate design that way would be wrong, and we do not make one. The
differentiation is in what is optimised and how.

CodonMPNN is a learned sampler, not an optimiser. It infers codon preference
from natural genes, so its only nucleic-acid signal is codon usage as it appears
in the training organisms. It carries no term for mRNA secondary structure, for
synthesisability, for cryptic regulatory elements, or for intracellular
degradation, and because those preferences are absorbed into network weights
there is no way for a user to add an objective or reweight an existing one
without retraining. Its stated evaluation is recovery of wild-type codons, a
metric that rewards imitating natural genes. For a de novo protein there is no
wild type to recover, which is the case where the target is least well defined
by imitation and most in need of an explicit objective.

JANUS states the objective, optimises it, and lets it be changed. The difference
is explicit multi-objective optimisation against learned imitation, not priority
of formulation.

## R2: MoMPNN and the preference-alignment genre

*Hou, Liu, Shi, Liu, Yang and Tang, arXiv:2603.06748, ProtAlign and MoMPNN.*

MoMPNN fine-tunes ProteinMPNN by semi-online multi-objective preference
alignment, using computational property predictors to build preference pairs
across developability objectives. It is the strongest current property-aware
inverse folder and it optimises protein biophysics: solubility, thermostability
and evolutionary perplexity. Nothing in its property panel is nucleic-acid level
and nothing is proteostatic.

The methodological point is narrower than "preference optimisation is the wrong
tool". Preference optimisation is the right tool for a black-box reward. Codon
adaptation, codon-pair bias and GC content are not black boxes: they are
additive over sequence positions or over adjacent codon pairs, which means the
optimum can be computed rather than approached. Spending a fine-tuning run to
approximate a quantity that decomposes exactly discards a solution that is
already available, and it fixes the objective into weights so that changing the
target organism or the vendor constraint set means retraining. We report the
measured consequence rather than asserting it: where the objective is separable,
per-position selection is exactly optimal and the parser buys nothing, and where
it couples, the parser leads by an amount that scales with the coupling.

## Relationship to LinearDesign and LinearDesign2

*Zhang et al., Nature 621:396 (2023); Liu, Gao, Zhang and Fang, arXiv:2410.20781.*

LinearDesign is exact over synonymous codon sequences for a fixed protein,
obtained by lattice parsing over a codon automaton. LinearDesign2 extends the
problem to a second axis, the 5' untranslated region, optimising translation
initiation efficiency, codon adaptation and folding energy together.

Two things about it matter here, and one corrects the assumption we began with.

First, it does not keep the parse. Its workflow alternates between evolutionary
optimisation of the 5' UTR conditioned on the CDS and optimisation of the CDS
conditioned on the 5' UTR, and its CDS step performs synonymous codon
substitutions. That is block-coordinate search with evolutionary operators, and
no global optimality is claimed. The precedent for extending this framework to a
second degree of freedom therefore comes with the precedent for giving up
exactness in doing so.

Second, the axis is different. LinearDesign2 varies the untranslated region and
holds the protein fixed; JANUS varies the protein and holds the untranslated
region fixed at a stated host sequence. The two are complementary rather than
competing, and combining all three axes is a later paper.

**Correction to our own framing.** We had treated LinearDesign2 as the
LinearDesign authors extending their own lattice, which sharpened the scooping
risk. It is not: LinearDesign2 comes from the PaddleHelix team at Baidu, and
LinearDesign from the group of Zhang and colleagues, though both carry Baidu
affiliations. The risk is better stated as a well-resourced industrial team
building on a published framework than as the original authors returning to it.

## Search record

Re-ran the section 4.3 term list against the arXiv API and Europe PMC on
19 August 2026, adding three terms aimed at the specific construction here:
"amino acid degenerate lattice", "joint protein codon dynamic programming",
"inverse folding lattice parsing".

Every term returned zero hits in both indexes, including the original eleven.
No published work implements an amino-acid-degenerate lattice for joint residue
and codon selection. The algorithmic claim stands.

Citation alerts to set: LinearDesign, LinearDesign2, CodonMPNN, MoMPNN,
CodonTransformer.
