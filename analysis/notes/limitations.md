# Limitations, manuscript-ready

Drafted as the measurements land, so the wording is what was measured rather
than what was hoped for. Section 12.3 requires that every use of "exact",
"optimal" or "globally optimal" names the objective it is exact with respect to.

## The surrogate objective

The dynamic program is exact with respect to ProteinMPNN's unconditional
single-pass marginals, not its autoregressive conditional posterior. The
conditional log-probability at a position depends on the residues chosen
elsewhere, so it does not decompose over lattice nodes and no dynamic program
can optimise it. We therefore optimise a stated surrogate exactly and rescore
candidates under the conditional model, and we measure the gap between the two
rankings rather than assuming it is small.

## Coupling in the lattice is codon-pair bias alone

Of the terms carried inside the parse, only codon-pair bias couples adjacent
positions. Fold compatibility, codon adaptation and GC content are position
local, so the objective they define is separable and per-position greedy attains
its optimum exactly. What the exact parse buys is therefore governed by the
strength of the codon-pair term, and we report that dependence directly.

## The folding term is realised in Tier 2, not in the lattice

mRNA base pairing forms at arbitrary range, so a path's folding energy is not a
sum over its positions and the term cannot enter the parse at this granularity.
Folding it in requires a nucleotide-level lattice, which is LinearDesign's own
construction. We chose not to build it: the amino-acid layer is the contribution
here, and the nucleotide lattice would move effort into the part of the method
that is not. The consequence is stated plainly rather than left implicit. The
term that motivates lattice parsing in the mRNA-design literature is, in this
work, a rescoring term.

## Generate-and-rank has a measured boundary

The two-tier design proposes candidates with the parser and ranks them with the
terms the parser cannot carry. That is sound while the out-of-parser terms are
secondary to the in-parser ones. We measured where it stops being sound: as the
folding weight rises, the Tier-1 prefix ceases to be the better candidate pool
and broad sampling of the same shell overtakes it. We report the crossover
rather than operating on one side of it silently.

## No wet-lab validation

Structural neutrality is assessed against measured folding free energies from
cDNA display proteolysis on the subset of backbones those measurements cover,
and the covered fraction is reported. No construct in this work has been
expressed.
