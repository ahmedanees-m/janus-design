# Which layer buys the initiation term

120 designed backbones, host *E. coli* BL21(DE3), delta 1.0 nats. Three arms
ranked under the same normalised objective, all measured against the Tier-1
optimum.

- **synonymous**: protein pinned to the Tier-1 optimum, only codons move. This is
  what a fixed-protein codon optimiser can reach.
- **residue**: residues drawn from the shell, each candidate's codons at their
  exact Tier-1 optimum, so any initiation change is attributable to the residue.
- **joint**: both move, which is the full method.

## Result

| weight | arm | initiation gained | Tier-1 cost | residue changes | codon changes |
|---|---|---|---|---|---|
| 0.125 | synonymous | 1.918 | 0.094 | 0.00 | 1.12 |
| 0.125 | residue | 2.543 | 8.234 | 17.35 | 19.88 |
| 0.125 | **joint** | **2.190** | **0.054** | **0.56** | **1.25** |
| 0.5 | synonymous | 1.976 | 0.123 | 0.00 | 1.27 |
| 0.5 | residue | 4.883 | 9.790 | 19.22 | 21.57 |
| 0.5 | **joint** | **2.258** | **0.113** | **0.71** | **1.47** |

Gains in kcal/mol, costs in nats, changes as counts per design.

**At the operating point, synonymous codon changes alone recover 87.6 percent of
the joint gain**, and 87.5 percent at a weight of 0.5. The joint arm changes no
residue at all on 55.8 percent of backbones, with a median of zero substitutions
and a maximum of three.

The identity figure of 0.987 quoted earlier is amino-acid identity, not
nucleotide identity. On a 43-residue design it is about half a substitution, and
that is the whole amino-acid contribution at this weight.

## What this means for the claim

The amino-acid axis is what distinguishes this work from fixed-protein codon
optimisation, and at the cheap operating point it contributes roughly an eighth
of the initiation gain. That has to be said plainly rather than left for a
reviewer to derive from an identity number.

The axis is not inert; it is expensive. The residue arm reaches 2.5 kcal/mol at a
weight of 0.125 and 4.9 at 0.5, more than the joint arm in both cases, but pays 8
to 10 nats of Tier-1 score to do it against the joint arm's 0.05 to 0.11. Roughly
a hundredfold more, for about twice the gain. The frontier is what prices that
trade, and the defensible operating point sits precisely where the amino-acid
axis matters least.

The correct summary is therefore: for the 5' initiation term specifically, most
of the cheap gain is synonymous, and the joint search matters when the gain has
to be larger than a fixed protein can supply. The claim that the joint space is
worth searching cannot rest on this term alone.

## The enumeration objection, answered rather than deferred

The initiation window covers codons 0 to 12, so exhaustive search over a single
substitution in that region is 13 positions times 19 alternatives times their
codons, on the order of a thousand candidates. Anyone can enumerate that. If the
intervention at the operating point is one codon change and half a substitution,
enumeration finds it.

That is the right objection and the answer is the one already in the methods: the
parser is a measurement instrument here. Enumeration finds a good candidate; it
does not tell you what the candidate cost, because the cost is measured against
the exact optimum of the fold and codon objective over roughly 3.5 to the power
43 residue assignments. Without that optimum there is no exchange rate to quote,
no attainable range to normalise shortfalls against, and no way to know whether
the enumerated candidate is near the frontier or far from it.
