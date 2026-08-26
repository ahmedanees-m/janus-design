# Measured stability

MegaScale folding free energies, 19 August 2026.

Mutant names carry a position whose numbering is offset against the synthesised
construct by a per-parent amount. The offset is recovered by agreement between
the stated wild-type residue and the translated construct rather than assumed,
and a parent is used only when at least 95 percent of its mutants agree. 710 of
the 713 backbones carrying measurements resolved; the remaining three did not
and were dropped.

## Coverage

| | |
|---|---|
| backbones with measurements | 713 |
| numbering resolved and backbone aligned | 710 |
| backbones contributing at least one covered substitution | 666 |
| covered substitutions | 58,456 |

This is far beyond what the plan budgeted for, which anticipated a minority of
designs being covered.

## Substitution cost

Each arm proposes one substitution per position from the same delta-shell, so
the contrast isolates which residue was chosen rather than how many were
changed.

| arm | n | median | mean | below -0.5 | below -1.0 |
|---|---|---|---|---|---|
| JANUS | 18,682 | +0.023 | +0.072 | 17.0% | 7.5% |
| ProteinMPNN, T = 0.1 | 18,579 | +0.018 | +0.062 | 17.4% | 7.7% |
| random from the shell | 21,195 | -0.023 | +0.013 | 19.3% | 8.4% |

The claim here is non-inferiority, which is what section 8.3 asked for, and it
holds with room. Effect sizes, which matter more than p values at this n:

| comparison | Cliff's delta | 95% interval | median shift |
|---|---|---|---|
| JANUS against random from the shell | +0.054 | [+0.043, +0.066] | +0.046 kcal/mol |
| JANUS against ProteinMPNN, T = 0.1 | +0.007 | [-0.005, +0.018] | +0.005 kcal/mol |

Both sit below the 0.15 effect-size floor fixed in the analysis plan, and that
floor applies here as much as it applied to the degron panel. The difference against
random shell draws is statistically unambiguous, one-sided p = 3.6e-21 on 18,682
and 21,195 substitutions, and negligible in size: the destabilisation rate below
-1 kcal/mol differs by 0.84 percentage points.

The correct statement is therefore that JANUS substitutions are **not more
destabilising** than random draws from the same shell, with the lower interval
bound excluding any disadvantage, and are indistinguishable from ProteinMPNN
sampling. An earlier draft said "measurably less destabilising"; that overstated
an effect of 0.054 and has been withdrawn.

The median substitution is neutral in all three arms. H4 holds as structural
neutrality, which is what it claimed.

## Additivity

230,992 double mutants, each compared against the sum of its two singles.

| set | n | median deviation | spread | r |
|---|---|---|---|---|
| all pairs | 230,992 | +0.787 | 1.785 | 0.775 |
| natural parents | 208,654 | +0.850 | 1.767 | 0.770 |
| designed parents | 22,338 | +0.323 | 1.737 | 0.628 |
| both sites inside the shell | 3,120 | +0.597 | 1.884 | 0.774 |
| designed, both inside the shell | 220 | +0.186 | 1.943 | 0.734 |

Deviation is observed minus additive, in kcal/mol.

The large positive bias across all pairs is very likely a measurement-range
effect rather than real positive epistasis. Proteolysis-derived free energies
saturate at the bottom of the dynamic range, so a pair of strongly
destabilising substitutions cannot register the sum of its parts and the
observed value is necessarily less negative than the additive prediction. The
bias falls as the sets narrow toward mild substitutions, from +0.85 on natural
parents to +0.19 on designed parents inside the shell, which is the pattern that
artefact predicts.

Inside the shell on designed parents, which is the case JANUS actually creates,
the mean deviation is -0.014 kcal/mol with a standard error of 0.131 across 220
pairs. **Additivity is unbiased to within 0.26 kcal/mol at 95 percent
confidence.** That bound is more useful than the p value of 0.48 it replaces,
because it states what the data can exclude rather than what they failed to
detect.

It is not precise. The per-pair spread is 1.94 kcal/mol and 50 percent of all
pairs deviate by more than 1 kcal/mol. Additive prediction of a multi-substitution
design is unbiased on average and unreliable for any individual design.

## What this supports and what it does not

The plan's decision rule asked whether to cap the number of simultaneous
substitutions. The median is unbiased inside the shell, so no cap is justified
on bias grounds, and the n of 220 is thin enough that this should be stated as a
limitation rather than as a general result about inverse folding.

The statement is narrower than the plan anticipated: additive scoring of
multiple shell substitutions introduces no systematic bias on designed
backbones, and carries a per-pair uncertainty of roughly 2 kcal/mol that
propagates into any multi-substitution prediction.
