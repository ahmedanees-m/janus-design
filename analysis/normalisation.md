# Term normalisation and the candidate pools

Two properties of the Tier-2 scoring path determine how a reported number should
be read: how a term with no spread is normalised, and which candidate pool an
exchange rate was measured on.

## A term with no spread over the pool

Each term is divided by its spread over the candidate pool being ranked, so that
a weight expresses a preference between terms rather than a conversion between
their units. A term taking the same value on every candidate carries no ordering
information. Its spread is therefore reported as infinite and the term drops out
of the ranking. Dividing by a small floor instead would turn a numerically zero
spread into a term that dominates every other, which inverts what a constant
term means.

This case is reached whenever a liability term is constant across a shell: a
design carrying no low-complexity window has none anywhere in its shell either.
Over a pool of several hundred distinct coding sequences the gene-level terms are
never constant, so the case arises on the protein-level panel alone.

`_spread` in `src/janus/rescore.py` implements this, and `tests/test_mrna.py`
covers it.

## Reported energies at zero weight

`Rescored.initiation_energy` is a reported field as well as a scoring term, so
the folding call is made whatever the initiation weight is. A zero weight means
the term does not steer the ranking, not that the energy is zero. A regression
test in `tests/test_mrna.py` asserts that a zero weight changes the ranking
rather than the reported energy.

## Which pool an exchange rate is measured on

The exchange rate between Tier-1 score and initiation-window energy is a property
of the pool being ranked as well as of the weight, so every reported rate names
its pool. Two are used:

- the **k-best prefix**, a 1000-path Tier-1 prefix, which is the pool the
  two-tier bound endorses across the operating region and the pool the reported
  rate is quoted on
- the **prefix plus shell draws**, which adds uniform draws from the same shell
  and reaches further into protein space at a higher Tier-1 cost

Medians over 120 designed backbones, at a budget of 92.58 nats:

| pool | weight | nats paid | kcal opened | fold change | budget |
|---|---|---|---|---|---|
| k-best prefix | 0.125 | 0.0078 | 1.450 | 1.92x | 0.008% |
| k-best prefix | **0.250** | **0.0362** | **2.000** | **2.46x** | **0.039%** |
| k-best prefix | 0.500 | 0.0735 | 2.500 | 3.08x | 0.079% |
| k-best prefix | 1.000 | 0.0985 | 2.500 | 3.08x | 0.106% |
| prefix and draws | 0.125 | 0.0492 | 2.000 | 2.46x | 0.053% |
| prefix and draws | 0.250 | 0.0566 | 2.000 | 2.46x | 0.061% |
| prefix and draws | 0.500 | 0.0606 | 2.000 | 2.46x | 0.065% |
| prefix and draws | 1.000 | 0.1016 | 3.350 | 4.52x | 0.110% |

Both pools reach 2.00 kcal/mol at the operating weight, the prefix for 0.0362
nats and the hybrid for 0.0566. Above a weight of 1.0 the hybrid pool's cost
rises to 10.1 nats, where the shell draws begin to win, which is the two-tier
boundary expressed as an exchange rate.

All values are medians. What the shell can open is right-skewed, so a mean sits
above the typical backbone and is not comparable with a median reported
elsewhere.

The pool is recorded in each output file. `analysis/expression.md` carries both
tables in full.

## Scope of the design-only inputs

The folding frontier, the layer decomposition and the two-tier bound are
statements about designs, so each filters its marginals directory to the 447
designed backbones rather than averaging in the 415 natural domains that share
it.

## Files

- `src/janus/rescore.py`, `_spread` and the initiation reporting
- `tests/test_mrna.py`, the two regression tests
- `analysis/scripts/folding_sweep.py`, `decompose.py`, `rescore_recall.py`,
  `expression.py`
