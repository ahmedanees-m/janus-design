"""Tier 2: rescoring k-best candidates under terms the parser cannot carry.

The parser proposes a small, high-quality candidate set from an astronomically
large space; terms that do not decompose over lattice nodes rank it. The
initiation window and the whole-transcript folding energy are scored here, and so
is the protein-level liability panel, which arrives as a caller-supplied function
so the panel can be swapped without touching the ranking.

The conditional ProteinMPNN likelihood is the third term of this kind. Scoring it
needs the model itself and a GPU, which the solver deliberately does not depend
on, so it is computed outside this package by the analysis scripts and folded in
there.

Terms arrive in different units: the Tier-1 score is in nats, folding energies
in kcal/mol. Adding them directly makes the weight an accidental unit
conversion, which is why an unnormalised weight does not transfer between
backbones of different length. Dividing each term by its spread across the
candidate pool first makes the weights dimensionless and comparable.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import inf
from statistics import pstdev

from .design import Design
from .hosts import Host
from .objectives import mrna


@dataclass(frozen=True)
class FoldingWeights:
    """Multipliers on the folding terms.

    ``initiation`` rewards an open 5' window, so a positive value pushes the
    window's folding energy toward zero. ``mfe`` is applied with the host's
    ``global_mfe_sign``, which is zero for *E. coli*: no prokaryotic result
    supports rewarding whole-transcript structure. ``liability`` penalises a
    protein-level burden supplied by the caller, which no codon choice can
    change.

    With ``Scales`` supplied the weights are in units of one pool standard
    deviation per term, so a value of one trades a standard deviation of Tier-1
    score for a standard deviation of the folding term.
    """

    initiation: float = 1.0
    mfe: float = 0.0
    liability: float = 0.0


@dataclass(frozen=True)
class Scales:
    """Spread of each term across a candidate pool, used to make weights unitless."""

    tier1: float = 1.0
    initiation: float = 1.0
    mfe: float = 1.0
    liability: float = 1.0


@dataclass(frozen=True)
class Rescored:
    design: Design
    tier1: float
    initiation_energy: float
    global_mfe: float
    total: float
    liability: float = 0.0

    @property
    def cds(self) -> str:
        return self.design.cds

    @property
    def protein(self) -> str:
        return self.design.protein


def _spread(values, floor=1e-9):
    """Population standard deviation of a term over the pool, or infinity if it has none.

    A term that takes the same value on every candidate cannot order them, so
    dividing it out is the correct treatment and infinity does exactly that. A
    small floor instead would amplify a numerically-zero spread into a term that
    dominates every other, which is the opposite of what a constant term means.
    """
    if len(values) < 2:
        return 1.0
    spread = pstdev(values)
    return spread if spread > floor else inf


def pool_scales(designs: list[Design], host: Host, use_mfe: bool = False,
                liability=None) -> Scales:
    """Measure each term's spread over the pool that will be ranked."""
    if not designs:
        return Scales()
    tier1 = [d.score for d in designs]
    initiation = [mrna.initiation_energy(d.cds, host) for d in designs]
    mfe = [mrna.global_mfe(d.cds, host) for d in designs] if use_mfe else [0.0]
    burden = [liability(d.protein) for d in designs] if liability else [0.0]
    return Scales(tier1=_spread(tier1), initiation=_spread(initiation),
                  mfe=_spread(mfe), liability=_spread(burden))


def rescore(
    designs: list[Design],
    host: Host,
    weights: FoldingWeights | None = None,
    scales: Scales | None = None,
    synthesisable_only: bool = False,
    liability=None,
) -> list[Rescored]:
    """Rank candidates under the folding terms, best first.

    Without ``scales`` the terms are combined in their native units, which is
    what the earlier unnormalised runs did and is kept so those are reproducible.
    """
    weights = weights or FoldingWeights()
    scales = scales or Scales()
    sign = host.initiation.global_mfe_sign

    scored = []
    for design in designs:
        if synthesisable_only and not design.synthesisable:
            continue
        # Always measured, never skipped when the weight is zero: this is a
        # reported field as well as a scoring term, and a zero weight means
        # the term does not steer the ranking, not that the energy is zero.
        initiation = mrna.initiation_energy(design.cds, host)
        mfe = mrna.global_mfe(design.cds, host) if (weights.mfe and sign) else 0.0
        burden = liability(design.protein) if liability else 0.0
        total = (
            design.score / scales.tier1
            + weights.initiation * initiation / scales.initiation
            - weights.mfe * sign * mfe / scales.mfe
            - weights.liability * burden / scales.liability
        )
        scored.append(
            Rescored(
                design=design,
                tier1=design.score,
                initiation_energy=initiation,
                global_mfe=mfe,
                liability=burden,
                total=total,
            )
        )

    scored.sort(key=lambda r: -r.total)
    return scored
