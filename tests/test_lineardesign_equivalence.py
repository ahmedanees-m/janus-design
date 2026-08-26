"""At zero budget the lattice collapses to a fixed-protein codon DFA.

LinearDesign's `lambda` weights CAI against folding energy and defaults to zero,
which is pure MFE; there is no exposed weight that removes the folding term, so
in general the two objectives never coincide. Correctness is therefore
established against exhaustive enumeration, which is ground truth rather than
agreement with another implementation.

One limit does coincide. At a large lambda the CAI term dominates and
LinearDesign returns a CAI-maximal assignment, which the codon layer must
reproduce exactly on the same codon table. That test runs when the binary and its
table are supplied through LINEARDESIGN_BIN and LINEARDESIGN_TABLE.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import itertools
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

from janus import Weights, design
from janus.genetic_code import AA_ALPHABET, SYNONYMOUS
from janus.objectives import assemble
from janus.lattice import build

LINEARDESIGN = importlib.util.find_spec("LinearDesign") or importlib.util.find_spec("lineardesign")


def argmax_protein(marginals: np.ndarray) -> str:
    return "".join(AA_ALPHABET[i] for i in marginals.argmax(axis=1))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_zero_budget_holds_the_protein_at_the_argmax(ecoli, marginals, seed):
    log_marginals = marginals(length=50, seed=seed)
    for weights in (Weights(mpnn=1.0), Weights(mpnn=1.0, cai=2.0, cpb=1.0, gc=0.5)):
        result = design(log_marginals, ecoli, weights=weights, delta=0.0, k=1)[0]
        assert result.protein == argmax_protein(log_marginals)


def test_zero_budget_lattice_admits_only_synonymous_codons(ecoli, marginals):
    log_marginals = marginals(length=30, seed=4)
    lattice = build(log_marginals, delta=0.0)
    expected = argmax_protein(log_marginals)
    for position, residue in enumerate(expected):
        assert lattice.amino_acids[position] == (residue,)
        assert lattice.codons[position] == SYNONYMOUS[residue]


@pytest.mark.parametrize(
    "weights",
    [
        Weights(mpnn=1.0, cai=1.0),
        Weights(mpnn=1.0, cpb=1.0),
        Weights(mpnn=1.0, cai=1.0, cpb=1.0, gc=-2.0),
    ],
)
def test_zero_budget_optimum_matches_exhaustive_search(ecoli, marginals, weights):
    """Six residues enumerates fully, which is the only way to tell an optimum
    from a strong local maximum."""
    log_marginals = marginals(length=6, seed=9)
    lattice = build(log_marginals, delta=0.0)
    node_scores, edge_scores = assemble(lattice, ecoli, log_marginals, weights)

    best = max(
        sum(node_scores[i][state] for i, state in enumerate(states))
        + sum(edge_scores[i][states[i], states[i + 1]] for i in range(len(states) - 1))
        for states in itertools.product(*(range(len(c)) for c in lattice.codons))
    )

    result = design(log_marginals, ecoli, weights=weights, delta=0.0, k=1)[0]
    assert result.score == pytest.approx(best)


def test_cai_weight_drives_the_sequence_toward_preferred_codons(ecoli, marginals):
    log_marginals = marginals(length=40, seed=13)
    neutral = design(log_marginals, ecoli, weights=Weights(mpnn=1.0), delta=0.0, k=1)[0]
    adapted = design(log_marginals, ecoli, weights=Weights(mpnn=1.0, cai=5.0), delta=0.0, k=1)[0]

    assert adapted.protein == neutral.protein
    assert adapted.terms["cai"] > neutral.terms["cai"]
    assert adapted.terms["cai"] == pytest.approx(1.0, abs=1e-9)


@pytest.mark.skip(
    reason="not applicable: LinearDesign's objective carries mRNA folding, which "
           "never enters this parse, so no weight makes the objectives coincide"
)
def test_zero_budget_reproduces_lineardesign():
    """Kept as a marker so the omission is visible rather than silent."""


LINEARDESIGN_BIN = os.environ.get("LINEARDESIGN_BIN")
LINEARDESIGN_TABLE = os.environ.get("LINEARDESIGN_TABLE")


def read_frequency_table(path):
    """Relative adaptiveness from LinearDesign's own codon usage CSV."""
    families = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        parts = line.strip().split(",")
        if len(parts) != 3 or parts[0].startswith("#"):
            continue
        codon, residue, frequency = parts[0].replace("U", "T"), parts[1], float(parts[2])
        families.setdefault(residue, {})[codon] = frequency
    adaptiveness = {}
    for residue, codons in families.items():
        if residue == "*":
            continue
        highest = max(codons.values())
        for codon, frequency in codons.items():
            adaptiveness[codon] = frequency / highest if highest else 0.0
    return adaptiveness


@pytest.mark.skipif(not (LINEARDESIGN_BIN and LINEARDESIGN_TABLE),
                    reason="LinearDesign binary and codon table not supplied")
def test_codon_layer_reproduces_lineardesign_in_the_cai_limit(ecoli, marginals):
    """LinearDesign's lambda weights CAI, so the limit is CAI-maximal, not MFE-free.

    There is no exposed weight that removes the folding term, so the two objectives
    still never coincide in general. At a large lambda, though, the CAI term
    dominates and LinearDesign returns a CAI-maximal assignment, which the codon
    layer must reproduce exactly on the same codon table. That is agreement with an
    independent implementation rather than with our own enumeration.
    """
    adaptiveness = read_frequency_table(LINEARDESIGN_TABLE)
    host = dataclasses.replace(ecoli, relative_adaptiveness=adaptiveness,
                               codon_pair_scores={})

    log_marginals = marginals(length=40, seed=7)
    ours = design(log_marginals, host, weights=Weights(mpnn=1.0, cai=1000.0),
                  delta=0.0, k=1)[0]

    # The binary resolves its shared library by a path relative to the working
    # directory, so it has to be run from beside its own src tree.
    result = subprocess.run(
        [str(Path(LINEARDESIGN_BIN).resolve()), "1000", "0",
         str(Path(LINEARDESIGN_TABLE).resolve())],
        cwd=Path(LINEARDESIGN_BIN).resolve().parent,
        input=ours.protein, capture_output=True, text=True, check=True,
    )
    match = re.search(r"mRNA sequence:\s+([ACGU]+)", result.stdout)
    assert match, result.stdout[-400:]
    theirs = match.group(1).replace("U", "T")

    assert theirs == ours.cds
