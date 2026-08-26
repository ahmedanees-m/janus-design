import numpy as np
import pytest

from janus import hosts


@pytest.fixture(scope="session")
def ecoli():
    return hosts.load("ecoli_bl21")


@pytest.fixture
def marginals():
    """Deterministic stand-in for ProteinMPNN unconditional marginals.

    ``concentration`` sets how peaked each position is, so a test can request a
    near-deterministic or near-flat posterior without a backbone or a GPU.
    """

    def build(length: int, seed: int = 0, concentration: float = 2.0) -> np.ndarray:
        rng = np.random.default_rng(seed)
        logits = rng.normal(scale=concentration, size=(length, 20))
        logits -= logits.max(axis=1, keepdims=True)
        return logits - np.log(np.exp(logits).sum(axis=1, keepdims=True))

    return build
