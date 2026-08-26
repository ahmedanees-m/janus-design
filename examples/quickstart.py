"""Smallest end-to-end run: marginals in, coding sequences out, no downloads."""
import time
import numpy as np
import janus

rng = np.random.default_rng(0)
logits = rng.normal(scale=1.2, size=(60, 20))
logits -= logits.max(axis=1, keepdims=True)
marginals = logits - np.log(np.exp(logits).sum(axis=1, keepdims=True))

started = time.perf_counter()
host = janus.load_host("ecoli_bl21")
designs = janus.design(marginals, host,
                       weights=janus.Weights(mpnn=1.0, cai=0.5, cpb=0.3),
                       delta=1.0, k=25)
elapsed = time.perf_counter() - started
print(f"{len(designs)} designs in {elapsed:.2f} s")
print("best", designs[0].cds[:60])
print("score", round(designs[0].score, 4), "violations", len(designs[0].violations))
