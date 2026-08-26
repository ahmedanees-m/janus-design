# CodonTransformer, kept in its own image because its pins conflict with ours.
#
# numpy is held below 2 and transformers below 4.50: BigBird's generate() path
# regressed in 4.50 and CodonTransformer depends on it. Torch is the CPU build,
# which is enough for a few hundred sequences of under a hundred residues and
# avoids carrying a CUDA runtime for a baseline.
FROM python:3.11-slim

RUN pip install --no-cache-dir \
        torch==2.2.2 --index-url https://download.pytorch.org/whl/cpu \
 && pip install --no-cache-dir \
        "numpy<2" \
        "transformers<4.50" \
        CodonTransformer

WORKDIR /w
