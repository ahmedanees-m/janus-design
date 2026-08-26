# Runtime for the solver and the analysis scripts. The image carries pinned
# dependencies only; the source tree is mounted at /work so a change does not
# need a rebuild.
#
#   docker build -t janus:0.3.0 .
#   docker run --rm --user "$(id -u):$(id -g)" \
#       -v "$PWD:/work" -v "/path/to/janus-data:/data" janus:0.3.0 pytest -q
#
# Pass --user on every run that writes to a mounted volume. Without it the
# container writes as root and the results come back undeletable from the host.

FROM python:3.11-slim

RUN pip install --no-cache-dir \
        numpy==2.1.3 \
        pyyaml==6.0.2 \
        pyarrow==18.1.0 \
        pandas==2.2.3 \
        scipy==1.14.1 \
        biopython==1.84 \
        biotite==1.6.0 \
        ViennaRNA==2.7.0 \
        matplotlib==3.9.2 \
        scikit-learn==1.5.2 \
        pytest==8.3.4

ENV PYTHONPATH=/work/src \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /work
