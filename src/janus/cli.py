"""Command line interface."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

from . import hosts
from .design import design
from .objectives import Weights
from .objectives.mpnn import load_unconditional


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="janus", description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)

    run = subcommands.add_parser("design", help="design coding sequences for one backbone")
    run.add_argument("--marginals", required=True, help="ProteinMPNN unconditional probability .npz")
    run.add_argument("--host", default="ecoli_bl21", help="shipped host name or path to a host YAML")
    run.add_argument("--delta", type=float, default=1.0, help="entropy budget in nats")
    run.add_argument("--k", type=int, default=1, help="number of paths to return")
    run.add_argument("--lambda-mpnn", type=float, default=1.0)
    run.add_argument("--lambda-cai", type=float, default=0.0)
    run.add_argument("--lambda-cpb", type=float, default=0.0)
    run.add_argument("--lambda-gc", type=float, default=0.0)
    run.add_argument("--fasta", help="write the coding sequences here instead of standard output")
    run.set_defaults(handler=_design)

    build = subcommands.add_parser("build-host", help="count codon usage for a host from a CDS set")
    build.add_argument("--cds", required=True, help="reference coding sequences, FASTA or FASTA.gz")
    build.add_argument("--host", required=True, help="host YAML whose count tables to regenerate")
    build.set_defaults(handler=_build_host)

    args = parser.parse_args(argv)
    return args.handler(args)


def _design(args: argparse.Namespace) -> int:
    host = hosts.load(args.host)
    marginals = load_unconditional(args.marginals)
    weights = Weights(
        mpnn=args.lambda_mpnn,
        cai=args.lambda_cai,
        cpb=args.lambda_cpb,
        gc=args.lambda_gc,
    )

    results = design(marginals, host, weights=weights, delta=args.delta, k=args.k)
    if not results:
        print("no path through the lattice", file=sys.stderr)
        return 1

    lines = []
    for rank, result in enumerate(results, start=1):
        terms = " ".join(f"{name}={value:.4f}" for name, value in result.terms.items())
        flags = "ok" if result.synthesisable else f"{len(result.violations)} violations"
        lines.append(f">rank{rank} score={result.score:.4f} {terms} synthesis={flags}")
        lines.append(result.cds)

    text = "\n".join(lines) + "\n"
    if args.fasta:
        Path(args.fasta).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


def _build_host(args: argparse.Namespace) -> int:
    host_path = Path(args.host)
    if not host_path.suffix:
        host_path = hosts.HOST_DIR / f"{host_path}.yaml"

    import yaml

    spec = yaml.safe_load(host_path.read_text(encoding="utf-8"))
    records = list(_read_fasta(Path(args.cds)))
    codon_counts, pair_counts = hosts.count_reference_cds(records)

    if not codon_counts:
        print(f"no usable coding sequences in {args.cds}", file=sys.stderr)
        return 1

    hosts.write_counts(
        host_path.parent / spec["codon_counts"],
        sorted(codon_counts.items(), key=lambda kv: (-kv[1], kv[0])),
    )
    hosts.write_counts(
        host_path.parent / spec["codon_pair_counts"],
        ((first, second, n) for (first, second), n in sorted(pair_counts.items())),
    )

    print(
        f"{len(records)} records read, {sum(codon_counts.values())} codons and "
        f"{len(pair_counts)} distinct pairs counted for {spec['name']}"
    )
    return 0


def _read_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if chunks:
                    yield "".join(chunks)
                    chunks = []
            else:
                chunks.append(line.strip())
        if chunks:
            yield "".join(chunks)


if __name__ == "__main__":
    raise SystemExit(main())
