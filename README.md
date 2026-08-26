# janus-design

[![tests](https://github.com/ahmedanees-m/janus-design/actions/workflows/tests.yml/badge.svg)](https://github.com/ahmedanees-m/janus-design/actions/workflows/tests.yml)
[![python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![licence](https://img.shields.io/badge/licence-MIT-blue)](LICENSE)

Joint optimisation of amino-acid and coding sequence for de novo designed proteins.

JANUS takes an inverse-folding posterior for a backbone and a host organism and
returns a coding sequence, choosing the protein sequence and the codons together
instead of one after the other.

A natural protein has a native coding sequence to fall back on. A de novo
designed protein does not, so its gene is built from scratch, usually by a codon
optimiser that treats the amino-acid sequence as fixed. But an inverse-folding
model does not produce a sequence; it produces a distribution, and at most
positions on a designed backbone several residues are nearly equally compatible
with the fold. JANUS keeps that degeneracy in the search space, which lets the
gene-design step reach properties that codon choice cannot touch: degron motifs,
low-complexity content, repeats and exposed hydrophobic surface are functions of
the residue sequence, and no synonymous substitution moves them.

## Installation

```bash
git clone https://github.com/ahmedanees-m/janus-design
cd janus-design
pip install -e .
```

The solver needs only numpy and pyyaml. The objective terms that read structure
and fold RNA need more, and come as an extra:

```bash
pip install -e ".[analysis]"
```

ProteinMPNN supplies the marginals and is run separately; LinearDesign is used by
one cross-implementation test. Both are optional, and the tests that need them
skip when they are absent.

A container definition is included:

```bash
docker build -t janus .
```

## Usage

JANUS scores residues by their inverse-folding marginals, so it starts from a
ProteinMPNN unconditional-probability archive rather than from the PDB file
directly:

```bash
python ProteinMPNN/protein_mpnn_run.py \
    --pdb_path backbone.pdb \
    --unconditional_probs_only 1 \
    --out_folder marginals/
```

That writes a `.npz` holding one probability vector per position. Pass it to
`janus design`:

```bash
janus design --marginals marginals/backbone.npz --host ecoli_bl21 --fasta gene.fasta
```

The command builds the lattice, parses it, and writes the best path as FASTA
with its score, per-term breakdown and synthesis-constraint status on the header
line. Without `--fasta` it writes to standard output.

To see alternatives rather than one answer, ask for more paths:

```bash
janus design --marginals marginals/backbone.npz --host ecoli_bl21 --k 50
```

`--delta` sets the entropy budget in nats, which controls how far from the
marginal optimum the search may wander. Half a nat captures most of what is
available and one nat captures nearly all of it. Zero pins the protein and
reduces JANUS to a codon optimiser, which is useful as a control:

```bash
janus design --marginals marginals/backbone.npz --host hek293 --delta 0.5
```

Term weights are set with `--lambda-mpnn`, `--lambda-cai`, `--lambda-cpb` and
`--lambda-gc`. Run `janus design --help` for the full set.

The Python interface carries the options the command line does not, including
pinning individual positions:

```python
import numpy as np
from janus import design, load_host
from janus.objectives.mpnn import load_unconditional

host = load_host("ecoli_bl21")
marginals = load_unconditional("marginals/backbone.npz")

results = design(marginals, host, delta=1.0, k=50, fixed={11: "H", 12: "E"})
best = results[0]
print(best.protein, best.cds, best.score, best.synthesisable)
```

`fixed` holds positions the design cannot afford to move: an active site, a
binding interface, a tag. `anchor` additionally admits a named residue sequence
wherever the shell does not already contain it, which is what lets a search start
from a given protein and then spend the budget from there.

For the non-decomposable terms, `janus.search.optimise` runs coordinate descent
over the shell and re-solves the codon layer exactly after every move.

## Objective

| Layer | Terms |
|---|---|
| Fold | inverse-folding marginal log-probability |
| Codon | codon adaptation, codon-pair bias, GC content |
| Transcript | translation-initiation window, global folding energy |
| Protein | degron load, low-complexity content, repeats, exposed hydrophobic surface |
| Constraint | restriction sites, homopolymers, repeats, local GC windows |

Weights are set per host and can be overridden. Terms are normalised by their
spread over the candidate pool, so a weight expresses a preference and not a unit
conversion, and the same weights transfer between backbones of different length.

Synthesis constraints are checked, not weighted. A six-base restriction site
spans up to three codons and a homopolymer limit of nine spans four, so neither
decomposes over adjacent codon pairs and neither can live inside the parse. They
are applied as a filter over the ranked list, which is why the interface returns
a list: about half of any given pool passes every constraint, so a method
returning a single sequence would violate one about half the time.

## Hosts

A host is a YAML policy file plus codon and codon-pair count tables. Two ship
with the package:

| | source | initiation model |
|---|---|---|
| `ecoli_bl21` | RefSeq GCF_000022665.1 | Shine-Dalgarno accessibility |
| `hek293` | Ensembl GRCh38 | cap-dependent scanning |

The two carry separate initiation models instead of one model with a sign
parameter. Initiation in bacteria depends on keeping the Shine-Dalgarno region
accessible; in eukaryotes there is no Shine-Dalgarno and initiation is by
cap-dependent scanning, with global transcript structure acting in the opposite
direction. Only the prokaryotic term rests on large synonymous libraries.

Adding an organism does not require touching the solver:

```bash
janus build-host --cds my_organism_cds.fasta.gz --host my_organism
```

That counts codons and codon pairs from the reference coding sequences and
rewrites the count tables named in the host YAML.

## How it works

The synonymous space of a fixed protein can be written as a deterministic finite
automaton over codons and searched exactly by lattice parsing. JANUS generalises
that automaton so each position admits the codons of every amino acid within a
log-probability tolerance of the marginal optimum, with each residue branch
carrying its marginal as an additive node weight. The automaton is then
amino-acid-degenerate as well as codon-degenerate, and parsing returns the best
(protein, gene) pair jointly. A 43-residue lattice solves in about 1.6 ms.

Optimisation runs in two tiers. The first parses the lattice exactly over the
terms that decompose over nodes and adjacent-codon edges and returns the *k* best
paths. The second scores those candidates under the terms that do not decompose,
being mRNA folding, the protein-level liability panel and the full autoregressive
inverse-folding likelihood, and refines within the shell.

ProteinMPNN decodes autoregressively in random order, so its conditional
per-position log-probability depends on residues chosen elsewhere and cannot
serve as an additive node weight. JANUS therefore uses the unconditional
single-pass marginals, which are additive by construction, and the parse is exact
with respect to those. It is not exact with respect to the conditional posterior.
The two differ by about 2.4 nats RMS, several times the width of the band the
parser's top candidates occupy, so the rescoring tier changes the answer.

Two search strategies are available for the residue layer. Coordinate descent is
the default because it reaches within about one percent of best-improvement
search using roughly a seventh of the evaluations.

## Tests

```bash
pytest
```

72 tests, 2 of which skip when ViennaRNA or LinearDesign is unavailable. On a
solver-only install a further 21 skip, since they need the analysis extra. Two
check the layers in isolation: an open shell scored on the marginal term alone
must return the marginal argmax, and a zero entropy budget must collapse the
lattice to a fixed-protein codon automaton whose optimum matches exhaustive
enumeration on short lattices. In the CAI limit the codon layer reproduces
LinearDesign path for path.

## Data

Derived data, run outputs and figure source data are deposited separately on
Zenodo. The deposit includes the inverse-folding marginals for all 862 backbones
used in the paper, so the analyses can be reproduced without a GPU.

## Citation

```
Mahaboob Ali AA, Delhibabu R, Nelson EJR.
Joint optimisation of amino-acid and coding sequence for de novo designed
proteins. (submitted)
```

`CITATION.cff` carries the machine-readable form. A version DOI is minted by
Zenodo on each tagged release and added there.

## Licence

MIT, in `LICENSE`. Host codon tables derived from RefSeq and Ensembl are
included. Motif definitions are not: they are fetched from ELM, whose terms of
use restrict redistribution and require a licence for non-academic use.
