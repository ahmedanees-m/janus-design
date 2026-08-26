"""What the amino-acid freedom is worth, as a function of how much of it there is.

Delta zero collapses the shell to the per-position marginal argmax, so the
lattice degenerates to a fixed-protein codon DP and the search has no residue
move available to it. Every delta above zero is the joint lattice. Holding the
whole objective fixed and sweeping delta therefore prices the amino-acid axis
directly, and decomposing the attained objective by term shows where the freedom
is being spent rather than only how much.

The objective here is the complete one. Tier 1 (MPNN marginal, CAI, codon-pair
bias) enters the parser and is exact. Initiation-window folding and the five
protein-level liabilities do not decompose over lattice nodes and are optimised
by hill climbing from the Tier-1 optimum of the same shell, one residue at a
time. Synthesis and cis-element constraints are checked rather than weighted, so
they are reported as violation counts at the attained point.

Every term is divided by its spread over one common candidate pool, drawn once
per backbone from the widest shell in the sweep. Sharing the scales across the
delta arms is what makes them comparable; recomputing them per arm would rescale
the axis being measured.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np

from janus import Weights, design, hosts
from janus.objectives import liability, mrna
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.sample import ShellSearch, shell_samples

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
UNBOUNDED = 1e9
DELTAS = [0.0, 0.5, 1.0, 2.0, 3.0, UNBOUNDED]
TARGETS = ["low_complexity", "protein_repeat", "exposed_hydrophobic",
           "exposed_hydrophobic_run", "degron"]
MAX_STEPS = 10
EPS = 1e-9

MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def accessibility(path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    area = struc.apply_residue_wise(atoms, struc.sasa(atoms, vdw_radii="Single"), np.nansum)
    residues = [THREE_TO_ONE.get(n, "X") for n in struc.get_residues(atoms)[1]]
    return [None if r not in MAX_ASA or np.isnan(area[i]) else float(area[i] / MAX_ASA[r])
            for i, r in enumerate(residues)]


def spread(values):
    """Spread of a term over the pool, or infinity when it is constant there.

    Infinity divides the term out, which is what a term that cannot order the
    candidates should contribute. Flooring a zero spread at some small number
    would instead let it dominate everything else.
    """
    value = float(np.std(values))
    return value if value > EPS else float("inf")


class Objective:
    """The full spread-normalised objective, and its decomposition."""

    def __init__(self, host, rsa, classes, scales, folding, penalty):
        self.host, self.rsa, self.classes = host, rsa, classes
        self.scales, self.folding, self.penalty = scales, folding, penalty
        self.cache: dict[str, dict] = {}

    def parts(self, candidate) -> dict:
        if candidate.cds not in self.cache:
            burdens = {name: liability.score(name, candidate.protein, self.rsa, self.classes)
                       for name in TARGETS}
            self.cache[candidate.cds] = {
                "tier1": candidate.score,
                "initiation": mrna.initiation_energy(candidate.cds, self.host),
                **burdens,
            }
        return self.cache[candidate.cds]

    def total(self, candidate) -> float:
        parts = self.parts(candidate)
        value = parts["tier1"] / self.scales["tier1"]
        value += self.folding * parts["initiation"] / self.scales["initiation"]
        for name in TARGETS:
            value -= self.penalty * parts[name] / self.scales[name]
        return value


def climb(search, start, objective, max_steps=MAX_STEPS):
    """Hill climb on the full objective by single admitted residue substitutions."""
    current, value = start, objective.total(start)
    steps = 0
    for _ in range(max_steps):
        best, best_value = None, value
        for position, residues in enumerate(search.admitted):
            for residue in residues:
                if residue == current.protein[position]:
                    continue
                candidate = search.best(
                    current.protein[:position] + residue + current.protein[position + 1:]
                )
                candidate_value = objective.total(candidate)
                if candidate_value > best_value + EPS:
                    best, best_value = candidate, candidate_value
        if best is None:
            break
        current, value, steps = best, best_value, steps + 1
    return current, value, steps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--folding", type=float, default=0.25,
                        help="initiation weight on the spread scale")
    parser.add_argument("--penalty", type=float, default=0.25,
                        help="weight on each liability, same scale")
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not mrna.available():
        raise SystemExit("ViennaRNA is required for this experiment")

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    rng = np.random.default_rng(args.seed)
    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    print(f"{len(files)} designed backbones")
    if args.limit:
        files = files[: args.limit]

    rows = []
    for number, path in enumerate(files, start=1):
        name = path.stem
        marginals = load_unconditional(path)
        try:
            rsa = accessibility(Path(args.backbones) / f"{name}.pdb")
        except Exception:
            continue
        if len(rsa) != len(marginals):
            continue

        # One pool from the widest shell fixes the scales for every arm.
        pool = shell_samples(marginals, host, weights=TIER1, delta=max(DELTAS),
                             count=args.pool, rng=rng)
        scales = {"tier1": spread([d.score for d in pool]),
                  "initiation": spread([mrna.initiation_energy(d.cds, host) for d in pool])}
        for term in TARGETS:
            scales[term] = spread([liability.score(term, d.protein, rsa, classes)
                                   for d in pool])

        objective = Objective(host, rsa, classes, scales, args.folding, args.penalty)
        baseline = None
        for delta in DELTAS:
            search = ShellSearch(marginals, host, TIER1, delta)
            start = design(marginals, host, weights=TIER1, delta=delta, k=1)[0]
            attained, value, steps = climb(search, start, objective)
            parts = objective.parts(attained)
            if baseline is None:
                baseline = value
            rows.append({
                "backbone": name,
                "length": len(marginals),
                "delta": None if delta >= UNBOUNDED else delta,
                "shell_mean": float(np.mean([len(r) for r in search.admitted])),
                "total": value,
                "gain_over_fixed_protein": value - baseline,
                "residue_changes": steps,
                "identity_to_start": 1.0 - steps / len(marginals),
                "violations": len(attained.violations),
                **{f"term_{k}": v for k, v in parts.items()},
            })
        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows, args.folding, args.penalty)
    print(f"\nwrote {out} ({len(rows)} rows)")


def summarise(rows, folding, penalty):
    backbones = {r["backbone"] for r in rows}
    print(f"\n=== delta sweep over {len(backbones)} designed backbones "
          f"(folding {folding}, liability {penalty}) ===")
    print("attained objective is spread-normalised, so it is unitless; the gain "
          "column is what the amino-acid freedom adds over the fixed-protein DP")

    labels = [None if d >= UNBOUNDED else d for d in DELTAS]
    header = ("delta".rjust(9) + "shell".rjust(8) + "gain".rjust(9) + "10th".rjust(8)
              + "90th".rjust(8) + "residues".rjust(10) + "wins".rjust(8)
              + "violations".rjust(12))
    print(header)
    for delta in labels:
        subset = [r for r in rows if r["delta"] == delta]
        if not subset:
            continue
        gain = np.array([r["gain_over_fixed_protein"] for r in subset])
        print(f"{('inf' if delta is None else f'{delta:.1f}'):>9}"
              f"{np.mean([r['shell_mean'] for r in subset]):>8.2f}"
              f"{np.median(gain):>9.4f}{np.percentile(gain, 10):>8.4f}"
              f"{np.percentile(gain, 90):>8.4f}"
              f"{np.mean([r['residue_changes'] for r in subset]):>10.2f}"
              f"{100 * np.mean(gain > EPS):>7.0f}%"
              f"{np.mean([r['violations'] for r in subset]):>12.2f}")

    print("\nwhere the freedom is spent, as change from the fixed-protein arm")
    fixed = {r["backbone"]: r for r in rows if r["delta"] == 0.0}
    terms = ["tier1", "initiation"] + TARGETS
    print("delta".rjust(9) + "".join(t[:11].rjust(13) for t in terms))
    for delta in labels[1:]:
        subset = [r for r in rows if r["delta"] == delta and r["backbone"] in fixed]
        if not subset:
            continue
        line = ("inf" if delta is None else f"{delta:.1f}").rjust(9)
        for term in terms:
            change = [r[f"term_{term}"] - fixed[r["backbone"]][f"term_{term}"] for r in subset]
            line += f"{np.median(change):>+13.4f}"
        print(line)
    print("tier1 and initiation are in native units, nats and kcal/mol; the "
          "liabilities are in their own counts and fractions")


if __name__ == "__main__":
    main()
