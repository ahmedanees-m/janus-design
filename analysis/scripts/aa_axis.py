"""Does the amino-acid axis buy anything a codon optimiser cannot?

The 5' initiation term is reachable by synonymous change, and the decomposition
showed 88 percent of the cheap gain is exactly that. The terms in
`janus.objectives.liability` are not reachable that way at all: they are
functions of the residue sequence and, where structure matters, of a fixed
accessibility. A fixed-protein codon designer's achievable improvement on them is
zero by construction, and the point of this script is to measure that zero beside
a nonzero rather than to assert it.

Two arms search the same shell under the same objective and the same weights:

- joint, in which every residue the shell admits is available
- codon only, in which each position is pinned to the residue the Tier-1 optimum
  chose, so only the codons can move

Searching rather than rescoring a sampled pool matters here. Ranking a fixed pool
of uniform shell draws under a rising weight gives a step function: the winner
stays at the Tier-1 optimum until some distant random draw overtakes it, and the
apparent exchange rate is then a property of what happened to be sampled. A
greedy exchange starting from the Tier-1 optimum traces a path along the
exchange, one residue at a time, and every weight in the grid reads off a point
on that trace. Being greedy it is a lower bound on what the shell holds rather
than the true Pareto frontier, which is enough here: the comparison is against an
arm whose achievable removal is zero, and a lower bound above zero settles it.

Two accounting conventions, stated rather than left implicit. Tier-1 cost is in
nats of the decomposable objective, positive meaning given up. Budget spent is
the ProteinMPNN marginal log-probability surrendered at the changed positions,
over the most the shell could surrender if every position took its worst
admitted residue.
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
from janus.genetic_code import AA_ALPHABET
from janus.objectives import liability
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import pool_scales
from janus.sample import ShellSearch, shell_samples

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
GRID = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]
TARGETS = ["low_complexity", "protein_repeat", "exposed_hydrophobic",
           "exposed_hydrophobic_run", "degron"]
MAX_STEPS = 12
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
AA_INDEX = {a: i for i, a in enumerate(AA_ALPHABET)}
TOPOLOGY = re.compile(r"^(HHH|EHEE|EEHEE|HEEH)(?:_[A-Z]{2})?_rd\d+")
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def accessibility(path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    area = struc.apply_residue_wise(atoms, struc.sasa(atoms, vdw_radii="Single"), np.nansum)
    residues = [THREE_TO_ONE.get(n, "X") for n in struc.get_residues(atoms)[1]]
    return [None if r not in MAX_ASA or np.isnan(area[i]) else float(area[i] / MAX_ASA[r])
            for i, r in enumerate(residues)]


def family(name):
    match = TOPOLOGY.match(name)
    if match:
        return match.group(1)
    return "hallucination" if ":" in name else "other"


def budget_denominator(marginals, admitted):
    """Most the shell could surrender: every position at its worst admitted residue."""
    total = 0.0
    for position, residues in enumerate(admitted):
        values = [marginals[position][AA_INDEX[r]] for r in residues]
        total += max(values) - min(values)
    return max(total, EPS)


def spent(marginals, reference, protein):
    """Marginal log-probability given up at the positions that moved."""
    return float(sum(
        marginals[i][AA_INDEX[a]] - marginals[i][AA_INDEX[b]]
        for i, (a, b) in enumerate(zip(reference, protein, strict=True)) if a != b
    ))


def descend(search, reference, burden, admitted, max_steps=MAX_STEPS):
    """Greedy single-residue exchange from the Tier-1 optimum, best rate first.

    Returns the trace including the starting point. Each step takes the admitted
    substitution with the largest burden reduction per nat of Tier-1 score given
    up, so the trace covers the exchange rather than one weight's answer. Restricting
    `admitted` to the reference residue gives the codon-only arm, which therefore
    runs this same code and is not simulated.
    """
    base = burden(reference.protein)
    trace = [{"design": reference, "burden": base, "cost": 0.0, "step": 0}]
    current, current_burden = reference, base

    for step in range(1, max_steps + 1):
        # Burden is cheap to evaluate and codon re-optimisation is not, so filter
        # on burden first and re-optimise only the substitutions that help.
        proposals = []
        for position, residues in enumerate(admitted):
            for residue in residues:
                if residue == current.protein[position]:
                    continue
                candidate = current.protein[:position] + residue + current.protein[position + 1:]
                value = burden(candidate)
                if value < current_burden - EPS:
                    proposals.append((value, candidate))
        if not proposals:
            break

        best = None
        for value, candidate in proposals:
            scored = search.best(candidate)
            cost = reference.score - scored.score
            rate = (base - value) / max(cost, EPS)
            if best is None or rate > best[0]:
                best = (rate, scored, value, cost)
        _, current, current_burden, cost = best
        trace.append({"design": current, "burden": current_burden, "cost": cost, "step": step})

    return trace


def read_off(trace, weight, scale_tier1, scale_burden):
    """Point on the trace maximising the penalised objective `rescore` uses."""
    base = trace[0]["burden"]
    best, value = trace[0], 0.0
    for point in trace[1:]:
        gain = -point["cost"] / scale_tier1 + weight * (base - point["burden"]) / scale_burden
        if gain > value + EPS:
            best, value = point, gain
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--pool", type=int, default=200,
                        help="shell draws used only to measure term spreads")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

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

        search = ShellSearch(marginals, host, TIER1, args.delta)
        reference = design(marginals, host, weights=TIER1, delta=args.delta, k=1)[0]
        pinned = [(r,) for r in reference.protein]
        denominator = budget_denominator(marginals, search.admitted)
        pool = shell_samples(marginals, host, weights=TIER1, delta=args.delta,
                             count=args.pool, rng=rng)

        for target in TARGETS:
            cache: dict[str, float] = {}

            def burden(protein, target=target, cache=cache):
                if protein not in cache:
                    cache[protein] = liability.score(target, protein, rsa, classes)
                return cache[protein]

            scales = pool_scales(pool, host, liability=burden)
            base = burden(reference.protein)
            traces = {
                "joint": descend(search, reference, burden, search.admitted),
                "codon_only": descend(search, reference, burden, pinned),
            }

            for arm, trace in traces.items():
                for weight in GRID:
                    point = read_off(trace, weight, scales.tier1, scales.liability)
                    protein = point["design"].protein
                    rows.append({
                        "backbone": name,
                        "family": family(name),
                        "length": len(marginals),
                        "target": target,
                        "arm": arm,
                        "weight": weight,
                        "reference_liability": base,
                        "liability": point["burden"],
                        "removed": base - point["burden"],
                        "tier1_cost": point["cost"],
                        "residue_changes": point["step"],
                        "codon_changes": sum(
                            a != b for a, b in zip(point["design"].cds, reference.cds, strict=True)
                        ),
                        "budget_spent": spent(marginals, reference.protein, protein) / denominator,
                        "identity": 1.0 - point["step"] / len(reference.protein),
                        "frontier_reach": base - trace[-1]["burden"],
                        "frontier_cost": trace[-1]["cost"],
                        "pool_spread": scales.liability,
                    })
        if number % 20 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows)
    print(f"\nwrote {args.out} ({len(rows)} rows)")


def summarise(rows):
    backbones = {r["backbone"] for r in rows}
    print(f"\n=== amino-acid axis over {len(backbones)} designed backbones ===")

    for target in TARGETS:
        subset = [r for r in rows if r["target"] == target and r["arm"] == "joint"]
        if not subset:
            continue
        base = {r["backbone"]: r["reference_liability"] for r in subset}
        reach = {r["backbone"]: r["frontier_reach"] for r in subset}
        carriers = [b for b, v in base.items() if v > 0]
        movable = [b for b in carriers if reach[b] > EPS]
        print(f"\n--- {target} ---")
        print(f"  carry a burden at the Tier-1 optimum: {len(carriers)}/{len(base)} "
              f"({100 * len(carriers) / max(len(base), 1):.0f}%); "
              f"of those, removable at all: {len(movable)}/{max(len(carriers), 1)} "
              f"({100 * len(movable) / max(len(carriers), 1):.0f}%)")

        header = ("weight".rjust(8) + "removed".rjust(11) + "codon arm".rjust(11)
                  + "nats".rjust(9) + "residues".rjust(10) + "budget".rjust(9)
                  + "moved".rjust(8))
        print(header)
        for weight in GRID:
            joint = [r for r in subset if r["weight"] == weight]
            codon = [r for r in rows if r["target"] == target and r["arm"] == "codon_only"
                     and r["weight"] == weight]
            print(f"{weight:>8.3f}{np.mean([r['removed'] for r in joint]):>11.4f}"
                  f"{np.mean([r['removed'] for r in codon]):>11.4f}"
                  f"{np.mean([r['tier1_cost'] for r in joint]):>9.3f}"
                  f"{np.mean([r['residue_changes'] for r in joint]):>10.2f}"
                  f"{100 * np.mean([r['budget_spent'] for r in joint]):>8.1f}%"
                  f"{100 * np.mean([r['residue_changes'] > 0 for r in joint]):>7.0f}%")

        # The tail the mean hides: what the exchange costs on designs that carry
        # something worth removing.
        loaded = [r for r in subset if r["weight"] == 1.0 and r["reference_liability"] > 0
                  and r["frontier_reach"] > EPS]
        if loaded:
            fraction = np.array([r["frontier_reach"] / r["reference_liability"] for r in loaded])
            cost = np.array([r["frontier_cost"] for r in loaded])
            print(f"  on the {len(loaded)} removable carriers, exhausting the frontier clears "
                  f"median {100 * np.median(fraction):.0f}% of the burden "
                  f"(10th {100 * np.percentile(fraction, 10):.0f}%, "
                  f"90th {100 * np.percentile(fraction, 90):.0f}%) "
                  f"for median {np.median(cost):.2f} nats "
                  f"(90th {np.percentile(cost, 90):.2f})")

    print("\n=== what the codon arm reached, every term and every weight ===")
    for target in TARGETS:
        codon = [r["removed"] for r in rows
                 if r["target"] == target and r["arm"] == "codon_only"]
        joint = [r["removed"] for r in rows
                 if r["target"] == target and r["arm"] == "joint" and r["weight"] == 1.0]
        if not codon:
            continue
        print(f"  {target:<26} codon-only largest removal {max(codon):.3e} over "
              f"{len(codon)} runs; joint mean at weight 1 {np.mean(joint):.4f}")


if __name__ == "__main__":
    main()
