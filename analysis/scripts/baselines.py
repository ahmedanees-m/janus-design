"""Every codon and joint design method, scored by one piece of code.

Six arms over the same backbones and the same objective:

- **vendor**, the design's residue sequence with the highest
  relative-adaptiveness codon at every position
- **CodonTransformer**, the same residue sequence, codons from the incumbent
  learned codon optimiser, read from a file its own container produced
- **rejection**, the best of N draws from the shell scored by the full objective,
  which is what a method without a lattice would have to do
- **codon DP**, the same residue sequence with the JANUS codon layer, pinned
- **CodonMPNN**, the published codon-level inverse folding model, which chooses
  residues and codons together from the backbone and so belongs with the joint
  arms rather than the codon ones
- **ProteinMPNN then codons**, residues sampled from the model whose marginals
  define the lattice, codons then solved exactly for them
- **SolubleMPNN then codons**, the same with the solubility-trained checkpoint
- **MoMPNN then codons**, the same with the developability-aligned checkpoint
- **coordinate descent**, a positionwise greedy pass that accepts the best
  improving substitution at each position before moving on, rather than the best
  move over the whole chain, with the codon layer re-solved exactly after each
- **joint**, the full search, taking the best move over the whole chain each step

The first three hold the protein fixed and differ only in how they choose codons.
The last four may move residues. Budgets are matched two ways and both are
reported: wall-clock, and the number of full-objective evaluations, since the
methods differ by orders of magnitude in cost per evaluation.

Two things this comparison is not.

CodonTransformer does not optimise this objective. It predicts the codons a host
would use, learned from natural genes, so scoring it here and finding it lower is
close to circular. Its per-term columns are reported for that reason: what is
informative is where it differs, not that it loses on a total it never saw.

The coordinate-descent arm is not the conditional-ProteinMPNN alternation that
the one comparable prior extension of this lattice used; both local searches below
move over the same neighbourhood and differ only in acceptance order.

Block alternation itself is included, in the only form this objective admits. A
block scheme alternates optimising one block given the other, which needs the
residue step to see what the codon step is optimising. An inverse folding model
has no channel for that: it scores residues against a backbone and knows nothing
of initiation windows, codon pairs or degron load. So alternation here collapses
to one round, propose residues with the model and then solve the codon layer
exactly for them, and running two arms of exactly that shape says what the
alternative is worth. `proteinmpnn` proposes with the same model whose marginals
define the lattice; `mompnn` proposes with the developability-aligned variant,
which is the closest published thing to optimising these liabilities directly.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from dataclasses import replace
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np

from janus import Weights, design, evaluate, hosts
from janus.genetic_code import SYNONYMOUS
from janus.objectives import liability, mrna
from janus.objectives.mpnn import load_unconditional
from janus.objectives.proteostasis import load_classes
from janus.rescore import FoldingWeights, pool_scales, rescore
from janus.sample import ShellSearch, shell_samples

TIER1 = Weights(mpnn=1.0, cai=0.5, cpb=0.3)
LIABILITIES = ["low_complexity", "protein_repeat", "exposed_hydrophobic",
               "exposed_hydrophobic_run", "degron"]
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
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


def read_backbone(path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    area = struc.apply_residue_wise(atoms, struc.sasa(atoms, vdw_radii="Single"), np.nansum)
    residues = [THREE_TO_ONE.get(name, "X") for name in struc.get_residues(atoms)[1]]
    rsa = [None if r not in MAX_ASA or np.isnan(area[i]) else float(area[i] / MAX_ASA[r])
           for i, r in enumerate(residues)]
    return "".join(residues), rsa


def vendor_cds(protein, host):
    adaptiveness = host.relative_adaptiveness
    best = {aa: max(codons, key=lambda c: adaptiveness.get(c, 0.0))
            for aa, codons in SYNONYMOUS.items()}
    return "".join(best[r] for r in protein if r in best)


def make_objective(marginals, host, rsa, classes, delta, native, pool_size, weights, rng):
    """One scoring function, and the pool its scales are measured on."""
    pool = design(marginals, host, weights=TIER1, delta=delta, k=pool_size, anchor=native)
    pool += shell_samples(marginals, host, weights=TIER1, delta=delta,
                          count=pool_size, rng=rng)

    spreads = {}
    for name in LIABILITIES:
        values = [liability.score(name, d.protein, rsa, classes) for d in pool]
        spread = float(np.std(values))
        spreads[name] = spread if spread > 1e-9 else float("inf")

    def burden(protein):
        return sum(liability.score(name, protein, rsa, classes) / spreads[name]
                   for name in LIABILITIES)

    scales = replace(pool_scales(pool, host), liability=1.0)

    def total(candidate):
        return rescore([candidate], host, weights, scales, liability=burden)[0].total

    return total, burden


def coordinate_descent(search, start, total, marginals, rounds):
    """Greedy positionwise descent, codons re-solved exactly after every move.

    Positions are visited in descending marginal entropy and the best admitted
    substitution at each is taken before moving on, rather than the best move
    over the whole chain being found first. Cheaper than the joint search by the
    evaluations it skips, and the question is whether that costs anything in the
    objective attained.
    """
    order = np.argsort(-np.array([
        -(np.exp(row) * row).sum() for row in marginals
    ]))
    current, value = start, total(start)
    evaluations = 1
    for _ in range(rounds):
        improved = False
        for position in order:
            best, best_value = None, value
            for residue in search.admitted[position]:
                if residue == current.protein[position]:
                    continue
                candidate = search.best(
                    current.protein[:position] + residue + current.protein[position + 1:]
                )
                evaluations += 1
                candidate_value = total(candidate)
                if candidate_value > best_value + 1e-12:
                    best, best_value = candidate, candidate_value
            if best is not None:
                current, value, improved = best, best_value, True
        if not improved:
            break
    return current, value, evaluations


def hill_climb(search, start, total, max_steps):
    current, value = start, total(start)
    evaluations, steps = 1, 0
    while steps < max_steps:
        best, best_value = None, value
        for position, residues in enumerate(search.admitted):
            for residue in residues:
                if residue == current.protein[position]:
                    continue
                candidate = search.best(
                    current.protein[:position] + residue + current.protein[position + 1:]
                )
                evaluations += 1
                candidate_value = total(candidate)
                if candidate_value > best_value + 1e-12:
                    best, best_value = candidate, candidate_value
        if best is None:
            break
        current, value, steps = best, best_value, steps + 1
    return current, value, evaluations, steps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--marginals", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--codontransformer", help="JSON of coding sequences from its container")
    parser.add_argument("--codonmpnn", help="JSON of coding sequences from its container")
    parser.add_argument("--proposals", nargs="*", default=[],
                        help="arm=path pairs of residue-only proposals, codons solved here")
    parser.add_argument("--host", default="ecoli_bl21")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--pool", type=int, default=200)
    parser.add_argument("--rejection", type=int, default=2000,
                        help="draws for the rejection-sampling arm")
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--max-steps", type=int, default=40)
    parser.add_argument("--folding", type=float, default=0.25)
    parser.add_argument("--penalty", type=float, default=0.25)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    host = hosts.load(args.host)
    classes = load_classes(args.elm)
    rng = np.random.default_rng(args.seed)
    weights = FoldingWeights(initiation=args.folding, liability=args.penalty)
    learned = {}
    if args.codontransformer:
        learned = json.loads(Path(args.codontransformer).read_text(encoding="utf-8"))
    inverse = {}
    if args.codonmpnn:
        inverse = json.loads(Path(args.codonmpnn).read_text(encoding="utf-8"))
    proposals = {}
    for item in args.proposals:
        arm, _, path = item.partition("=")
        proposals[arm] = json.loads(Path(path).read_text(encoding="utf-8"))

    files = [p for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        files = files[: args.limit]
    print(f"{len(files)} designed backbones, {len(learned)} CodonTransformer and "
          f"{len(inverse)} CodonMPNN sequences")

    rows = []
    for number, path in enumerate(files, start=1):
        name = path.stem
        marginals = load_unconditional(path)
        try:
            native, rsa = read_backbone(Path(args.backbones) / f"{name}.pdb")
        except Exception:
            continue
        if len(rsa) != len(marginals) or "X" in native:
            continue

        total, burden = make_objective(marginals, host, rsa, classes, args.delta,
                                       native, args.pool, weights, rng)
        search = ShellSearch(marginals, host, TIER1, args.delta, anchor=native)

        def record(arm, candidate, seconds, evaluations):
            rows.append({
                "backbone": name, "arm": arm, "length": len(marginals),
                "objective": total(candidate), "tier1": candidate.score,
                "burden": burden(candidate.protein),
                "violations": len(candidate.violations),
                "seconds": seconds, "evaluations": evaluations,
                "cai": candidate.terms.get("cai"), "cpb": candidate.terms.get("cpb"),
                "gc": candidate.terms.get("gc"),
                "initiation": mrna.initiation_energy(candidate.cds, host),
                "identity_to_native": sum(
                    a == b for a, b in zip(candidate.protein, native, strict=True)
                ) / len(native),
            })

        started = time.perf_counter()
        pinned = design(marginals, host, weights=TIER1, delta=args.delta, k=1,
                        fixed=dict(enumerate(native)), anchor=native)[0]
        codon_seconds = time.perf_counter() - started

        started = time.perf_counter()
        record("vendor", evaluate(vendor_cds(native, host), host, marginals, TIER1),
               time.perf_counter() - started, 1)

        if name in learned:
            try:
                scored = evaluate(learned[name]["cds"], host, marginals, TIER1)
            except ValueError:
                scored = None
            if scored is not None and scored.protein == native:
                record("codontransformer", scored, learned[name]["seconds"], 1)

        if name in inverse:
            try:
                scored = evaluate(inverse[name]["cds"], host, marginals, TIER1)
            except ValueError:
                scored = None
            if scored is not None:
                record("codonmpnn", scored, inverse[name]["seconds"], 1)

        for arm, table in proposals.items():
            entry = table.get(name)
            if entry is None or len(entry["protein"]) != len(native):
                continue
            if any(r not in search.slots[i] for i, r in enumerate(entry["protein"])):
                # Outside the shell, so the lattice cannot hold it; solve its
                # codons on a lattice pinned to that protein instead.
                started = time.perf_counter()
                solved = design(marginals, host, weights=TIER1, delta=args.delta, k=1,
                                fixed=dict(enumerate(entry["protein"])),
                                anchor=entry["protein"])[0]
                seconds = entry.get("seconds", 0.0) + time.perf_counter() - started
            else:
                started = time.perf_counter()
                solved = search.best(entry["protein"])
                seconds = entry.get("seconds", 0.0) + time.perf_counter() - started
            record(arm, solved, seconds, 1)

        record("codon_dp", pinned, codon_seconds, 1)

        started = time.perf_counter()
        drawn = shell_samples(marginals, host, weights=TIER1, delta=args.delta,
                              count=args.rejection, rng=rng)
        best = max(drawn, key=total)
        record("rejection", best, time.perf_counter() - started, args.rejection)

        started = time.perf_counter()
        descended, _, evaluations = coordinate_descent(
            search, pinned, total, marginals, args.rounds)
        record("coordinate_descent", descended, time.perf_counter() - started, evaluations)

        started = time.perf_counter()
        joint, _, evaluations, _ = hill_climb(search, pinned, total, args.max_steps)
        record("joint", joint, time.perf_counter() - started, evaluations)

        if number % 10 == 0:
            print(f"  {number}/{len(files)} backbones", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows), encoding="utf-8")
    summarise(rows)
    print(f"\nwrote {out} ({len(rows)} rows)")


def summarise(rows):
    arms = ["vendor", "codontransformer", "codon_dp", "codonmpnn",
            "proteinmpnn", "solublempnn", "mompnn", "rejection",
            "coordinate_descent", "joint"]
    present = [a for a in arms if any(r["arm"] == a for r in rows)]
    backbones = {r["backbone"] for r in rows}
    print(f"\n=== {len(backbones)} backbones, {len(present)} arms ===")

    reference = {r["backbone"]: r["objective"] for r in rows if r["arm"] == "codon_dp"}
    header = ("arm".ljust(20) + "objective".rjust(11) + "vs codon DP".rjust(13)
              + "Tier-1".rjust(10) + "burden".rjust(9) + "violations".rjust(12)
              + "identity".rjust(10) + "seconds".rjust(10) + "evals".rjust(9))
    print(header)
    for arm in present:
        subset = [r for r in rows if r["arm"] == arm]
        gap = [r["objective"] - reference[r["backbone"]] for r in subset
               if r["backbone"] in reference]
        print(f"{arm:<20}{np.median([r['objective'] for r in subset]):>11.3f}"
              f"{np.median(gap):>+13.3f}"
              f"{np.median([r['tier1'] for r in subset]):>10.3f}"
              f"{np.median([r['burden'] for r in subset]):>9.3f}"
              f"{np.mean([r['violations'] for r in subset]):>12.2f}"
              f"{np.mean([r['identity_to_native'] for r in subset]):>10.3f}"
              f"{np.median([r['seconds'] for r in subset]):>10.3f}"
              f"{np.median([r['evaluations'] for r in subset]):>9.0f}")

    print("\nper term, so a method that optimises something else is legible")
    print("arm".ljust(20) + "CAI".rjust(9) + "codon pair".rjust(12)
          + "GC".rjust(8) + "initiation".rjust(12))
    for arm in present:
        subset = [r for r in rows if r["arm"] == arm]
        print(f"{arm:<20}{np.median([r['cai'] for r in subset]):>9.3f}"
              f"{np.median([r['cpb'] for r in subset]):>12.3f}"
              f"{np.median([r['gc'] for r in subset]):>8.3f}"
              f"{np.median([r['initiation'] for r in subset]):>12.3f}")

    print("\nshare of backbones on which each arm beats the codon DP")
    for arm in present:
        subset = [r for r in rows if r["arm"] == arm and r["backbone"] in reference]
        wins = np.mean([r["objective"] > reference[r["backbone"]] + 1e-9 for r in subset])
        print(f"  {arm:<20}{100 * wins:>6.0f}%")

    joint = {r["backbone"]: r["objective"] for r in rows if r["arm"] == "joint"}
    if joint:
        print("\nshare of the joint arm's gain over the codon DP that each arm recovers")
        for arm in present:
            if arm in ("joint", "codon_dp"):
                continue
            subset = [r for r in rows if r["arm"] == arm]
            shares = [(r["objective"] - reference[r["backbone"]])
                      / (joint[r["backbone"]] - reference[r["backbone"]])
                      for r in subset
                      if r["backbone"] in reference and r["backbone"] in joint
                      and abs(joint[r["backbone"]] - reference[r["backbone"]]) > 1e-9]
            if shares:
                print(f"  {arm:<20}{100 * np.median(shares):>6.0f}%")


if __name__ == "__main__":
    main()
