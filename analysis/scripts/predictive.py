"""H5: does a sequence-level liability score add to ESMFold pLDDT?

The benchmark supplies sequences and labels but no backbones. Without a structure
the accessibility weighting that carries the atlas cannot be applied and every
motif feature is a raw count, which is the weaker form of the score. With
`--structures` pointing at ESMFold models of the same sequences, the weighted
form is available too and both are reported, so the narrowing is measured rather
than assumed.

A weighting taken from a predicted structure is only as good as the prediction,
so models below a pLDDT floor do not supply accessibility.

That creates a circularity the comparison has to control. The structures come
from ESMFold and the comparator is ESMFold pLDDT, so any feature read off the
model carries some of the confidence signal: a design the model folds badly gets
an extended, highly accessible structure, and its accessibility-weighted
liabilities are then partly a restatement of a low pLDDT. Two controls are
reported. The weighted arm is repeated on only the designs with usable models, so
that no part of the gain comes from which designs had a model at all; and the
features read off the model are reported as their own family in the attribution,
so their contribution is visible rather than folded into the liability total.

Comparison is against ESMFold pLDDT, which the source paper found to be the best
single discriminator, by DeLong's test on correlated ROC curves.

A single AUC for a bag of forty features says nothing about which liability is
doing the work, and there is a tension that needs it: the atlas finds designs are
not degron-enriched against natural proteins, while this score adds AUC over
pLDDT. Those are compatible, one being a between-group comparison and the other
within-design, but which features carry the discrimination decides whether the
two results tell one story. So each family of features is reported alone and with
itself removed.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import biotite.structure as struc
import biotite.structure.io.pdb as pdb

from features import protein_features
from janus.objectives import liability
from janus.objectives.proteostasis import load_classes

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


def read_model(path, expected):
    """Per-residue accessibility, secondary structure and mean pLDDT from a model.

    Returns nothing when the model does not match the sequence it is supposed to
    be, which is the only case where a silent misalignment could put one
    design's accessibility on another design's residues.
    """
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1, extra_fields=["b_factor"])
    atoms = atoms[struc.filter_amino_acids(atoms)]
    residues = [THREE_TO_ONE.get(name, "X") for name in struc.get_residues(atoms)[1]]
    if "".join(residues) != expected:
        return None
    area = struc.apply_residue_wise(atoms, struc.sasa(atoms, vdw_radii="Single"), np.nansum)
    rsa = [None if r not in MAX_ASA or np.isnan(area[i]) else float(area[i] / MAX_ASA[r])
           for i, r in enumerate(residues)]
    sse = "".join(struc.annotate_sse(atoms))
    plddt = float(np.mean(atoms.b_factor)) if atoms.array_length() else float("nan")
    return rsa, sse, plddt

DROP = {"n_end_class", "c_end_residue", "length"}


def bootstrap_auc(labels, scores, rng, resamples=2000):
    """Percentile interval for one AUC, resampling designs with replacement."""
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    draws = []
    for _ in range(resamples):
        index = rng.integers(0, len(labels), len(labels))
        if len(set(labels[index])) < 2:
            continue
        draws.append(roc_auc_score(labels[index], scores[index]))
    if not draws:
        return float("nan"), float("nan")
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def delong(labels, first, second):
    """DeLong's test for two correlated ROC curves."""
    labels = np.asarray(labels, dtype=bool)
    positives, negatives = np.sum(labels), np.sum(~labels)

    def structural(scores):
        pos, neg = scores[labels], scores[~labels]
        v10 = np.array([(np.sum(neg < p) + 0.5 * np.sum(neg == p)) / negatives for p in pos])
        v01 = np.array([(np.sum(pos > n) + 0.5 * np.sum(pos == n)) / positives for n in neg])
        return v10, v01

    a10, a01 = structural(first)
    b10, b01 = structural(second)
    auc_a, auc_b = a10.mean(), b10.mean()
    s10 = np.cov(np.vstack([a10, b10]))
    s01 = np.cov(np.vstack([a01, b01]))
    covariance = s10 / positives + s01 / negatives
    variance = covariance[0, 0] + covariance[1, 1] - 2 * covariance[0, 1]
    if variance <= 0:
        return auc_a, auc_b, np.nan
    z = (auc_a - auc_b) / np.sqrt(variance)
    return auc_a, auc_b, float(2 * (1 - stats.norm.cdf(abs(z))))


GROUPS = {
    "degron motifs": lambda name: "degron" in name,
    "other ELM motifs": lambda name: any(
        stem in name for stem in ("protease_site", "targeting", "modification")),
    "complexity and repeats": lambda name: name in ("low_complexity_fraction",
                                                    "protein_repeat"),
    "bulk composition": lambda name: name in ("gravy", "isoelectric_point",
                                              "net_charge", "free_cysteines"),
    "length": lambda name: name == "seq_len",
    "model geometry": lambda name: name in (
        "mean_rsa", "n_term_rsa", "c_term_rsa", "helix_fraction", "strand_fraction",
        "coil_fraction", "exposed_hydrophobic", "exposed_hydrophobic_run",
        "exposed_hydrophobic_area", "exposed_hydrophobic_patches",
        "longest_exposed_hydrophobic_run"),
}


def cross_validated(matrix, labels, folds, seed):
    """Out-of-fold probabilities from the same model the headline score uses."""
    if matrix.shape[1] == 0:
        return np.full(len(labels), 0.5)
    out = np.zeros(len(labels))
    for train, test in folds.split(matrix, labels):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1))
        model.fit(matrix[train], labels[train])
        out[test] = model.predict_proba(matrix[test])[:, 1]
    return out


def attribution(matrix, names, labels, folds, seed, full_auc):
    """Each family of features alone, and with itself removed."""
    print()
    print("=== which liabilities carry the discrimination ===")
    print(f"{'family':<26}{'n':>4}{'alone':>9}{'without':>10}{'drop':>8}")
    results = {}
    for label, belongs in GROUPS.items():
        indices = [i for i, name in enumerate(names) if belongs(name)]
        if not indices:
            continue
        rest = [i for i in range(len(names)) if i not in indices]
        alone = roc_auc_score(labels, cross_validated(matrix[:, indices], labels, folds, seed))
        without = roc_auc_score(labels, cross_validated(matrix[:, rest], labels, folds, seed))
        results[label] = {"n": len(indices), "alone": float(alone),
                          "without": float(without), "drop": float(full_auc - without)}
        print(f"{label:<26}{len(indices):>4}{alone:>9.3f}{without:>10.3f}"
              f"{full_auc - without:>+8.3f}")
    unassigned = [n for n in names if not any(belongs(n) for belongs in GROUPS.values())]
    if unassigned:
        print(f"  {len(unassigned)} features in no family: {', '.join(sorted(unassigned)[:6])}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--structures", help="directory of ESMFold models, named by design")
    parser.add_argument("--plddt-floor", type=float, default=70.0,
                        help="models below this mean pLDDT do not supply accessibility")
    parser.add_argument("--resamples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    classes = load_classes(args.elm)
    with Path(args.benchmark).open(newline="", encoding="utf-8", errors="replace") as fh:
        raw = list(csv.DictReader(fh))

    records = []
    for row in raw:
        sequence = (row.get("sequence") or "").strip().upper()
        label = (row.get("Exp_success") or "").strip().upper()
        plddt = row.get("ESMFold_pLDDT")
        if not sequence or label not in ("TRUE", "FALSE") or not plddt:
            continue
        if set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        records.append({
            "name": row["Name"], "sequence": sequence, "fold": row.get("Fold", ""),
            "success": label == "TRUE", "plddt": float(plddt),
            "mpnn": float(row["MPNN_score"]) if row.get("MPNN_score") else np.nan,
        })

    print(f"{len(records)} usable designs of {len(raw)}, "
          f"{sum(r['success'] for r in records)} successful "
          f"({100 * np.mean([r['success'] for r in records]):.1f}%)")

    structures = {}
    if args.structures:
        weak = 0
        for record in records:
            path = Path(args.structures) / f"{record['name']}.pdb"
            if not path.exists():
                continue
            model = read_model(path, record["sequence"])
            if model is None:
                continue
            rsa, sse, plddt = model
            if plddt < args.plddt_floor:
                weak += 1
                continue
            structures[record["name"]] = (rsa, sse)
        print(f"{len(structures)} models usable of {len(records)} designs, "
              f"{weak} below the pLDDT floor of {args.plddt_floor}")

    def build(record, use_structure):
        rsa, sse = structures.get(record["name"], (None, None)) if use_structure else (None, None)
        computed = protein_features(record["sequence"], classes, rsa, sse)
        computed = {k: v for k, v in computed.items()
                    if k not in DROP and isinstance(v, (int, float))}
        length = len(record["sequence"])
        for key in list(computed):
            if key.endswith("_raw"):
                computed[f"{key}_density"] = 100.0 * computed[key] / length
        computed["seq_len"] = length
        # Repeat content needs no structure, so it belongs in both arms.
        computed["protein_repeat"] = liability.score("protein_repeat", record["sequence"])
        if rsa is not None:
            for name in ("exposed_hydrophobic", "exposed_hydrophobic_run"):
                computed[name] = liability.score(name, record["sequence"], rsa)
        return {**record, **computed}

    rows = [build(record, False) for record in records]
    weighted_rows = [build(record, True) for record in records] if structures else None

    names = sorted({k for r in rows for k in r
                    if isinstance(r.get(k), (int, float)) and k not in ("success", "plddt", "mpnn")})
    matrix = np.array([[float(r.get(k, 0.0)) for k in names] for r in rows])
    matrix = np.nan_to_num(matrix)
    labels = np.array([r["success"] for r in rows])
    plddt = np.array([r["plddt"] for r in rows])

    print(f"{len(names)} sequence-level features")

    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.seed)
    predictions = {"liability": np.zeros(len(rows)), "combined": np.zeros(len(rows))}
    for train, test in folds.split(matrix, labels):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1))
        model.fit(matrix[train], labels[train])
        predictions["liability"][test] = model.predict_proba(matrix[test])[:, 1]

        joint = np.column_stack([matrix, plddt])
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1))
        model.fit(joint[train], labels[train])
        predictions["combined"][test] = model.predict_proba(joint[test])[:, 1]

    print("\n=== discrimination, five-fold cross-validated ===")
    print(f"{'score':<26}{'ROC AUC':>10}{'PR AUC':>10}")
    scores = {"ESMFold pLDDT": plddt,
              "liability score": predictions["liability"],
              "pLDDT and liability": predictions["combined"]}
    results = {}
    for label, values in scores.items():
        roc = roc_auc_score(labels, values)
        pr = average_precision_score(labels, values)
        results[label] = {"roc": float(roc), "pr": float(pr)}
        print(f"{label:<26}{roc:>10.3f}{pr:>10.3f}")

    results["attribution"] = attribution(
        matrix, names, labels, folds, args.seed, results["liability score"]["roc"])

    if weighted_rows:
        weighted_names = sorted({k for r in weighted_rows for k in r
                                 if isinstance(r.get(k), (int, float))
                                 and k not in ("success", "plddt", "mpnn")})
        weighted = np.nan_to_num(np.array(
            [[float(r.get(k, 0.0)) for k in weighted_names] for r in weighted_rows]))

        # Only the designs with a usable model, so no part of the comparison
        # comes from which designs had one.
        keep = np.array([r["name"] in structures for r in rows])
        print()
        print("=== accessibility weighting, and what of it is the folding model ===")
        print(f"{len(weighted_names)} features against {len(names)} sequence-only")
        print(f"{'cohort':<16}{'n':>6}{'pLDDT':>9}{'sequence':>10}{'weighted':>10}"
              f"{'gain':>8}{'combined':>10}")
        for label, mask in (("all designs", np.ones(len(rows), dtype=bool)),
                            ("usable models", keep)):
            subset_labels = labels[mask]
            subset_plddt = plddt[mask]
            sequence_only = cross_validated(matrix[mask], subset_labels, folds, args.seed)
            with_model = cross_validated(weighted[mask], subset_labels, folds, args.seed)
            together = cross_validated(
                np.column_stack([weighted[mask], subset_plddt]), subset_labels,
                folds, args.seed)
            entry = {
                "n": int(mask.sum()),
                "plddt": float(roc_auc_score(subset_labels, subset_plddt)),
                "sequence": float(roc_auc_score(subset_labels, sequence_only)),
                "weighted": float(roc_auc_score(subset_labels, with_model)),
                "combined": float(roc_auc_score(subset_labels, together)),
            }
            entry["gain"] = entry["weighted"] - entry["sequence"]
            # Whether the liability score beats pLDDT is a different question from
            # whether weighting helps, and it reverses between the two cohorts, so
            # it is tested inside each rather than only on the pooled set.
            for name, values in (("sequence", sequence_only), ("weighted", with_model),
                                 ("combined", together)):
                _, _, p_value = delong(subset_labels, values, subset_plddt)
                entry[f"delong_{name}_vs_plddt"] = None if np.isnan(p_value) else float(p_value)
            results[f"weighted, {label}"] = entry
            print(f"{label:<16}{entry['n']:>6}{entry['plddt']:>9.3f}"
                  f"{entry['sequence']:>10.3f}{entry['weighted']:>10.3f}"
                  f"{entry['gain']:>+8.3f}{entry['combined']:>10.3f}")

        print()
        print("against pLDDT alone, within each cohort, DeLong two-sided")
        print(f"{'cohort':<16}{'sequence':>12}{'weighted':>12}{'combined':>12}")
        for label in ("all designs", "usable models"):
            entry = results[f"weighted, {label}"]
            cells = "".join(
                (f"{entry[f'delong_{name}_vs_plddt']:>12.2g}"
                 if entry.get(f"delong_{name}_vs_plddt") is not None else f"{'-':>12}")
                for name in ("sequence", "weighted", "combined"))
            print(f"{label:<16}{cells}")

        results["attribution_weighted"] = attribution(
            weighted[keep], weighted_names, labels[keep], folds, args.seed,
            results["weighted, usable models"]["weighted"])

    print("\n=== DeLong against pLDDT alone ===")
    for label in ("liability score", "pLDDT and liability"):
        a, b, p = delong(labels, scores[label], plddt)
        results[label]["delong_p"] = None if np.isnan(p) else float(p)
        print(f"  {label:<24} AUC {a:.3f} against {b:.3f}, p = {p:.3g}")

    print("\n=== stratified by fold class, classes with at least 15 designs ===")
    strata = {}
    for record in rows:
        strata.setdefault(record["fold"], []).append(record)
    rng = np.random.default_rng(args.seed)
    by_fold = {}
    for fold, subset in sorted(strata.items(), key=lambda kv: -len(kv[1]))[:8]:
        if len(subset) < 15:
            continue
        index = [rows.index(r) for r in subset]
        y = labels[index]
        if len(set(y)) < 2:
            print(f"  {fold:<18} n={len(subset):<4} single class, skipped")
            continue
        entry = {
            "n": len(subset),
            "success_rate": float(y.mean()),
            "plddt": float(roc_auc_score(y, plddt[index])),
            "liability": float(roc_auc_score(y, predictions["liability"][index])),
        }
        for key, values in (("plddt", plddt[index]),
                            ("liability", predictions["liability"][index])):
            entry[f"{key}_ci"] = list(bootstrap_auc(y, values, rng, args.resamples))
        by_fold[fold] = entry
        low, high = entry["liability_ci"]
        print(f"  {fold:<18} n={len(subset):<4} success {100 * y.mean():>5.1f}%  "
              f"pLDDT {entry['plddt']:.3f}  "
              f"liability {entry['liability']:.3f} [{low:.2f}, {high:.2f}]")
    results["fold_classes"] = by_fold

    Path(args.out).write_text(json.dumps({
        "results": results,
        "n": len(rows),
        "positives": int(labels.sum()),
        "features": names,
    }), encoding="utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
