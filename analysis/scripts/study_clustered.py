"""Discrimination on the outcome benchmark when designs are grouped by campaign.

The benchmark pools eleven design studies, and success rates differ between them.
Designs from one study share a design method, a laboratory and an expression
protocol, so treating the 614 as independent draws can inflate an AUC: a score
that happens to track campaign identity is rewarded for it. This is the same
non-independence the atlas handles by resampling source proteins.

The released benchmark carries no study identifier. It carries a fold class, and
fold classes are nested within studies, so two grouping levels are reported:

``fold``      the 53 fold classes as shipped. No inference, but finer than a
              study, so it under-corrects: two topologies from one laboratory
              land in different groups.
``campaign``  fold classes gathered by name stem, giving 13 groups. Closer to
              study level and the more conservative of the two, but the mapping
              is read off the naming and is not the authors' own label.

Neither is exactly leave-one-study-out. Both replace the pooled estimate with one
that never scores a design using a model trained on its own group, and both
replace the record bootstrap with a bootstrap over groups.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from predictive import DROP, delong, load_classes, read_model
from features import protein_features
from janus.objectives import liability

CAMPAIGNS = [
    (re.compile(r"^Foldit$"), "Foldit"),
    (re.compile(r"^NF"), "NF series"),
    (re.compile(r"^Di-"), "Di series"),
    (re.compile(r"^Fd_"), "Fd series"),
    (re.compile(r"^dcs_"), "dcs series"),
    (re.compile(r"^H\d_fold"), "helical bundles"),
    (re.compile(r"^(R\dx\d|Pl\dx\d|Rsmn)"), "repeat and propeller"),
    (re.compile(r"^BH"), "BH series"),
    (re.compile(r"^HBI"), "HBI series"),
    (re.compile(r"^DA?(05R)?$"), "D series"),
    (re.compile(r"^CA$"), "CA"),
    (re.compile(r"^BB$"), "BB"),
    (re.compile(r"^K$"), "K"),
]


def campaign_of(fold: str) -> str:
    for pattern, label in CAMPAIGNS:
        if pattern.match(fold):
            return label
    return f"other ({fold})"


def grouped_predictions(matrix, labels, groups, seed=0):
    """Out-of-fold probabilities where no design is scored by its own group."""
    if matrix.shape[1] == 0:
        return np.full(len(labels), 0.5)
    out = np.zeros(len(labels))
    splitter = LeaveOneGroupOut()
    for train, test in splitter.split(matrix, labels, groups):
        if len(set(labels[train])) < 2:
            out[test] = labels[train].mean() if len(train) else 0.5
            continue
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, C=0.1, random_state=seed))
        model.fit(matrix[train], labels[train])
        out[test] = model.predict_proba(matrix[test])[:, 1]
    return out


def cluster_auc_interval(labels, scores, groups, rng, resamples):
    """Percentile interval for one AUC, resampling groups with replacement."""
    blocks = {}
    for index, key in enumerate(groups):
        blocks.setdefault(key, []).append(index)
    keys = list(blocks)
    draws = []
    for _ in range(resamples):
        picked = np.concatenate([blocks[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
        if len(set(labels[picked])) < 2:
            continue
        draws.append(roc_auc_score(labels[picked], scores[picked]))
    if not draws:
        return float("nan"), float("nan")
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def cluster_difference(labels, first, second, groups, rng, resamples):
    """Interval and two-sided p for AUC(first) - AUC(second), resampling groups."""
    blocks = {}
    for index, key in enumerate(groups):
        blocks.setdefault(key, []).append(index)
    keys = list(blocks)
    draws = []
    for _ in range(resamples):
        picked = np.concatenate([blocks[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
        if len(set(labels[picked])) < 2:
            continue
        draws.append(roc_auc_score(labels[picked], first[picked])
                     - roc_auc_score(labels[picked], second[picked]))
    if not draws:
        return float("nan"), float("nan"), float("nan")
    draws = np.asarray(draws)
    share = float(np.mean(draws > 0))
    p = 2 * min(share, 1 - share)
    return (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)),
            float(min(1.0, max(p, 1 / len(draws)))))


def per_group(labels, scores, groups):
    """AUC on each held-out group, for the groups that carry both outcomes."""
    out = []
    for key in sorted(set(groups)):
        mask = np.array([g == key for g in groups])
        if len(set(labels[mask])) < 2:
            out.append({"group": key, "n": int(mask.sum()), "auc": None,
                        "success_rate": float(labels[mask].mean())})
            continue
        out.append({"group": key, "n": int(mask.sum()),
                    "auc": float(roc_auc_score(labels[mask], scores[mask])),
                    "success_rate": float(labels[mask].mean())})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--elm", required=True)
    parser.add_argument("--structures", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--floors", type=float, nargs="+", default=[60.0, 70.0, 80.0])
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
        records.append({"name": row["Name"], "sequence": sequence,
                        "fold": row.get("Fold", ""), "success": label == "TRUE"})

    models = {}
    for record in records:
        path = Path(args.structures) / f"{record['name']}.pdb"
        if not path.exists():
            continue
        read = read_model(path, record["sequence"])
        if read is not None:
            models[record["name"]] = read

    print(f"{len(records)} designs, {len(models)} with an ESMFold model")
    folds = {r["fold"] for r in records}
    campaigns = {campaign_of(r["fold"]) for r in records}
    print(f"{len(folds)} fold classes, {len(campaigns)} campaign groups")
    counts = {}
    for r in records:
        counts[campaign_of(r["fold"])] = counts.get(campaign_of(r["fold"]), 0) + 1
    for key, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"    {key:<24}{n:>5}")

    def build(record, rsa, sse):
        computed = protein_features(record["sequence"], classes, rsa, sse)
        computed = {k: v for k, v in computed.items()
                    if k not in DROP and isinstance(v, (int, float))}
        length = len(record["sequence"])
        for key in list(computed):
            if key.endswith("_raw"):
                computed[f"{key}_density"] = 100.0 * computed[key] / length
        computed["seq_len"] = length
        computed["protein_repeat"] = liability.score("protein_repeat", record["sequence"])
        if rsa is not None:
            for name in ("exposed_hydrophobic", "exposed_hydrophobic_run"):
                computed[name] = liability.score(name, record["sequence"], rsa)
        return computed

    plain = [build(r, None, None) for r in records]
    names = sorted({k for row in plain for k in row})
    sequence_matrix = np.nan_to_num(
        np.array([[float(row.get(k, 0.0)) for k in names] for row in plain]))
    labels = np.array([r["success"] for r in records])
    plddt = np.array([models[r["name"]][2] if r["name"] in models else np.nan
                      for r in records])

    rng = np.random.default_rng(args.seed)
    results = {"campaign_counts": counts, "cohorts": {}}

    for floor in args.floors:
        usable = np.array([r["name"] in models and models[r["name"]][2] >= floor
                           for r in records])
        for cohort, mask in (("all designs", np.ones(len(records), dtype=bool)),
                             (f"pLDDT at least {floor:g}", usable)):
            if cohort in results["cohorts"] or mask.sum() < 40:
                continue
            weighted = []
            for keep, record in zip(mask, records, strict=True):
                if not keep:
                    weighted.append(None)
                    continue
                rsa, sse, _ = models.get(record["name"], (None, None, None))
                weighted.append(build(record, rsa, sse))
            present = [w for w in weighted if w is not None]
            wnames = sorted({k for row in present for k in row})
            weighted_matrix = np.nan_to_num(
                np.array([[float(row.get(k, 0.0)) for k in wnames] for row in present]))

            subset = {
                "labels": labels[mask],
                "plddt": plddt[mask],
                "sequence_matrix": sequence_matrix[mask],
                "weighted_matrix": weighted_matrix,
                "fold": [r["fold"] for r, m in zip(records, mask, strict=True) if m],
            }
            subset["campaign"] = [campaign_of(f) for f in subset["fold"]]

            entry = {"n": int(mask.sum()),
                     "positives": int(subset["labels"].sum()),
                     "levels": {}}

            for level in ("fold", "campaign"):
                groups = subset[level]
                sequence = grouped_predictions(subset["sequence_matrix"],
                                               subset["labels"], groups, args.seed)
                weighted_pred = grouped_predictions(subset["weighted_matrix"],
                                                    subset["labels"], groups, args.seed)
                joint = np.column_stack([subset["weighted_matrix"], subset["plddt"]])
                combined = grouped_predictions(np.nan_to_num(joint),
                                               subset["labels"], groups, args.seed)
                scores = {"plddt": subset["plddt"], "sequence": sequence,
                          "weighted": weighted_pred, "combined": combined}

                block = {"n_groups": len(set(groups)), "scores": {}}
                for name, values in scores.items():
                    values = np.nan_to_num(values, nan=np.nanmin(values))
                    auc = float(roc_auc_score(subset["labels"], values))
                    low, high = cluster_auc_interval(subset["labels"], values, groups,
                                                     rng, args.resamples)
                    record = {"auc": auc, "ci": [low, high]}
                    if name != "plddt":
                        d_low, d_high, p = cluster_difference(
                            subset["labels"], values, np.nan_to_num(subset["plddt"]),
                            groups, rng, args.resamples)
                        record["vs_plddt"] = {"difference": auc - float(
                            roc_auc_score(subset["labels"], np.nan_to_num(subset["plddt"]))),
                            "ci": [d_low, d_high], "p": p}
                        # delong returns (auc_a, auc_b, p)
                        record["delong_p"] = float(delong(
                            subset["labels"], values, np.nan_to_num(subset["plddt"]))[2])
                    block["scores"][name] = record
                block["per_group"] = per_group(subset["labels"], sequence, groups)
                entry["levels"][level] = block

            results["cohorts"][cohort] = entry

            print()
            print(f"=== {cohort}: n = {entry['n']}, {entry['positives']} successful ===")
            for level, block in entry["levels"].items():
                print(f"  grouped by {level} ({block['n_groups']} groups), "
                      f"leave-one-group-out")
                print(f"    {'score':<12}{'AUC':>8}{'95% CI':>20}"
                      f"{'vs pLDDT':>11}{'clustered p':>14}{'DeLong p':>11}")
                for name, record in block["scores"].items():
                    low, high = record["ci"]
                    if "vs_plddt" in record:
                        v = record["vs_plddt"]
                        print(f"    {name:<12}{record['auc']:>8.3f}"
                              f"{f'[{low:.2f}, {high:.2f}]':>20}"
                              f"{v['difference']:>+11.3f}{v['p']:>14.3g}"
                              f"{record['delong_p']:>11.3g}")
                    else:
                        print(f"    {name:<12}{record['auc']:>8.3f}"
                              f"{f'[{low:.2f}, {high:.2f}]':>20}")

    print()
    print("=== held-out group performance, sequence-only score, campaign level ===")
    for cohort, entry in results["cohorts"].items():
        block = entry["levels"]["campaign"]
        print(f"  {cohort}")
        for row in sorted(block["per_group"], key=lambda r: -r["n"]):
            auc = f"{row['auc']:.3f}" if row["auc"] is not None else "one class"
            print(f"    {row['group']:<24}{row['n']:>5}  success {100*row['success_rate']:>5.1f}%"
                  f"  AUC {auc}")

    Path(args.out).write_text(json.dumps(results, indent=1), encoding="utf-8")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
