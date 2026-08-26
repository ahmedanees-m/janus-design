"""What the initiation opening is worth in units a wet-lab reader can price.

The exchange rate is measured in kcal/mol of initiation-window folding energy
bought per nat of unconditional fold log-likelihood. Neither unit means anything
to someone deciding whether to order a gene, so this converts the first one into
a predicted fold-change in protein output.

The conversion is the apparent Boltzmann factor from the Ribosome Binding Site
Calculator, which relates translation initiation rate to the total free energy of
initiation as rate proportional to exp(-beta * dG_total). In an ideal system beta
would be 1/RT, about 1.62 mol/kcal at 37 C; fitted against measured expression in
the crowded interior of E. coli it is 0.45 plus or minus 0.05 mol/kcal (Salis,
Mirsky and Voigt 2009).

This is an extrapolation and not a measurement, in three specific ways.

The RBS Calculator's dG_total is a composite of ribosome-mRNA hybridisation,
spacing, standby site and the cost of unfolding mRNA structure. The term measured
here is the folding energy of the initiation window alone. Applying beta to it
assumes that window's energy enters dG_total roughly one for one, which is what
the model intends but is not an equivalence anyone has measured for this term.

beta was fitted in E. coli, so nothing here transfers to the HEK293 host, whose
initiation is cap-scanning rather than Shine-Dalgarno.

And a fold-change in initiation rate is not a fold-change in yield. Initiation is
rate-limiting for expression over a wide range but not everywhere, and the
prediction saturates where something else becomes limiting.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BETA = 0.45
BETA_ERROR = 0.05


def fold_change(energy, beta=BETA):
    """Predicted change in initiation rate for an opening of ``energy`` kcal/mol.

    Opening the window means a folding energy closer to zero, so a positive
    opening lowers the cost of unfolding it and raises the rate.
    """
    return float(np.exp(beta * energy))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recall", help="rescore_recall output, normalised scale")
    parser.add_argument("--sweep", help="folding_sweep output, the hybrid pool")
    parser.add_argument("--entropy", help="entropy budget output, for the budget share")
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if bool(args.recall) == bool(args.sweep):
        raise SystemExit("give exactly one of --recall or --sweep")

    # The two sources rank different pools, and the exchange rate is a property of
    # the pool as much as of the weight, so which one produced a number has to
    # travel with it. --recall is the k-best prefix alone, which is the pool the
    # two-tier bound endorses at the operating weight. --sweep is that prefix plus
    # shell draws, which reaches further into protein space for a higher price.
    if args.recall:
        raw = json.loads(Path(args.recall).read_text(encoding="utf-8"))
        if not raw[0].get("normalised"):
            raise SystemExit("this expects the normalised run of rescore_recall")
        pool = "k-best prefix"
        rows = [{"lambda_initiation": r["lambda_initiation"],
                 "opened": r["initiation_winner"] - r["initiation_top"],
                 "paid": r["tier1_top"] - r["tier1_of_winner"]} for r in raw]
    else:
        raw = json.loads(Path(args.sweep).read_text(encoding="utf-8"))
        pool = "k-best prefix plus shell draws"
        rows = [{"lambda_initiation": r["lam"],
                 "opened": r["initiation"] - r["initiation_reference"],
                 "paid": r["tier1_reference"] - r["tier1"]} for r in raw]

    budget = None
    if args.entropy:
        entropy = json.loads(Path(args.entropy).read_text(encoding="utf-8"))
        designed = [r for r in entropy if r.get("origin") == "designed"]
        if designed:
            # Total marginal entropy over the chain, which is the freedom the
            # backbone carries before any shell is drawn around it.
            budget = float(np.median([sum(r["entropy"]) for r in designed]))

    weights = sorted({r["lambda_initiation"] for r in rows})
    results = []
    print(f"beta = {BETA} plus or minus {BETA_ERROR} mol/kcal, E. coli, "
          f"Salis Mirsky and Voigt 2009")
    print(f"pool: {pool}, {len({r['lambda_initiation'] for r in rows})} weights over "
          f"{len(rows) // max(len({r['lambda_initiation'] for r in rows}), 1)} backbones")
    print()
    header = ("weight".rjust(8) + "opened".rjust(10) + "paid".rjust(9)
              + "fold change".rjust(14) + "range".rjust(16))
    if budget:
        header += "budget".rjust(10)
    print(header)
    for weight in weights:
        subset = [r for r in rows if r["lambda_initiation"] == weight]
        opened = float(np.median([r["opened"] for r in subset]))
        paid = float(np.median([r["paid"] for r in subset]))
        low = fold_change(opened, BETA - BETA_ERROR)
        high = fold_change(opened, BETA + BETA_ERROR)
        record = {"weight": weight, "opened_kcal": opened, "paid_nats": paid,
                  "fold_change": fold_change(opened), "low": low, "high": high}
        line = (f"{weight:>8.3f}{opened:>+10.3f}{paid:>9.4f}"
                f"{fold_change(opened):>13.2f}x{f'{low:.2f} to {high:.2f}':>16}")
        if budget:
            record["budget_share"] = paid / budget
            line += f"{100 * paid / budget:>9.3f}%"
        results.append(record)
        print(line)

    print()
    print("the 90th percentile of what the shell can open, at each weight")
    print("weight".rjust(8) + "opened".rjust(10) + "paid".rjust(9) + "fold change".rjust(14))
    for weight in weights:
        subset = [r for r in rows if r["lambda_initiation"] == weight]
        opened = float(np.percentile([r["opened"] for r in subset], 90))
        paid = float(np.percentile([r["paid"] for r in subset], 90))
        print(f"{weight:>8.3f}{opened:>+10.3f}{paid:>9.4f}{fold_change(opened):>13.2f}x")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"beta": BETA, "beta_error": BETA_ERROR, "pool": pool,
                               "budget_nats": budget, "rows": results}),
                   encoding="utf-8")
    print()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
