"""Are the natural controls whole proteins, or domains cut out of larger ones?

The atlas reports that designs bury their termini more than natural domains do,
and that this drives a lower accessibility-weighted degron burden. MegaScale's
natural arm is drawn from the PDB at 40 to 80 residues, which is the size range
where excised domains live. A domain cut out of a larger chain has termini that
were never solvent exposed in the parent, so the excision protocol alone could
manufacture the effect.

This classifies each natural entry by how much of its UniProt reference sequence
the deposited entity covers, using RCSB's own alignment, and re-runs the
comparison on the whole-protein subset.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.request
from pathlib import Path

import numpy as np

PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
UA = {"User-Agent": "janus-excision-check", "Content-Type": "application/json"}

QUERY = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities{
      entity_poly{ rcsb_sample_sequence_length }
      rcsb_polymer_entity_align{
        reference_database_name
        reference_database_accession
        aligned_regions{ length ref_beg_seq_id }
      }
      rcsb_polymer_entity_container_identifiers{
        reference_sequence_identifiers{ database_accession database_name }
      }
    }
  }
}"""


def fetch(ids):
    body = json.dumps({"query": QUERY, "variables": {"ids": ids}}).encode()
    request = urllib.request.Request("https://data.rcsb.org/graphql", data=body, headers=UA)
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.load(response)["data"]["entries"]


def uniprot_length(accession, cache):
    if accession in cache:
        return cache[accession]
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json?fields=length"
    try:
        request = urllib.request.Request(url, headers={"User-Agent": UA["User-Agent"]})
        with urllib.request.urlopen(request, timeout=60) as response:
            cache[accession] = int(json.load(response)["sequence"]["length"])
    except Exception:
        cache[accession] = None
    return cache[accession]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--coverage", type=float, default=0.90,
                        help="fraction of the reference a whole protein must cover")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    atlas = json.loads(Path(args.atlas).read_text(encoding="utf-8"))
    natural = [r for r in atlas if r["group"] == "natural"]
    codes = sorted({r["name"].split("_for_")[0] for r in natural})
    codes = [c for c in codes if PDBLIKE.match(c)]
    print(f"{len(natural)} natural records, {len(codes)} distinct PDB entries")

    records = []
    for start in range(0, len(codes), 150):
        chunk = codes[start : start + 150]
        try:
            records.extend(fetch(chunk))
        except Exception as exc:
            print(f"  chunk at {start} failed: {type(exc).__name__}")
        time.sleep(1)
    print(f"retrieved {len(records)} entries")

    cache: dict[str, int | None] = {}
    coverage = {}
    for entry in records:
        if not entry:
            continue
        best = None
        for entity in entry.get("polymer_entities") or []:
            length = (entity.get("entity_poly") or {}).get("rcsb_sample_sequence_length")
            aligns = entity.get("rcsb_polymer_entity_align") or []
            for align in aligns:
                if align.get("reference_database_name") != "UniProt":
                    continue
                accession = align.get("reference_database_accession")
                reference = uniprot_length(accession, cache) if accession else None
                if not reference or not length:
                    continue
                aligned = sum(r.get("length") or 0 for r in align.get("aligned_regions") or [])
                fraction = aligned / reference
                best = fraction if best is None else max(best, fraction)
        if best is not None:
            coverage[entry["rcsb_id"]] = best

    print(f"coverage resolved for {len(coverage)} of {len(codes)} entries")
    whole = {k for k, v in coverage.items() if v >= args.coverage}
    excised = {k for k, v in coverage.items() if v < args.coverage}
    print(f"  whole proteins (covering at least {args.coverage:.0%}): {len(whole)}")
    print(f"  excised fragments: {len(excised)}")
    if coverage:
        values = np.array(list(coverage.values()))
        print(f"  coverage: median {np.median(values):.2f}, "
              f"10th {np.percentile(values, 10):.2f}, 90th {np.percentile(values, 90):.2f}")

    Path(args.out).write_text(json.dumps(coverage), encoding="utf-8")
    compare(atlas, whole, excised)
    print(f"\nwrote {args.out}")


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    ordered = np.sort(b)
    greater = np.searchsorted(ordered, a, side="left").sum()
    less = (len(b) - np.searchsorted(ordered, a, side="right")).sum()
    return (greater - less) / (len(a) * len(b))


def bootstrap_delta(a, b, rng, resamples=2000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return (np.nan, np.nan)
    draws = np.empty(resamples)
    for i in range(resamples):
        draws[i] = cliffs_delta(rng.choice(a, len(a)), rng.choice(b, len(b)))
    return tuple(np.percentile(draws, [2.5, 97.5]))


def compare(atlas, whole, excised):
    designs = [r for r in atlas if r["group"] == "design"]
    natural = [r for r in atlas if r["group"] == "natural"]

    def subset(codes):
        return [r for r in natural if r["name"].split("_for_")[0] in codes]

    features = ["n_term_rsa", "c_term_rsa", "mean_rsa",
                "n_degron_density_weighted", "n_degron_density_raw",
                "all_degron_density_weighted", "all_degron_density_raw"]

    rng = np.random.default_rng(1)
    print("\n=== designs against natural, split by whether the entry is a whole protein ===")
    print(f"{'feature':<30}{'all':>8}{'whole only':>24}{'excised only':>24}")
    for feature in features:
        d = [r.get(feature, np.nan) for r in designs]
        row = f"{feature:<30}"
        row += f"{cliffs_delta(d, [r.get(feature, np.nan) for r in natural]):>8.3f}"
        for group in (subset(whole), subset(excised)):
            values = [r.get(feature, np.nan) for r in group]
            if len(values) > 5:
                delta = cliffs_delta(d, values)
                low, high = bootstrap_delta(d, values, rng)
                row += f"{delta:>+10.3f} [{low:+.2f},{high:+.2f}]"
            else:
                row += f"{'n/a':>24}"
        print(row)

    print(f"\n  n: all natural {len(natural)}, whole {len(subset(whole))}, "
          f"excised {len(subset(excised))}")
    print("\n  If the terminal-accessibility and degron effects hold on the")
    print("  whole-protein subset, the excision protocol does not explain them.")


if __name__ == "__main__":
    main()
