"""Assemble a natural control arm of whole proteins, not excised domains.

Selection is by annotation rather than by hand. Reviewed Swiss-Prot entries of 26
to 74 residues, protein-level evidence, not fragments, with an AlphaFold model,
and carrying no signal peptide, propeptide, transit peptide or peptide feature.
The cleavage filters are what make them whole: an entry with a cleaved region has
a precursor sequence and a different mature protein, and slicing one out of the
other would reintroduce the excision problem at a smaller scale.

Transmembrane and intramembrane annotations are excluded as well. Without that
the band fills with membrane microproteins, whose hydrophobic helices are buried
in a bilayer the monomeric model does not contain, so they would read as exposed
hydrophobic surface.

Structures come from AlphaFold DB so that accessibility is computed the same way
as for the designed backbones, which are themselves AlphaFold models."""

from __future__ import annotations

import argparse
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

FILTERS = (
    "(reviewed:true) AND (existence:1) AND (fragment:false)"
    " AND (database:alphafolddb)"
    " NOT (ft_signal:*) NOT (ft_propep:*) NOT (ft_transit:*) NOT (ft_peptide:*)"
    " NOT (ft_transmem:*) NOT (ft_intramem:*)"
)
QUERY = f"(length:[26 TO 74]) AND {FILTERS}"
FIELDS = "accession,length,sequence,organism_name,protein_name,cc_subcellular_location"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
MODEL = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def get(url, retries=4, pause=2.0):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read(), response.headers
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            if attempt == retries - 1:
                raise
            reason = getattr(error, "code", error)
            print(f"    retry {attempt + 1} after {reason}", flush=True)
            time.sleep(pause * (attempt + 1))
    raise RuntimeError("unreachable")


def search(query, fields, limit):
    """Walk the cursor pagination and return one row per entry."""
    url = f"{SEARCH}?" + urllib.parse.urlencode(
        {"query": query, "fields": fields, "format": "tsv", "size": 500}
    )
    rows, header = [], None
    while url and (limit is None or len(rows) < limit):
        body, headers = get(url)
        lines = body.decode("utf-8").splitlines()
        if header is None:
            header = lines[0].split("\t")
        rows.extend(dict(zip(header, line.split("\t"), strict=False)) for line in lines[1:])
        print(f"  {len(rows)} entries", flush=True)
        link = headers.get("Link", "")
        url = link.split(";")[0].strip("<> ") if 'rel="next"' in link else None
    return rows[:limit] if limit else rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="directory for the models")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--version", type=int, default=6, help="AlphaFold DB file version")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--pause", type=float, default=0.05)
    parser.add_argument("--length", nargs=2, type=int, metavar=("LOW", "HIGH"),
                        help="residue range, in place of the default 26 to 74")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    query = QUERY
    if args.length:
        query = f"(length:[{args.length[0]} TO {args.length[1]}]) AND {FILTERS}"
    print("querying UniProt")
    entries = search(query, FIELDS, args.limit)

    kept, dropped = [], 0
    for entry in entries:
        sequence = entry.get("Sequence", "")
        if not sequence or set(sequence) - STANDARD:
            dropped += 1
            continue
        kept.append(entry)
    if args.limit:
        kept = kept[: args.limit]
    print(f"{len(kept)} entries with standard residues only, {dropped} dropped")

    written, cached, failed = 0, 0, []
    for number, entry in enumerate(kept, start=1):
        accession = entry["Entry"]
        target = out / f"{accession}.pdb"
        if target.exists() and target.stat().st_size > 0:
            cached += 1
            continue
        try:
            body, _ = get(MODEL.format(accession=accession, version=args.version))
        except Exception as error:
            failed.append((accession, str(error)))
            continue
        target.write_bytes(body)
        written += 1
        if args.pause:
            time.sleep(args.pause)
        if number % 200 == 0:
            print(f"  {number}/{len(kept)} models, {written} new, {len(failed)} failed",
                  flush=True)

    have = {p.stem for p in out.glob("*.pdb")}
    lines = ["accession\tlength\tsequence\torganism\tprotein"]
    for entry in kept:
        if entry["Entry"] not in have:
            continue
        lines.append("\t".join([
            entry["Entry"], entry.get("Length", ""), entry["Sequence"],
            entry.get("Organism", "").replace("\t", " "),
            entry.get("Protein names", "").replace("\t", " "),
        ]))
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\n{len(lines) - 1} whole natural proteins with models "
          f"({written} downloaded, {cached} already present)")
    if failed:
        print(f"{len(failed)} models could not be fetched, first few:")
        for accession, reason in failed[:5]:
            print(f"  {accession} {reason}")
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
