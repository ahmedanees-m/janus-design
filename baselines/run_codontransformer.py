"""Generate CodonTransformer coding sequences for a set of proteins.

Runs in its own image, since CodonTransformer pins numpy below 2 and
transformers below 4.50. It writes coding sequences and nothing else; scoring
happens in the main container so every baseline is scored by the same code.

CodonTransformer chooses codons for a given protein, so it is a codon-layer
method and is given the design's own residue sequence, exactly as the vendor
baseline is. Comparing it against the joint arm without saying that would be
comparing two different problems.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import torch
from CodonTransformer.CodonPrediction import predict_dna_sequence
from transformers import AutoTokenizer, BigBirdForMaskedLM

CHECKPOINT = "adibvafa/CodonTransformer"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True, help="JSON mapping name to sequence")
    parser.add_argument("--organism", default="Escherichia coli general")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    proteins = json.loads(Path(args.proteins).read_text(encoding="utf-8"))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(CHECKPOINT)
    model = BigBirdForMaskedLM.from_pretrained(CHECKPOINT).to(device)
    print(f"{len(proteins)} proteins on {device}", flush=True)

    results = {}
    started = time.perf_counter()
    for number, (name, protein) in enumerate(sorted(proteins.items()), start=1):
        elapsed = time.perf_counter()
        output = predict_dna_sequence(
            protein=protein, organism=args.organism, device=device,
            tokenizer=tokenizer, model=model, deterministic=True,
            match_protein=True,
        )
        results[name] = {"cds": output.predicted_dna,
                         "seconds": time.perf_counter() - elapsed}
        if number % 25 == 0:
            print(f"  {number}/{len(proteins)}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results), encoding="utf-8")
    total = time.perf_counter() - started
    print(f"\nwrote {out}, {total:.1f} s total, "
          f"{total / max(len(results), 1):.2f} s per sequence")


if __name__ == "__main__":
    main()
