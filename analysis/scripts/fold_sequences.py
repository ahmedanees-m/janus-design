"""Fold sequences with ESMFold and record per-model confidence.

Input is a JSONL file with one object per line carrying "name" and "seq".
Each model is written as <out>/<name>.pdb and skipped if it already exists, so
an interrupted run resumes where it stopped. Sequences of equal length are
folded together, which avoids padding and keeps the language model's forward
pass from dominating the run.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import time
from pathlib import Path

import torch
from transformers import AutoTokenizer, EsmForProteinFolding

MODEL = "facebook/esmfold_v1"


def read_jsonl(path):
    out = []
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                record = json.loads(line)
                out.append((record["name"], record["seq"].strip().upper()))
    return out


def mean_plddt(pdb_text):
    values = [float(line[60:66]) for line in pdb_text.splitlines()
              if line.startswith("ATOM")]
    return sum(values) / len(values) if values else None


def batches(pending, size):
    by_length = collections.defaultdict(list)
    for name, sequence in pending:
        by_length[len(sequence)].append((name, sequence))
    for length in sorted(by_length):
        group = by_length[length]
        for start in range(0, len(group), size):
            yield group[start:start + size]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--recycles", type=int, default=3)
    parser.add_argument("--batch", type=int, default=8)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    designs = read_jsonl(args.sequences)
    pending = [(name, seq) for name, seq in designs
               if not (out / f"{name}.pdb").exists()]
    print(f"{len(designs)} sequences, {len(pending)} still to fold", flush=True)

    if pending:
        tokeniser = AutoTokenizer.from_pretrained(MODEL)
        model = EsmForProteinFolding.from_pretrained(MODEL, low_cpu_mem_usage=True)
        model = model.eval()
        # Halve the language model before it reaches the card: moving the whole
        # model across in single precision needs more memory than the card has.
        model.esm = model.esm.half()
        model = model.cuda()
        model.trunk.set_chunk_size(args.chunk_size)
        torch.backends.cuda.matmul.allow_tf32 = True
        print("model loaded", flush=True)

        done, started = 0, time.time()
        for group in batches(pending, args.batch):
            names = [name for name, _ in group]
            tokens = tokeniser([seq for _, seq in group], return_tensors="pt",
                               add_special_tokens=False)["input_ids"].cuda()
            with torch.no_grad():
                output = model(tokens, num_recycles=args.recycles)
            for name, text in zip(names, model.output_to_pdb(output)):
                (out / f"{name}.pdb").write_text(text, encoding="utf-8")
            done += len(group)
            if done % 200 < args.batch:
                rate = (time.time() - started) / done
                print(f"  {done}/{len(pending)}  {rate:.2f} s each", flush=True)

    records = []
    for name, sequence in designs:
        path = out / f"{name}.pdb"
        if not path.exists():
            continue
        records.append({"name": name, "length": len(sequence),
                        "mean_plddt": mean_plddt(path.read_text(encoding="utf-8"))})

    summary = Path(args.summary)
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["name", "length", "mean_plddt"])
        writer.writeheader()
        writer.writerows(records)

    folded = [r["mean_plddt"] for r in records if r["mean_plddt"] is not None]
    print(f"\n{len(records)} models, mean pLDDT {sum(folded) / len(folded):.1f}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
