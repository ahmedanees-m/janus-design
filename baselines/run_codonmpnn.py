"""Generate CodonMPNN coding sequences for a set of backbones.

CodonMPNN is the baseline a reviewer is certain to ask about: a codon-level
inverse folding model that chooses residues and codons together from structure,
conditioned on an organism. It is the closest published thing to what this
project does, and the comparison is the point.

The published repository has weights and a training script but no inference
entry point: `likelihood_eval.py` carries absolute paths to the authors' cluster
and evaluates likelihood on one specific dataset class. Their `ProteinMPNN`
module, however, depends only on torch and numpy and exposes the usual
autoregressive `sample`, so this drives that directly and skips the Lightning
wrapper, openfold and wandb that the rest of the repository pulls in.

Two facts had to be recovered rather than read from documentation. The codon
vocabulary is 65 wide, 64 codons and an unknown, ordered by `codon_order`. And
the organism conditioning is not an NCBI taxon id but an index into a 20000-way
grouping built for training; *E. coli*, NCBI 562, sits at index 1893, read off
the grouping table the authors publish alongside the weights.

Sampling is at low temperature rather than argmax, because the decoder is
autoregressive over a random decoding order and a greedy pass is not the model's
own estimate of its best sequence. The decoding order is seeded so the run
reproduces.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch

from codon.utils.codon_const import codon_order
from codon.utils.pmpnn import ProteinMPNN, get_weird_pmpnn_stuff

BACKBONE_ATOMS = ("N", "CA", "C", "O")
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def read_backbone(path):
    """Backbone coordinates in ProteinMPNN's N, CA, C, O order."""
    residues: dict[int, dict[str, list[float]]] = {}
    order: list[int] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom not in BACKBONE_ATOMS:
            continue
        number = int(line[22:26])
        if number not in residues:
            residues[number] = {}
            order.append(number)
        residues[number][atom] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

    coordinates = []
    for number in order:
        atoms = residues[number]
        if any(name not in atoms for name in BACKBONE_ATOMS):
            return None
        coordinates.append([atoms[name] for name in BACKBONE_ATOMS])
    return np.array(coordinates, dtype=np.float32)


def load_model(checkpoint, device):
    blob = torch.load(checkpoint, map_location="cpu")
    hyper = blob.get("hyper_parameters", {})
    args = SimpleNamespace(
        taxon_condition=hyper.get("taxon_condition", True),
        num_taxon_ids=hyper.get("num_taxon_ids", 20000),
    )
    weights = {k[len("model."):]: v for k, v in blob["state_dict"].items()
               if k.startswith("model.")}
    model = ProteinMPNN(
        args,
        hidden_dim=hyper.get("hidden_dim", 128),
        num_encoder_layers=hyper.get("num_encoder_layers", 3),
        num_decoder_layers=hyper.get("num_decoder_layers", 3),
        vocab=weights["W_out.weight"].shape[0],
        k_neighbors=hyper.get("num_neighbors", 48),
        augment_eps=0.0,
        dropout=hyper.get("dropout", 0.1),
    )
    missing, unexpected = model.load_state_dict(weights, strict=False)
    if missing or unexpected:
        print(f"  {len(missing)} missing and {len(unexpected)} unexpected parameters")
        for name in list(missing)[:5] + list(unexpected)[:5]:
            print(f"    {name}")
    return model.to(device).eval(), args


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--marginals", required=True,
                        help="only used to pick the same backbones the other arms use")
    parser.add_argument("--taxon", type=int, default=1893,
                        help="grouping index for the organism, 1893 is E. coli")
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, model_args = load_model(args.checkpoint, device)
    codons = list(codon_order)
    print(f"vocabulary of {len(codons)} codons, taxon index {args.taxon}, on {device}")

    names = [p.stem for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        names = names[: args.limit]

    generator = torch.Generator(device="cpu").manual_seed(args.seed)
    results, skipped = {}, 0
    started = time.perf_counter()
    for number, name in enumerate(names, start=1):
        path = Path(args.backbones) / f"{name}.pdb"
        if not path.exists():
            skipped += 1
            continue
        coordinates = read_backbone(path)
        if coordinates is None:
            skipped += 1
            continue

        length = len(coordinates)
        X = torch.from_numpy(coordinates)[None].to(device)
        mask = torch.ones(1, length, device=device)
        chain_mask = torch.ones(1, length, device=device)
        chain_M_pos = torch.ones(1, length, device=device)
        residue_idx, chain_encoding = get_weird_pmpnn_stuff(torch.zeros(length, dtype=torch.long))
        randn = torch.randn(1, length, generator=generator).to(device)
        taxon = torch.tensor([args.taxon], dtype=torch.long, device=device)

        elapsed = time.perf_counter()
        with torch.no_grad():
            output = model.sample(
                X, randn, torch.zeros(1, length, dtype=torch.long, device=device),
                taxon, chain_mask, chain_encoding[None].to(device),
                residue_idx[None].to(device), mask, args.temperature,
                np.zeros(model.vocab, dtype=np.float32),
                np.zeros(model.vocab, dtype=np.float32),
                chain_M_pos, torch.zeros(1, length, model.vocab, device=device),
            )
        sampled = output["S"] if isinstance(output, dict) else output[0]
        cds = "".join(codons[i] for i in sampled[0].tolist())
        results[name] = {"cds": cds, "seconds": time.perf_counter() - elapsed}
        if number % 50 == 0:
            print(f"  {number}/{len(names)}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results), encoding="utf-8")
    total = time.perf_counter() - started
    print(f"\nwrote {out}: {len(results)} sequences, {skipped} skipped, "
          f"{total:.1f} s total, {total / max(len(results), 1):.3f} s each")


if __name__ == "__main__":
    main()
