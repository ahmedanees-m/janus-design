"""Sample residue sequences from a ProteinMPNN-format checkpoint.

The method comparison gives every stochastic model one draw per backbone, which
is the cheapest point on its curve and not the one a user would operate at. A
search method that spends thousands of objective evaluations should be compared
against a sampler allowed to spend the same, so this draws a batch per backbone
at a stated temperature and writes every draw. The scoring step then reads them
back and reports the best of the first n for each n.

The vanilla, soluble and developability-aligned checkpoints share a format and
load through the same path, so the arm is selected by checkpoint rather than by
code.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
BACKBONE_ATOMS = ("N", "CA", "C", "O")
PDBLIKE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def read_backbone(path):
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
        residues[number][atom] = [float(line[30:38]), float(line[38:46]),
                                  float(line[46:54])]
    coordinates = []
    for number in order:
        atoms = residues[number]
        if any(name not in atoms for name in BACKBONE_ATOMS):
            return None
        coordinates.append([atoms[name] for name in BACKBONE_ATOMS])
    return np.array(coordinates, dtype=np.float32)


def load_model(checkpoint, device):
    from protein_mpnn_utils import ProteinMPNN

    blob = torch.load(checkpoint, map_location="cpu")
    hidden = 128
    model = ProteinMPNN(num_letters=21, node_features=hidden,
                        edge_features=hidden, hidden_dim=hidden,
                        num_encoder_layers=3, num_decoder_layers=3,
                        augment_eps=0.0, k_neighbors=blob["num_edges"])
    state = blob.get("model_state_dict", blob.get("state_dict", blob))
    state = {k[len("model."):] if k.startswith("model.") else k: v
             for k, v in state.items()}
    model.load_state_dict(state)
    return model.to(device).eval()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--code", required=True,
                        help="directory holding protein_mpnn_utils.py")
    parser.add_argument("--backbones", required=True)
    parser.add_argument("--marginals", required=True,
                        help="only used to pick the same backbones the other arms use")
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--num-samples", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    sys.path.insert(0, args.code)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = load_model(args.checkpoint, device)
    print(f"{Path(args.checkpoint).name} on {device}, "
          f"{args.num_samples} draws at temperature {args.temperature}")

    names = [p.stem for p in sorted(Path(args.marginals).glob("*.npz"))
             if not PDBLIKE.match(p.stem)]
    if args.limit:
        names = names[: args.limit]

    omit = np.array([aa == "X" for aa in ALPHABET], dtype=np.float32)
    bias = np.zeros(len(ALPHABET), dtype=np.float32)

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
        batch = args.num_samples
        X = torch.from_numpy(coordinates)[None].repeat(batch, 1, 1, 1).to(device)
        ones = torch.ones(batch, length, device=device)
        residue_idx = torch.arange(length, device=device)[None].repeat(batch, 1)
        generator = torch.Generator(device="cpu").manual_seed(
            args.seed * 1000003 + number)
        randn = torch.randn(batch, length, generator=generator).to(device)

        elapsed = time.perf_counter()
        with torch.no_grad():
            output = model.sample(
                X, randn, torch.zeros(batch, length, dtype=torch.long, device=device),
                ones, ones, residue_idx, mask=ones,
                temperature=args.temperature, omit_AAs_np=omit, bias_AAs_np=bias,
                chain_M_pos=ones,
                bias_by_res=torch.zeros(batch, length, len(ALPHABET), device=device))
        seconds = (time.perf_counter() - elapsed) / batch
        sampled = output["S"].cpu().numpy()
        proteins = ["".join(ALPHABET[i] for i in row) for row in sampled]
        results[name] = {"proteins": proteins, "seconds": seconds,
                         "temperature": args.temperature}
        if number % 50 == 0:
            print(f"  {number}/{len(names)}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results), encoding="utf-8")
    total = time.perf_counter() - started
    print(f"\nwrote {out}: {len(results)} backbones, {skipped} skipped, "
          f"{total:.1f} s total")


if __name__ == "__main__":
    main()
