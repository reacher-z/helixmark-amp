"""Train the disclosed second-order Markov ensemble checkpoint."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
INDEX = {residue: index for index, residue in enumerate(ALPHABET)}


def read_fasta(path: Path) -> list[str]:
    sequences: list[str] = []
    parts: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if parts:
                sequences.append("".join(parts).upper())
            parts = []
        else:
            parts.append(line)
    if parts:
        sequences.append("".join(parts).upper())
    return sequences


def normalize_last_axis(values: np.ndarray) -> np.ndarray:
    return values / values.sum(axis=-1, keepdims=True)


def train(input_fasta: Path, output: Path, pseudocount: float = 0.5) -> int:
    size = len(ALPHABET)
    length_counts = np.full(43, pseudocount, dtype=float)
    first_counts = np.full(size, pseudocount, dtype=float)
    second_counts = np.full((size, size), pseudocount, dtype=float)
    transition_counts = np.full((size, size, size), pseudocount, dtype=float)
    valid = 0

    for sequence in read_fasta(input_fasta):
        if not 8 <= len(sequence) <= 50 or set(sequence) - set(ALPHABET):
            continue
        indices = [INDEX[residue] for residue in sequence]
        length_counts[len(sequence) - 8] += 1
        first_counts[indices[0]] += 1
        second_counts[indices[0], indices[1]] += 1
        for first, second, third in zip(indices, indices[1:], indices[2:]):
            transition_counts[first, second, third] += 1
        valid += 1

    if not valid:
        raise ValueError("training FASTA contains no valid sequences")

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        alphabet=np.array(list(ALPHABET)),
        length_probabilities=length_counts / length_counts.sum(),
        first_probabilities=first_counts / first_counts.sum(),
        second_probabilities=normalize_last_axis(second_counts),
        transition_probabilities=normalize_last_axis(transition_counts),
        training_sequence_count=np.array(valid),
        pseudocount=np.array(pseudocount),
    )
    return valid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/antibacterial.fasta"))
    parser.add_argument(
        "--output", type=Path, default=Path("checkpoint/markov_model.npz")
    )
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()
    count = train(args.input, args.output, args.pseudocount)
    print(f"Trained on {count} sequences -> {args.output}")


if __name__ == "__main__":
    main()
