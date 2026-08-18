"""Generate and rank HelixMark AMP candidates."""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import Levenshtein
import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AMINO_ACIDS = set(ALPHABET)
HYDROPHOBIC = set("AVILMFWY")
HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CHECKPOINT = PROJECT_ROOT / "checkpoint" / "markov_model.npz"
DEFAULT_REFERENCE = PROJECT_ROOT / "data" / "antibacterial.fasta"


@dataclass(frozen=True)
class SequenceFeatures:
    length: int
    net_charge: float
    hydrophobic_fraction: float
    hydrophobic_moment: float
    entropy: float
    proline_fraction: float
    glycine_fraction: float
    cysteines: int
    longest_run: int
    longest_hydrophobic_run: int


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


def _write_fasta(sequences: list[str], path: Path, prefix: str) -> None:
    with path.open("w") as handle:
        for rank, sequence in enumerate(sequences, start=1):
            handle.write(f">{prefix}_{rank:05d}\n{sequence}\n")


def _normalize_last_axis(values: np.ndarray) -> np.ndarray:
    return values / values.sum(axis=-1, keepdims=True)


def _profiled_probabilities(
    probabilities: np.ndarray, bias: np.ndarray, temperature: float
) -> np.ndarray:
    adjusted = np.power(probabilities, 1.0 / temperature) * bias
    return _normalize_last_axis(adjusted)


def _load_profiled_model(
    checkpoint: Path,
) -> tuple[np.ndarray, list[dict[str, np.ndarray]]]:
    with np.load(checkpoint) as model:
        length_probabilities = model["length_probabilities"]
        first = model["first_probabilities"]
        second = model["second_probabilities"]
        transition = model["transition_probabilities"]

    neutral = np.ones(len(ALPHABET))
    cationic = np.ones(len(ALPHABET))
    hydrophobic = np.ones(len(ALPHABET))
    for residue in "KR":
        cationic[ALPHABET.index(residue)] = 1.45
    cationic[ALPHABET.index("H")] = 1.12
    for residue in "DE":
        cationic[ALPHABET.index(residue)] = 0.68
    for residue in HYDROPHOBIC:
        hydrophobic[ALPHABET.index(residue)] = 1.28
    for residue in "KR":
        hydrophobic[ALPHABET.index(residue)] = 1.08

    profiles = []
    for bias, temperature in (
        (neutral, 0.92),
        (neutral, 1.08),
        (cationic, 1.0),
        (hydrophobic, 1.0),
    ):
        profiles.append(
            {
                "first": _profiled_probabilities(first, bias, temperature),
                "second": _profiled_probabilities(second, bias, temperature),
                "transition": _profiled_probabilities(transition, bias, temperature),
            }
        )
    return length_probabilities, profiles


def _sample_sequence(
    rng: np.random.Generator,
    length: int,
    profile: dict[str, np.ndarray],
) -> str:
    indices = [int(rng.choice(len(ALPHABET), p=profile["first"]))]
    if length > 1:
        indices.append(int(rng.choice(len(ALPHABET), p=profile["second"][indices[0]])))
    while len(indices) < length:
        indices.append(
            int(
                rng.choice(
                    len(ALPHABET),
                    p=profile["transition"][indices[-2], indices[-1]],
                )
            )
        )
    return "".join(ALPHABET[index] for index in indices)


def generate(
    n_sequences: int,
    *,
    seed: int = 42,
    length: int | None = None,
    checkpoint: Path = DEFAULT_CHECKPOINT,
    reference_fasta: Path = DEFAULT_REFERENCE,
) -> list[str]:
    """Generate unique, novel sequences from the disclosed Markov ensemble."""
    if n_sequences <= 0:
        raise ValueError("n_sequences must be positive")
    if length is not None and not 8 <= length <= 50:
        raise ValueError("length must be between 8 and 50")

    length_probabilities, profiles = _load_profiled_model(checkpoint)
    reference_sequences = set(read_fasta(reference_fasta))
    rng = np.random.default_rng(seed)
    lengths = np.arange(8, 51)
    sequences: list[str] = []
    seen: set[str] = set()
    max_attempts = n_sequences * 30

    for _ in range(max_attempts):
        sampled_length = (
            length
            if length is not None
            else int(rng.choice(lengths, p=length_probabilities))
        )
        profile = profiles[int(rng.choice(len(profiles), p=[0.36, 0.24, 0.24, 0.16]))]
        sequence = _sample_sequence(rng, sampled_length, profile)
        if sequence in seen or sequence in reference_sequences:
            continue
        if _longest_run(sequence, lambda left, right: left == right) > 4:
            continue
        seen.add(sequence)
        sequences.append(sequence)
        if len(sequences) == n_sequences:
            break

    if len(sequences) != n_sequences:
        raise RuntimeError(
            f"generated only {len(sequences)} unique sequences after {max_attempts} attempts"
        )

    return sequences


def _longest_run(sequence: str, predicate) -> int:
    longest = current = 0
    previous = ""
    for residue in sequence:
        if previous and predicate(previous, residue):
            current += 1
        else:
            current = 1
        longest = max(longest, current)
        previous = residue
    return longest


def _longest_hydrophobic_run(sequence: str) -> int:
    longest = current = 0
    for residue in sequence:
        current = current + 1 if residue in HYDROPHOBIC else 0
        longest = max(longest, current)
    return longest


def sequence_features(sequence: str) -> SequenceFeatures:
    length = len(sequence)
    charge = (
        sequence.count("K")
        + sequence.count("R")
        + 0.1 * sequence.count("H")
        - sequence.count("D")
        - sequence.count("E")
    )
    hydrophobic_fraction = sum(r in HYDROPHOBIC for r in sequence) / length

    x_component = 0.0
    y_component = 0.0
    for index, residue in enumerate(sequence):
        angle = math.radians(index * 100.0)
        scaled_hydropathy = (HYDROPATHY[residue] + 4.5) / 9.0
        x_component += scaled_hydropathy * math.cos(angle)
        y_component += scaled_hydropathy * math.sin(angle)
    hydrophobic_moment = math.hypot(x_component, y_component) / length

    counts = np.array([sequence.count(residue) for residue in ALPHABET], dtype=float)
    nonzero = counts[counts > 0] / length
    entropy = float(-(nonzero * np.log(nonzero)).sum() / math.log(len(ALPHABET)))

    return SequenceFeatures(
        length=length,
        net_charge=charge,
        hydrophobic_fraction=hydrophobic_fraction,
        hydrophobic_moment=hydrophobic_moment,
        entropy=entropy,
        proline_fraction=sequence.count("P") / length,
        glycine_fraction=sequence.count("G") / length,
        cysteines=sequence.count("C"),
        longest_run=_longest_run(sequence, lambda left, right: left == right),
        longest_hydrophobic_run=_longest_hydrophobic_run(sequence),
    )


def score(sequence: str) -> float:
    """Rank for cationicity, amphipathicity, diversity, and synthesizability."""
    features = sequence_features(sequence)
    value = 0.0
    value += 2.1 * math.exp(-(((features.net_charge - 4.8) / 3.1) ** 2))
    value += 1.8 * math.exp(-(((features.hydrophobic_fraction - 0.48) / 0.17) ** 2))
    value += 1.5 * min(features.hydrophobic_moment / 0.22, 1.25)
    value += 1.0 * math.exp(-(((features.length - 22.0) / 11.0) ** 2))
    value += 0.9 * min(features.entropy / 0.78, 1.0)
    value -= 1.2 * max(features.proline_fraction - 0.10, 0.0)
    value -= 0.8 * max(features.glycine_fraction - 0.18, 0.0)
    value -= 0.35 * max(features.cysteines - 1, 0)
    value -= 0.45 * max(features.longest_run - 2, 0)
    value -= 0.35 * max(features.longest_hydrophobic_run - 4, 0)
    if not 2.0 <= features.net_charge <= 9.0:
        value -= 1.0
    return value


def passes_synthesizability_filter(sequence: str) -> bool:
    """Apply the stringent residue filters published with HydrAMP."""
    if "C" in sequence:
        return False
    if any(
        window.count("K") + window.count("R") >= 3
        for window in (
            sequence[index : index + 5] for index in range(len(sequence) - 4)
        )
    ):
        return False
    return not any(
        first == second == third
        for first, second, third in zip(sequence, sequence[1:], sequence[2:])
    )


def _references_by_length(reference_fasta: Path) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for sequence in read_fasta(reference_fasta):
        grouped[len(sequence)].append(sequence)
    return grouped


def max_reference_similarity(
    sequence: str, references: dict[int, list[str]], threshold: float = 0.8
) -> float:
    """Return max ratio, pruning lengths that cannot exceed the threshold."""
    maximum = 0.0
    candidate_length = len(sequence)
    for reference_length, candidates in references.items():
        upper_bound = (
            2.0
            * min(candidate_length, reference_length)
            / (candidate_length + reference_length)
        )
        if upper_bound <= threshold:
            continue
        for reference in candidates:
            ratio = Levenshtein.ratio(sequence, reference)
            maximum = max(maximum, ratio)
            if maximum > threshold:
                return maximum
    return maximum


def select_top(
    sequences: list[str],
    *,
    top_k: int,
    reference_fasta: Path = DEFAULT_REFERENCE,
    novelty_threshold: float = 0.8,
    diversity_threshold: float = 0.75,
) -> list[str]:
    if not 0 < top_k <= len(sequences):
        raise ValueError("top_k must be positive and no larger than the library")

    references = _references_by_length(reference_fasta)
    ranked = sorted(sequences, key=lambda sequence: (-score(sequence), sequence))
    selected: list[str] = []
    for sequence in ranked:
        if not passes_synthesizability_filter(sequence):
            continue
        if (
            max_reference_similarity(sequence, references, novelty_threshold)
            > novelty_threshold
        ):
            continue
        if any(
            Levenshtein.ratio(sequence, existing) > diversity_threshold
            for existing in selected
        ):
            continue
        selected.append(sequence)
        if len(selected) == top_k:
            return selected

    raise RuntimeError(
        f"only {len(selected)} candidates passed novelty and diversity filters"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-sequences", type=int, default=50_000)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument(
        "--length",
        type=int,
        default=None,
        help="fixed sequence length; by default sample the learned length distribution",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REFERENCE)
    args = parser.parse_args()

    out_dir = Path(Path(sys.argv[0]).stem)
    out_dir.mkdir(parents=True, exist_ok=True)
    sequences = generate(
        args.n_sequences,
        seed=args.seed,
        length=args.length,
        checkpoint=args.checkpoint,
        reference_fasta=args.reference_fasta,
    )
    top_sequences = select_top(
        sequences,
        top_k=args.top_k,
        reference_fasta=args.reference_fasta,
    )

    library_path = out_dir / "library.fasta"
    top_path = out_dir / "top.fasta"
    _write_fasta(sequences, library_path, "helixmark_library")
    _write_fasta(top_sequences, top_path, "helixmark_rank")
    print(f"Generated {len(sequences)} sequences -> {library_path}")
    print(f"Selected {len(top_sequences)} candidates -> {top_path}")


if __name__ == "__main__":
    main()
