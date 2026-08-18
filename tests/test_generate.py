from pathlib import Path

import Levenshtein

from amp_challenge_2027.generate import (
    ALPHABET,
    DEFAULT_REFERENCE,
    _references_by_length,
    generate,
    max_reference_similarity,
    passes_synthesizability_filter,
    score,
    select_top,
    sequence_features,
)


def test_generation_is_reproducible_and_valid():
    first = generate(128, seed=7)
    second = generate(128, seed=7)

    assert first == second
    assert len(first) == len(set(first)) == 128
    assert all(8 <= len(sequence) <= 50 for sequence in first)
    assert all(set(sequence) <= set(ALPHABET) for sequence in first)


def test_fixed_length_is_supported():
    sequences = generate(32, seed=11, length=19)
    assert {len(sequence) for sequence in sequences} == {19}


def test_feature_score_prefers_balanced_candidate():
    balanced = "KWKLLKALKKLAGWAL"
    acidic = "DEDEDEDEDEDEDEDE"
    features = sequence_features(balanced)

    assert features.net_charge > 2
    assert 0.3 < features.hydrophobic_fraction < 0.8
    assert score(balanced) > score(acidic)


def test_synthesizability_filter_rejects_published_failure_patterns():
    assert passes_synthesizability_filter("KWLALQRLAMG")
    assert not passes_synthesizability_filter("KWLACQRLAMG")
    assert not passes_synthesizability_filter("KKRLAQGLAMG")
    assert not passes_synthesizability_filter("KWLAAAQLAMG")


def test_max_similarity_prunes_impossible_lengths():
    references = {50: ["A" * 50], 10: ["ACDEFGHIKL"]}
    assert max_reference_similarity("ACDEFGHIKL", references) == 1.0
    assert max_reference_similarity("W" * 8, {50: ["W" * 50]}) == 0.0


def test_top_candidates_meet_novelty_and_diversity_rules():
    sequences = generate(600, seed=23)
    top = select_top(sequences, top_k=8)
    references = _references_by_length(DEFAULT_REFERENCE)

    assert len(top) == 8
    assert all(passes_synthesizability_filter(sequence) for sequence in top)
    assert all(
        max_reference_similarity(sequence, references) <= 0.8 for sequence in top
    )
    assert all(
        Levenshtein.ratio(left, right) <= 0.75
        for index, left in enumerate(top)
        for right in top[index + 1 :]
    )


def test_checkpoint_and_training_data_are_disclosed():
    root = Path(__file__).resolve().parents[1]
    assert (root / "checkpoint" / "markov_model.npz").is_file()
    assert (root / "data" / "antibacterial.fasta").is_file()
