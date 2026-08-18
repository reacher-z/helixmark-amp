# HelixMark-AMP

HelixMark-AMP is a CPU-friendly, fully reproducible submission for the
[AMP Challenge 2027](https://szczurek-lab.github.io/amp-challenge-website/).
It learns local sequence structure and the length distribution of known
antibacterial peptides, samples a deliberately mixed library, and ranks novel
candidates using transparent biophysical objectives.

## Abstract

We train a smoothed second-order Markov model on the challenge-provided MarLys
antibacterial reference snapshot. Four sampling profiles combine a conservative
temperature, a diversity-increasing temperature, a cationic bias, and a
hydrophobic bias. The mixture preserves broad coverage while enriching motifs
associated with antimicrobial activity. Candidate ranking combines approximate
net charge, hydrophobic fraction, alpha-helical hydrophobic moment, sequence
entropy, and length, with penalties for excessive proline, glycine, cysteine,
homopolymers, and hydrophobic runs. The final list applies the published HydrAMP
synthesizability filter, excludes candidates above the challenge's 0.80
Levenshtein-ratio novelty ceiling, and limits pairwise similarity within the top
100 to 0.75.

## Reproduce the submission

Python 3.10+ and [`uv`](https://docs.astral.sh/uv/) are required.

```bash
uv sync --locked
uv run generate
```

This creates:

```text
generate/
├── library.fasta  # 50,000 unique sequences
└── top.fasta      # ranked top 100
```

The default seed is 42. Repeating `uv run generate` produces byte-identical
files. Useful development commands are:

```bash
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run python scripts/train_model.py
```

Retraining is not needed to reproduce the submitted library; the fitted
`checkpoint/markov_model.npz` is committed. The training command is provided so
that every weight can also be regenerated from the disclosed input.

## Method

### Training

The model uses maximum-likelihood residue, first-order, and second-order
transition counts with a 0.5 pseudocount. It also learns the empirical length
distribution over the challenge's 8-50-residue interval. All probability arrays
are stored in `checkpoint/markov_model.npz`.

### Library generation

Each sequence first receives one of four fixed sampling profiles:

| Profile | Weight | Purpose |
| --- | ---: | --- |
| conservative (`T=0.92`) | 0.36 | preserve frequent AMP motifs |
| exploratory (`T=1.08`) | 0.24 | expand sequence-space coverage |
| cationic | 0.24 | enrich K/R and suppress D/E |
| hydrophobic | 0.16 | enrich membrane-interacting residues |

Exact training/reference matches, duplicate sequences, and homopolymer runs
longer than four are rejected. No manual sequence curation is performed.

### Candidate ranking

Ranking is deterministic and favors the following soft targets:

- approximate net charge near +4.8;
- hydrophobic fraction near 0.48;
- high alpha-helical hydrophobic moment (100-degree residue spacing);
- length near 22 residues and normalized residue entropy near 0.78;
- low proline/glycine burden, at most one unpenalized cysteine, and short
  homopolymer/hydrophobic runs.

Candidates are visited in score order. Before novelty checks, the stringent
HydrAMP residue filter from Szymczak et al. (Nature Communications 14, 1453,
2023) rejects cysteine-containing sequences, any sequence with three or more K/R
residues in a five-residue window, and any three-residue homopolymer. A candidate
is then admitted only when its Levenshtein ratio is at most 0.80 against every
challenge reference and at most 0.75 against every already selected top
candidate. This makes all 100 entries eligible for uniform experimental sampling
while retaining internal diversity.

## Training data and provenance

The only training data are `data/antibacterial.fasta`, copied unchanged from the
official AMP Challenge starter repository. It is the challenge's MarLys-AMP
reference snapshot (39,448 FASTA records, all valid in the required
8-50-residue standard alphabet). The upstream record is:

> Bogdan Marczak, Aleksandra Bocian, and Andrzej Lyskowski. *MarLys AMP
> database - MLAMP_db* (2026). <https://doi.org/10.17632/w4hb5grjwb.3>

No proprietary data, external database, pretrained language model, wet-lab
label, or manually selected peptide is used. The bundled snapshot and template
were redistributed under the starter repository's BSD-3-Clause license; its
original copyright notice is retained in `LICENSE`.

## AI-assistance disclosure

OpenAI Codex was used to help design and implement the generator, tests, and
documentation. The human participant remains responsible for the submission,
scientific claims, provenance, and competition compliance. Generated peptide
sequences are produced only by the disclosed Markov model and deterministic
ranking code; an LLM does not generate or manually select any submitted peptide.

## Limitations

This is a transparent statistical baseline, not evidence that any generated
peptide is safe or clinically effective. The property score is heuristic and
cannot replace toxicity, stability, synthesis, or antimicrobial assays. The
sequences must not be used in humans, animals, or the environment outside an
appropriately reviewed research protocol.

## License

BSD-3-Clause. See `LICENSE`.
