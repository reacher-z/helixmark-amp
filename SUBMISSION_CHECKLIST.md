# AMP Challenge submission checklist

HelixMark-AMP is technically ready, but joining and submitting on Kaggle binds
the participant to the competition-specific and Kaggle Foundational Rules.
Complete the identity, eligibility, rights, and consent items personally.

## Participant decisions and attestations

- [ ] Sign in to the intended single Kaggle account and verify that the account
  identity is accurate.
- [ ] Confirm that the participant is at least 18 or the age of majority in the
  relevant jurisdiction.
- [ ] Confirm residence/export-control/sanctions eligibility under the current
  Kaggle Foundational Rules.
- [ ] If entering for an employer, university, or other entity, obtain any
  required institutional knowledge and consent and verify that participation
  complies with its policies.
- [ ] Review and personally accept the competition-specific and Kaggle
  Foundational Rules by selecting **Join Hackathon**.
- [ ] Confirm the participant owns or has the rights needed to submit the code,
  model weights, training snapshot, generated outputs, and documentation.
- [ ] Decide whether the entry is individual or a team entry. Each team member
  must join through a separate Kaggle account before becoming an official team.

## Public-code compliance

- [ ] After joining, post the public
  [`reacher-z/helixmark-amp`](https://github.com/reacher-z/helixmark-amp)
  repository link in the competition Discussion or an associated Kaggle
  Notebook so it is available to all participants, as required by the
  Foundational Rules for public competition-code sharing.
- [ ] Keep the repository public under BSD-3-Clause and retain both copyright
  notices.
- [ ] Keep the AI-assistance disclosure, data provenance, training script,
  fitted checkpoint, `uv.lock`, and fixed default seed intact.

## Final technical preflight

Run from a clean clone:

```bash
uv sync --locked
uv run pytest -q
uv run ruff check .
uv run ruff format --check .
uv run generate
validation_root=$(mktemp -d)
uv run python scripts/verify_submission.py \
  https://github.com/reacher-z/helixmark-amp \
  --dir "$validation_root/submission"
```

- [ ] Confirm all GitHub Actions jobs are green on the exact final commit.
- [ ] Confirm `generate/library.fasta` has 50,000 unique records and
  `generate/top.fasta` has 100 records drawn from the library.
- [ ] Confirm the generated SHA-256 values match the preflight table in the
  README, or update the documented hashes after any intentional model change.
- [ ] Re-read the live Kaggle rules and deadline because the sponsor may update
  the timeline or requirements.

## Kaggle entry

- [ ] Upload the 50,000-sequence library and ranked top-100 in the format shown
  by the live competition submission interface.
- [ ] Provide the method abstract and the documented ranking procedure.
- [ ] Disclose `data/antibacterial.fasta` as the only training source and state
  that there was no proprietary data or manual peptide selection.
- [ ] Provide the public repository URL and any requested Kaggle access details
  for the organizers.
- [ ] Select the intended entry as the final submission and save the resulting
  Kaggle receipt/confirmation outside the generated-output directory.

Creating this repository does not constitute Kaggle entry, experimental-phase
qualification, acceptance of any rules, or a present co-authorship claim.
