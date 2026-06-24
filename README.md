# targetscan-python

A Python rewrite of the original TargetScan miRNA target-prediction
toolkit -- same algorithm, same output, no Perl.

[TargetScan](https://www.targetscan.org) predicts microRNA targets by
scanning 3' UTR sequences for miRNA seed matches and scoring how
conserved and effective each site is. It was originally released as a
set of Perl scripts; this project reimplements all of it natively in
Python, with no Perl interpreter and no BioPerl dependency anywhere in
the chain.

**Validated for exact compatibility**: every stage's output has been
checked byte-for-byte against the original Perl scripts on sample data
(see `tests/test_pipeline.py`). A few quirky edge-case behaviors from the
original Perl are intentionally kept as-is (documented inline in the
code), since the goal is a drop-in replacement, not a "corrected"
reimplementation.

## Pipeline stages

| Stage | Module | Replaces | Purpose |
|---|---|---|---|
| 1 | `targetscan/site_prediction.py` | `targetscan_70.pl` | Find 6mer/7mer/8mer miRNA seed-match sites in aligned 3' UTRs, group overlapping cross-species sites |
| 2 | `targetscan/bl_bins.py` | `targetscan_70_BL_bins.pl` | Assign each UTR to a branch-length conservation bin (1-10) |
| 3 | `targetscan/bl_pct.py` | `targetscan_70_BL_PCT.pl` | Compute branch length (BL) and probability of conserved targeting (PCT) per site |
| 4 | `targetscan/count_8mers.py` | `targetscan_count_8mers.pl` | Count 8mer sites in ORFs (input to context scoring) |
| 5 | `targetscan/context_scores.py` | `targetscan_70_context_scores.pl` | Compute context++ scores (Agarwal et al. 2015 model), including RNAplfold-based site accessibility |

`targetscan/phylo.py` is a shared helper (used by stages 2 and 3) that
replaces BioPerl's tree handling with Biopython.

Run all five stages with one command via `run_pipeline.py`, or run each
stage individually with the scripts in `scripts/`.

## What changed going from Perl to Python

- No Perl runtime, no BioPerl. The only non-stdlib dependency is
  `biopython`, used for reading phylogenetic trees.
- One package (`targetscan/`) instead of five standalone scripts.
- A single end-to-end orchestrator (`run_pipeline.py`).
- A test suite (`tests/test_pipeline.py`) that checks every stage's
  output against the original tool's results, so changes can't quietly
  drift from the original behavior.
- The algorithm and scoring model themselves are unchanged -- same
  seed-matching rules, same branch-length/PCT math, same context++
  coefficients.

## Running

### Full pipeline

```bash
python3 run_pipeline.py \
    --mirna-family miR_Family_info.txt \
    --utr UTR_Sequences_clean.txt \
    --orf ORF_Sequences_clean.txt \
    --mirna-context miR_for_context_scores.txt \
    --out-dir results/
```

This writes `predicted_targets.txt`, `utr_bl_bins.txt`,
`predicted_targets.bl_pct.txt`, `orf_lengths.txt`, `orf_8mer_counts.txt`,
and the final `context_scores.txt` into `results/`.

### Individual stages

```bash
python3 scripts/targetscan_70.py miR_Family_info.txt UTR_Sequences_clean.txt predicted_targets.txt
python3 scripts/targetscan_70_BL_bins.py UTR_Sequences_clean.txt data/PCT_parameters/Tree.generic.txt utr_bl_bins.txt
python3 scripts/targetscan_70_BL_PCT.py miR_Family_info.txt predicted_targets.txt utr_bl_bins.txt data/PCT_parameters predicted_targets.bl_pct.txt
python3 scripts/targetscan_count_8mers.py miR_Family_info.txt ORF_Sequences_clean.txt orf_8mer_counts.txt
python3 scripts/targetscan_70_context_scores.py miR_for_context_scores.txt UTR_Sequences_clean.txt predicted_targets.bl_pct.txt orf_lengths.txt orf_8mer_counts.txt context_scores.txt data RNAplfold_in_out
```

### Input data

The pipeline expects the standard TargetScan input files: a miRNA family
table (family, seed, species list), aligned 3' UTR sequences (gene ID,
species ID, aligned sequence), and aligned ORF sequences in the same
format (used for context scoring). Small sample files are included in
`samples/` so you can try the pipeline right away.

## Requirements

- Python 3.9+
- `biopython`
- `RNAplfold` (from the ViennaRNA Package 2) on `$PATH`, for stage 5's
  site-accessibility contribution. If it's missing, the pipeline still
  runs, just scoring that contribution as 0.

```bash
pip install -r requirements.txt
```

## Testing

```bash
pip install pytest
pytest tests/test_pipeline.py -v
```

These tests run every stage on the bundled sample data and check the
output. Stage 5's test uses cached RNAplfold output to get exact
byte-for-byte equality; running it with a freshly-invoked RNAplfold may
produce tiny numeric differences in the SA contribution column if your
installed ViennaRNA version differs from the one used to generate the
cached results -- that's expected, not a bug.

## Repo layout

```
targetscan_py/
  targetscan/             # the library (one module per pipeline stage)
  scripts/                # CLI wrappers, one per pipeline stage
  data/                   # model parameters + PCT trees (sample AIRs included)
  samples/                # small sample input files
  tests/                  # regression tests
  run_pipeline.py         # end-to-end orchestrator
```

## Roadmap

- Adapt and validate the pipeline for hg38-based data (next step).
