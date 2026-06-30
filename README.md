# targetscan-python

A Python rewrite of the TargetScan miRNA target-prediction toolkit, plus
some tooling for moving its hg19-based data over to hg38.

[TargetScan](https://www.targetscan.org) predicts microRNA targets by
scanning 3' UTR sequences for miRNA seed matches and scoring how
conserved and effective each site looks. It was originally a set of Perl
scripts. This repo has two parts: a full Python port of that pipeline
(Chapter 1), and a toolkit for re-anchoring its hg19 data to hg38
(Chapter 2).

---

# Chapter 1: Perl to Python

Same algorithm, same output, no Perl. Every stage is checked byte-for-byte
against the original Perl output (see `tests/test_pipeline.py`).

## Pipeline stages

| Stage | Module | Replaces | What it does |
|---|---|---|---|
| 1 | `targetscan/site_prediction.py` | `targetscan_70.pl` | Finds 6mer/7mer/8mer miRNA seed-match sites in aligned 3' UTRs, groups overlapping cross-species sites |
| 2 | `targetscan/bl_bins.py` | `targetscan_70_BL_bins.pl` | Assigns each UTR to a branch-length conservation bin (1-10) |
| 3 | `targetscan/bl_pct.py` | `targetscan_70_BL_PCT.pl` | Computes branch length and probability of conserved targeting per site |
| 4 | `targetscan/count_8mers.py` | `targetscan_count_8mers.pl` | Counts 8mer sites in ORFs (feeds into context scoring) |
| 5 | `targetscan/context_scores.py` | `targetscan_70_context_scores.pl` | Computes context++ scores (Agarwal et al. 2015), including RNAplfold-based site accessibility |

`targetscan/phylo.py` is a shared helper (stages 2 and 3) that replaces
BioPerl's tree handling with Biopython.

Run all five stages at once with `run_pipeline.py`, or run them one at a
time with the scripts in `scripts/`.

## What changed going Perl to Python

- No Perl, no BioPerl. The only real dependency is `biopython`, used
  just for reading phylogenetic trees.
- One package instead of five loose scripts.
- A single orchestrator (`run_pipeline.py`) instead of chaining five
  Perl calls by hand.
- A test suite that pins every stage's output against the original
  tool's results, so future changes can't quietly drift.
- The algorithm itself didn't change: same seed-matching rules, same
  branch-length/PCT math, same context++ coefficients.

## Running it

### Full pipeline

```bash
python3 run_pipeline.py \
    --mirna-family miR_Family_info.txt \
    --utr UTR_Sequences_clean.txt \
    --orf ORF_Sequences_clean.txt \
    --mirna-context miR_for_context_scores.txt \
    --out-dir results/
```

Writes `predicted_targets.txt`, `utr_bl_bins.txt`,
`predicted_targets.bl_pct.txt`, `orf_lengths.txt`, `orf_8mer_counts.txt`,
and the final `context_scores.txt` into `results/`.

### One stage at a time

```bash
python3 scripts/targetscan_70.py miR_Family_info.txt UTR_Sequences_clean.txt predicted_targets.txt
python3 scripts/targetscan_70_BL_bins.py UTR_Sequences_clean.txt data/PCT_parameters/Tree.generic.txt utr_bl_bins.txt
python3 scripts/targetscan_70_BL_PCT.py miR_Family_info.txt predicted_targets.txt utr_bl_bins.txt data/PCT_parameters predicted_targets.bl_pct.txt
python3 scripts/targetscan_count_8mers.py miR_Family_info.txt ORF_Sequences_clean.txt orf_8mer_counts.txt
python3 scripts/targetscan_70_context_scores.py miR_for_context_scores.txt UTR_Sequences_clean.txt predicted_targets.bl_pct.txt orf_lengths.txt orf_8mer_counts.txt context_scores.txt data RNAplfold_in_out
```

### Input data

Standard TargetScan files: a miRNA family table (family, seed, species
list), aligned 3' UTR sequences (gene ID, species ID, aligned sequence),
and aligned ORF sequences in the same format for context scoring. Small
samples are in `samples/` if you just want to try it.

## Requirements

- Python 3.9+
- `biopython`
- `RNAplfold` (ViennaRNA Package 2) on `$PATH`, for the site-accessibility
  term in stage 5. If it's missing the pipeline still runs, that term
  just scores as 0.

```bash
pip install -r requirements.txt
```

## Testing

```bash
pip install pytest
pytest tests/ -v
```

Runs every stage on the bundled samples and checks the output. The
context-scores test uses cached RNAplfold output to get exact equality;
running with a freshly-invoked RNAplfold can shift the "SA contribution"
column slightly if your ViennaRNA version differs from whatever generated
the cached files. That's just a version difference in the external tool,
not a bug here.

---

# Chapter 2: hg19 to hg38

TargetScan's data is all hg19. This re-anchors it to hg38, first the
underlying alignment, then TargetScan's own published per-site BED files.
Every region gets tagged `ok`, `shifted`, `split`, or `failed` depending
on how confidently it mapped over. Only `ok` regions get touched.

## `hg38_liftover.py`

Re-anchors the human row of TargetScan's UTR alignment to hg38, leaving
the other ~80 species untouched. Uses the Ensembl REST API for coordinate
mapping and sequence lookup, so no genome FASTA or chain file needed.

```bash
python3 scripts/hg38_liftover.py hg19_3utr_regions.bed liftover_report.tsv \
    --utr-file UTR_Sequences_clean.txt --utr-out UTR_Sequences_hg38.txt --applied-out liftover_applied.tsv

# or against TargetScan's full coordinate file (resumable, ~1-1.5h for all ~28k transcripts)
python3 scripts/hg38_liftover.py --gff TSHuman_7_hg19_3UTRs.gff liftover_report.tsv --workers 4
```

Ran it against the real vert80 GFF: of 28,347 transcripts, 27,455 (97%)
came back clean. Multi-exon UTRs are handled by lifting each block and
concatenating in the right order.

## `hg38_annotate_sites.py`

The pipeline in Chapter 1 only reports a site's position as an offset
within its UTR. This converts that offset into a real `chrom:start-end`
using the liftover above.

```bash
python3 scripts/hg38_annotate_sites.py predicted_targets.txt liftover_report.tsv \
    TSHuman_7_hg19_3UTRs.gff predicted_targets.hg38_annotated.txt
```

## `hg38_liftover_sites.py`

TargetScan also publishes genomic coordinates for every predicted site
directly (e.g. `All_Target_Locations.hg19.bed.zip`, 8 files split by
conservation). Since those coordinates are already genomic and the
liftover above already worked out how each gene's region shifts between
assemblies, lifting them is just arithmetic with no new API calls.

```bash
python3 scripts/hg38_liftover_sites.py Gene_info.txt liftover_report.tsv in.bed out.bed
```

Ran this against all 8 real files: 12.3 million sites total, finished in
under a minute. About 86% came back `ok`. Most of the rest trace to genes
whose model shifted slightly between TargetScan's 2016 GFF and the 2021
BED export, which gets flagged rather than papered over.

`examples/real_hg38_demo/` ties all of this together on 6 real genes end
to end, if you want something small to look at first.

---

## Repo layout

```
targetscan_py/
  targetscan/             # the library: pipeline stages (Ch.1) + liftover tools (Ch.2)
  scripts/                # CLI wrappers, one per tool
  data/                   # model parameters + PCT trees (Ch.1; sample AIRs included)
  samples/                # small sample input files (Ch.1)
  examples/real_hg38_demo/  # real 6-gene example, hg19->hg38->pipeline end to end (Ch.2)
  tests/                  # regression tests for both chapters
  run_pipeline.py         # Ch.1 end-to-end orchestrator
```
