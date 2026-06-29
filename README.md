# targetscan-python

A Python rewrite of the original TargetScan miRNA target-prediction
toolkit, plus a toolkit for re-anchoring its hg19-based data to hg38.

[TargetScan](https://www.targetscan.org) predicts microRNA targets by
scanning 3' UTR sequences for miRNA seed matches and scoring how
conserved and effective each site is. This project has two parts:

1. **Perl to Python**: a faithful, byte-for-byte-verified rewrite of
   TargetScan's five-stage Perl pipeline, with no Perl/BioPerl dependency.
2. **hg19 to hg38**: a toolkit that re-anchors TargetScan's hg19-based
   data (both the underlying alignments and its published per-site BED
   coordinate files) to hg38, with explicit safety tagging so nothing is
   silently guessed.

---

# Chapter 1: Perl to Python

A faithful rewrite of TargetScan's Perl pipeline -- same algorithm, same
output, no Perl. Originally released as a set of Perl scripts; this
reimplements all of it natively in Python, with no Perl interpreter and
no BioPerl dependency anywhere in the chain.

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
pytest tests/ -v
```

These tests run every stage on the bundled sample data and check the
output. The context-scores test uses cached RNAplfold output to get exact
byte-for-byte equality; running it with a freshly-invoked RNAplfold may
produce tiny numeric differences in the SA contribution column if your
installed ViennaRNA version differs from the one used to generate the
cached results -- that's expected, not a bug.

---

# Chapter 2: hg19 to hg38

TargetScan's data is hg19-based. This chapter re-anchors it to hg38 in
two layers: first the underlying multi-species **alignment** (so the
prediction pipeline above can run on hg38-verified sequence), then
TargetScan's own published **per-site genomic coordinate files** (so
existing prediction results can be annotated with real hg38 locations
without re-running the pipeline at all).

Every tool here follows the same rule: only act where the mapping is
*verified safe* (identical sequence, unambiguous coordinates). Anything
else is tagged and left alone rather than guessed.

## Tags

| Tag | Meaning | Auto-applied? |
|---|---|---|
| `ok` | Clean 1:1 coordinate mapping, hg38 sequence byte-identical to hg19 | Yes |
| `shifted` | Mapped cleanly but the sequence differs (small edits/indels between assemblies) | Only if length still matches exactly; otherwise flagged |
| `split` | Region maps to more than one block in hg38 (e.g. a segmental-duplication/rearrangement area) | No -- flagged |
| `failed` | No mapping at all (e.g. centromeric/assembly-gap region) | No -- flagged |

## 2.1 Re-anchoring the alignment (`hg38_liftover.py`)

Re-anchors the **human row** of TargetScan's existing multi-species UTR
alignment to hg38, leaving the other ~80 species exactly as they were
(they stay on whatever assembly they were originally aligned to). Uses
the Ensembl REST API for both coordinate mapping and sequence lookup --
no local genome FASTA or UCSC chain file needed.

```bash
python3 scripts/hg38_liftover.py hg19_3utr_regions.bed liftover_report.tsv \
    --utr-file UTR_Sequences_clean.txt \
    --utr-out UTR_Sequences_hg38.txt \
    --applied-out liftover_applied.tsv
```

`hg19_3utr_regions.bed` is `gene_id, chrom, start, end, strand` (1-based
inclusive). About a quarter of human transcripts have a 3' UTR split
across multiple non-adjacent genomic blocks (multi-exon UTRs); each block
is lifted independently and a transcript's tag is the worst of its
blocks' tags, with hg38 sequence concatenated in correct 5'->3'
transcript order when all blocks are `ok` -- verified against real
production data in `tests/test_hg38_liftover.py`.

### Running against TargetScan's full coordinate file

```bash
python3 scripts/hg38_liftover.py --gff TSHuman_7_hg19_3UTRs.gff liftover_report.tsv --workers 4
```

Processes every transcript in TargetScan's hg19 3' UTR GFF (28,347
transcripts for vert80). Coordinate mapping runs concurrently, sequence
verification is batched, and results are written incrementally so the
run is resumable -- if interrupted, re-running the same command skips
transcripts already in the output file. A full run takes roughly
1-1.5 hours.

**Actual result on the real vert80 GFF (28,347 transcripts):**

| Tag | Count | % |
|---|---|---|
| `ok` | 27,455 | 96.9% |
| `shifted` | 475 | 1.7% |
| `split` | 148 | 0.5% |
| `failed` | 269 | 0.9% |

### Annotating predicted sites with real hg38 coordinates (`hg38_annotate_sites.py`)

`site_prediction.py` (Chapter 1) reports each site's position as a
UTR-relative offset, not a genome coordinate. This converts that offset
into a real `chrom:start-end` hg38 location, using the per-transcript
hg38 blocks the liftover above worked out. Sites spanning a splice
junction in a multi-exon UTR get multiple semicolon-separated blocks.
Only works for transcripts tagged `ok`; everything else gets `NA`.

```bash
python3 scripts/hg38_annotate_sites.py predicted_targets.txt liftover_report.tsv \
    TSHuman_7_hg19_3UTRs.gff predicted_targets.hg38_annotated.txt
```

Verified by hand against real predicted sites on both strands (e.g. a
NLRP1 site at UTR_start=2226 on this `-` strand gene's hg38 region
17:5499430-5501813 correctly lands at 17:5499583-5499588) -- see
`tests/test_hg38_annotate_sites.py`.

### Real-data example

`examples/real_hg38_demo/` runs the full chain against real TargetScan
vert80 data (not synthetic examples) for 6 genes: their actual hg19 3'
UTR coordinates, their actual 84-species alignment rows, liftover to
hg38, and the re-anchored result fed through the Chapter 1 pipeline end
to end. See `examples/real_hg38_demo/README.md` for exact commands.

## 2.2 Lifting TargetScan's own per-site BED files (`hg38_liftover_sites.py`)

TargetScan also publishes the genomic locations of every predicted site
directly, as hg19 BED files (`score` = context++ score percentile) --
e.g. `All_Target_Locations.hg19.bed.zip`, split into 8 files by miRNA
family conservation and site conservation. These give exact per-site
genomic coordinates already, so rather than re-deriving them from
UTR-relative offsets (2.1's approach), this lifts them directly.

The key trick: BED coordinates are plain genomic coordinates (not
transcript-relative), and 2.1 already verified, per transcript, exactly
how its hg19 UTR block(s) map to hg38. Within a block tagged `ok` that
mapping is a constant additive offset (no strand bookkeeping needed --
BED coordinates are always plus-strand on both assemblies). So every site
under an `ok`-tagged gene gets its hg38 coordinate by **pure arithmetic**,
reusing 2.1's already Ensembl-verified mapping -- no new API calls. Sites
under genes that weren't tagged `ok` get that same tag and no coordinate.

```bash
python3 scripts/hg38_liftover_sites.py Gene_info.txt liftover_report.tsv in.bed out.bed
```

A site can legitimately start in one GFF-recorded block and end in the
next (TargetScan's GFF splits annotations by feature type -- e.g. "UTR"
vs "ORF" -- even where they're genomically contiguous), so blocks aren't
required to fully contain a site; the start and end positions are
resolved independently and only fail if they disagree (a real
`site_spans_inconsistent_blocks`, meaning the two ends sit in genuinely
differently-shifted regions) -- verified against a real multi-block case
(MCM10) in `tests/test_hg38_liftover_sites.py`.

**Actual result, all 8 real vert80 files (12,344,655 total predicted sites):**

| File | Total sites | `ok` | % ok |
|---|---|---|---|
| broadConsFam.consSite | 81,886 | 75,735 | 92.5% |
| broadConsFam.nonConsSite | 531,508 | 453,364 | 85.3% |
| consFam.consSite | 44,813 | 41,063 | 91.6% |
| consFam.nonConsSite | 624,224 | 534,836 | 85.7% |
| nonConsFam.consSite | 1,060 | 958 | 90.4% |
| nonConsFam.nonConsSite | 1,282,120 | 1,102,298 | 86.0% |
| otherFam.consSite | 397 | 372 | 93.7% |
| otherFam.nonConsSite | 9,778,647 | 8,391,695 | 85.8% |
| **Total** | **12,344,655** | **10,600,321** | **85.9%** |

The ~14% not tagged `ok` break down as: site outside any block TargetScan's
GFF actually covers (10.3%, mostly a small number of sites whose
underlying gene model differs between this BED file's 2021 vintage and
the 2016 GFF used for the transcript-level liftover -- a genuine data
versioning mismatch, not a bug, hence flagged rather than guessed),
`shifted` (2.1%), `split` (1.0%), `failed` (0.7%), and a handful of
`no_gene_mapping`/`no_liftover_data` (<0.02%).

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

## Roadmap

- The `site_outside_lifted_blocks` residual in 2.2 traces to gene-model
  version drift between TargetScan's 2016 GFF and 2021 BED exports;
  refreshing the transcript-level liftover (2.1) against a matching-vintage
  GFF would close most of that gap.
