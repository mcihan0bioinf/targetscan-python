# Real-data hg38 demo

Real TargetScan vert80 data for 6 genes (84-species alignment rows each),
extracted from the official downloads at
https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
(`UTR_Sequences.txt.zip`, `ORF_Sequences.txt.zip`, `Gene_info.txt.zip`,
`TSHuman_7_hg19_3UTRs.gff.zip`), used here to demonstrate the hg19->hg38
liftover against genuine production data rather than synthetic examples.

Genes (old symbol -> current symbol, since some were renamed since the
data was generated): CDC2L6 -> **CDK19**, FNDC3A, LIN28B, LPHN1, NLRP1,
ZNF197.

## Files

- `hg19_3utr_regions.bed` -- each gene's hg19 3' UTR genomic coordinates
  (from `TSHuman_7_hg19_3UTRs.gff`, cross-referenced to gene symbols via
  `Gene_info.txt`'s representative transcript)
- `UTR_Sequences_6genes.txt` -- the real 84-species aligned 3' UTRs for
  these 6 genes (504 rows = 6 genes x ~84 species), straight from
  `UTR_Sequences.txt`
- `ORF_Sequences_6genes.txt` -- matching ORF alignments, from
  `ORF_Sequences.txt`
- `liftover_report.tsv` -- the actual liftover result: all 6 regions
  tag `ok` (clean 1:1 hg19->hg38 mapping, byte-identical sequence)

## Reproducing

```bash
# 1. Lift hg19 coordinates to hg38 and verify they're safe to re-anchor
python3 ../../scripts/hg38_liftover.py hg19_3utr_regions.bed liftover_report.tsv \
    --utr-file UTR_Sequences_6genes.txt \
    --utr-out UTR_Sequences_6genes.hg38.txt \
    --applied-out liftover_applied.tsv

# 2. Run the full prediction pipeline on the hg38-verified alignment
python3 ../../run_pipeline.py \
    --mirna-family ../../samples/miR_Family_info_sample.txt \
    --utr UTR_Sequences_6genes.hg38.txt \
    --orf ORF_Sequences_6genes.txt \
    --mirna-context ../../samples/miR_for_context_scores.sample.txt \
    --out-dir out/
```

Step 1 confirms all 6 genes' human rows are safe to re-anchor (sequence
unchanged between hg19 and hg38 at these loci -- `liftover_applied.tsv`
will say "applied (sequence unchanged...)" for every gene). Step 2 then
runs site prediction, conservation scoring, and context++ scoring on
that data end to end, producing real predicted targets and scores for
these 6 genes in `out/context_scores.txt`.

Note: the bundled `samples/miR_Family_info_sample.txt` and AIRs reference
data are generic/demo files, not specific to these particular genes, so
AIR lookups will fall back to a default of 1 for them (a documented
fallback, not an error) -- this only matters for the AIR-weighted context
score terms, not the liftover step itself.
