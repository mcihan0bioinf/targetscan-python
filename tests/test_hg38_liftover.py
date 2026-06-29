"""Tests for the hg19->hg38 liftover/tagging helper.

These hit the live Ensembl REST API with a handful of real genomic
regions chosen to exercise each tag:
  - chr1:1,000,000-1,000,100   -> "ok" (clean 1:1 mapping, identical sequence)
  - chr1:121,500,000-121,500,100 -> "failed" (centromeric, no mapping in hg38)
  - chr1:144,080,000-144,082,000 -> "split" (segmental-duplication region,
    maps to 2 separate blocks in hg38)

Skipped automatically if there's no network access.
"""

import os
import sys
import urllib.request

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import tempfile

from targetscan.hg38_liftover import (  # noqa: E402
    LiftoverResult,
    Region,
    Transcript,
    apply_liftover_to_utr_file,
    liftover_region,
    liftover_transcripts,
    splice_into_alignment,
)


def _network_available() -> bool:
    try:
        urllib.request.urlopen("https://rest.ensembl.org/info/ping", timeout=5)
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _network_available(), reason="no network access to Ensembl REST API")


def test_clean_region_is_tagged_ok():
    region = Region("TEST_OK", "1", 1000000, 1000100, 1)
    result = liftover_region(region)
    assert result.tag == "ok"
    assert result.hg38_region == "1:1064620-1064720"


def test_centromeric_region_is_tagged_failed():
    region = Region("TEST_FAILED", "1", 121500000, 121500100, 1)
    result = liftover_region(region)
    assert result.tag == "failed"


def test_segdup_region_is_tagged_split():
    region = Region("TEST_SPLIT", "1", 144080000, 144082000, 1)
    result = liftover_region(region)
    assert result.tag == "split"
    assert "blocks" in result.note


def test_splice_into_alignment_preserves_gaps():
    aligned = "AC--GT"
    new_seq = "ACGA"  # same ungapped length (4) as aligned's non-gap count
    spliced = splice_into_alignment(aligned, new_seq)
    assert spliced == "AC--GA"


def test_splice_into_alignment_rejects_length_mismatch():
    aligned = "AC--GT"
    new_seq = "ACGAT"  # 5 chars, but aligned has only 4 non-gap positions
    assert splice_into_alignment(aligned, new_seq) is None


def test_apply_ok_tag_normalizes_rna_vs_dna_before_comparing(tmp_path):
    """Regression test: found while validating against real TargetScan data --
    Ensembl returns DNA (T), but TargetScan's UTR_Sequences.txt is RNA (U),
    so the "does this row match the hg19 region" check must normalize
    T<->U before comparing, or every real "ok" region would be wrongly
    rejected as not matching."""
    utr_file = tmp_path / "utr.txt"
    utr_file.write_text("GENE1\t9606\tAC-GU\nGENE1\t9615\tAC-GA\n")

    result = LiftoverResult(
        gene_id="GENE1",
        hg19_region="1:1-4",
        tag="ok",
        hg38_region="1:101-104",
        hg19_sequence="ACGT",  # DNA, as Ensembl would return it
    )

    out_file = tmp_path / "utr_out.txt"
    applied_file = tmp_path / "applied.tsv"
    apply_liftover_to_utr_file(str(utr_file), [result], str(out_file), str(applied_file))

    applied_text = applied_file.read_text()
    assert "applied (sequence unchanged" in applied_text
    assert "NOT applied" not in applied_text


def test_multiexon_transcript_orders_and_concatenates_segments():
    """ENST00000423372.3 is a real multi-exon-UTR transcript (- strand,
    2 genomic blocks). Confirmed by direct comparison against the real
    UTR_Sequences.txt row for this transcript: the correctly-ordered,
    concatenated hg19 sequence is exactly 1811 nt and starts with
    "CUGUGAGGCCAUUUCCAGGCC..." -- this test pins that down without
    requiring the (multi-GB) production file to be present."""
    t = Transcript(
        "ENST00000423372.3",
        [
            Region("ENST00000423372.3", "1", 134901, 135802, -1),
            Region("ENST00000423372.3", "1", 137621, 138529, -1),
        ],
    )
    # On the "-" strand, the higher-coordinate block is transcribed first.
    ordered = t.ordered_segments()
    assert [(r.start, r.end) for r in ordered] == [(137621, 138529), (134901, 135802)]

    results = liftover_transcripts([t])
    r = results["ENST00000423372.3"]
    assert r.tag == "ok"
    assert r.n_segments == 2
    assert len(r.hg19_sequence) == 1811
    assert r.hg19_sequence.upper().startswith("CTGTGAGGCCATTTCCAGGCC")


def test_real_demo_genes_lift_cleanly(tmp_path):
    """All 6 genes in examples/real_hg38_demo/ are real TargetScan vert80
    data and should lift hg19->hg38 cleanly and apply without any manual
    review needed."""
    demo_dir = os.path.join(ROOT, "examples", "real_hg38_demo")
    regions_file = os.path.join(demo_dir, "hg19_3utr_regions.bed")
    utr_file = os.path.join(demo_dir, "UTR_Sequences_6genes.txt")
    if not (os.path.exists(regions_file) and os.path.exists(utr_file)):
        pytest.skip("examples/real_hg38_demo/ not present")

    from targetscan.hg38_liftover import apply_liftover_to_utr_file, liftover_regions, read_regions_bed

    regions = read_regions_bed(regions_file)
    results = liftover_regions(regions)
    assert {r.tag for r in results} == {"ok"}

    out_file = tmp_path / "utr_out.txt"
    applied_file = tmp_path / "applied.tsv"
    apply_liftover_to_utr_file(utr_file, results, str(out_file), str(applied_file))

    applied_text = applied_file.read_text()
    assert applied_text.count("applied (sequence unchanged") == len(regions)
