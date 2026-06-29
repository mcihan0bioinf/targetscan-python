"""Tests for lifting TargetScan's per-site hg19 BED files to hg38 by
reusing the already-verified transcript-level liftover report (pure
arithmetic, no new API calls).

The AJAP1 case is pinned to a value cross-checked against a live Ensembl
/map call (verified manually): hg19 chr1:4834591-4834598 (1-based) maps
to hg38 chr1:4774531-4774538.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from targetscan.hg38_liftover_sites import (  # noqa: E402
    TranscriptLiftInfo,
    liftover_site,
)


def test_ok_single_block():
    gene_to_tx = {"AJAP1": "ENST00000378191.4"}
    transcript_info = {
        "ENST00000378191.4": TranscriptLiftInfo(
            tag="ok",
            blocks=[("1", 4834560, 4834618, 4774500, 4774558)],
        )
    }
    chrom, start, end, tag = liftover_site("1", 4834591, 4834598, "AJAP1", gene_to_tx, transcript_info)
    assert tag == "ok"
    assert (chrom, start, end) == ("1", 4774531, 4774538)


def test_site_spanning_two_contiguous_blocks_with_same_offset():
    """Mirrors the real MCM10 case: a site starting in one GFF-recorded
    block (UTR) and ending in the next (ORF), where both blocks shift by
    the same offset since they're genomically contiguous."""
    gene_to_tx = {"MCM10": "ENST00000378694.1"}
    transcript_info = {
        "ENST00000378694.1": TranscriptLiftInfo(
            tag="ok",
            blocks=[
                ("10", 13251161, 13251226, 13209161, 13209226),
                ("10", 13251227, 13251310, 13209227, 13209310),
            ],
        )
    }
    # Site spans 13251222-13251230, crossing the 13251226/13251227 boundary.
    chrom, start, end, tag = liftover_site("10", 13251222, 13251230, "MCM10", gene_to_tx, transcript_info)
    assert tag == "ok"
    assert (chrom, start, end) == ("10", 13209222, 13209230)


def test_site_spanning_blocks_with_different_offsets_is_flagged():
    gene_to_tx = {"GENE1": "TX1"}
    transcript_info = {
        "TX1": TranscriptLiftInfo(
            tag="ok",
            blocks=[
                ("1", 100, 200, 1000, 1100),  # offset +900
                ("1", 201, 300, 1102, 1201),  # offset +901 (inconsistent)
            ],
        )
    }
    _, _, _, tag = liftover_site("1", 195, 205, "GENE1", gene_to_tx, transcript_info)
    assert tag == "site_spans_inconsistent_blocks"


def test_no_gene_mapping():
    _, _, _, tag = liftover_site("1", 100, 110, "UNKNOWN_GENE", {}, {})
    assert tag == "no_gene_mapping"


def test_no_liftover_data_for_transcript():
    gene_to_tx = {"GENE1": "TX1"}
    _, _, _, tag = liftover_site("1", 100, 110, "GENE1", gene_to_tx, {})
    assert tag == "no_liftover_data"


def test_propagates_non_ok_transcript_tag():
    gene_to_tx = {"GENE1": "TX1"}
    transcript_info = {"TX1": TranscriptLiftInfo(tag="split", blocks=[])}
    chrom, start, end, tag = liftover_site("1", 100, 110, "GENE1", gene_to_tx, transcript_info)
    assert tag == "split"
    assert chrom is None and start is None and end is None


def test_site_outside_any_block():
    gene_to_tx = {"GENE1": "TX1"}
    transcript_info = {
        "TX1": TranscriptLiftInfo(tag="ok", blocks=[("1", 100, 200, 1000, 1100)])
    }
    _, _, _, tag = liftover_site("1", 5000, 5010, "GENE1", gene_to_tx, transcript_info)
    assert tag == "site_outside_lifted_blocks"
