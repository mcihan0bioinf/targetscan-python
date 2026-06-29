"""Tests for converting predicted sites' UTR-relative positions into real
hg38 genomic coordinates.

The plus/minus-strand single-block cases here are pinned to values
verified manually against real TargetScan data (FNDC3A and NLRP1 from
examples/real_hg38_demo/): e.g. NLRP1 is "-" strand with hg38 region
17:5499430-5501813 (length 2384); a site at UTR_start=2226, UTR_end=2231
should land at genomic 17:5499583-5499588 (since position 1 of the UTR,
5' end, is the *highest* genomic coordinate on "-" strand).
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from targetscan.hg38_annotate_sites import (  # noqa: E402
    TranscriptGenomicMap,
    annotate_predicted_targets,
    annotate_utr_range,
    parse_region_blocks,
    utr_pos_to_genomic,
)


def test_parse_region_blocks_single():
    assert parse_region_blocks("17:5499430-5501813") == [("17", 5499430, 5501813)]


def test_parse_region_blocks_multi():
    assert parse_region_blocks("1:137621-138529; 1:134901-135802") == [
        ("1", 137621, 138529),
        ("1", 134901, 135802),
    ]


def test_plus_strand_single_block():
    # FNDC3A, real data: hg38 region 13:49207396-49209779, "+" strand
    tmap = TranscriptGenomicMap("13", 1, [(49207396, 49209779)])
    assert utr_pos_to_genomic(tmap, 465) == ("13", 49207860)
    assert annotate_utr_range(tmap, 465, 471) == "13:49207860-49207866"


def test_minus_strand_single_block():
    # NLRP1, real data: hg38 region 17:5499430-5501813, "-" strand
    tmap = TranscriptGenomicMap("17", -1, [(5499430, 5501813)])
    # Position 1 (UTR 5' end) is the highest genomic coordinate on "-" strand.
    assert utr_pos_to_genomic(tmap, 1) == ("17", 5501813)
    assert annotate_utr_range(tmap, 2226, 2231) == "17:5499583-5499588"


def test_multiblock_minus_strand_continuous_site():
    # Mirrors ENST00000423372.3's real ordered blocks (137621-138529 then
    # 134901-135802, "-" strand); a site entirely within the first block.
    tmap = TranscriptGenomicMap("1", -1, [(137621, 138529), (134901, 135802)])
    # Block 1 has length 909 (138529-137621+1); position 905-909 is near its end.
    assert annotate_utr_range(tmap, 905, 909) == "1:137621-137625"


def test_multiblock_site_spans_splice_junction():
    """A site that straddles the boundary between two exonic blocks
    should report two separate genomic blocks, not one contiguous range."""
    tmap = TranscriptGenomicMap("1", -1, [(137621, 138529), (134901, 135802)])
    # Block 1 is positions 1-909; this site spans positions 907-911, i.e.
    # the last 3 nt of block 1 and the first 2 nt of block 2.
    assert annotate_utr_range(tmap, 907, 911) == "1:137621-137623; 1:135801-135802"


def test_utr_pos_beyond_transcript_returns_none():
    tmap = TranscriptGenomicMap("1", 1, [(100, 200)])  # length 101
    assert utr_pos_to_genomic(tmap, 102) is None


def test_annotate_predicted_targets_only_annotates_human_ok_genes(tmp_path):
    predicted = tmp_path / "predicted.txt"
    predicted.write_text(
        "a_Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\n"
        "GENE1\tfam\t9606\t1\t10\t5\t10\n"  # human, has a map -> annotated
        "GENE1\tfam\t9615\t1\t10\t5\t10\n"  # not human -> NA
        "GENE2\tfam\t9606\t1\t10\t5\t10\n"  # human, no map for GENE2 -> NA
    )

    genomic_maps = {"GENE1": TranscriptGenomicMap("1", 1, [(1000, 1100)])}

    out = tmp_path / "out.txt"
    annotate_predicted_targets(str(predicted), genomic_maps, str(out))

    lines = out.read_text().splitlines()
    assert lines[0].endswith("\thg38_location")
    assert lines[1].endswith("1:1004-1009")  # GENE1, 9606
    assert lines[2].endswith("\tNA")  # GENE1, not human
    assert lines[3].endswith("\tNA")  # GENE2, no map
