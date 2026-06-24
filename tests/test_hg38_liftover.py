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

from targetscan.hg38_liftover import Region, liftover_region, splice_into_alignment  # noqa: E402


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
