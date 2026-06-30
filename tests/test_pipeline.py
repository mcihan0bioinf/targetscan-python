"""Regression tests: run each pipeline stage on the bundled samples and
compare against the original Perl tool's golden output files.

Run with: pytest tests/test_pipeline.py -v
"""

import os
import subprocess
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(THIS_DIR)
SAMPLES = os.path.join(ROOT, "samples")
DATA = os.path.join(ROOT, "data")
GOLDEN = os.path.join(THIS_DIR, "golden")


def _read(path):
    with open(path) as fh:
        return fh.read()


def test_site_prediction(tmp_path):
    sys.path.insert(0, ROOT)
    from targetscan import site_prediction

    out = tmp_path / "out.txt"
    site_prediction.run(
        os.path.join(SAMPLES, "miR_Family_info_sample.txt"),
        os.path.join(SAMPLES, "UTR_Sequences_sample.txt"),
        str(out),
        verbose=False,
    )
    assert _read(str(out)) == _read(os.path.join(GOLDEN, "targetscan_70_output.txt"))


def test_bl_bins(tmp_path):
    sys.path.insert(0, ROOT)
    from targetscan import bl_bins

    out = tmp_path / "out.txt"
    bl_bins.run(
        os.path.join(SAMPLES, "UTR_Sequences_sample.txt"),
        os.path.join(DATA, "PCT_parameters", "Tree.generic.txt"),
        str(out),
    )
    assert _read(str(out)) == _read(os.path.join(GOLDEN, "UTRs_median_BLs_bins.txt"))


def test_bl_pct(tmp_path):
    sys.path.insert(0, ROOT)
    from targetscan import bl_pct

    out = tmp_path / "out.txt"
    bl_pct.run(
        os.path.join(SAMPLES, "miR_Family_info_sample.txt"),
        os.path.join(GOLDEN, "targetscan_70_output.txt"),
        os.path.join(GOLDEN, "UTRs_median_BLs_bins.txt"),
        os.path.join(DATA, "PCT_parameters"),
        str(out),
    )
    assert _read(str(out)) == _read(os.path.join(GOLDEN, "targetscan_70_output.BL_PCT.txt"))


def test_count_8mers(tmp_path):
    sys.path.insert(0, ROOT)
    from targetscan import count_8mers

    out = tmp_path / "out.txt"
    count_8mers.run(
        os.path.join(SAMPLES, "miR_Family_info_sample.txt"),
        os.path.join(SAMPLES, "ORF_Sequences_sample.txt"),
        str(out),
    )
    assert sorted(_read(str(out)).splitlines()) == sorted(
        _read(os.path.join(GOLDEN, "ORF_8mer_counts_sample.txt")).splitlines()
    )


def test_context_scores_with_cached_rnaplfold(tmp_path):
    sys.path.insert(0, ROOT)
    from targetscan import context_scores

    rnaplfold_dir = tmp_path / "RNAplfold_in_out"
    subprocess.run(
        ["cp", "-r", os.path.join(GOLDEN, "RNAplfold_in_out"), str(rnaplfold_dir)],
        check=True,
    )

    out = tmp_path / "out.txt"
    context_scores.run(
        os.path.join(SAMPLES, "miR_for_context_scores.sample.txt"),
        os.path.join(SAMPLES, "UTR_Sequences_sample.txt"),
        os.path.join(GOLDEN, "targetscan_70_output.BL_PCT.txt"),
        os.path.join(SAMPLES, "ORF_Sequences_sample.lengths.txt"),
        os.path.join(SAMPLES, "ORF_8mer_counts_sample.txt"),
        str(out),
        DATA,
        str(rnaplfold_dir),
    )
    assert _read(str(out)) == _read(os.path.join(GOLDEN, "Targets.BL_PCT.context_scores.txt"))
