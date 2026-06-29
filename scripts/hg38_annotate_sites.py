#!/usr/bin/env python3
"""CLI: annotate predicted sites' UTR positions with real hg38 genomic coordinates.

USAGE: hg38_annotate_sites.py predicted_targets.txt liftover_report.tsv hg19_3utr.gff out.txt
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.hg38_annotate_sites import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
