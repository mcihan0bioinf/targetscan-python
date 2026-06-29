#!/usr/bin/env python3
"""CLI: lift human 3' UTR regions from hg19 to hg38 and tag each one.

USAGE: hg38_liftover.py hg19_regions.bed report_out.tsv [--no-sequence-check]
       [--utr-file UTR_Sequences.txt --utr-out UTR_Sequences.hg38.txt --applied-out applied.tsv]
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.hg38_liftover import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
