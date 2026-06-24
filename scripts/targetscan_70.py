#!/usr/bin/env python3
"""CLI: predict miRNA target sites (port of targetscan_70.pl).

USAGE: targetscan_70.py miRNA_file UTR_file PredictedTargetsOutputFile
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.site_prediction import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
