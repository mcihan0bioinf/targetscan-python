#!/usr/bin/env python3
"""CLI: calculate branch length and PCT (port of targetscan_70_BL_PCT.pl).

USAGE: targetscan_70_BL_PCT.py miRNA_file predicted_targets UTR_bin_info pct_params_dir [out_file]
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.bl_pct import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
