#!/usr/bin/env python3
"""CLI: calculate context++ scores (port of targetscan_70_context_scores.pl).

USAGE: targetscan_70_context_scores.py miRNA_file UTR_file PredictedTargetsBL_PCT_file \\
           ORF_lengths_file ORF_8mer_counts_file ContextScoresOutput_file [data_dir] [rnaplfold_dir]
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.context_scores import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
