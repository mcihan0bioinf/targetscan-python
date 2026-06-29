#!/usr/bin/env python3
"""CLI: assign UTRs to branch-length bins (port of targetscan_70_BL_bins.pl).

USAGE: targetscan_70_BL_bins.py tabbedAlignmentFile [tree_file] [out_file]
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.bl_bins import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
