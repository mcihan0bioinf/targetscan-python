#!/usr/bin/env python3
"""CLI: count 8mer sites in ORFs (port of targetscan_count_8mers.pl).

USAGE: targetscan_count_8mers.py miRNA_seeds_file UTRs out_file
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.count_8mers import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
