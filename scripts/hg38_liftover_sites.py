#!/usr/bin/env python3
"""CLI: lift a TargetScan per-site hg19 BED file to hg38, reusing the
already-verified transcript-level liftover report (pure arithmetic, no
new Ensembl API calls).

USAGE: hg38_liftover_sites.py gene_info.txt transcript_liftover_report.tsv in.bed out.bed
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from targetscan.hg38_liftover_sites import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
