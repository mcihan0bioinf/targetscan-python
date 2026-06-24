#!/usr/bin/env python3
"""Run the full TargetScan prediction pipeline end to end.

Chains all five stages (site prediction -> BL binning -> BL/PCT
conservation -> ORF 8mer counting -> context++ scoring). Works on
TargetScan vert70 (hg19) or vert80 (hg38) data downloads -- only the
input files differ; see README.md for where to get vert80/hg38 files
(https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi).

USAGE:
    python run_pipeline.py \\
        --mirna-family miR_Family_info.txt \\
        --utr UTR_Sequences.txt \\
        --orf ORF_Sequences.txt \\
        --mirna-context miR_for_context_scores.txt \\
        --out-dir results/

Outputs (written to --out-dir):
    predicted_targets.txt              (site prediction)
    utr_bl_bins.txt                    (branch-length bins)
    predicted_targets.bl_pct.txt       (branch length + PCT)
    orf_lengths.txt, orf_8mer_counts.txt
    context_scores.txt                 (final context++ scores)
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from targetscan import bl_bins, bl_pct, context_scores, count_8mers, site_prediction  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.join(THIS_DIR, "data")
DEFAULT_TREE_FILE = os.path.join(DEFAULT_DATA_DIR, "PCT_parameters", "Tree.generic.txt")
DEFAULT_PCT_PARAMS_DIR = os.path.join(DEFAULT_DATA_DIR, "PCT_parameters")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mirna-family", required=True, help="miR_Family_info.txt (family, seed, species list)")
    p.add_argument("--utr", required=True, help="UTR_Sequences.txt (aligned 3' UTRs)")
    p.add_argument("--orf", required=True, help="ORF_Sequences.txt (matching ORFs, for context scores)")
    p.add_argument("--mirna-context", required=True, help="mature miRNA file for context scores (family, species, mirbase ID, sequence)")
    p.add_argument("--out-dir", required=True, help="directory to write all pipeline output into")
    p.add_argument("--data-dir", default=DEFAULT_DATA_DIR, help="dir with Agarwal/TA_SPS/AIRs reference files")
    p.add_argument("--tree-file", default=DEFAULT_TREE_FILE, help="generic phylogenetic tree for BL binning")
    p.add_argument("--pct-params-dir", default=DEFAULT_PCT_PARAMS_DIR, help="dir with PCT parameter + per-bin tree files")
    p.add_argument("--rnaplfold-dir", default=None, help="dir for RNAplfold input/output (default: <out-dir>/RNAplfold_in_out)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rnaplfold_dir = args.rnaplfold_dir or os.path.join(args.out_dir, "RNAplfold_in_out")

    predicted_targets = os.path.join(args.out_dir, "predicted_targets.txt")
    bl_bins_file = os.path.join(args.out_dir, "utr_bl_bins.txt")
    bl_pct_file = os.path.join(args.out_dir, "predicted_targets.bl_pct.txt")
    orf_lengths_file = os.path.join(args.out_dir, "orf_lengths.txt")
    orf_8mers_file = os.path.join(args.out_dir, "orf_8mer_counts.txt")
    context_scores_file = os.path.join(args.out_dir, "context_scores.txt")

    print("[1/5] Predicting miRNA target sites...", file=sys.stderr)
    site_prediction.run(args.mirna_family, args.utr, predicted_targets)

    print("[2/5] Assigning UTRs to branch-length bins...", file=sys.stderr)
    bl_bins.run(args.utr, args.tree_file, bl_bins_file)

    print("[3/5] Calculating branch length + PCT conservation...", file=sys.stderr)
    bl_pct.run(args.mirna_family, predicted_targets, bl_bins_file, args.pct_params_dir, bl_pct_file)

    print("[4/5] Counting 8mer sites in ORFs...", file=sys.stderr)
    count_8mers.run(args.mirna_family, args.orf, orf_8mers_file)
    # count_8mers.run() also writes "<orf>.lengths.txt" next to --orf; move/copy it here.
    auto_lengths = (
        args.orf[: -len(".txt")] + ".lengths.txt" if args.orf.endswith(".txt") else args.orf + ".lengths.txt"
    )
    if os.path.abspath(auto_lengths) != os.path.abspath(orf_lengths_file):
        os.replace(auto_lengths, orf_lengths_file)

    print("[5/5] Calculating context++ scores...", file=sys.stderr)
    context_scores.run(
        args.mirna_context,
        args.utr,
        bl_pct_file,
        orf_lengths_file,
        orf_8mers_file,
        context_scores_file,
        args.data_dir,
        rnaplfold_dir,
    )

    print(f"\nDone. Final output: {context_scores_file}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
