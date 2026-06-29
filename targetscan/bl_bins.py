"""Assign each UTR (gene) to a branch-length conservation bin (1-10).

Python port of ``tsh_orig/TargetScan7_BL_PCT/targetscan_70_BL_bins.pl``.
For every alignment column that has a non-gap nucleotide in the reference
genome (human, taxon 9606 by default), compute the branch length spanning
all species sharing that nucleotide; the median of those values across the
whole UTR alignment determines its bin.
"""

from __future__ import annotations

import statistics
import sys
from typing import Dict, Iterator, List, Optional, Tuple

from .phylo import Tree

REF_GENOME = "9606"

# Minimum BL for each bin (1-10), from TargetScan 7's 3' UTR partitioning.
BL_THRESHOLDS = [
    0,
    1.21207417,
    2.17396073,
    2.80215414,
    3.26272822,
    3.65499277,
    4.01461968,
    4.40729032,
    4.90457274,
    5.78196252,
]


def _format_number(value: float) -> str:
    """Format like Perl's default number stringification (up to 15 sig figs)."""
    if value == int(value):
        return str(int(value))
    return f"{value:.15g}"


def assign_bin(bl: float) -> int:
    bin_this_bl = 1
    for i, threshold in enumerate(BL_THRESHOLDS, start=1):
        if bl > threshold:
            bin_this_bl = i
    return bin_this_bl


def _read_alignment_blocks(path: str) -> Iterator[Tuple[str, Dict[str, str]]]:
    """Yield (gene_id, {species_id: raw_aligned_sequence}) gene by gene."""
    last_gene: Optional[str] = None
    species_to_alignment: Dict[str, str] = {}

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            gene_id, species_id, seq = fields[0], fields[1], fields[2]

            if gene_id != last_gene and last_gene is not None:
                yield last_gene, species_to_alignment
                species_to_alignment = {}

            species_to_alignment[species_id] = seq
            last_gene = gene_id

    if last_gene is not None:
        yield last_gene, species_to_alignment


def median_branch_length(
    tree: Tree, species_to_alignment: Dict[str, str], ref_genome: str = REF_GENOME
) -> Optional[float]:
    """Median branch length across alignment columns sharing the ref nt."""
    if ref_genome not in species_to_alignment:
        return None

    alignment_length = len(species_to_alignment[ref_genome])
    species_list_to_bl: Dict[Tuple[str, ...], float] = {}
    branch_lengths: List[float] = []

    species = list(species_to_alignment.keys())

    for i in range(alignment_length):
        consensus = species_to_alignment[ref_genome][i]
        if consensus == "-":
            continue

        matching = sorted(
            sp
            for sp in species
            if species_to_alignment[sp][i].upper() == consensus.upper()
        )
        if not matching:
            continue

        key = tuple(matching)
        if len(matching) == 1:
            bl = 0.0
        elif key in species_list_to_bl:
            bl = species_list_to_bl[key]
        else:
            bl = tree.branch_length(ref_genome, matching)
            species_list_to_bl[key] = bl

        branch_lengths.append(bl)

    if not branch_lengths:
        return None
    return statistics.median(branch_lengths)


def run(utr_file: str, tree_file: str, out_path: str) -> None:
    tree = Tree(tree_file)
    with open(out_path, "w") as out:
        for gene_id, species_to_alignment in _read_alignment_blocks(utr_file):
            median_bl = median_branch_length(tree, species_to_alignment)
            if median_bl is None:
                continue
            bin_this_utr = assign_bin(median_bl)
            out.write(f"{gene_id}\t{_format_number(median_bl)}\t{bin_this_utr}\n")


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 1:
        print(
            "USAGE: targetscan_70_BL_bins.py tabbedAlignmentFile [tree_file] > Gene_BL_bin_file",
            file=sys.stderr,
        )
        return 0
    utr_file = argv[0]
    tree_file = argv[1] if len(argv) > 1 else "data/PCT_parameters/Tree.generic.txt"
    out_path = argv[2] if len(argv) > 2 else "/dev/stdout"
    run(utr_file, tree_file, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
