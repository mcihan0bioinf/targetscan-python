"""Lift TargetScan's per-site hg19 BED files (predicted target site genomic
locations) to hg38, at the scale of millions of rows -- without making a
new Ensembl API call per site.

Each BED row already has plain genomic (not transcript-relative) hg19
coordinates for one predicted site. Critically, that's the *same kind* of
coordinate hg38_liftover.py already verified at the transcript level
(```full_liftover_report.tsv```): for any transcript tagged "ok" there,
its hg19 3' UTR block(s) map to hg38 block(s) of identical length with
byte-identical sequence -- which means the hg19->hg38 mapping *within*
that block is a constant additive offset (BED coordinates on both
assemblies are genomic plus-strand, so no strand bookkeeping is needed
here, unlike hg38_annotate_sites.py's UTR-relative-offset case). So each
site's new coordinate is `site_pos + (hg38_block_start - hg19_block_start)`
for whichever block contains it -- pure arithmetic, reusing the already
Ensembl-verified mapping, not a fresh per-site liftover call.

Sites under genes whose transcript wasn't tagged "ok" (shifted/split/
failed/not found at all) get that same tag and no hg38 coordinate --
consistent with the same conservative "only act on verified-safe regions"
policy as the rest of this project.
"""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

NULL_TAGS = {"no_gene_mapping", "no_liftover_data", "site_outside_lifted_blocks"}


@dataclass
class TranscriptLiftInfo:
    tag: str
    # (chrom, hg19_start, hg19_end, hg38_start, hg38_end) per block, 1-based inclusive
    blocks: List[Tuple[str, int, int, int, int]]


def _parse_region_field(region_field: str) -> List[Tuple[str, int, int]]:
    blocks = []
    for part in region_field.split(";"):
        part = part.strip()
        if not part:
            continue
        chrom, coords = part.rsplit(":", 1)
        start, end = coords.split("-")
        blocks.append((chrom, int(start), int(end)))
    return blocks


def load_gene_to_transcript(gene_info_path: str, species: str = "9606") -> Dict[str, str]:
    """Gene symbol -> representative transcript ID, for one species."""
    mapping = {}
    with open(gene_info_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["Species ID"] == species and row["Representative transcript?"] == "1":
                mapping[row["Gene symbol"]] = row["Transcript ID"]
    return mapping


def load_transcript_lift_info(report_path: str) -> Dict[str, TranscriptLiftInfo]:
    """Parse the transcript-level liftover report into per-transcript
    block pairs (only meaningful/populated for tag == "ok")."""
    info: Dict[str, TranscriptLiftInfo] = {}
    with open(report_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            blocks: List[Tuple[str, int, int, int, int]] = []
            if row["tag"] == "ok" and row.get("hg38_region"):
                hg19_blocks = _parse_region_field(row["hg19_region"])
                hg38_blocks = _parse_region_field(row["hg38_region"])
                if len(hg19_blocks) == len(hg38_blocks):
                    for (chrom19, s19, e19), (chrom38, s38, e38) in zip(hg19_blocks, hg38_blocks):
                        blocks.append((chrom19, s19, e19, s38, e38))
            info[row["gene_id"]] = TranscriptLiftInfo(tag=row["tag"], blocks=blocks)
    return info


def _offset_for_position(chrom: str, pos: int, blocks: List[Tuple[str, int, int, int, int]]) -> Optional[int]:
    for block_chrom, hg19_start, hg19_end, hg38_start, hg38_end in blocks:
        if block_chrom == chrom and hg19_start <= pos <= hg19_end:
            return hg38_start - hg19_start
    return None


def liftover_site(
    chrom: str,
    start_1based: int,
    end_1based: int,
    gene_symbol: str,
    gene_to_tx: Dict[str, str],
    transcript_info: Dict[str, TranscriptLiftInfo],
) -> Tuple[Optional[str], Optional[int], Optional[int], str]:
    """Returns (hg38_chrom, hg38_start_1based, hg38_end_1based, tag).

    A transcript's recorded blocks come straight from TargetScan's GFF,
    which splits annotations by feature type (e.g. "UTR" vs "ORF") even
    when they're genomically contiguous -- so a short site can legitimately
    start in one recorded block and end in the very next one. Rather than
    requiring strict single-block containment, we look up the offset for
    the site's start and end positions independently and only fail if
    either is unresolvable or they disagree (which would mean the site
    spans an actual structural difference, not just an annotation split).
    """
    tx = gene_to_tx.get(gene_symbol)
    if tx is None:
        return None, None, None, "no_gene_mapping"

    info = transcript_info.get(tx)
    if info is None:
        return None, None, None, "no_liftover_data"

    if info.tag != "ok":
        return None, None, None, info.tag

    start_offset = _offset_for_position(chrom, start_1based, info.blocks)
    end_offset = _offset_for_position(chrom, end_1based, info.blocks)

    if start_offset is None or end_offset is None:
        return None, None, None, "site_outside_lifted_blocks"
    if start_offset != end_offset:
        return None, None, None, "site_spans_inconsistent_blocks"

    return chrom, start_1based + start_offset, end_1based + end_offset, "ok"


def liftover_bed_file(
    in_path: str,
    out_path: str,
    gene_to_tx: Dict[str, str],
    transcript_info: Dict[str, TranscriptLiftInfo],
) -> Dict[str, int]:
    """Stream a TargetScan site BED file, append hg38 coordinates + tag
    columns, write the result. Returns tag counts."""
    counts: Dict[str, int] = {}

    with open(in_path) as fh, open(out_path, "w") as out:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom, start0, end0, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            gene_symbol = name.split(":", 1)[0]

            # The transcript-level report (derived from TargetScan's GFF)
            # stores chromosomes without the "chr" prefix; BED files have it.
            bare_chrom = chrom[3:] if chrom.startswith("chr") else chrom

            start1, end1 = start0 + 1, end0  # BED half-open 0-based -> 1-based inclusive

            hg38_chrom, hg38_start1, hg38_end1, tag = liftover_site(
                bare_chrom, start1, end1, gene_symbol, gene_to_tx, transcript_info
            )
            counts[tag] = counts.get(tag, 0) + 1

            if tag == "ok":
                hg38_start0, hg38_end0 = hg38_start1 - 1, hg38_end1
                extra = [f"chr{hg38_chrom}", str(hg38_start0), str(hg38_end0), tag]
            else:
                extra = ["", "", "", tag]

            out.write("\t".join(fields + extra) + "\n")

    return counts


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 4:
        print(
            "USAGE: hg38_liftover_sites.py gene_info.txt transcript_liftover_report.tsv "
            "in.bed out.bed",
            file=sys.stderr,
        )
        return 0

    gene_info_path, report_path, in_path, out_path = argv[:4]

    gene_to_tx = load_gene_to_transcript(gene_info_path)
    transcript_info = load_transcript_lift_info(report_path)
    print(f"Loaded {len(gene_to_tx)} gene->transcript mappings, {len(transcript_info)} transcript liftover records.", file=sys.stderr)

    counts = liftover_bed_file(in_path, out_path, gene_to_tx, transcript_info)
    print(f"{in_path} -> {out_path}: {counts}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
