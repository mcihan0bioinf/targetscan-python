"""Annotate predicted miRNA target sites with their actual hg38 genomic
coordinates.

``site_prediction.py`` (and the rest of the pipeline) reports each site's
position as a UTR-relative offset (``UTR_start``/``UTR_end``), since the
algorithm only ever looks at sequence, never genome coordinates. This
module closes that gap for the human (9606) rows: using the per-transcript
hg38 block(s) recorded by ``hg38_liftover.py``'s report, it converts a
UTR-relative position back into a real ``chrom:start-end`` hg38 location.

This only works for transcripts tagged ``ok`` by the liftover (clean,
single-or-multi-block mapping with verified-identical sequence) -- those
are the only ones with an unambiguous hg38 coordinate to report.

A site can itself span a splice junction (rare, but possible for
multi-exon 3' UTRs if a junction falls inside a 6-8nt seed match), so the
output is a ``chrom:start-end[; chrom:start-end ...]`` block list, same
style as the liftover report's own ``hg38_region`` field, rather than a
single start/end pair.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .hg38_liftover import read_transcripts_gff

HUMAN_SPECIES_ID = "9606"


@dataclass
class TranscriptGenomicMap:
    chrom: str
    strand: int  # 1 or -1
    blocks: List[Tuple[int, int]]  # hg38 (start, end), in 5'->3' transcript order


def parse_region_blocks(region_field: str) -> List[Tuple[str, int, int]]:
    """Parse a "chrom:start-end; chrom:start-end" field (as written by
    hg38_liftover's report) into a list of (chrom, start, end)."""
    blocks = []
    for part in region_field.split(";"):
        part = part.strip()
        if not part:
            continue
        chrom, coords = part.rsplit(":", 1)
        start, end = coords.split("-")
        blocks.append((chrom, int(start), int(end)))
    return blocks


def build_genomic_maps(report_path: str, gff_path: str) -> Dict[str, TranscriptGenomicMap]:
    """Build {transcript_id: TranscriptGenomicMap} for every "ok"-tagged
    transcript in a liftover report, using strand from the original GFF
    (the report itself doesn't carry strand)."""
    transcripts = read_transcripts_gff(gff_path)
    strand_by_id = {t.transcript_id: t.segments[0].strand for t in transcripts}

    maps: Dict[str, TranscriptGenomicMap] = {}
    with open(report_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["tag"] != "ok" or not row.get("hg38_region"):
                continue
            blocks = parse_region_blocks(row["hg38_region"])
            chrom = blocks[0][0]
            strand = strand_by_id.get(row["gene_id"], 1)
            maps[row["gene_id"]] = TranscriptGenomicMap(chrom, strand, [(s, e) for _c, s, e in blocks])
    return maps


def utr_pos_to_genomic(tmap: TranscriptGenomicMap, utr_pos: int) -> Optional[Tuple[str, int]]:
    """Convert a 1-based UTR-relative position (5'->3') to (chrom, genomic_pos)."""
    remaining = utr_pos
    for start, end in tmap.blocks:
        length = end - start + 1
        if remaining <= length:
            genomic = start + remaining - 1 if tmap.strand == 1 else end - remaining + 1
            return tmap.chrom, genomic
        remaining -= length
    return None  # utr_pos beyond the end of the transcript's recorded UTR length


def annotate_utr_range(tmap: TranscriptGenomicMap, utr_start: int, utr_end: int) -> str:
    """Convert a UTR-relative [utr_start, utr_end] range to hg38 genomic
    block(s). Usually one block; more than one if the site spans a splice
    junction in a multi-exon UTR."""
    positions = []
    for p in range(utr_start, utr_end + 1):
        res = utr_pos_to_genomic(tmap, p)
        if res is None:
            return "NA"
        positions.append(res)

    step = 1 if tmap.strand == 1 else -1
    blocks = []
    chrom0, g0 = positions[0]
    cur_start = cur_end = g0
    for chrom, g in positions[1:]:
        if chrom == chrom0 and g == cur_end + step:
            cur_end = g
        else:
            lo, hi = (cur_start, cur_end) if cur_start <= cur_end else (cur_end, cur_start)
            blocks.append(f"{chrom0}:{lo}-{hi}")
            chrom0, cur_start, cur_end = chrom, g, g
    lo, hi = (cur_start, cur_end) if cur_start <= cur_end else (cur_end, cur_start)
    blocks.append(f"{chrom0}:{lo}-{hi}")
    return "; ".join(blocks)


def annotate_predicted_targets(
    predicted_targets_path: str,
    genomic_maps: Dict[str, TranscriptGenomicMap],
    out_path: str,
    utr_start_col: int = 5,
    utr_end_col: int = 6,
    gene_col: int = 0,
    species_col: int = 2,
) -> None:
    """Add an "hg38_location" column to a predicted-targets-style file
    (works for both ``predicted_targets.txt`` and ``*.bl_pct.txt``, which
    share the same first 7 columns) for human rows of transcripts with a
    clean ("ok") liftover; everything else gets "NA"."""
    with open(predicted_targets_path) as fh, open(out_path, "w") as out:
        header = fh.readline().rstrip("\n")
        out.write(header + "\thg38_location\n")

        for line in fh:
            raw = line.rstrip("\n")
            if not raw:
                continue
            fields = raw.split("\t")
            gene_id = fields[gene_col]
            species_id = fields[species_col]

            location = "NA"
            if species_id == HUMAN_SPECIES_ID and gene_id in genomic_maps:
                utr_start, utr_end = int(fields[utr_start_col]), int(fields[utr_end_col])
                location = annotate_utr_range(genomic_maps[gene_id], utr_start, utr_end)

            out.write(raw + "\t" + location + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 4:
        print(
            "USAGE: hg38_annotate_sites.py predicted_targets.txt liftover_report.tsv "
            "hg19_3utr.gff out.txt",
            file=sys.stderr,
        )
        return 0
    predicted_targets_path, report_path, gff_path, out_path = argv[:4]
    genomic_maps = build_genomic_maps(report_path, gff_path)
    print(f"Loaded hg38 genomic maps for {len(genomic_maps)} transcripts.", file=sys.stderr)
    annotate_predicted_targets(predicted_targets_path, genomic_maps, out_path)
    print(f"Annotated file written to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
