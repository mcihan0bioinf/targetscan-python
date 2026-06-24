"""Lift human 3' UTR regions from hg19 (GRCh37) to hg38 (GRCh38) and tag
each one with whether it's safe to splice the new hg38 sequence into an
existing TargetScan multi-species alignment as-is.

This only re-anchors the *human* side of TargetScan's existing
multi-species alignment -- the other ~80 species' sequences are untouched
(they stay on whatever assembly TargetScan originally aligned them to).
That's why each region needs a clear "is this actually safe to swap in"
tag rather than just blindly trusting the coordinate mapping: a clean
1-to-1 coordinate mapping that also returns byte-identical sequence is
safe to substitute directly (no realignment needed, since the alignment
columns for the human row don't change). Anything else needs a human to
look at it before it's used.

Uses the Ensembl REST API (https://rest.ensembl.org) for both coordinate
mapping and sequence retrieval, so no local genome FASTA or UCSC chain
file download is required.

Tags
----
ok       - exactly one mapped block, same length as the input region, and
           the hg38 sequence is byte-identical to the hg19 sequence.
           Safe to substitute the human row's coordinates (sequence is
           unchanged, so the existing alignment columns still apply).
shifted  - exactly one mapped block, but the sequence differs from hg19
           (small edits/indels between assemblies) or the mapped length
           differs from the input length. NOT safe to auto-splice --
           the human row would need realignment against the other
           species in that UTR.
split    - more than one mapped block returned (the region is broken
           across multiple hg38 locations, e.g. an insertion or
           rearrangement between assemblies). NOT safe to auto-splice.
failed   - no mapping returned at all (e.g. region overlaps an assembly
           gap, or was removed in hg38). Cannot be lifted.
"""

from __future__ import annotations

import csv
import sys
import time
import urllib.error
import urllib.request
import json
from dataclasses import dataclass, field
from typing import List, Optional

ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES = "human"
HG19_ASSEMBLY = "GRCh37"
HG38_ASSEMBLY = "GRCh38"

# Be polite to the public Ensembl REST API (it rate-limits at ~15 req/s).
REQUEST_DELAY_SECONDS = 0.08
MAX_RETRIES = 3


def _get_json(url: str) -> dict:
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
            time.sleep(REQUEST_DELAY_SECONDS)
            return data
        except urllib.error.HTTPError as exc:
            if exc.code == 429:  # rate limited
                time.sleep(1.0)
                last_error = exc
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(0.5)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


@dataclass
class Region:
    gene_id: str
    chrom: str
    start: int  # 1-based, inclusive
    end: int  # 1-based, inclusive
    strand: int = 1  # 1 or -1


@dataclass
class LiftoverResult:
    gene_id: str
    hg19_region: str
    tag: str
    hg38_region: Optional[str] = None
    note: str = ""
    hg38_sequence: Optional[str] = None
    hg19_sequence: Optional[str] = None


def map_region_to_hg38(region: Region) -> List[dict]:
    """Return Ensembl's list of mapped blocks for this hg19 region in hg38."""
    loc = f"{region.chrom}:{region.start}..{region.end}:{region.strand}"
    url = f"{ENSEMBL_REST}/map/{SPECIES}/{HG19_ASSEMBLY}/{loc}/{HG38_ASSEMBLY}?content-type=application/json"
    data = _get_json(url)
    return data.get("mappings", [])


def fetch_sequence(chrom: str, start: int, end: int, strand: int, assembly: str) -> str:
    loc = f"{chrom}:{start}..{end}:{strand}"
    url = (
        f"{ENSEMBL_REST}/sequence/region/{SPECIES}/{loc}"
        f"?content-type=application/json;coord_system_version={assembly}"
    )
    data = _get_json(url)
    return data["seq"]


def liftover_region(region: Region, fetch_sequences: bool = True) -> LiftoverResult:
    hg19_region_str = f"{region.chrom}:{region.start}-{region.end}"
    input_length = region.end - region.start + 1

    try:
        mappings = map_region_to_hg38(region)
    except Exception as exc:  # network/API failure -- not a liftover failure per se
        return LiftoverResult(region.gene_id, hg19_region_str, "failed", note=f"API error: {exc}")

    if not mappings:
        return LiftoverResult(region.gene_id, hg19_region_str, "failed", note="No mapping returned (assembly gap or removed region)")

    if len(mappings) > 1:
        blocks = "; ".join(
            f"{m['mapped']['seq_region_name']}:{m['mapped']['start']}-{m['mapped']['end']}" for m in mappings
        )
        return LiftoverResult(
            region.gene_id, hg19_region_str, "split",
            note=f"Region split across {len(mappings)} blocks in hg38: {blocks}",
        )

    mapped = mappings[0]["mapped"]
    hg38_chrom = mapped["seq_region_name"]
    hg38_start, hg38_end = mapped["start"], mapped["end"]
    hg38_strand = mapped["strand"]
    hg38_region_str = f"{hg38_chrom}:{hg38_start}-{hg38_end}"
    mapped_length = hg38_end - hg38_start + 1

    if mapped_length != input_length:
        return LiftoverResult(
            region.gene_id, hg19_region_str, "shifted", hg38_region=hg38_region_str,
            note=f"Mapped length {mapped_length} != input length {input_length} (indel between assemblies)",
        )

    if not fetch_sequences:
        return LiftoverResult(region.gene_id, hg19_region_str, "ok", hg38_region=hg38_region_str, note="Length matches (sequence not checked)")

    try:
        hg19_seq = fetch_sequence(region.chrom, region.start, region.end, region.strand, HG19_ASSEMBLY)
        hg38_seq = fetch_sequence(hg38_chrom, hg38_start, hg38_end, hg38_strand, HG38_ASSEMBLY)
    except Exception as exc:
        return LiftoverResult(
            region.gene_id, hg19_region_str, "shifted", hg38_region=hg38_region_str,
            note=f"Could not fetch sequence to verify: {exc}",
        )

    if hg19_seq.upper() == hg38_seq.upper():
        return LiftoverResult(
            region.gene_id, hg19_region_str, "ok", hg38_region=hg38_region_str,
            note="Sequence identical; safe to re-anchor coordinates only",
            hg38_sequence=hg38_seq, hg19_sequence=hg19_seq,
        )

    return LiftoverResult(
        region.gene_id, hg19_region_str, "shifted", hg38_region=hg38_region_str,
        note="Coordinates map cleanly but sequence content differs between assemblies",
        hg38_sequence=hg38_seq, hg19_sequence=hg19_seq,
    )


def liftover_regions(regions: List[Region], fetch_sequences: bool = True) -> List[LiftoverResult]:
    results = []
    for region in regions:
        results.append(liftover_region(region, fetch_sequences=fetch_sequences))
    return results


def read_regions_bed(path: str) -> List[Region]:
    """Read a BED-like file: gene_id, chrom, start, end, strand (+/- or 1/-1).

    1-based inclusive coordinates are expected (not 0-based BED-style half
    open), matching how TargetScan's own 3' UTR coordinate files are laid
    out. Convert first if your source file is 0-based half-open BED.
    """
    regions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            gene_id, chrom, start, end = fields[0], fields[1], int(fields[2]), int(fields[3])
            strand_field = fields[4] if len(fields) > 4 else "+"
            strand = -1 if strand_field in ("-", "-1") else 1
            regions.append(Region(gene_id, chrom, start, end, strand))
    return regions


HUMAN_SPECIES_ID = "9606"


def splice_into_alignment(aligned_seq: str, new_ungapped_seq: str) -> Optional[str]:
    """Substitute a new ungapped human sequence into an existing gapped
    alignment row, preserving the row's gap positions.

    Only safe when the new sequence has exactly as many characters as the
    old row had non-gap positions -- i.e. same alignment column count, so
    no realignment against the other species is needed. Returns None if
    that doesn't hold (caller should leave the row untouched and flag it).
    """
    non_gap_count = sum(1 for c in aligned_seq if c != "-")
    if non_gap_count != len(new_ungapped_seq):
        return None

    out = []
    i = 0
    for c in aligned_seq:
        if c == "-":
            out.append("-")
        else:
            out.append(new_ungapped_seq[i])
            i += 1
    return "".join(out)


def apply_liftover_to_utr_file(
    utr_file: str, results: List[LiftoverResult], out_path: str, applied_report_path: str
) -> None:
    """Re-anchor the human (9606) row of each gene in a TargetScan
    UTR_Sequences.txt-style file using its liftover result.

    - tag "ok": only applied if the existing alignment row's sequence
      actually matches the hg19 sequence the liftover was computed from
      (verified here, not assumed) -- otherwise the file's row doesn't
      correspond to the genomic region we looked up, and re-anchoring it
      would be wrong. If it matches, the row is left as-is (only its
      real-world coordinates changed); recorded as "applied (no change
      needed)".
    - tag "shifted" with a hg38 sequence available: spliced in *only* if
      it has the same ungapped length as the existing row (so it slots
      into the same alignment columns without realigning the other
      species); otherwise left untouched and flagged "NOT applied".
    - tags "split"/"failed": always left untouched and flagged
      "NOT applied" -- there is no safe hg38 sequence to use.
    """
    result_by_gene = {r.gene_id: r for r in results}
    applied_log = []

    with open(utr_file) as fh, open(out_path, "w") as out:
        for line in fh:
            raw = line.rstrip("\r\n")
            if not raw:
                continue
            fields = raw.split("\t")
            gene_id, species_id, seq = fields[0], fields[1], fields[2]

            if species_id == HUMAN_SPECIES_ID and gene_id in result_by_gene:
                r = result_by_gene[gene_id]
                if r.tag == "ok":
                    # TargetScan's UTR sequences are RNA (U); Ensembl returns DNA (T).
                    existing_ungapped = seq.replace("-", "").upper().replace("T", "U")
                    hg19_as_rna = r.hg19_sequence.upper().replace("T", "U") if r.hg19_sequence else None
                    if hg19_as_rna and existing_ungapped == hg19_as_rna:
                        applied_log.append((gene_id, "applied (sequence unchanged, coordinates re-anchored)"))
                    else:
                        applied_log.append(
                            (gene_id, "NOT applied: alignment row doesn't match the hg19 region that was lifted")
                        )
                elif r.tag == "shifted" and r.hg38_sequence:
                    spliced = splice_into_alignment(seq, r.hg38_sequence)
                    if spliced is not None:
                        seq = spliced
                        applied_log.append((gene_id, "applied (spliced new hg38 sequence into existing columns)"))
                    else:
                        applied_log.append((gene_id, "NOT applied: hg38 sequence length differs from alignment row"))
                else:
                    applied_log.append((gene_id, f"NOT applied: tag={r.tag}"))

            out.write(f"{gene_id}\t{species_id}\t{seq}\n")

    with open(applied_report_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id", "status"])
        for gene_id, status in applied_log:
            writer.writerow([gene_id, status])


def write_report(results: List[LiftoverResult], out_path: str) -> None:
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id", "hg19_region", "tag", "hg38_region", "note"])
        for r in results:
            writer.writerow([r.gene_id, r.hg19_region, r.tag, r.hg38_region or "", r.note])


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        print(
            "USAGE: hg38_liftover.py hg19_regions.bed report_out.tsv [--no-sequence-check]\n"
            "       [--utr-file UTR_Sequences.txt --utr-out UTR_Sequences.hg38.txt --applied-out applied.tsv]",
            file=sys.stderr,
        )
        return 0
    regions_file, out_path = argv[0], argv[1]
    fetch_sequences = "--no-sequence-check" not in argv

    regions = read_regions_bed(regions_file)
    print(f"Lifting {len(regions)} region(s) from hg19 to hg38...", file=sys.stderr)
    results = liftover_regions(regions, fetch_sequences=fetch_sequences)
    write_report(results, out_path)

    counts: dict = {}
    for r in results:
        counts[r.tag] = counts.get(r.tag, 0) + 1
    print(f"Done. {counts}. See {out_path}", file=sys.stderr)

    if "--utr-file" in argv:
        utr_file = argv[argv.index("--utr-file") + 1]
        utr_out = argv[argv.index("--utr-out") + 1] if "--utr-out" in argv else utr_file + ".hg38.txt"
        applied_out = argv[argv.index("--applied-out") + 1] if "--applied-out" in argv else out_path + ".applied.tsv"
        apply_liftover_to_utr_file(utr_file, results, utr_out, applied_out)
        print(f"Re-anchored UTR file written to {utr_out}; per-gene status in {applied_out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
