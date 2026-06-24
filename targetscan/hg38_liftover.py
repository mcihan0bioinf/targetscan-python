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
file download is required. Coordinate mapping is done concurrently (a
small thread pool); sequence verification is batched via Ensembl's POST
endpoint (up to 50 regions per request) since that's what actually makes
a ~40,000-transcript run finish in well under an hour instead of several.

A 3' UTR can span multiple, non-adjacent genomic blocks (multi-exon UTRs
-- about a quarter of human transcripts have this). Each block is lifted
independently; a transcript's overall tag is the worst of its blocks'
tags (failed > split > shifted > ok), and its hg38 sequence (when safe)
is the per-block hg38 sequences concatenated in 5'->3' transcript order
(which is descending genomic order for "-" strand transcripts).

Tags (per genomic block, and aggregated per transcript)
---------------------------------------------------------
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
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES = "human"
HG19_ASSEMBLY = "GRCh37"
HG38_ASSEMBLY = "GRCh38"

# Ensembl's public REST API documents a ~55,000 requests/hour quota
# (roughly 15 req/s); a small concurrent pool gets close to that without
# triggering its short-window burst limiter (HTTP 429).
DEFAULT_MAX_WORKERS = 4
SEQUENCE_BATCH_SIZE = 50
MAX_RETRIES = 5

HUMAN_SPECIES_ID = "9606"

# Aggregate-tag priority when combining multiple blocks of one transcript
# (worst wins): higher number = takes precedence.
_TAG_SEVERITY = {"ok": 0, "shifted": 1, "split": 2, "failed": 3}


def _request(url: str, data: Optional[bytes] = None) -> dict:
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(
                url, data=data, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                retry_after = exc.headers.get("Retry-After")
                time.sleep(float(retry_after) if retry_after else 1.0 + attempt * 0.5)
                last_error = exc
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(0.5 + attempt * 0.5)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


@dataclass
class Region:
    gene_id: str
    chrom: str
    start: int  # 1-based, inclusive
    end: int  # 1-based, inclusive
    strand: int = 1  # 1 or -1


@dataclass
class Transcript:
    """One transcript's 3' UTR, possibly spanning multiple genomic blocks."""

    transcript_id: str
    segments: List[Region] = field(default_factory=list)

    def ordered_segments(self) -> List[Region]:
        """Segments in 5'->3' transcript order."""
        strand = self.segments[0].strand
        return sorted(self.segments, key=lambda r: r.start, reverse=(strand == -1))


@dataclass
class LiftoverResult:
    gene_id: str
    hg19_region: str
    tag: str
    hg38_region: Optional[str] = None
    note: str = ""
    hg38_sequence: Optional[str] = None
    hg19_sequence: Optional[str] = None
    n_segments: int = 1


def map_region_to_hg38(region: Region) -> List[dict]:
    """Return Ensembl's list of mapped blocks for this hg19 region in hg38."""
    loc = f"{region.chrom}:{region.start}..{region.end}:{region.strand}"
    url = f"{ENSEMBL_REST}/map/{SPECIES}/{HG19_ASSEMBLY}/{loc}/{HG38_ASSEMBLY}?content-type=application/json"
    data = _request(url)
    return data.get("mappings", [])


def fetch_sequence(chrom: str, start: int, end: int, strand: int, assembly: str) -> str:
    loc = f"{chrom}:{start}..{end}:{strand}"
    url = (
        f"{ENSEMBL_REST}/sequence/region/{SPECIES}/{loc}"
        f"?content-type=application/json;coord_system_version={assembly}"
    )
    data = _request(url)
    return data["seq"]


def fetch_sequences_batch(locations: List[str], assembly: str) -> List[str]:
    """Fetch many regions' sequences in one POST request (<=50 at a time).

    ``locations`` are Ensembl-style ``chrom:start..end:strand`` strings.
    Returns sequences in the same order as the input (Ensembl's batch
    response is ordered the same as the request).
    """
    results: List[str] = []
    for i in range(0, len(locations), SEQUENCE_BATCH_SIZE):
        batch = locations[i : i + SEQUENCE_BATCH_SIZE]
        url = f"{ENSEMBL_REST}/sequence/region/{SPECIES}?coord_system_version={assembly};content-type=application/json"
        body = json.dumps({"regions": batch}).encode()
        data = _request(url, data=body)
        results.extend(item["seq"] for item in data)
    return results


# --------------------------------------------------------------------------
# Phase 1: coordinate mapping (concurrent, one call per genomic block)
# --------------------------------------------------------------------------


def _map_one_block(region: Region) -> dict:
    """Map a single genomic block; returns a dict describing the outcome
    (still needs sequence verification if status == "pending_ok")."""
    hg19_region_str = f"{region.chrom}:{region.start}-{region.end}"
    input_length = region.end - region.start + 1

    try:
        mappings = map_region_to_hg38(region)
    except Exception as exc:  # network/API failure -- treat like "no mapping"
        return {"region": region, "status": "failed", "note": f"API error: {exc}", "hg19_region": hg19_region_str}

    if not mappings:
        return {
            "region": region, "status": "failed", "hg19_region": hg19_region_str,
            "note": "No mapping returned (assembly gap or removed region)",
        }

    if len(mappings) > 1:
        blocks = "; ".join(
            f"{m['mapped']['seq_region_name']}:{m['mapped']['start']}-{m['mapped']['end']}" for m in mappings
        )
        return {
            "region": region, "status": "split", "hg19_region": hg19_region_str,
            "note": f"Region split across {len(mappings)} blocks in hg38: {blocks}",
        }

    mapped = mappings[0]["mapped"]
    hg38_chrom = mapped["seq_region_name"]
    hg38_start, hg38_end = mapped["start"], mapped["end"]
    hg38_strand = mapped["strand"]
    hg38_region_str = f"{hg38_chrom}:{hg38_start}-{hg38_end}"
    mapped_length = hg38_end - hg38_start + 1

    if mapped_length != input_length:
        return {
            "region": region, "status": "shifted", "hg19_region": hg19_region_str, "hg38_region": hg38_region_str,
            "note": f"Mapped length {mapped_length} != input length {input_length} (indel between assemblies)",
        }

    return {
        "region": region, "status": "pending_ok", "hg19_region": hg19_region_str, "hg38_region": hg38_region_str,
        "hg38_chrom": hg38_chrom, "hg38_start": hg38_start, "hg38_end": hg38_end, "hg38_strand": hg38_strand,
    }


def map_blocks_concurrent(regions: List[Region], max_workers: int = DEFAULT_MAX_WORKERS) -> List[dict]:
    if not regions:
        return []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(_map_one_block, regions))


# --------------------------------------------------------------------------
# Phase 2: batched sequence verification (only for "pending_ok" blocks)
# --------------------------------------------------------------------------


def verify_sequences(map_results: List[dict], fetch_sequences: bool = True) -> List[dict]:
    """For every "pending_ok" mapping, fetch+compare hg19 vs hg38 sequence
    (batched), finalizing each as "ok" or "shifted"."""
    pending = [r for r in map_results if r["status"] == "pending_ok"]

    if not pending:
        return map_results

    if not fetch_sequences:
        for r in pending:
            r["status"] = "ok"
            r["note"] = "Length matches (sequence not checked)"
        return map_results

    hg19_locs = [f"{r['region'].chrom}:{r['region'].start}..{r['region'].end}:{r['region'].strand}" for r in pending]
    hg38_locs = [f"{r['hg38_chrom']}:{r['hg38_start']}..{r['hg38_end']}:{r['hg38_strand']}" for r in pending]

    try:
        hg19_seqs = fetch_sequences_batch(hg19_locs, HG19_ASSEMBLY)
        hg38_seqs = fetch_sequences_batch(hg38_locs, HG38_ASSEMBLY)
    except Exception as exc:
        for r in pending:
            r["status"] = "shifted"
            r["note"] = f"Could not fetch sequence to verify: {exc}"
        return map_results

    for r, hg19_seq, hg38_seq in zip(pending, hg19_seqs, hg38_seqs):
        r["hg19_sequence"] = hg19_seq
        r["hg38_sequence"] = hg38_seq
        if hg19_seq.upper() == hg38_seq.upper():
            r["status"] = "ok"
            r["note"] = "Sequence identical; safe to re-anchor coordinates only"
        else:
            r["status"] = "shifted"
            r["note"] = "Coordinates map cleanly but sequence content differs between assemblies"

    return map_results


def _block_result_to_liftover_result(r: dict) -> LiftoverResult:
    return LiftoverResult(
        gene_id=r["region"].gene_id,
        hg19_region=r["hg19_region"],
        tag=r["status"],
        hg38_region=r.get("hg38_region"),
        note=r.get("note", ""),
        hg38_sequence=r.get("hg38_sequence"),
        hg19_sequence=r.get("hg19_sequence"),
    )


def liftover_region(region: Region, fetch_sequences: bool = True) -> LiftoverResult:
    """Lift a single genomic block. Convenience wrapper for callers that
    don't need the concurrent/batched machinery (tests, small scripts)."""
    [mapped] = map_blocks_concurrent([region], max_workers=1)
    [verified] = verify_sequences([mapped], fetch_sequences=fetch_sequences)
    return _block_result_to_liftover_result(verified)


def liftover_regions(regions: List[Region], fetch_sequences: bool = True) -> List[LiftoverResult]:
    """Lift a flat list of independent (single-block) regions."""
    mapped = map_blocks_concurrent(regions, max_workers=1)
    verified = verify_sequences(mapped, fetch_sequences=fetch_sequences)
    return [_block_result_to_liftover_result(r) for r in verified]


# --------------------------------------------------------------------------
# Multi-segment transcripts
# --------------------------------------------------------------------------


def liftover_transcripts(
    transcripts: List[Transcript],
    max_workers: int = DEFAULT_MAX_WORKERS,
    fetch_sequences: bool = True,
) -> Dict[str, LiftoverResult]:
    """Lift every block of every transcript (concurrently + batched), then
    aggregate each transcript's blocks into one LiftoverResult.

    A transcript's tag is the worst of its blocks' tags. If "ok", its
    hg38_sequence is the per-block hg38 sequences concatenated in 5'->3'
    transcript order (so it's a drop-in ungapped replacement sequence).
    """
    all_regions: List[Region] = []
    for t in transcripts:
        all_regions.extend(t.ordered_segments())

    mapped = map_blocks_concurrent(all_regions, max_workers=max_workers)
    verified = verify_sequences(mapped, fetch_sequences=fetch_sequences)

    result_by_region_id = {id(r["region"]): r for r in verified}

    results: Dict[str, LiftoverResult] = {}
    for t in transcripts:
        ordered = t.ordered_segments()
        block_results = [result_by_region_id[id(r)] for r in ordered]
        tags = [b["status"] for b in block_results]
        worst_tag = max(tags, key=lambda tag: _TAG_SEVERITY[tag])

        hg19_region_summary = "; ".join(b["hg19_region"] for b in block_results)

        if worst_tag == "ok":
            hg38_sequence = "".join(b["hg38_sequence"] for b in block_results)
            hg19_sequence = "".join(b["hg19_sequence"] for b in block_results)
            hg38_region_summary = "; ".join(b["hg38_region"] for b in block_results)
            note = "Sequence identical; safe to re-anchor coordinates only" + (
                f" ({len(ordered)} exonic blocks)" if len(ordered) > 1 else ""
            )
            results[t.transcript_id] = LiftoverResult(
                t.transcript_id, hg19_region_summary, "ok", hg38_region_summary, note,
                hg38_sequence=hg38_sequence, hg19_sequence=hg19_sequence, n_segments=len(ordered),
            )
        else:
            bad_blocks = [b for b in block_results if b["status"] == worst_tag]
            note = f"{len(bad_blocks)}/{len(ordered)} block(s) tagged {worst_tag}: " + " | ".join(
                b.get("note", "") for b in bad_blocks
            )
            results[t.transcript_id] = LiftoverResult(
                t.transcript_id, hg19_region_summary, worst_tag, note=note, n_segments=len(ordered),
            )

    return results


def read_transcripts_gff(path: str) -> List[Transcript]:
    """Read TargetScan's hg19 3' UTR GFF (e.g. ``TSHuman_7_hg19_3UTRs.gff``),
    grouping rows by transcript ID into multi-block Transcripts."""
    by_id: Dict[str, Transcript] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("browser") or line.startswith("track") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, _feature, start, end, _score, strand_field, _frame, transcript_id = fields[:9]
            strand = -1 if strand_field == "-" else 1
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            region = Region(transcript_id, chrom, int(start), int(end), strand)
            by_id.setdefault(transcript_id, Transcript(transcript_id)).segments.append(region)
    return list(by_id.values())


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


# --------------------------------------------------------------------------
# Splicing results into an existing UTR_Sequences.txt-style alignment file
# --------------------------------------------------------------------------


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
    utr_file: str, results: Iterable[LiftoverResult], out_path: str, applied_report_path: str
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


def write_report(results: Iterable[LiftoverResult], out_path: str) -> None:
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id", "hg19_region", "tag", "hg38_region", "n_segments", "note"])
        for r in results:
            writer.writerow([r.gene_id, r.hg19_region, r.tag, r.hg38_region or "", r.n_segments, r.note])


# --------------------------------------------------------------------------
# Resumable full-scale runner (designed for tens of thousands of transcripts)
# --------------------------------------------------------------------------


def run_gff_liftover(
    gff_path: str,
    out_path: str,
    batch_size: int = 200,
    max_workers: int = DEFAULT_MAX_WORKERS,
    fetch_sequences: bool = True,
    progress_every: int = 1,
) -> None:
    """Lift every transcript in a TargetScan hg19 3' UTR GFF file, writing
    results incrementally (one row per transcript, flushed after every
    batch) so the run can be safely interrupted and resumed: transcripts
    already present in ``out_path`` are skipped on restart.
    """
    all_transcripts = read_transcripts_gff(gff_path)

    done_ids: set = set()
    write_header = True
    try:
        with open(out_path) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            write_header = header is None
            for row in reader:
                if row:
                    done_ids.add(row[0])
    except FileNotFoundError:
        pass

    todo = [t for t in all_transcripts if t.transcript_id not in done_ids]
    print(
        f"{len(all_transcripts)} transcripts total, {len(done_ids)} already done, {len(todo)} remaining.",
        file=sys.stderr,
    )

    with open(out_path, "a", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if write_header:
            writer.writerow(["gene_id", "hg19_region", "tag", "hg38_region", "n_segments", "note"])
            fh.flush()

        n_done = 0
        for i in range(0, len(todo), batch_size):
            batch = todo[i : i + batch_size]
            results = liftover_transcripts(batch, max_workers=max_workers, fetch_sequences=fetch_sequences)
            for t in batch:
                r = results[t.transcript_id]
                writer.writerow([r.gene_id, r.hg19_region, r.tag, r.hg38_region or "", r.n_segments, r.note])
            fh.flush()
            n_done += len(batch)
            if (i // batch_size) % progress_every == 0:
                tags_so_far = ", ".join(f"{tag}={sum(1 for t in batch if results[t.transcript_id].tag == tag)}" for tag in _TAG_SEVERITY)
                print(f"[{n_done}/{len(todo)}] batch tags: {tags_so_far}", file=sys.stderr)

    print(f"Done. Results in {out_path}", file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        print(
            "USAGE: hg38_liftover.py hg19_regions.bed report_out.tsv [--no-sequence-check]\n"
            "       [--utr-file UTR_Sequences.txt --utr-out UTR_Sequences.hg38.txt --applied-out applied.tsv]\n"
            "   or: hg38_liftover.py --gff TSHuman_7_hg19_3UTRs.gff report_out.tsv [--workers N] [--batch-size N]",
            file=sys.stderr,
        )
        return 0

    if argv[0] == "--gff":
        gff_path, out_path = argv[1], argv[2]
        max_workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else DEFAULT_MAX_WORKERS
        batch_size = int(argv[argv.index("--batch-size") + 1]) if "--batch-size" in argv else 200
        run_gff_liftover(gff_path, out_path, batch_size=batch_size, max_workers=max_workers)
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
