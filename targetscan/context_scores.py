"""Calculate TargetScan context++ scores (Agarwal et al. 2015 model).

Python port of
``tsh_orig/TargetScan7_context_scores/targetscan_70_context_scores.pl``.
This is the TargetScan 8 (context++) scoring model already embedded in
that script (3' UTR length, ORF length, ORF 8mer count, offset 6mer count,
min-distance-to-end, and site-accessibility -- via RNAplfold -- terms were
all added "for TargetScan 8"), so no algorithm changes are needed to score
vert80/hg38 data; only the downloaded data files differ.

Several quirks/edge-case behaviors of the original Perl (including ones
that look like bugs) are intentionally preserved -- see inline comments --
because the goal is byte-for-byte compatible output, not a "corrected"
reimplementation.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Site type numbering used *only* within this module (distinct from the
# numbering used in site_prediction.py) -- matches the Perl script exactly.
SITE_TYPE_NUM = {"7mer-1a": 1, "7mer-m8": 2, "8mer-1a": 3, "6mer": 4}
SITE_TYPE_NAME = {v: k for k, v in SITE_TYPE_NUM.items()}

MIN_DIST_TO_CDS = 15
TOO_CLOSE_TO_CDS = "too_close"
DESIRED_UTR_ALIGNMENT_LENGTH = 23
DIGITS_AFTER_DECIMAL = 3

MAX_CONTEXT_SCORE = {1: -0.01, 2: -0.02, 3: -0.03, 4: 0}

SET_AIRS_TO_1 = False
MIN_AIR = 0.0001

REF_SPECIES = "9606"
SPECIES = ["10090", "10116", "13616", "8364", "9031", "9544", "9598", "9606", "9615", "9913"]

_NUMBER_RE = re.compile(r"^[+-]?(?=\d|\.\d)\d*(\.\d*)?([Ee][+-]?\d+)?$")


def _looks_like_number(value) -> bool:
    if value is None:
        return False
    return bool(_NUMBER_RE.match(str(value)))


def is_nonzero_number(value) -> bool:
    return _looks_like_number(value) and float(value) != 0


def _format_number(value: float) -> str:
    """Format like Perl's default number stringification (up to 15 sig figs)."""
    if value == int(value):
        return str(int(value))
    return f"{value:.15g}"


def _safe_float(value) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


# --------------------------------------------------------------------------
# Reference-data loading
# --------------------------------------------------------------------------


@dataclass
class AgarwalParams:
    intercept: Dict[int, float] = field(default_factory=dict)
    coeff: Dict[int, Dict[str, float]] = field(default_factory=dict)
    minv: Dict[int, Dict[str, float]] = field(default_factory=dict)
    maxv: Dict[int, Dict[str, float]] = field(default_factory=dict)

    @classmethod
    def read(cls, path: str) -> "AgarwalParams":
        obj = cls()
        with open(path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                fields = line.rstrip("\r\n").split("\t")
                feature = fields[0]
                if not feature:
                    continue
                vals = fields[1:13] + [""] * max(0, 13 - len(fields))
                coeffs = vals[0:4]
                mins = vals[4:8]
                maxs = vals[8:12]
                # Column order is 8mer(3), 7mer-m8(2), 7mer-A1(1), 6mer(4)
                site_types = [3, 2, 1, 4]
                if i == 1:  # Intercept row
                    for st, c in zip(site_types, coeffs):
                        c = _safe_float(c)
                        if c is not None:
                            obj.intercept[st] = c
                else:
                    for st, c, mn, mx in zip(site_types, coeffs, mins, maxs):
                        c = _safe_float(c)
                        mn = _safe_float(mn)
                        mx = _safe_float(mx)
                        if c is not None:
                            obj.coeff.setdefault(st, {})[feature] = c
                        if mn is not None:
                            obj.minv.setdefault(st, {})[feature] = mn
                        if mx is not None:
                            obj.maxv.setdefault(st, {})[feature] = mx
        return obj


def get_agarwal_contribution(agarwal: AgarwalParams, site_type: int, contribution_type: str, raw_score) -> str:
    if not _looks_like_number(raw_score):
        print(
            f"In getAgarwalContribution({site_type}, {contribution_type}, {raw_score}), "
            "score (last argument) is not a number",
            file=sys.stderr,
        )
        return "0.0000000"

    raw = float(raw_score)
    scale_these = {
        "TA_3UTR",
        "SPS",
        "Local_AU",
        "3P_score",
        "SA",
        "Len_ORF",
        "Len_3UTR",
        "Min_dist",
        "PCT",
    }
    if contribution_type in scale_these:
        mn = agarwal.minv[site_type][contribution_type]
        mx = agarwal.maxv[site_type][contribution_type]
        scaled = (raw - mn) / (mx - mn)
    else:
        scaled = raw

    coeff = agarwal.coeff.get(site_type, {}).get(contribution_type)
    if coeff is None:
        print(f"PROBLEM: No coeff for {contribution_type} (site type = {site_type})", file=sys.stderr)
        coeff = 0.0

    contribution = coeff * scaled
    return f"{contribution:.{DIGITS_AFTER_DECIMAL}f}"


def read_utrs(utr_file: str, use_species: set) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, float]], set]:
    utr_seq: Dict[str, Dict[str, str]] = {}
    have_utrs_this_species: set = set()

    with open(utr_file) as fh:
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            transcript_id, species_id, seq = fields[0], fields[1], fields[2]
            have_utrs_this_species.add(species_id)
            if species_id not in use_species:
                continue
            seq = seq.replace("-", "")
            if seq:
                seq = re.sub(r"[Tt]", "U", seq)
                utr_seq.setdefault(transcript_id, {})[species_id] = seq

    utr_length_scaling_factor: Dict[str, Dict[str, float]] = {}
    for transcript_id, by_species in utr_seq.items():
        ref_len = len(by_species.get(REF_SPECIES, ""))
        for species_id, seq in by_species.items():
            if len(seq):
                utr_length_scaling_factor.setdefault(transcript_id, {})[species_id] = ref_len / len(seq)

    return utr_seq, utr_length_scaling_factor, have_utrs_this_species


def read_mirnas_context(mirna_file: str) -> Tuple[Dict[Tuple[str, str], List[Tuple[str, str]]], Dict[str, str]]:
    mature_seq: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}
    family_to_seed: Dict[str, str] = {}

    with open(mirna_file) as fh:
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            fam_id, species_id, mature_id, mature_seq_str = fields[0], fields[1], fields[2], fields[3]
            mature_seq_str = mature_seq_str.upper()
            seed_region = mature_seq_str[1:8]

            if fam_id not in family_to_seed:
                family_to_seed[fam_id] = seed_region
            elif seed_region != family_to_seed[fam_id]:
                print(
                    f"ERROR: miRNA family {fam_id} seems to have more than 1 seed region: "
                    f"{family_to_seed[fam_id]}, {seed_region}",
                    file=sys.stderr,
                )
                print("Please correct family definitions and re-run analysis", file=sys.stderr)
                sys.exit(1)

            mature_seq.setdefault((fam_id, species_id), []).append((mature_id, mature_seq_str))

    return mature_seq, family_to_seed


def read_isoform_ratios(
    air_file: str,
    utr_length_scaling_factor: Dict[str, Dict[str, float]],
    have_utrs_this_species: set,
) -> Tuple[Dict[str, Dict[str, Dict[int, float]]], Dict[str, Dict[str, Dict[int, int]]]]:
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]] = {}
    utr_scaled_end_to_start: Dict[str, Dict[str, Dict[int, int]]] = {}

    with open(air_file) as fh:
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            utr, start, end, air = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
            if SET_AIRS_TO_1:
                air = 100.0

            for species_id in have_utrs_this_species:
                scale = utr_length_scaling_factor.get(utr, {}).get(species_id)
                if scale is None:
                    continue
                scaled_end = round(end / scale)
                scaled_start = round(start / scale) if start > 1 else 1
                utr_scaled_end_to_air.setdefault(utr, {}).setdefault(species_id, {})[scaled_end] = air
                utr_scaled_end_to_start.setdefault(utr, {}).setdefault(species_id, {})[scaled_end] = scaled_start

    return utr_scaled_end_to_air, utr_scaled_end_to_start


def read_orf_data(
    orf_lengths_file: str, orf_8mers_file: str, use_species: set
) -> Tuple[Dict[str, Dict[str, int]], Dict[str, Dict[str, Dict[str, int]]]]:
    orf2length: Dict[str, Dict[str, int]] = {}
    orf_to_8mer_counts: Dict[str, Dict[str, Dict[str, int]]] = {}

    with open(orf_lengths_file) as fh:
        for line in fh:
            orf, species, length = line.rstrip("\r\n").split("\t")
            if species in use_species:
                orf2length.setdefault(orf, {})[species] = int(length)

    with open(orf_8mers_file) as fh:
        for line in fh:
            orf, species, family, count = line.rstrip("\r\n").split("\t")
            if species in use_species:
                orf_to_8mer_counts.setdefault(orf, {}).setdefault(species, {})[family] = int(count)

    return orf2length, orf_to_8mer_counts


def read_ta_sps(path: str) -> Tuple[Dict[int, Dict[str, float]], Dict[str, float]]:
    garcia_sps: Dict[int, Dict[str, float]] = {1: {}, 2: {}, 3: {}, 4: {}}
    garcia_ta: Dict[str, float] = {}

    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            fields = line.rstrip("\r\n").split("\t")
            seed_region, sps_1, sps_2, ta = fields[0], float(fields[1]), float(fields[2]), float(fields[3])
            garcia_sps[3][seed_region] = sps_1
            garcia_sps[2][seed_region] = sps_1
            garcia_sps[1][seed_region] = sps_2
            garcia_sps[4][seed_region] = sps_2
            garcia_ta[seed_region] = ta

    return garcia_sps, garcia_ta


def get_replace_length(longest: int) -> Dict[int, str]:
    return {i: " " * i for i in range(1, longest + 1)}


# --------------------------------------------------------------------------
# RNAplfold (site accessibility)
# --------------------------------------------------------------------------


def rnaplfold_installed() -> bool:
    try:
        subprocess.run(
            ["RNAplfold", "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return True
    except FileNotFoundError:
        return False


def run_rnaplfold_all_utrs(utr_seq: Dict[str, Dict[str, str]], rnaplfold_dir: str) -> bool:
    if not rnaplfold_installed():
        print(
            "Failed to run RNAplfold; We'll keep going but will not calculate the SA contribution.",
            file=sys.stderr,
        )
        return False

    os.makedirs(rnaplfold_dir, exist_ok=True)
    print("Running RNAplfold on UTRs (if we didn't do so before)....", file=sys.stderr)

    cwd = os.getcwd()
    os.chdir(rnaplfold_dir)
    try:
        for transcript_id in sorted(utr_seq.keys()):
            for species in sorted(utr_seq[transcript_id].keys()):
                lunp_file = f"{transcript_id}.{species}_lunp"
                if os.path.exists(lunp_file):
                    continue
                fasta_file = f"{transcript_id}.{species}.fa"
                with open(fasta_file, "a") as fh:
                    fh.write(f">{transcript_id}.{species}\n{utr_seq[transcript_id][species]}\n")
                with open(fasta_file) as fasta_in:
                    subprocess.run(
                        ["RNAplfold", "-L", "40", "-W", "80", "-u", "20"],
                        stdin=fasta_in,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                        check=False,
                    )
                ps_file = f"{transcript_id}.{species}_dp.ps"
                if os.path.exists(ps_file):
                    os.remove(ps_file)
    finally:
        os.chdir(cwd)
    return True


_LUNP_CACHE: Dict[str, Dict[int, List[str]]] = {}


def _read_lunp(path: str) -> Dict[int, List[str]]:
    if path in _LUNP_CACHE:
        return _LUNP_CACHE[path]
    rows: Dict[int, List[str]] = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            if not fields or not fields[0].lstrip("-").isdigit():
                continue
            rows[int(fields[0])] = fields[1:]
    _LUNP_CACHE[path] = rows
    return rows


def get_sa_contribution(
    agarwal: AgarwalParams,
    rnaplfold_dir: str,
    transcript_id: str,
    species_id: str,
    utr_start: int,
    site_type: int,
    missing_files: set,
) -> str:
    if site_type in (1, 5):
        utr_start -= 1

    lunp_path = os.path.join(rnaplfold_dir, f"{transcript_id}.{species_id}_lunp")

    if not os.path.exists(lunp_path):
        missing_files.add(lunp_path)
        return "0"

    rows = _read_lunp(lunp_path)
    if not rows:
        return "0"
    max_row = max(rows.keys())
    if utr_start not in rows:
        return "0"
    target_row = min(utr_start + 7, max_row)
    values = rows[target_row]
    plfold = values[13] if len(values) > 13 else "NA"

    if not plfold or plfold == "NA":
        return "0"

    plfold_f = _safe_float(plfold)
    if plfold_f and is_nonzero_number(plfold_f):
        log10_plfold = math.log10(plfold_f)
    else:
        log10_plfold = 0.0

    return get_agarwal_contribution(agarwal, site_type, "SA", log10_plfold)


# --------------------------------------------------------------------------
# Per-site contributions
# --------------------------------------------------------------------------


def get_air_this_site(
    transcript_id: str,
    species_id: str,
    site_end: int,
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]],
    utr_scaled_end_to_start: Dict[str, Dict[str, Dict[int, int]]],
) -> float:
    air_by_end = utr_scaled_end_to_air.get(transcript_id, {}).get(species_id)
    if not air_by_end:
        print(
            f"No isoform ratios (AIRs) for UTR {transcript_id} of species {species_id}.  Setting AIR to 1.",
            file=sys.stderr,
        )
        return 1.0

    start_by_end = utr_scaled_end_to_start[transcript_id][species_id]
    closest_air_region = 10_000_000
    air_closest_region = None

    for air_end, air in air_by_end.items():
        air_start = start_by_end[air_end]
        if air_end >= site_end and air_start <= site_end:
            return 1.0 * air / 100
        dist = min(abs(air_end - site_end), abs(air_start - site_end))
        if dist < closest_air_region:
            closest_air_region = dist
            air_closest_region = air

    if not air_closest_region:
        air_closest_region = 100
    air_this_region = 1.0 * air_closest_region / 100
    if air_this_region < MIN_AIR:
        air_this_region = MIN_AIR
    return air_this_region


def extract_subseq_for_alignment(
    utr_seq: Dict[str, Dict[str, str]], transcript_id: str, species_id: str, utr_start: int, utr_end: int, site_type: int
) -> str:
    if site_type < 5:
        real_start = utr_start - 16
    else:
        real_start = utr_start - 19
    if real_start < 0:
        real_start = 0

    real_end = utr_end + 1 if site_type in (1, 2, 4) else utr_end
    if real_start >= real_end:
        real_start = 0

    length = real_end - real_start
    seq = utr_seq.get(transcript_id, {}).get(species_id, "")
    subseq = seq[real_start : real_start + length]

    subseq_len = len(subseq)
    spacer = ""
    if DESIRED_UTR_ALIGNMENT_LENGTH > subseq_len:
        spacer = "N" * (DESIRED_UTR_ALIGNMENT_LENGTH - subseq_len + 1)
        subseq = spacer + subseq

    if len(seq) < real_start + length:
        length_diff = real_start + length - len(seq)
        for _ in range(length_diff):
            subseq += "N"
            subseq = re.sub(r"^NN", "  ", subseq)

    return subseq


def modify_subseq_for_alignment(mature_mirna: str, subseq_for_alignment: str, site_type: int) -> Tuple[str, str]:
    subseq_len = len(subseq_for_alignment)
    spacer1_len = len(mature_mirna) - DESIRED_UTR_ALIGNMENT_LENGTH
    spacer1 = " " * max(spacer1_len, 0)
    spacer2_len = -spacer1_len if spacer1_len < 0 else 0
    spacer2 = " " * spacer2_len

    if DESIRED_UTR_ALIGNMENT_LENGTH > subseq_len:
        spacer2 += " " * (DESIRED_UTR_ALIGNMENT_LENGTH - subseq_len + 1)

    if site_type in (1, 3):
        spacer2 = spacer2[:-1]

    final_subseq = spacer1 + subseq_for_alignment
    mature_for_alignment = spacer2 + mature_mirna[::-1]

    return final_subseq, mature_for_alignment


def get_local_au_contribution(
    agarwal: AgarwalParams,
    utr_seq: Dict[str, Dict[str, str]],
    transcript_id: str,
    species_id: str,
    site_type: int,
    utr_start: int,
    utr_end: int,
) -> str:
    seq = utr_seq.get(transcript_id, {}).get(species_id, "")
    seq_len = len(seq)

    utr_subseq_len = 30
    utr_up_start = utr_start - 31
    if utr_up_start < 0:
        utr_up_start = 0
        utr_subseq_len = utr_start - 1

    utr_down_start = utr_end

    utr_up = seq[utr_up_start : utr_up_start + max(utr_subseq_len, 0)]
    utr_down = seq[utr_down_start : utr_down_start + 30]

    total_up = 0.0
    total_down = 0.0
    max_raw = 0.0

    up_3to5 = utr_up.upper()[::-1]
    for i, ch in enumerate(up_3to5):
        score = 1 / (i + 1) if site_type in (2, 3) else 1 / (i + 2)
        if ch in "UA":
            total_up += score
        max_raw += score

    down_5to3 = utr_down.upper()
    for i, ch in enumerate(down_5to3):
        score = 1 / (i + 2) if site_type in (3, 1) else 1 / (i + 1)
        if ch in "UA":
            total_down += score
        max_raw += score

    raw = total_up + total_down
    fraction = round(raw / max_raw, DIGITS_AFTER_DECIMAL) if max_raw != 0 else 0.0

    return get_agarwal_contribution(agarwal, site_type, "Local_AU", fraction)


def get_min_dist_weighted_contribution(
    agarwal: AgarwalParams,
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]],
    transcript_id: str,
    species_id: str,
    site_type: int,
    site_start: int,
    site_end: int,
) -> str:
    nearest_ends: List[float] = []
    nearest_end_airs: List[float] = []

    for air_end, air in utr_scaled_end_to_air.get(transcript_id, {}).get(species_id, {}).items():
        if air_end >= site_end:
            dist5 = site_start - 1
            dist3 = air_end - site_end
            nearest_ends.append(dist5 if dist5 <= dist3 else dist3)
            nearest_end_airs.append(air)

    airs_sum = sum(nearest_end_airs)
    product_sum = sum(d * a for d, a in zip(nearest_ends, nearest_end_airs))
    weighted_mean = round(product_sum / airs_sum) if airs_sum != 0 else 0

    log10_dist = math.log10(weighted_mean) if is_nonzero_number(weighted_mean) else 0.0
    return get_agarwal_contribution(agarwal, site_type, "Min_dist", log10_dist)


def get_len3utr_weighted_contribution(
    agarwal: AgarwalParams,
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]],
    transcript_id: str,
    species_id: str,
    site_type: int,
    site_end: int,
) -> str:
    utr_lengths: List[int] = []
    utr_airs: List[float] = []

    for air_end, air in utr_scaled_end_to_air.get(transcript_id, {}).get(species_id, {}).items():
        if air_end >= site_end:
            utr_lengths.append(air_end)
            # Preserve original behavior: an AIR of 0 is "falsy" in Perl, so
            # it is *not* appended here, which can misalign the two lists.
            if air:
                utr_airs.append(air)

    utr_length_sum = 0.0
    airs_sum = 0.0
    for i in range(len(utr_lengths)):
        air_i = utr_airs[i] if i < len(utr_airs) else 0
        utr_length_sum += utr_lengths[i] * air_i
        airs_sum += air_i

    utr_length = round(utr_length_sum / airs_sum) if airs_sum != 0 else 0

    log10_len = math.log10(utr_length) if is_nonzero_number(utr_length) else 0.0
    return get_agarwal_contribution(agarwal, site_type, "Len_3UTR", log10_len)


def get_orflength_contribution(
    agarwal: AgarwalParams, orf2length: Dict[str, Dict[str, int]], transcript_id: str, species_id: str, site_type: int
) -> str:
    orf_length = orf2length.get(transcript_id, {}).get(species_id)
    log10_len = math.log10(orf_length) if is_nonzero_number(orf_length) else 0.0
    return get_agarwal_contribution(agarwal, site_type, "Len_ORF", log10_len)


def get_orf8mer_contribution(
    agarwal: AgarwalParams,
    orf_to_8mer_counts: Dict[str, Dict[str, Dict[str, int]]],
    family_to_seed: Dict[str, str],
    transcript_id: str,
    mirna_family_id: str,
    species_id: str,
    site_type: int,
) -> str:
    seed = family_to_seed.get(mirna_family_id)
    count = orf_to_8mer_counts.get(transcript_id, {}).get(species_id, {}).get(seed, 0)
    return get_agarwal_contribution(agarwal, site_type, "ORF8m", count)


_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def _reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_COMPLEMENT)


def get_offset6mer_sites(
    utr_seq: Dict[str, Dict[str, str]],
    family_to_seed: Dict[str, str],
    cache: Dict[Tuple[str, str, str], List[int]],
    transcript_id: str,
    species_id: str,
    mirna_family_id: str,
) -> None:
    """Find offset-6mer sites for a (transcript, species, miRNA family).

    Preserved bug: the original Perl searches ``$utr_seq{$transcriptID}``
    (missing the ``{$speciesID}`` subscript), so it actually matches
    against the *stringified hash reference* (something like
    "HASH(0x55b2...)") rather than a real sequence. Since real seed-derived
    site strings always contain a G or U (never present in that
    stringified-hashref text), the search effectively never finds a real
    site -- but Perl's resulting sentinel value (``-1``) is then treated as
    truthy and, due to ``-1 <= $AIR_end`` always being true, contributes a
    constant offset-6mer count of 1 to every site. We replicate that
    observed, constant behavior directly (verified against the bundled
    sample context-score output) rather than trying to reproduce an
    unreproducible memory address.
    """
    cache[(transcript_id, species_id, mirna_family_id)] = [-1]


def get_offset6mer_weighted_contribution(
    agarwal: AgarwalParams,
    utr_seq: Dict[str, Dict[str, str]],
    family_to_seed: Dict[str, str],
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]],
    cache: Dict[Tuple[str, str, str], List[int]],
    transcript_id: str,
    mirna_family_id: str,
    species_id: str,
    site_type: int,
    site_end: int,
) -> str:
    key = (transcript_id, species_id, mirna_family_id)
    if key not in cache:
        get_offset6mer_sites(utr_seq, family_to_seed, cache, transcript_id, species_id, mirna_family_id)

    offset6mer_sites = cache.get(key, [])

    if offset6mer_sites:
        counts: List[int] = []
        airs: List[float] = []
        for air_end, air in utr_scaled_end_to_air.get(transcript_id, {}).get(species_id, {}).items():
            if air_end >= site_end:
                count_this_utr = sum(1 for pos in offset6mer_sites if pos <= air_end)
                counts.append(count_this_utr)
                airs.append(air)

        airs_sum = sum(airs)
        product_sum = sum(c * a for c, a in zip(counts, airs))
        offset6mer_count = round(product_sum / airs_sum) if airs_sum != 0 else 0
    else:
        offset6mer_count = 0

    return get_agarwal_contribution(agarwal, site_type, "Off6m", offset6mer_count)


def get_sRNA1_8_contributions(agarwal: AgarwalParams, mirna: str, site_type: int) -> Tuple[str, str, str, str, str, str]:
    s1 = mirna[0:1]
    s8 = mirna[7:8]

    is_pos1 = {"A": 0, "C": 0, "G": 0}
    if s1 in is_pos1:
        is_pos1[s1] = 1
    is_pos8 = {"A": 0, "C": 0, "G": 0}
    if s8 in is_pos8:
        is_pos8[s8] = 1

    s1a = s1c = s1g = s8a = s8c = s8g = "0"

    if s1 != "U":
        s1a = get_agarwal_contribution(agarwal, site_type, "sRNA1A", is_pos1["A"])
        s1c = get_agarwal_contribution(agarwal, site_type, "sRNA1C", is_pos1["C"])
        s1g = get_agarwal_contribution(agarwal, site_type, "sRNA1G", is_pos1["G"])
    if s8 != "U":
        s8a = get_agarwal_contribution(agarwal, site_type, "sRNA8A", is_pos8["A"])
        s8c = get_agarwal_contribution(agarwal, site_type, "sRNA8C", is_pos8["C"])
        s8g = get_agarwal_contribution(agarwal, site_type, "sRNA8G", is_pos8["G"])

    return s1a, s1c, s1g, s8a, s8c, s8g


def get_site8_contribution(agarwal: AgarwalParams, subseq_for_alignment: str, site_type: int) -> Tuple[str, str, str]:
    site_pos8 = subseq_for_alignment[14:15].upper()

    is_site_pos8 = {"A": 0, "C": 0, "G": 0}
    if site_pos8 in is_site_pos8:
        is_site_pos8[site_pos8] = 1

    if site_pos8 == "U":
        return "0", "0", "0"

    a = get_agarwal_contribution(agarwal, site_type, "Site8A", is_site_pos8["A"])
    c = get_agarwal_contribution(agarwal, site_type, "Site8C", is_site_pos8["C"])
    g = get_agarwal_contribution(agarwal, site_type, "Site8G", is_site_pos8["G"])
    return a, c, g


def get_pct_contribution(agarwal: AgarwalParams, site_type: int, pct) -> str:
    return get_agarwal_contribution(agarwal, site_type, "PCT", pct)


def get_weighted_context_score(context_score, air: float) -> str:
    score = _safe_float(context_score)
    if score is None:
        score = 0.0
    weighted = math.log((2 ** score - 1) * air + 1) / math.log(2)
    return f"{weighted:.4f}"


# --------------------------------------------------------------------------
# 3' pairing (consequential pairing) -- direct, literal port
# --------------------------------------------------------------------------

_SEED_INFO = {
    "utrstart": {1: 8, 2: 8, 3: 8, 4: 8},
    "mirnastart": {1: 7, 2: 8, 3: 8, 4: 8},
    "offset": {1: 1, 2: 0, 3: 1, 4: 2},
    "overhang": {1: 1, 2: 0, 3: 0, 4: 0},
    "seedspan": {1: 6, 2: 7, 3: 7, 4: 6},
}

_NT_TO_NUM = str.maketrans("AUCGN", "12345")


def _match_score(a: str, b: str) -> bool:
    product = int(a) * int(b) if a.isdigit() and b.isdigit() else 0
    return product in (2, 12)


def _scan_alignment(utr_nums: str, mirna_nums: str, offset: int, type_: int, on_top: bool) -> Tuple[float, str]:
    score = 0.0
    string = ""
    tempstring = ""
    bestmatch = 0
    tempscore = 0.0
    prevmatch = 0

    if on_top:
        max_i = min(len(mirna_nums) - 1 - offset, len(utr_nums) - 1)
        pairs = [(utr_nums[i], mirna_nums[i + offset], i + offset) for i in range(max_i + 1)] if max_i >= 0 else []
    else:
        max_i = min(len(utr_nums) - 1 - offset, len(mirna_nums) - 1)
        pairs = [(utr_nums[i + offset], mirna_nums[i], i) for i in range(max_i + 1)] if max_i >= 0 else []

    for u_ch, m_ch, pos_for_overhang in pairs:
        if _match_score(u_ch, m_ch):
            in_seed = 4 <= (pos_for_overhang - _SEED_INFO["overhang"][type_]) <= 7
            tempstring += "|"
            if prevmatch == 0:
                tempscore = 0
            tempscore += 1 if in_seed else 0.5
            prevmatch += 1
        elif prevmatch >= 2:
            if tempscore == score:
                string += tempstring
            elif tempscore > score:
                bestmatch = prevmatch
                string = re.sub(r"[|X]", " ", string)
                string += tempstring
                score = tempscore
            else:
                tempstring = re.sub(r"[|X]", " ", tempstring)
                string += tempstring
            string += " "
            tempstring = ""
            tempscore = 0
            prevmatch = 0
        else:
            tempstring = re.sub(r"[|X]", " ", tempstring)
            string += tempstring
            string += " "
            tempstring = ""
            tempscore = 0
            prevmatch = 0

    if prevmatch >= 2:
        if tempscore == score:
            string += tempstring
        elif tempscore > score:
            bestmatch = prevmatch
            string = re.sub(r"[|X]", " ", string)
            string += tempstring
            score = tempscore

    score = score - max(0.0, (offset - 2) / 2)
    string = re.sub(r"\s([|X])\s", "   ", string)
    string = re.sub(r"^([|X])\s", "  ", string)
    string = re.sub(r"\s([|X])$", "  ", string)
    return score, string


def get3_prime_pairing_contribution(
    agarwal: AgarwalParams, site_type: int, utr: str, mirna: str
) -> Tuple[str, str, str, str, float]:
    utr = re.sub(r"[ \n]", "", utr).upper().replace("T", "U")
    mirna = re.sub(r"[ \n]", "", mirna).upper().replace("T", "U")

    seedinfo = _SEED_INFO
    utr_rev = utr[::-1]
    mirna_rev = mirna[::-1]

    utr_num_str = utr_rev[seedinfo["utrstart"][site_type] :]
    mirna_num_str = mirna_rev[seedinfo["mirnastart"][site_type] :]
    maxscore = max(len(utr_num_str), len(mirna_num_str))

    utr_nums = utr_num_str.translate(_NT_TO_NUM)
    mirna_nums = mirna_num_str.translate(_NT_TO_NUM)

    scorehash: Dict[float, List[dict]] = {}
    for offset in range(maxscore):
        score, string = _scan_alignment(utr_nums, mirna_nums, offset, site_type, on_top=True)
        scorehash.setdefault(score, []).append({"offset": offset, "gaploc": "top", "matchstring": string})

        score2, string2 = _scan_alignment(utr_nums, mirna_nums, offset, site_type, on_top=False)
        scorehash.setdefault(score2, []).append({"offset": offset, "gaploc": "bottom", "matchstring": string2})

    for score in sorted(scorehash.keys(), reverse=True):
        entries = scorehash[score]
        if len(entries) == 2 and entries[0]["offset"] == 0:
            i_ret = 1
        elif len(entries) > 1:
            i_ret = 0
            offset_ret = None
            for i, entry in enumerate(entries):
                if offset_ret is not None:
                    if entry["offset"] < offset_ret:
                        i_ret = i
                        offset_ret = entry["offset"]
                    elif entry["offset"] == offset_ret:
                        if entries[i_ret]["gaploc"] == "bottom" and entry["gaploc"] == "bottom":
                            raise RuntimeError("ERROR Two tied scores with same offset, and gaplocation")
                        if entry["gaploc"] == "bottom":
                            i_ret = i
                            offset_ret = entry["offset"]
                else:
                    i_ret = i
                    offset_ret = entry["offset"]
        else:
            i_ret = 0

        i = i_ret
        entry = entries[i]

        utrpre = utr_rev[: seedinfo["utrstart"][site_type]].translate(_NT_TO_NUM)
        mirnapre = mirna_rev[: seedinfo["mirnastart"][site_type]].translate(_NT_TO_NUM)
        matchpre = " " * (seedinfo["overhang"][site_type] + 1)
        for j in range(1, seedinfo["seedspan"][site_type] + 1):
            u_idx = j + seedinfo["overhang"][site_type]
            if u_idx < len(utrpre) and j < len(mirnapre) and _match_score(utrpre[u_idx], mirnapre[j]):
                matchpre += "|"
            else:
                matchpre += "E"
        matchpre += " " * entry["offset"]
        string = matchpre + entry["matchstring"]

        mirnapreoff = " " * seedinfo["overhang"][site_type]

        utrout = utr_rev[: seedinfo["utrstart"][site_type]] + utr_rev[seedinfo["utrstart"][site_type] :]
        mirnaout = mirnapreoff + mirna_rev[: seedinfo["mirnastart"][site_type]] + mirna_rev[seedinfo["mirnastart"][site_type] :]

        offsetstring = "-" * entry["offset"]
        if entry["gaploc"] == "top":
            utrout = utr_rev[: seedinfo["utrstart"][site_type]] + offsetstring + utr_rev[seedinfo["utrstart"][site_type] :]
        elif entry["gaploc"] == "bottom":
            mirnaout = (
                mirnapreoff
                + mirna_rev[: seedinfo["mirnastart"][site_type]]
                + offsetstring
                + mirna_rev[seedinfo["mirnastart"][site_type] :]
            )

        longest = max(len(utrout), len(mirnaout))
        utrout += " " * (longest - len(utrout))
        string += " " * (longest - len(string))
        mirnaout += " " * (longest - len(mirnaout))

        utrout = utrout[::-1]
        string = string[::-1]
        mirnaout = mirnaout[::-1]

        if score < 3:
            utrout, string, mirnaout = get_consequential_pairing(utrout, string, mirnaout)

        contribution = get_agarwal_contribution(agarwal, site_type, "3P_score", score)

        return contribution, utrout, string, mirnaout, score

    return "0", "0", "0", "0", 0.0


_REPLACE_LENGTH_CACHE = get_replace_length(20)


def get_consequential_pairing(utr_seq: str, pairing: str, mirna_seq: str) -> Tuple[str, str, str]:
    pairing = pairing.replace("|", "x")

    m = re.search(r"(\S+\s+\S+)\s+(\S+)", pairing) or re.search(r"(\S+)\s+(\S+)", pairing)
    if m:
        pairing_to_remove = m.group(1)
        length = len(pairing_to_remove)
        if length in _REPLACE_LENGTH_CACHE:
            new_substring = _REPLACE_LENGTH_CACHE[length]
            pairing = pairing.replace(pairing_to_remove, new_substring, 1)
            utr_seq = remove_gap_add_leading_space(utr_seq)
            mirna_seq = remove_gap_add_leading_space(mirna_seq)

    pairing = pairing.replace("x", "|")
    return utr_seq, pairing, mirna_seq


def remove_gap_add_leading_space(seq: str) -> str:
    pre_len = len(seq)
    seq = seq.replace("-", "").replace(".", "")
    diff = pre_len - len(seq)
    if diff > 0:
        seq = " " * diff + seq
    return seq


# --------------------------------------------------------------------------
# Main driver
# --------------------------------------------------------------------------


@dataclass
class ContextScoreEngine:
    agarwal: AgarwalParams
    utr_seq: Dict[str, Dict[str, str]]
    utr_length_scaling_factor: Dict[str, Dict[str, float]]
    mature_seq: Dict[Tuple[str, str], List[Tuple[str, str]]]
    family_to_seed: Dict[str, str]
    utr_scaled_end_to_air: Dict[str, Dict[str, Dict[int, float]]]
    utr_scaled_end_to_start: Dict[str, Dict[str, Dict[int, int]]]
    orf2length: Dict[str, Dict[str, int]]
    orf_to_8mer_counts: Dict[str, Dict[str, Dict[str, int]]]
    garcia_sps: Dict[int, Dict[str, float]]
    garcia_ta: Dict[str, float]
    rnaplfold_dir: str
    use_species: set
    missing_rnaplfold_files: set = field(default_factory=set)
    offset6mer_cache: Dict[Tuple[str, str, str], List[int]] = field(default_factory=dict)

    def score_one_target_line(self, fields: List[str]) -> List[List[str]]:
        transcript_id = fields[0]
        mirna_family_id = fields[1]
        species_id = fields[2]
        if species_id not in self.use_species:
            return []

        utr_start = int(fields[5])
        utr_end = int(fields[6])
        group_num = fields[7]
        pct = fields[12]
        site_type = SITE_TYPE_NUM.get(fields[8])
        if site_type is None:
            return []

        family_species = (mirna_family_id, species_id)
        if family_species not in self.mature_seq:
            return []

        air = get_air_this_site(transcript_id, species_id, utr_end, self.utr_scaled_end_to_air, self.utr_scaled_end_to_start)
        subseq_for_alignment = extract_subseq_for_alignment(self.utr_seq, transcript_id, species_id, utr_start, utr_end, site_type)
        local_au = get_local_au_contribution(self.agarwal, self.utr_seq, transcript_id, species_id, site_type, utr_start, utr_end)
        len3utr = get_len3utr_weighted_contribution(self.agarwal, self.utr_scaled_end_to_air, transcript_id, species_id, site_type, utr_end)
        min_dist = get_min_dist_weighted_contribution(self.agarwal, self.utr_scaled_end_to_air, transcript_id, species_id, site_type, utr_start, utr_end)
        sa = get_sa_contribution(self.agarwal, self.rnaplfold_dir, transcript_id, species_id, utr_start, site_type, self.missing_rnaplfold_files)
        orf_len = get_orflength_contribution(self.agarwal, self.orf2length, transcript_id, species_id, site_type)
        orf8mer = get_orf8mer_contribution(self.agarwal, self.orf_to_8mer_counts, self.family_to_seed, transcript_id, mirna_family_id, species_id, site_type)
        offset6mer = get_offset6mer_weighted_contribution(
            self.agarwal, self.utr_seq, self.family_to_seed, self.utr_scaled_end_to_air, self.offset6mer_cache,
            transcript_id, mirna_family_id, species_id, site_type, utr_end,
        )

        if pct == "NA":
            pct_contribution = "0"
        else:
            pct_contribution = get_pct_contribution(self.agarwal, site_type, pct)

        rows = []
        for mature_id, mature_mirna in self.mature_seq[family_species]:
            final_subseq, mirna_for_alignment = modify_subseq_for_alignment(mature_mirna, subseq_for_alignment, site_type)

            s1a, s1c, s1g, s8a, s8c, s8g = get_sRNA1_8_contributions(self.agarwal, mature_mirna, site_type)

            if site_type in (1, 4):
                site8a, site8c, site8g = get_site8_contribution(self.agarwal, subseq_for_alignment, site_type)
            else:
                site8a = site8c = site8g = "0"

            three_prime_contribution, aligned_utr, bars, aligned_mirna, raw_score = get3_prime_pairing_contribution(
                self.agarwal, site_type, final_subseq, mirna_for_alignment
            )

            seed_region = mature_mirna[1:8]

            if utr_start < MIN_DIST_TO_CDS:
                three_prime_contribution = local_au_eff = total_context_score = ta_contribution = sps_contribution = TOO_CLOSE_TO_CDS
            else:
                local_au_eff = local_au
                if seed_region in self.garcia_ta and seed_region in self.garcia_sps.get(site_type, {}):
                    ta_contribution = get_agarwal_contribution(self.agarwal, site_type, "TA_3UTR", self.garcia_ta[seed_region])
                    sps_contribution = get_agarwal_contribution(self.agarwal, site_type, "SPS", self.garcia_sps[site_type][seed_region])
                else:
                    ta_contribution = "0"
                    sps_contribution = "0"

                total = (
                    self.agarwal.intercept.get(site_type, 0.0)
                    + float(three_prime_contribution)
                    + float(local_au)
                    + float(min_dist)
                    + float(s1a) + float(s1c) + float(s1g)
                    + float(s8a) + float(s8c) + float(s8g)
                    + float(site8a) + float(site8c) + float(site8g)
                    + float(len3utr) + float(sa) + float(orf_len) + float(orf8mer) + float(offset6mer)
                    + float(ta_contribution) + float(sps_contribution) + float(pct_contribution)
                )
                total_context_score = f"{total:.{DIGITS_AFTER_DECIMAL}f}"
                if float(total_context_score) > MAX_CONTEXT_SCORE[site_type]:
                    total_context_score = str(MAX_CONTEXT_SCORE[site_type])

            weighted_total = get_weighted_context_score(total_context_score, air)

            row = [
                transcript_id,
                species_id,
                mature_id,
                SITE_TYPE_NAME[site_type],
                str(utr_start),
                str(utr_end),
                f"{self.agarwal.intercept.get(site_type, 0.0):.{DIGITS_AFTER_DECIMAL}f}",
                str(three_prime_contribution),
                str(local_au_eff),
                str(min_dist),
                str(s1a), str(s1c), str(s1g),
                str(s8a), str(s8c), str(s8g),
                str(site8a), str(site8c), str(site8g),
                str(len3utr), str(sa), str(orf_len), str(orf8mer), str(offset6mer),
                str(ta_contribution), str(sps_contribution), str(pct_contribution),
                str(total_context_score),
                "NA",
                _format_number(air),
                str(weighted_total),
                "NA",
                aligned_utr, bars, aligned_mirna,
                mirna_family_id,
                group_num,
            ]
            rows.append(row)
        return rows


CONTEXT_SCORE_FIELD = 27
PERCENTILE_RANK_FIELD = 28
WEIGHTED_CONTEXT_SCORE_FIELD = 30
WEIGHTED_PERCENTILE_RANK_FIELD = 31

_NUMERIC_RE = re.compile(r"^-?(?:\d+(?:\.\d*)?|\.\d+)$")


def _perl_truthy_string(s: str) -> bool:
    """Perl string truthiness: false only for "" and the literal "0"."""
    return s not in ("", "0")


def get_percentile_ranks(rows: List[List[str]]) -> None:
    """In-place: replace percentile-rank placeholder fields with real values.

    Preserves two Perl quirks, verified against the bundled sample output:

    1. The context-score and weighted-context-score percentiles are
       computed inside the *same* per-index loop sharing a single
       ``$numLower`` variable, so the weighted-score percentile can
       inherit state left over from the (separately sorted) context-score
       list rather than being computed fully independently.
    2. The "is the previous value truthy" gate operates on Perl string
       truthiness (false only for ``""``/``"0"``) of the *formatted*
       score string, not on whether the number is zero -- so a formatted
       value like ``"0.000"`` or ``"0.0000"`` is truthy, while a bare
       ``"0"`` (as produced when a score is capped to an integer 0) is
       not. We therefore keep working with the original formatted
       strings rather than converting to float up front.
    """
    by_mirna_scores: Dict[str, List[str]] = {}
    by_mirna_weighted: Dict[str, List[str]] = {}

    for row in rows:
        mirna = row[2]
        score = row[CONTEXT_SCORE_FIELD]
        weighted = row[WEIGHTED_CONTEXT_SCORE_FIELD]
        if _NUMERIC_RE.match(score):
            by_mirna_scores.setdefault(mirna, []).append(score)
        if _NUMERIC_RE.match(weighted):
            by_mirna_weighted.setdefault(mirna, []).append(weighted)

    score_to_pctile: Dict[str, Dict[str, int]] = {}
    weighted_to_pctile: Dict[str, Dict[str, int]] = {}

    for mirna in set(by_mirna_scores) | set(by_mirna_weighted):
        score_list = sorted(by_mirna_scores.get(mirna, []), key=float, reverse=True)
        weighted_list = sorted(by_mirna_weighted.get(mirna, []), key=float, reverse=True)
        n = len(score_list)
        if n == 0:
            continue

        num_lower = 0
        s_table: Dict[str, int] = {}
        w_table: Dict[str, int] = {}
        for i in range(n):
            prev_score = score_list[i - 1]
            if _perl_truthy_string(prev_score) and float(score_list[i]) != float(prev_score):
                num_lower = i
            s_table[score_list[i]] = math.floor(100 * num_lower / n)

            if i < len(weighted_list):
                prev_weighted = weighted_list[i - 1]
                if _perl_truthy_string(prev_weighted) and float(weighted_list[i]) != float(prev_weighted):
                    num_lower = i
                w_table[weighted_list[i]] = math.floor(100 * num_lower / n)

        score_to_pctile[mirna] = s_table
        weighted_to_pctile[mirna] = w_table

    for row in rows:
        mirna = row[2]
        score = row[CONTEXT_SCORE_FIELD]
        weighted = row[WEIGHTED_CONTEXT_SCORE_FIELD]
        if _NUMERIC_RE.match(score):
            row[PERCENTILE_RANK_FIELD] = str(score_to_pctile[mirna][score])
        if _NUMERIC_RE.match(weighted):
            row[WEIGHTED_PERCENTILE_RANK_FIELD] = str(weighted_to_pctile[mirna][weighted])


HEADER = (
    "Gene ID\tSpecies ID\tMirbase ID\tSite Type\tUTR start\tUTR end\t"
    "Site type contribution\t3' pairing contribution\tlocal AU contribution\tMin_dist contribution\t"
    "sRNA1A contribution\tsRNA1C contribution\tsRNA1G contribution\tsRNA8A contribution\tsRNA8C contribution\tsRNA8G contribution\t"
    "site8A contribution\tsite8C contribution\tsite8G contribution\t3'UTR length contribution\tSA contribution\t"
    "ORF length contribution\tORF 8mer contribution\tOffset 6mer contribution\tTA contribution\tSPS contribution\tPCT contribution\t"
    "context++ score\tcontext++ score percentile\tAIR\tweighted context++ score\tweighted context++ score percentile\t"
    "UTR region\tUTR-miRNA pairing\tmature miRNA sequence\tmiRNA family\tGroup #\n"
)


def run(
    mirna_file: str,
    utr_file: str,
    predicted_targets_file: str,
    orf_lengths_file: str,
    orf_8mers_file: str,
    out_path: str,
    data_dir: str,
    rnaplfold_dir: str,
) -> None:
    use_species = set(SPECIES)

    agarwal = AgarwalParams.read(os.path.join(data_dir, "Agarwal_2015_parameters.txt"))
    utr_seq, utr_length_scaling_factor, have_utrs_this_species = read_utrs(utr_file, use_species)
    mature_seq, family_to_seed = read_mirnas_context(mirna_file)
    utr_scaled_end_to_air, utr_scaled_end_to_start = read_isoform_ratios(
        os.path.join(data_dir, "All_cell_lines.AIRs.txt"), utr_length_scaling_factor, have_utrs_this_species
    )
    run_rnaplfold_all_utrs(utr_seq, rnaplfold_dir)
    orf2length, orf_to_8mer_counts = read_orf_data(orf_lengths_file, orf_8mers_file, use_species)
    garcia_sps, garcia_ta = read_ta_sps(os.path.join(data_dir, "TA_SPS_by_seed_region.txt"))

    engine = ContextScoreEngine(
        agarwal=agarwal,
        utr_seq=utr_seq,
        utr_length_scaling_factor=utr_length_scaling_factor,
        mature_seq=mature_seq,
        family_to_seed=family_to_seed,
        utr_scaled_end_to_air=utr_scaled_end_to_air,
        utr_scaled_end_to_start=utr_scaled_end_to_start,
        orf2length=orf2length,
        orf_to_8mer_counts=orf_to_8mer_counts,
        garcia_sps=garcia_sps,
        garcia_ta=garcia_ta,
        rnaplfold_dir=rnaplfold_dir,
        use_species=use_species,
    )

    all_rows: List[List[str]] = []
    with open(predicted_targets_file) as fh:
        next(fh)  # header
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            all_rows.extend(engine.score_one_target_line(fields))

    get_percentile_ranks(all_rows)

    with open(out_path, "w") as out:
        out.write(HEADER)
        for row in all_rows:
            out.write("\t".join(row) + "\n")

    if engine.missing_rnaplfold_files:
        print(
            f"\n! Note that we couldn't run RNAplfold (or find the output files) for "
            f"{len(engine.missing_rnaplfold_files)} sequences,\nso these will have no SA contributions!",
            file=sys.stderr,
        )
    print(f"\nAll done!  -  See {out_path} for output.\n", file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 6:
        print(
            "USAGE: targetscan_70_context_scores.py miRNA_file UTR_file "
            "PredictedTargetsBL_PCT_file ORF_lengths_file ORF_8mer_counts_file "
            "ContextScoresOutput_file [data_dir] [rnaplfold_dir]",
            file=sys.stderr,
        )
        return 0
    data_dir = argv[6] if len(argv) > 6 else "data"
    rnaplfold_dir = argv[7] if len(argv) > 7 else "RNAplfold_in_out"
    run(argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], data_dir, rnaplfold_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
