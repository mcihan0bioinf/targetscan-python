"""Predict miRNA target sites in aligned 3' UTRs.

Python port of ``tsh_orig/targetscan_70.pl`` (TargetScan miRNA target site
prediction, "TargetScanS" algorithm). Behavior, including the exact set of
seed-match site types searched for and the grouping/merging rules, is kept
identical to the original Perl so output is line-for-line compatible.

Input/output file formats are unchanged from the original tool, so this
module works on TargetScan vert70 (hg19) data as well as the newer vert80
(hg38) data downloaded from
https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, TextIO, Tuple

# Find sites in species even when the miRNA hasn't been annotated there.
FIND_SITES_ALL_SPECIES = False
# Minimum overlap (nt) required to group sites found in different species.
REQUIRED_OVERLAP = 2
# Ribosome-shadow length at the start of the UTR to mask out.
BEG_UTR_MASK_LENGTH = 14

# Which site types to search for (matches $GET_MATCH in the Perl script).
GET_MATCH = {1: True, 2: True, 3: True, 6: True, 4: False, 5: False}

SITE_ID_TO_TYPE = {
    1: "7mer-1a",
    2: "7mer-m8",
    3: "8mer-1a",
    4: "8mer-1u",
    5: "6mer-1a",
    6: "6mer",
}
SITE_ID_TO_LENGTH = {1: 7, 2: 7, 3: 8, 4: 8, 5: 6, 6: 6}

_COMPLEMENT = str.maketrans("AUCG", "UAGC")


def make_seed_match_regex(seed_match: str) -> "re.Pattern[str]":
    """Turn a seed-match sequence into a regex tolerant of alignment gaps."""
    pieces = []
    n = len(seed_match)
    for i, nt in enumerate(seed_match):
        pieces.append(nt if i == n - 1 else f"{nt}-{{0,}}")
    return re.compile("".join(pieces))


def get_seeds(seed_region: str) -> Dict[int, "re.Pattern[str]"]:
    """Build the regexes for all site types from a 7nt miRNA seed region."""
    rseed2 = seed_region[::-1].translate(_COMPLEMENT)  # 7mer-m8
    rseed6 = rseed2[1:]  # 6mer
    rseed5 = (rseed6 + "A")[1:]  # 6mer-1a
    rseed1 = rseed6 + "A"  # 7mer-1a
    rseed3 = rseed2 + "A"  # 8mer-1a
    rseed4 = rseed2 + "U"  # 8mer-1u

    return {
        1: make_seed_match_regex(rseed1),
        2: make_seed_match_regex(rseed2),
        3: make_seed_match_regex(rseed3),
        4: make_seed_match_regex(rseed4),
        5: make_seed_match_regex(rseed5),
        6: make_seed_match_regex(rseed6),
    }


@dataclass
class MiRNAFamilies:
    seed: Dict[str, str] = field(default_factory=dict)
    species: Dict[Tuple[str, str], bool] = field(default_factory=dict)
    patterns: Dict[str, Dict[int, "re.Pattern[str]"]] = field(default_factory=dict)

    @classmethod
    def read(cls, path: str) -> "MiRNAFamilies":
        obj = cls()
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                fam_id, seed, species_list = line.split("\t")[:3]
                seed = seed.upper().replace("T", "U")
                obj.seed[fam_id] = seed
                for sp in species_list.split(";"):
                    if sp:
                        obj.species[(fam_id, sp)] = True
        for fam_id, seed in obj.seed.items():
            obj.patterns[fam_id] = get_seeds(seed)
        return obj


def mask_5prime(seq: str, mask_length: int = BEG_UTR_MASK_LENGTH) -> str:
    """Mask the first ``mask_length`` nucleotides (gaps don't count)."""
    out = list(seq)
    masked = 0
    for i, ch in enumerate(out):
        if masked >= mask_length:
            break
        if ch.upper() in "ACGTU":
            out[i] = "N"
            masked += 1
    return "".join(out)


def get_utr_coords(alignment: str, msa_start: int, site_type: int) -> Tuple[int, int]:
    """Convert a 1-based MSA *start* coordinate to 1-based UTR (gapless) coords."""
    utr_beg = alignment[:msa_start].replace("-", "")
    start = len(utr_beg)
    end = start + SITE_ID_TO_LENGTH[site_type] - 1
    return start, end


def _find_matches(alignment_upper: str, pattern: "re.Pattern[str]") -> Iterator[Tuple[int, int, str]]:
    """Yield (1-based start, 1-based end, matched substring) for all,
    possibly overlapping, matches -- mirroring the Perl backtracking loop."""
    pos = 0
    n = len(alignment_upper)
    while pos <= n:
        m = pattern.search(alignment_upper, pos)
        if not m:
            return
        start1 = m.start() + 1
        end1 = m.end()
        yield start1, end1, m.group(0)
        pos = m.start() + 1


def _drop_subset_coords(gapped_match: str) -> Tuple[int, int, int]:
    """Offsets of the 2nd, 3rd, and second-to-last non-gap characters."""
    m1 = re.match(r"^[^-]-*([^-])", gapped_match)
    start_plus_one = m1.end() - 1 if m1 else 0
    m2 = re.match(r"^[^-]-*[^-]-*([^-])", gapped_match)
    start_plus_two = m2.end() - 1 if m2 else 0
    m3 = re.search(r"([^-])-*[^-]$", gapped_match)
    end_minus_one = (len(m3.group(0)) - 1) if m3 else 0
    return start_plus_one, start_plus_two, end_minus_one


class _GeneMirState:
    """Per-(gene, miRNA family) working state, mirroring Perl globals."""

    def __init__(self) -> None:
        self.species_start_end: Dict[Tuple[str, int, int], int] = {}
        self.species_start_end_match: Dict[Tuple[str, int, int], str] = {}
        self.species_start_end_removed: Dict[Tuple[str, int, int], int] = {}
        self.species_start_end_match_removed: Dict[Tuple[str, int, int], str] = {}
        self.species_start_end_masked: Dict[Tuple[str, int, int], int] = {}

    def get_matches(self, species_id: str, alignment: str, pattern: "re.Pattern[str]", match_type: int) -> None:
        upper = alignment.upper()
        for start, end, matched in _find_matches(upper, pattern):
            key = (species_id, start, end)
            self.species_start_end[key] = match_type
            self.species_start_end_match[key] = matched
            length = end - start + 1
            original_region = alignment[start - 1 : start - 1 + length]
            uc_region = upper[start - 1 : start - 1 + length]
            if uc_region != original_region:
                self.species_start_end_masked[key] = 1

    def drop_site(self, species: str, start: int, end: int) -> bool:
        key = (species, start, end)
        if key in self.species_start_end:
            self.species_start_end_removed[key] = self.species_start_end[key]
            self.species_start_end_match_removed[key] = self.species_start_end_match[key]
            del self.species_start_end[key]
            del self.species_start_end_match[key]
            return True
        if key in self.species_start_end_removed:
            return True
        return False

    def find_remove_match_subsets(self) -> None:
        for key in list(self.species_start_end.keys()):
            species, start, end = key
            gapped_match = self.species_start_end_match.get(key)
            if gapped_match and "-" not in gapped_match:
                start_plus_one = start + 1
                start_plus_two = start + 2
                end_minus_one = end - 1
            elif gapped_match:
                off1, off2, off3 = _drop_subset_coords(gapped_match)
                start_plus_one = start + off1
                start_plus_two = start + off2
                end_minus_one = end - off3
            else:
                continue

            site_type = self.species_start_end.get(key)
            if site_type == 1:  # 7mer-1a
                if GET_MATCH[6]:
                    self.drop_site(species, start, end_minus_one)
                if GET_MATCH[5]:
                    self.drop_site(species, start_plus_one, end)
            elif site_type == 2:  # 7mer-m8
                if GET_MATCH[6]:
                    self.drop_site(species, start_plus_one, end)
            elif site_type == 3:  # 8mer-1a
                if GET_MATCH[2]:
                    self.drop_site(species, start, end_minus_one)
                if GET_MATCH[1]:
                    self.drop_site(species, start_plus_one, end)
                if GET_MATCH[5]:
                    self.drop_site(species, start_plus_two, end)
                if GET_MATCH[6]:
                    self.drop_site(species, start_plus_one, end_minus_one)
            elif site_type == 4:  # 8mer-1u
                if GET_MATCH[2]:
                    self.drop_site(species, start, end_minus_one)
                if GET_MATCH[6]:
                    self.drop_site(species, start_plus_one, end_minus_one)


def _site_sort_key(site: Tuple[str, int, int]) -> str:
    """Replicate Perl's default (lexicographic string) sort of the
    ``"species::start::end"`` hash key, so group numbers come out identical."""
    return f"{site[0]}::{site[1]}::{site[2]}"


def _sites_overlap(s1: int, e1: int, s2: int, e2: int) -> Tuple[bool, int]:
    """Return (overlap?, num_overlap_nt) following the Perl branch logic."""
    if s1 == s2 and e1 == e2:
        return True, e1 - s1 + 1
    if s1 == s2:
        return True, max(e1, e2) - s1 + 1
    if e1 == e2:
        return True, e1 - min(s1, s2) + 1
    if s1 > s2 and s1 <= e2:
        n = e2 - s1 + 1
        return n >= REQUIRED_OVERLAP, n
    if e1 >= s2 and e1 < e2:
        n = e1 - s2 + 1
        return n >= REQUIRED_OVERLAP, n
    if (s1 > s2 and e1 < e2) or (s2 > s1 and e2 < e1):
        return True, 0
    return False, 0


def _group_sites(state: _GeneMirState) -> Tuple[Dict[Tuple[str, int, int], int], int]:
    """Group overlapping cross-species sites; returns (site->group, last group#)."""
    sites = sorted(state.species_start_end.keys(), key=_site_sort_key)
    site_to_group: Dict[Tuple[str, int, int], int] = {}
    group_num = 0

    for site1 in sites:
        sp1, s1, e1 = site1
        for site2 in sites:
            sp2, s2, e2 = site2
            if sp1 == sp2:
                continue
            overlap, _ = _sites_overlap(s1, e1, s2, e2)
            if overlap:
                g1 = site_to_group.get(site1)
                g2 = site_to_group.get(site2)
                if not g1 or not g2:
                    if g1:
                        site_to_group[site2] = g1
                    elif g2:
                        site_to_group[site1] = g2
                    else:
                        group_num += 1
                        site_to_group[site1] = group_num
                        site_to_group[site2] = group_num

    for site in sites:
        if site not in site_to_group:
            group_num += 1
            site_to_group[site] = group_num

    return site_to_group, group_num


def _make_nonredundant_sorted(values: Iterable[str]) -> List[str]:
    seen = sorted(set(values), key=lambda v: float(v))
    return seen


class TargetScanPredictor:
    """Runs miRNA target-site prediction for a UTR alignment file."""

    def __init__(self, mirnas: MiRNAFamilies) -> None:
        self.mirnas = mirnas
        self.group_num = 0

    def predict_gene(self, gene_id: str, species_to_utr: Dict[str, str], out: TextIO) -> None:
        for fam_id in sorted(self.mirnas.seed.keys()):
            state = _GeneMirState()

            for species_id, alignment in sorted(species_to_utr.items()):
                if FIND_SITES_ALL_SPECIES or self.mirnas.species.get((fam_id, species_id)):
                    for match_type, enabled in GET_MATCH.items():
                        if enabled:
                            pattern = self.mirnas.patterns[fam_id][match_type]
                            state.get_matches(species_id, alignment, pattern, match_type)
                    state.find_remove_match_subsets()

            if not state.species_start_end:
                continue

            self._group_and_print(gene_id, fam_id, species_to_utr, state, out)

    def _group_and_print(
        self,
        gene_id: str,
        fam_id: str,
        species_to_utr: Dict[str, str],
        state: _GeneMirState,
        out: TextIO,
    ) -> None:
        sites = sorted(state.species_start_end.keys(), key=_site_sort_key)

        # Pairwise grouping (offset group numbers so they're unique across the file)
        local_site_to_group, max_local_group = _group_sites(state)
        site_to_group = {
            site: self.group_num + g for site, g in local_site_to_group.items()
        }
        self.group_num += max_local_group

        group_to_site_types: Dict[int, List[int]] = {}
        group_to_species: Dict[int, List[str]] = {}
        group_type_species: Dict[Tuple[int, int], List[str]] = {}

        rows = []  # (group_num, site_type_id, line_fields...)

        for site in sites:
            species_this_site, msa_start, msa_end = site
            g = site_to_group[site]
            site_type = state.species_start_end[site]

            group_to_site_types.setdefault(g, []).append(site_type)
            group_to_species.setdefault(g, []).append(species_this_site)
            group_type_species.setdefault((g, site_type), []).append(species_this_site)

            # If a wider site is present, its subset site types are implied present too.
            if site_type == 1:
                if GET_MATCH[6]:
                    group_type_species.setdefault((g, 6), []).append(species_this_site)
                if GET_MATCH[5]:
                    group_type_species.setdefault((g, 5), []).append(species_this_site)
            elif site_type == 2:
                if GET_MATCH[6]:
                    group_type_species.setdefault((g, 6), []).append(species_this_site)
            elif site_type == 3:
                if GET_MATCH[1]:
                    group_type_species.setdefault((g, 1), []).append(species_this_site)
                if GET_MATCH[2]:
                    group_type_species.setdefault((g, 2), []).append(species_this_site)
                if GET_MATCH[5]:
                    group_type_species.setdefault((g, 5), []).append(species_this_site)
                if GET_MATCH[6]:
                    group_type_species.setdefault((g, 6), []).append(species_this_site)
            elif site_type == 4:
                if GET_MATCH[2]:
                    group_type_species.setdefault((g, 2), []).append(species_this_site)
                if GET_MATCH[6]:
                    group_type_species.setdefault((g, 6), []).append(species_this_site)

            utr_start, utr_end = get_utr_coords(species_to_utr[species_this_site], msa_start, site_type)

            annotated = "x" if self.mirnas.species.get((fam_id, species_this_site)) else " "

            is_masked = 1 if site in state.species_start_end_masked else 0

            rows.append(
                {
                    "group": g,
                    "site_type": site_type,
                    "species": species_this_site,
                    "msa_start": msa_start,
                    "msa_end": msa_end,
                    "utr_start": utr_start,
                    "utr_end": utr_end,
                    "annotated": annotated,
                    "is_masked": is_masked,
                }
            )

        # Determine group type (unique, sorted-by-type-id site-type-name list) per group.
        group_type_name: Dict[int, str] = {}
        for g, types in group_to_site_types.items():
            unique_types = sorted(set(types))
            group_type_name[g] = "+".join(SITE_ID_TO_TYPE[t] for t in unique_types)

        rows.sort(key=lambda r: r["group"])

        for r in rows:
            g = r["group"]
            site_type = r["site_type"]
            group_type = group_type_name[g]
            species_in_group = _make_nonredundant_sorted(group_to_species[g])

            site_type_name = SITE_ID_TO_TYPE[site_type]

            species_with_this_type = ""
            if group_type != site_type_name:
                species_with_this_type = " ".join(
                    _make_nonredundant_sorted(group_type_species.get((g, site_type), []))
                )

            out.write(
                "\t".join(
                    str(x)
                    for x in (
                        gene_id,
                        fam_id,
                        r["species"],
                        r["msa_start"],
                        r["msa_end"],
                        r["utr_start"],
                        r["utr_end"],
                        g,
                        site_type_name,
                        r["annotated"],
                        group_type,
                        " ".join(species_in_group),
                        species_with_this_type,
                        r["is_masked"],
                    )
                )
                + "\n"
            )


def read_utr_blocks(path: str) -> Iterator[Tuple[str, Dict[str, str]]]:
    """Yield (gene_id, {species_id: masked_aligned_utr}) gene by gene."""
    last_gene = None
    species_to_utr: Dict[str, str] = {}

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            gene_id, species_id, utr = fields[0], fields[1], fields[2]
            if not species_id:
                continue
            utr = re.sub(r"[Tt]", "U", utr)
            utr = mask_5prime(utr)

            if gene_id and gene_id != last_gene:
                if last_gene is not None:
                    yield last_gene, species_to_utr
                species_to_utr = {}

            species_to_utr[species_id] = utr
            last_gene = gene_id

    if last_gene is not None:
        yield last_gene, species_to_utr


def run(mirna_file: str, utr_file: str, out_path: str, verbose: bool = True) -> None:
    mirnas = MiRNAFamilies.read(mirna_file)
    predictor = TargetScanPredictor(mirnas)

    with open(out_path, "w") as out:
        out.write(
            "a_Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\t"
            "Group_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\t"
            "Species_in_this_group_with_this_site_type\tORF_overlap\n"
        )
        for gene_id, species_to_utr in read_utr_blocks(utr_file):
            if verbose:
                print(f"Processing {gene_id}")
            predictor.predict_gene(gene_id, species_to_utr, out)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 3:
        print(
            "USAGE: targetscan_70.py miRNA_file UTR_file PredictedTargetsOutputFile",
            file=sys.stderr,
        )
        return 0
    run(argv[0], argv[1], argv[2])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
