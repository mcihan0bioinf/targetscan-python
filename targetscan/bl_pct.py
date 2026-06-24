"""Calculate branch length (BL) and probability of conserved targeting (PCT).

Python port of ``tsh_orig/TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl``.

Note on a quirk preserved from the original: the Perl script's per-site
branch is chosen with::

    if ($siteType eq $groupType || $siteType ne "8mer") { ... simple case ... }
    elsif ($siteType eq "8mer-1a" && ...) { ... combination case ... }

Since the predicted-target site type is always named ``"8mer-1a"`` (never
the literal string ``"8mer"``), ``$siteType ne "8mer"`` is always true, so
the first ("simple") branch is *always* taken and the combination branch is
dead code. This was verified against the bundled sample BL/PCT output file
(an "8mer-1a" site in a "7mer-m8+8mer-1a" group is scored using the simple
8mer-1a branch length, not the combination logic). We replicate that
observed behavior -- i.e. only the simple branch is implemented -- so output
matches the original tool exactly.
"""

from __future__ import annotations

import itertools
import math
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from .phylo import Tree

E_CONSTANT = math.e

SITE_TYPE_NAMES = ("8mer-1a", "7mer-m8", "7mer-1a", "6mer")

# Conservation thresholds (BL) for TargetScanHuman, TS7 trees.
SITE_TYPE_TO_CONS_THRESHOLD = {
    "8mer-1a": 1.8,
    "7mer-m8": 2.8,
    "7mer-1a": 3.6,
    "6mer": 100,  # 6mers should never be scored as "conserved"
}

NULL_VALUE = "NA"

_NUMBER_RE = re.compile(r"^-?(?:\d+(?:\.\d*)?|\.\d+)$")


def _looks_like_number(value: object) -> bool:
    return isinstance(value, (int, float)) or bool(_NUMBER_RE.match(str(value)))


def read_mirnas(path: str) -> Tuple[Dict[str, str], set]:
    """Return (family_id -> seed region, set of families/seeds to score)."""
    mir_id_to_seed: Dict[str, str] = {}
    get_bls_this_family: set = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            fam_id, seed = fields[0], fields[1]
            seed = re.sub(r"[Tt]", "U", seed)
            mir_id_to_seed[fam_id] = seed
            get_bls_this_family.add(fam_id)
            get_bls_this_family.add(seed)
    return mir_id_to_seed, get_bls_this_family


def read_bl_bins(path: str) -> Dict[str, str]:
    refseq_to_bin: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\r\n").split("\t")
            refseq_to_bin[fields[0]] = fields[2]
    return refseq_to_bin


@dataclass
class PCTData:
    family_type_to_coeff: Dict[str, Dict[str, List[float]]] = field(default_factory=dict)
    has_pct: set = field(default_factory=set)

    @classmethod
    def read(cls, pct_dir: str) -> "PCTData":
        obj = cls()
        files = {
            "8mer-1a": f"{pct_dir}/8mer_PCT_parameters.txt",
            "7mer-m8": f"{pct_dir}/7mer_m8_PCT_parameters.txt",
            "7mer-1a": f"{pct_dir}/7mer_1a_PCT_parameters.txt",
        }
        for site_type, path in files.items():
            with open(path) as fh:
                for i, line in enumerate(fh):
                    if i == 0:
                        continue
                    fields = line.rstrip("\r\n").split("\t")
                    family = fields[0]
                    obj.has_pct.add(family)
                    coeffs = [float(x) for x in fields[1:5]]
                    obj.family_type_to_coeff.setdefault(family, {})[site_type] = coeffs
        return obj


def calculate_pct_this_bl(pct_data: PCTData, family: str, site_type: str, bl: float) -> str:
    b0, b1, b2, b3 = pct_data.family_type_to_coeff[family][site_type]
    pct = b0 + (b1 / (1 + E_CONSTANT ** ((0 - b2) * bl + b3)))
    formatted = f"{pct:.4f}"
    # Match Perl: round to 4 decimals *first*, then check sign of the
    # rounded value (so e.g. raw -0.00001 -> "-0.0000", which is not < 0).
    if float(formatted) < 0:
        return "0.0"
    return formatted


def get_pct(pct_data: PCTData, site_type: str, family: str, branch_length) -> str:
    if not _looks_like_number(branch_length):
        return "EMPTY"
    bl = float(branch_length)
    if family in pct_data.has_pct and bl > 0:
        return "0.0" if site_type == "6mer" else calculate_pct_this_bl(pct_data, family, site_type, bl)
    if family in pct_data.has_pct and bl == 0:
        return "0.0"
    return NULL_VALUE


def select_conservation_for_this_site(bl, site_type: str) -> str:
    if site_type == "6mer":
        return ""
    threshold = SITE_TYPE_TO_CONS_THRESHOLD.get(site_type)
    if threshold is None or not _looks_like_number(bl):
        return ""
    return "x" if float(bl) >= threshold else ""


def get_group_type_another_way(site_types: Iterable[str]) -> str:
    return "+".join(sorted(set(site_types)))


class BranchLengthCache:
    """Loads bin-specific trees lazily and caches BL by (species list, bin)."""

    def __init__(self, pct_params_dir: str) -> None:
        self._pct_params_dir = pct_params_dir
        self._trees: Dict[str, Tree] = {}
        self._cache: Dict[Tuple[str, str], float] = {}

    def _tree_for_bin(self, bin_num: str) -> Tree:
        bin_padded = bin_num.zfill(2)
        if bin_padded not in self._trees:
            path = f"{self._pct_params_dir}/Tree.bin_{bin_padded}.txt"
            self._trees[bin_padded] = Tree(path)
        return self._trees[bin_padded]

    def look_at_branch_length(
        self, refseq: str, species_list_str: Optional[str], refseq_to_bin: Dict[str, str]
    ) -> object:
        if not species_list_str:
            return NULL_VALUE

        refseq_bin = refseq_to_bin.get(refseq)
        if refseq_bin is None:
            print(f"ERROR: No bin assignment for {refseq} ; assigning to bin 1", file=sys.stderr)
            refseq_bin = "1"

        cache_key = (species_list_str, refseq_bin)
        if cache_key in self._cache:
            return self._cache[cache_key]

        species = species_list_str.split(" ")
        if len(species) == 1:
            bl: object = 0
        else:
            tree = self._tree_for_bin(refseq_bin)
            bl = tree.branch_length(species[0], species)

        self._cache[cache_key] = bl
        return bl


def _format_bl(bl: object) -> str:
    """Match the Perl getBranchLength() formatting: bare "0" for a
    single-species (no-tree-lookup) site, sprintf("%.4f", ...) otherwise."""
    if bl == NULL_VALUE:
        return NULL_VALUE
    if isinstance(bl, int):
        return str(bl)
    return f"{bl:.4f}"


def process_group(
    lines: List[List[str]],
    mir_id_to_seed: Dict[str, str],
    pct_data: PCTData,
    bl_cache: BranchLengthCache,
    refseq_to_bin: Dict[str, str],
    out,
) -> None:
    site_types_this_group = [f[8] for f in lines]
    rechecked_group_type = get_group_type_another_way(site_types_this_group)

    # Map siteType -> species list (using the most specific list available),
    # mirroring the read loop in the Perl script.
    site_type_to_species_list: Dict[str, str] = {}
    for f in lines:
        site_type = f[8]
        group_species_list = f[11]
        group_site_type_species_list = f[12]
        site_type_to_species_list[site_type] = group_site_type_species_list or group_species_list

    for f in lines:
        f = list(f)
        refseq = f[0]
        fam_id = f[1]
        site_type = f[8]
        group_type = f[10]
        is_masked = f[13]

        mirna_family = mir_id_to_seed.get(fam_id)
        if mirna_family is None:
            print(f"No seed for {fam_id}", file=sys.stderr)
            mirna_family = "UNKNOWN"

        if group_type != rechecked_group_type:
            group_type = rechecked_group_type
            f[10] = rechecked_group_type

        # See module docstring: the "combination" branch in the original
        # Perl is unreachable, so only the simple per-site-type branch
        # is implemented here.
        species_list = site_type_to_species_list.get(site_type)
        if species_list is None:
            print(f"No list for {f[7]} site type of {f}", file=sys.stderr)
            bl_selected: object = NULL_VALUE
        else:
            bl_selected = bl_cache.look_at_branch_length(refseq, species_list, refseq_to_bin)

        group_type_for_this_site = site_type
        cons_score = select_conservation_for_this_site(bl_selected, group_type_for_this_site)

        if mirna_family in pct_data.has_pct and is_masked != "1":
            pct = get_pct(pct_data, group_type_for_this_site, mirna_family, bl_selected)
        else:
            pct = NULL_VALUE

        f[10] = group_type_for_this_site
        f[11] = _format_bl(bl_selected)
        f[12] = pct if pct else ""
        f[13] = cons_score

        out.write("\t".join(f) + "\n")


def run(
    mirna_file: str,
    predicted_targets_file: str,
    gene_bin_file: str,
    pct_params_dir: str,
    out_path: str,
) -> None:
    mir_id_to_seed, get_bls_this_family = read_mirnas(mirna_file)
    refseq_to_bin = read_bl_bins(gene_bin_file)
    pct_data = PCTData.read(pct_params_dir)
    bl_cache = BranchLengthCache(pct_params_dir)

    with open(predicted_targets_file) as fh, open(out_path, "w") as out:
        out.write(
            "Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\t"
            "Group_num\tSite_type\tmiRNA in this species\tGroup_type\tBranch length score\t"
            "Pct\tConserved\n"
        )
        next(fh)  # skip header

        data_lines = [line.rstrip("\r\n").split("\t") for line in fh if line.strip()]
        data_lines.sort(key=lambda f: int(f[7]))

        for group_num, group_iter in itertools.groupby(data_lines, key=lambda f: f[7]):
            group_lines = list(group_iter)
            qualifying = [f for f in group_lines if f[1] in get_bls_this_family]
            non_qualifying = [f for f in group_lines if f[1] not in get_bls_this_family]

            for f in non_qualifying:
                f = list(f)
                f[11] = "NA"
                f[12] = "NA"
                out.write("\t".join(f) + "\n")

            if qualifying:
                process_group(qualifying, mir_id_to_seed, pct_data, bl_cache, refseq_to_bin, out)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 4:
        print(
            "USAGE: targetscan_70_BL_PCT.py miRNA_file predicted_targets "
            "UTR_bin_info pct_params_dir [out_file]",
            file=sys.stderr,
        )
        return 0
    out_path = argv[4] if len(argv) > 4 else "/dev/stdout"
    run(argv[0], argv[1], argv[2], argv[3], out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
