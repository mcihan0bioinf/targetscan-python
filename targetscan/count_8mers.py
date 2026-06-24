"""Count 8mer sites (defined by miRNA seed regions) in a set of sequences.

Python port of ``tsh_orig/TargetScan7_context_scores/targetscan_count_8mers.pl``.
Used to generate the "ORF 8mer counts" and "ORF lengths" files that
``context_scores.py`` needs, by running this against the ORF_Sequences
file for the matching TargetScan release (vert70/hg19 or vert80/hg38):

    python -m targetscan.count_8mers miR_Family_info.txt ORF_Sequences.txt \\
        ORF_8mer_counts.txt

Note: despite the historical name (this also doubled as an offset-6mer-site
finder in some uses), this implementation -- like the original -- builds an
8mer-1a-equivalent site string (reverse complement of the seed, plus a
trailing "A") from each miRNA family's seed region, keyed by the seed region
itself (matching how ``context_scores.py``/``getORF8mer_contribution`` looks
counts up by seed region, not by family name).
"""

from __future__ import annotations

import re
import sys
from typing import Dict, List, Optional, Tuple

_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_COMPLEMENT)


def read_seed_sites(mirna_file: str) -> Dict[str, str]:
    """Return {seed_region: site_to_find} (site = revcomp(seed) + "A")."""
    seed_to_site: Dict[str, str] = {}
    with open(mirna_file) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            seed = fields[1]
            site = reverse_complement(seed)
            site = re.sub(r"[Tt]", "U", site)
            site = site + "A"
            seed_to_site[seed] = site
    return seed_to_site


def run(mirna_file: str, sequence_file: str, out_path: str) -> None:
    seed_to_site = read_seed_sites(mirna_file)
    sites_to_find = sorted(set(seed_to_site.values()))

    if sequence_file.endswith(".txt"):
        lengths_path = sequence_file[: -len(".txt")] + ".lengths.txt"
    else:
        lengths_path = sequence_file + ".lengths.txt"

    with open(sequence_file) as fh, open(out_path, "w") as out, open(lengths_path, "w") as lengths_out:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            sequence_id, species, sequence = line.split("\t")[:3]
            sequence = sequence.replace("-", "").replace(".", "")
            sequence = re.sub(r"[Tt]", "U", sequence)

            lengths_out.write(f"{sequence_id}\t{species}\t{len(sequence)}\n")

            upper_seq = sequence.upper()
            site_to_count: Dict[str, int] = {}
            for site in sites_to_find:
                site_to_count[site] = upper_seq.count(site.upper())

            for seed in sorted(seed_to_site.keys()):
                site = seed_to_site[seed]
                count = site_to_count.get(site, 0)
                if count:
                    out.write(f"{sequence_id}\t{species}\t{seed}\t{count}\n")

    print(f"\nAll done -- Also see {lengths_path} for sequence lengths.\n", file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        print(
            "USAGE: targetscan_count_8mers.py miRNA_seeds_file UTRs out_file",
            file=sys.stderr,
        )
        return 0
    out_path = argv[2] if len(argv) > 2 else "/dev/stdout"
    run(argv[0], argv[1], out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
