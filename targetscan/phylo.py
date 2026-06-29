"""Phylogenetic branch-length calculations shared by BL_bins and BL_PCT.

Ports the "Robin Friedman" branch-length algorithm (shared, near-identical
code in ``targetscan_70_BL_bins.pl`` and ``targetscan_70_BL_PCT.pl``) from
BioPerl's ``Bio::TreeIO`` to Biopython's ``Bio.Phylo``.

The algorithm computes the total branch length of the minimal subtree
(Steiner tree) connecting a reference leaf to a set of other leaves, by
walking each leaf up to the root and stopping as soon as a previously
visited (or reference) node is reached.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

from Bio import Phylo


class Tree:
    """A parsed Newick tree with parent pointers and taxon-ID lookup."""

    def __init__(self, newick_path: str) -> None:
        self.tree = Phylo.read(newick_path, "newick")
        self._parents: Dict[int, object] = {}
        self._by_name: Dict[str, object] = {}
        for clade in self.tree.find_clades(order="level"):
            if clade.name:
                self._by_name[clade.name] = clade
            for child in clade.clades:
                self._parents[id(child)] = clade

    def node(self, taxon_id: str) -> Optional[object]:
        return self._by_name.get(str(taxon_id))

    def _parent(self, clade: object) -> Optional[object]:
        return self._parents.get(id(clade))

    def branch_length(self, ref_id: str, org_ids: Iterable[str]) -> float:
        """Total branch length of the subtree spanning ``ref_id`` and ``org_ids``."""
        ref_node = self.node(ref_id)
        if ref_node is None:
            raise KeyError(f"Taxon {ref_id} not found in tree")

        ref_ancestors: Dict[object, List[object]] = {}
        ref_cumul_dist: Dict[object, float] = {}

        place = ref_node
        cumul = 0.0
        path: List[object] = []
        while place is not None:
            path.append(place)
            ref_ancestors[place] = list(path)
            ref_cumul_dist[place] = cumul
            if place.branch_length:
                cumul += place.branch_length
            place = self._parent(place)

        included: set = set()
        total = 0.0

        for org_id in org_ids:
            org_node = self.node(org_id)
            if org_node is None:
                raise KeyError(f"Taxon {org_id} not found in tree")

            place = org_node
            cumul = 0.0
            while place is not None:
                if place in included:
                    total += cumul
                    break
                if place in ref_ancestors:
                    total += ref_cumul_dist[place] + cumul
                    for cur in reversed(ref_ancestors[place]):
                        if cur in included:
                            total -= ref_cumul_dist[cur]
                            break
                        included.add(cur)
                    included.add(place)
                    break
                included.add(place)
                cumul += place.branch_length or 0.0
                place = self._parent(place)

        return total
