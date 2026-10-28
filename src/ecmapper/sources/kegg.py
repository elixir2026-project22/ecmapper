"""KEGG REACTION -> EC numbers via the KEGG REST API (free for academic use only)."""

from __future__ import annotations

from pathlib import Path

from .base import CachedResource, ECSource


class KEGGSource(ECSource):
    """EC numbers of KEGG reactions from ``rest.kegg.jp/link/ec/rn``.

    The file has lines like ``rn:R00004<TAB>ec:3.6.1.1``.
    """

    name = "kegg"
    namespace = "kegg.reaction"
    files = {"kegg_link_ec_rn.tsv": "https://rest.kegg.jp/link/ec/rn"}

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        table: dict[str, set[str]] = {}
        with open(paths["kegg_link_ec_rn.tsv"], encoding="utf-8") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) != 2:
                    continue
                reaction, ec = fields
                self._add(table, reaction.removeprefix("rn:"), ec.removeprefix("ec:"))
        return table


class KEGGReactionCompounds(CachedResource):
    """Compounds of KEGG reactions from ``rest.kegg.jp/link/cpd/rn`` (``rn:R00004<TAB>cpd:C00001``).

    The links do not say on which side of the reaction a compound is.
    """

    name = "kegg-compounds"
    files = {"kegg_link_cpd_rn.tsv": "https://rest.kegg.jp/link/cpd/rn"}

    def _parse(self, paths: dict[str, Path]) -> dict[str, frozenset[str]]:
        compounds: dict[str, set[str]] = {}
        with open(paths["kegg_link_cpd_rn.tsv"], encoding="utf-8") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) == 2 and fields[1].startswith("cpd:"):
                    compounds.setdefault(fields[0].removeprefix("rn:"), set()).add(fields[1][4:])
        return {reaction: frozenset(c) for reaction, c in compounds.items()}

    def compounds(self, reaction: str) -> frozenset[str]:
        """KEGG compound IDs of the KEGG reaction `reaction`."""
        return self.data.get(reaction, frozenset())
