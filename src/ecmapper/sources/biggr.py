"""BiGG universal reactions -> EC numbers and cross-references from BiGGr (bigg.bio)."""

from __future__ import annotations

import json
from pathlib import Path

from .base import ECSource


class BiGGrSource(ECSource):
    """EC numbers and cross-references of BiGG reactions from the BiGGr reactions download.

    The download is a JSON list with one entry per (charge variant of a) reaction.
    The EC numbers are the ``cross_references`` with namespace ``ec-code``; in
    BiGGr these come from the Rhea reference reaction of the universal reaction.
    The other cross-references (KEGG, ModelSEED, MetaCyc, Rhea, ...) are kept
    and are available through :meth:`xrefs`. Entries are indexed by the
    universal reaction ID and by its superseded legacy IDs
    (``universalreaction__old_bigg_ids``).

    Parameters
    ----------
    include_obsolete : bool
        Also keep cross-references that BiGGr took from ModelSEED reactions
        that are marked as obsolete.
    """

    name = "biggr"
    namespace = "bigg.reaction"
    files = {"biggr_reactions.json": "https://bigg.bio/api/v3/download/reactions"}

    def __init__(self, *args, include_obsolete: bool = False, **kwargs):
        super().__init__(*args, **kwargs)
        self.include_obsolete = include_obsolete

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        with open(paths["biggr_reactions.json"], encoding="utf-8") as f:
            entries = json.load(f)
        if not isinstance(entries, list):
            raise ValueError("unexpected format of the BiGGr reactions download")
        table: dict[str, set[str]] = {}
        xrefs: dict[str, dict[str, set[str]]] = {}
        for entry in entries:
            universal_id = entry["universalreaction__bigg_id"]
            ids = [universal_id, *(entry.get("universalreaction__old_bigg_ids") or [])]
            for xref in entry.get("cross_references") or []:
                if xref.get("source_obsolete") and not self.include_obsolete:
                    continue
                namespace, identifier = xref["namespace"], str(xref["identifier"])
                for bigg_id in ids:
                    if namespace == "ec-code":
                        self._add(table, bigg_id, identifier)
                    else:
                        xrefs.setdefault(bigg_id, {}).setdefault(namespace, set()).add(identifier)
        self._xrefs = xrefs
        return table
