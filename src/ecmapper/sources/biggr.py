"""BiGG universal reactions -> EC numbers and cross-references from BiGGr (bigg.bio)."""

from __future__ import annotations

import json
import re
from pathlib import Path

from .base import CachedResource, ECSource


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


def bigg_key(bigg_id: str) -> str:
    """Key for comparing BiGG IDs: legacy IDs (``asn_L``) and current ones (``asn__L``) differ in underscores."""
    return re.sub(r"_+", "_", bigg_id)


class BiGGrMetabolites(CachedResource):
    """Cross-references of BiGG metabolites from the BiGGr metabolites download.

    The download is a JSON list with one entry per (charge variant of a)
    metabolite with ``cross_references`` (KEGG, ModelSEED, MetaCyc, MetaNetX,
    ChEBI, ...) and ``references`` (ChEBI). They are indexed by universal
    metabolite ID and its superseded IDs, compared by :func:`bigg_key`. The
    names of the metabolite are kept as pseudo-namespace ``name``.
    """

    name = "biggr-metabolites"
    files = {"biggr_metabolites.json": "https://bigg.bio/api/v3/download/metabolites"}

    def _parse(self, paths: dict[str, Path]) -> dict[str, dict[str, set[str]]]:
        with open(paths["biggr_metabolites.json"], encoding="utf-8") as f:
            entries = json.load(f)
        if not isinstance(entries, list):
            raise ValueError("unexpected format of the BiGGr metabolites download")
        xrefs: dict[str, dict[str, set[str]]] = {}
        for entry in entries:
            ids = [entry["universalcomponent__bigg_id"], *(entry.get("universalcomponent__old_bigg_ids") or [])]
            found = [(x["namespace"].lower(), str(x["identifier"])) for x in entry.get("cross_references") or []]
            found += [("chebi", r.split(":", 1)[1]) for r in entry.get("references") or [] if r.upper().startswith("CHEBI:")]
            found += [("name", n) for n in (entry.get("name"), entry.get("universalcomponent__name")) if n]
            for bigg_id in ids:
                target = xrefs.setdefault(bigg_key(bigg_id), {})
                for namespace, identifier in found:
                    target.setdefault(namespace, set()).add(identifier)
        return xrefs

    def xrefs(self, bigg_id: str) -> dict[str, set[str]]:
        """Cross-references {namespace: identifiers} of the universal BiGG metabolite `bigg_id`."""
        return {ns: set(ids) for ns, ids in self.data.get(bigg_key(bigg_id), {}).items()}
