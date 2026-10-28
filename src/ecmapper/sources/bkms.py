"""MetaCyc reactions -> EC numbers from BKMS-react (bkms.brenda-enzymes.org)."""

from __future__ import annotations

import csv
import io
import tarfile
from pathlib import Path

from .base import ECSource


class BKMSSource(ECSource):
    """EC numbers of MetaCyc reactions from the BKMS-react download.

    ``Reactions_BKMS.tar.gz`` contains ``Reactions_BKMS.csv``, which is actually
    tab separated. Every row is one non-redundant reaction with one EC number in
    ``EC_Number`` and the aligned MetaCyc reactions, comma separated, in
    ``Reaction_ID_MetaCyc``. Rows without valid EC number (empty, ``SPONTANEOUS``)
    are skipped. Note that the MetaCyc data in BKMS-react is from MetaCyc 20.0,
    reactions added to MetaCyc later are not covered.

    The aligned KEGG reactions (``Reaction_ID_KEGG``) are indexed as well and
    available through :meth:`lookup_kegg`. They are not used as EC source of
    their own because BKMS aligns KEGG reactions with generic BRENDA reactions
    (e.g. R00754 ethanol dehydrogenase with 24 EC numbers of alcohol
    dehydrogenases), but they are used for :meth:`consensus`.
    """

    name = "bkms"
    namespace = "metacyc.reaction"
    files = {"Reactions_BKMS.tar.gz": "https://bkms.brenda-enzymes.org/download/Reactions_BKMS.tar.gz"}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._kegg: dict[str, set[str]] = {}

    def lookup_kegg(self, identifier: str) -> set[str]:
        """EC numbers of the BKMS reactions aligned with the KEGG reaction `identifier`."""
        self.data  # noqa: B018, loads the files
        return set(self._kegg.get(identifier, ()))

    def reset(self) -> None:
        super().reset()
        self._kegg = {}

    def consensus(self, identifiers: dict[str, list[str]]) -> set[str]:
        """EC numbers of all BKMS reactions that match the KEGG or MetaCyc reactions in `identifiers`.

        `identifiers` is {namespace: [identifiers]}, only ``kegg.reaction`` and
        ``metacyc.reaction`` are used.
        """
        ecs = set()
        for identifier in identifiers.get("metacyc.reaction", []):
            ecs |= self.lookup(identifier)
        for identifier in identifiers.get("kegg.reaction", []):
            ecs |= self.lookup_kegg(identifier)
        return ecs

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        table: dict[str, set[str]] = {}
        kegg: dict[str, set[str]] = {}
        with tarfile.open(paths["Reactions_BKMS.tar.gz"], "r:gz") as tar:
            member = next((m for m in tar.getmembers() if m.isfile() and m.name.endswith((".csv", ".tsv"))), None)
            if member is None:
                raise ValueError("no reaction table found in the BKMS-react download")
            raw = tar.extractfile(member).read()
        text = raw.decode("utf-8", errors="replace").replace("\r", "")
        for row in csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE):
            ec = row.get("EC_Number") or ""
            for reaction in (row.get("Reaction_ID_MetaCyc") or "").split(","):
                if reaction.strip():
                    self._add(table, reaction.strip(), ec)
            for reaction in (row.get("Reaction_ID_KEGG") or "").split(","):
                if reaction.strip():
                    self._add(kegg, reaction.strip(), ec)
        self._kegg = kegg
        return table
