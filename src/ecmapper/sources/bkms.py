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
    """

    name = "bkms"
    namespace = "metacyc.reaction"
    files = {"Reactions_BKMS.tar.gz": "https://bkms.brenda-enzymes.org/download/Reactions_BKMS.tar.gz"}

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        table: dict[str, set[str]] = {}
        with tarfile.open(paths["Reactions_BKMS.tar.gz"], "r:gz") as tar:
            member = next((m for m in tar.getmembers() if m.isfile() and m.name.endswith((".csv", ".tsv"))), None)
            if member is None:
                raise ValueError("no reaction table found in the BKMS-react download")
            raw = tar.extractfile(member).read()
        text = raw.decode("utf-8", errors="replace").replace("\r", "")
        for row in csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE):
            for reaction in (row.get("Reaction_ID_MetaCyc") or "").split(","):
                if reaction.strip():
                    self._add(table, reaction.strip(), row.get("EC_Number") or "")
        return table
