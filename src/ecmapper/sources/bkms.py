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

    The reaction equations (``Reaction``) are indexed by BRENDA reaction ID
    (:meth:`equations`) and by EC number (:meth:`ec_equations`). BRENDA
    equations (``A + B = C + D``) are written in the direction of the natural
    substrates; equations of KEGG reactions without BRENDA reaction (``<=>``)
    are not indexed.
    """

    name = "bkms"
    namespace = "metacyc.reaction"
    files = {"Reactions_BKMS.tar.gz": "https://bkms.brenda-enzymes.org/download/Reactions_BKMS.tar.gz"}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._kegg: dict[str, set[str]] = {}
        self._brenda: dict[str, dict[str, set[str]]] = {}
        self._equations: dict[str, list[tuple[str, str]]] = {}
        self._ec_equations: dict[str, list[str]] = {}

    def brenda_reaction(self, identifier: str) -> dict[str, set[str]]:
        """KEGG and MetaCyc reactions aligned with the BRENDA reaction `identifier` (e.g. ``BR344``)."""
        self.data  # noqa: B018, loads the files
        return {ns: set(ids) for ns, ids in self._brenda.get(identifier, {}).items()}

    def equations(self, identifier: str) -> list[tuple[str, str]]:
        """(EC number, equation) of the BRENDA reaction `identifier` (e.g. ``BR344``), one per BKMS row."""
        self.data  # noqa: B018, loads the files
        return list(self._equations.get(identifier, ()))

    def ec_equations(self, ec: str) -> list[str]:
        """The BRENDA equations of all BKMS rows with EC number `ec`."""
        self.data  # noqa: B018, loads the files
        return list(self._ec_equations.get(ec, ()))

    def lookup_kegg(self, identifier: str) -> set[str]:
        """EC numbers of the BKMS reactions aligned with the KEGG reaction `identifier`."""
        self.data  # noqa: B018, loads the files
        return set(self._kegg.get(identifier, ()))

    def reset(self) -> None:
        super().reset()
        self._kegg, self._brenda, self._equations, self._ec_equations = {}, {}, {}, {}

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
        brenda: dict[str, dict[str, set[str]]] = {}
        equations: dict[str, list[tuple[str, str]]] = {}
        ec_equations: dict[str, list[str]] = {}
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
            equation = (row.get("Reaction") or "").strip()
            brenda_ids = [r.strip() for r in (row.get("Reaction_ID_BRENDA") or "").split(",") if r.strip()]
            if brenda_ids and " = " in equation:
                for reaction in brenda_ids:
                    equations.setdefault(reaction, []).append((ec.strip(), equation))
                if ec.strip():
                    ec_equations.setdefault(ec.strip(), []).append(equation)
            for reaction in (row.get("Reaction_ID_BRENDA") or "").split(","):
                if reaction.strip():
                    aligned = brenda.setdefault(reaction.strip(), {"kegg.reaction": set(), "metacyc.reaction": set()})
                    for namespace, column in (("kegg.reaction", "Reaction_ID_KEGG"), ("metacyc.reaction", "Reaction_ID_MetaCyc")):
                        aligned[namespace] |= {i.strip() for i in (row.get(column) or "").split(",") if i.strip()}
        self._kegg, self._brenda, self._equations, self._ec_equations = kegg, brenda, equations, ec_equations
        return table
