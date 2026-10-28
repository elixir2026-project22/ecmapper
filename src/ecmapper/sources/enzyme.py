"""Current EC numbers from the ExPASy ENZYME nomenclature database."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .base import CachedResource

_EC = re.compile(r"\d+\.\d+\.\d+\.n?\d+")


@dataclass
class EnzymeEntry:
    """An entry of the ENZYME database.

    Attributes
    ----------
    ec : str
    name : str
        Recommended name, empty for transferred and deleted entries.
    transferred_to : tuple of str
        The EC numbers the entry was transferred to.
    deleted : bool
    organisms : frozenset of str
        Swiss-Prot organism mnemonics (e.g. ``ECOLI``) of the proteins with this EC number.
    """

    ec: str
    name: str = ""
    transferred_to: tuple[str, ...] = ()
    deleted: bool = False
    organisms: frozenset[str] = field(default_factory=frozenset)


class EnzymeNomenclature(CachedResource):
    """Transferred and deleted EC numbers and Swiss-Prot organisms per EC number from ExPASy ENZYME.

    ``enzyme.dat`` has one entry per EC number. Transferred entries have a
    description like ``Transferred entry: 1.7.1.1, 1.7.1.2 and 1.9.6.1.``,
    deleted ones ``Deleted entry.``; the ``DR`` lines list the Swiss-Prot entries
    (``Q9RHN2, ACDH2_COMTE``) whose organism mnemonic is the part after ``_``.
    """

    name = "enzyme"
    files = {"enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat"}

    def _parse(self, paths: dict[str, Path]) -> dict[str, EnzymeEntry]:
        entries: dict[str, EnzymeEntry] = {}
        ec, description, organisms = None, [], set()

        def store():
            if ec is None:
                return
            text = " ".join(description)
            entry = EnzymeEntry(ec, organisms=frozenset(organisms))
            if text.startswith("Transferred entry:"):
                entry.transferred_to = tuple(_EC.findall(text))
            elif text.startswith("Deleted entry"):
                entry.deleted = True
            else:
                entry.name = text.rstrip(".")
            entries[ec] = entry

        with open(paths["enzyme.dat"], encoding="utf-8") as f:
            for line in f:
                code, value = line[:2], line[5:].strip()
                if code == "ID":
                    ec, description, organisms = value, [], set()
                elif code == "DE":
                    description.append(value)
                elif code == "DR":
                    organisms |= set(re.findall(r"_([A-Z0-9]+)\s*;", line))
                elif code == "//":
                    store()
                    ec = None
        store()
        if not entries:
            raise ValueError("no entries found in enzyme.dat")
        return entries

    def entry(self, ec: str) -> EnzymeEntry | None:
        return self.data.get(ec)

    def current(self, ec: str) -> set[str]:
        """The current EC number(s) for `ec`.

        Transferred EC numbers are replaced by the ones they were transferred
        to (recursively), deleted ones give an empty set. EC numbers that are not
        in ENZYME (incomplete or BRENDA provisional numbers, new or unknown EC
        numbers) are returned unchanged.
        """
        result, todo, seen = set(), [ec], set()
        while todo:
            ec = todo.pop()
            if ec in seen:
                continue
            seen.add(ec)
            entry = self.data.get(ec)
            if entry is None or not (entry.transferred_to or entry.deleted):
                result.add(ec)
            elif entry.transferred_to:
                todo.extend(entry.transferred_to)
        return result

    def organisms(self, ec: str) -> frozenset[str]:
        """Swiss-Prot organism mnemonics of the proteins with EC number `ec`."""
        entry = self.data.get(ec)
        return entry.organisms if entry else frozenset()
