"""Rhea reactions -> EC numbers from the Rhea TSV downloads."""

from __future__ import annotations

import csv
from pathlib import Path

from .base import ECSource

_RHEA_TSV = "https://ftp.expasy.org/databases/rhea/tsv/"


class RheaSource(ECSource):
    """EC numbers of Rhea reactions from ``rhea2ec.tsv``.

    ``rhea2ec.tsv`` only lists the undirected master reactions. Models are also
    annotated with the directed variants (left-to-right, right-to-left and
    bidirectional), these are mapped to their master reaction with
    ``rhea-directions.tsv``.
    """

    name = "rhea"
    namespace = "rhea"
    files = {
        "rhea2ec.tsv": _RHEA_TSV + "rhea2ec.tsv",
        "rhea-directions.tsv": _RHEA_TSV + "rhea-directions.tsv",
    }

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        table: dict[str, set[str]] = {}
        with open(paths["rhea2ec.tsv"], encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                self._add(table, row["MASTER_ID"].strip(), row["ID"])
        with open(paths["rhea-directions.tsv"], encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                ecs = table.get(row["RHEA_ID_MASTER"].strip())
                if ecs:
                    for column in ("RHEA_ID_LR", "RHEA_ID_RL", "RHEA_ID_BI"):
                        table.setdefault(row[column].strip(), set()).update(ecs)
        return table
