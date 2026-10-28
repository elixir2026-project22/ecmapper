"""ModelSEED reactions -> EC numbers from the ModelSEEDDatabase on GitHub."""

from __future__ import annotations

import csv
from pathlib import Path

from .base import ECSource


class ModelSEEDSource(ECSource):
    """EC numbers of ModelSEED reactions from ``Unique_ModelSEED_Reaction_ECs.txt``.

    The file is a TSV with the columns ``ModelSEED ID``, ``External ID`` (the EC
    number) and ``Source``. It also contains entries that are no EC numbers
    (e.g. transporter classification numbers), these are skipped.
    """

    name = "modelseed"
    namespace = "seed.reaction"
    files = {
        "Unique_ModelSEED_Reaction_ECs.txt": "https://raw.githubusercontent.com/ModelSEED/"
        "ModelSEEDDatabase/master/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt"
    }

    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        table: dict[str, set[str]] = {}
        with open(paths["Unique_ModelSEED_Reaction_ECs.txt"], encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                self._add(table, row["ModelSEED ID"].strip(), row["External ID"])
        return table
