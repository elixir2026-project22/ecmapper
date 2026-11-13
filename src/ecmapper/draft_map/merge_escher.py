"""One Escher map of a map's subsystems, merged by MergingEscher.

Every subsystem of the map is exported as an Escher map of its own, in one
of the :data:`STYLES`: with smooth curves (:func:`.escher.to_escher` with
`smooth`) or curved as the Escher maps of the BRENDA map
(:func:`ecmapper.brenda_escher.from_draft`); MergingEscher's
`escher_merge_auto.py` (a separate tool, run as a process: it has no
Python API) joins them on their common primary metabolites and places them
with few crossings.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from .core import DraftMap, ModelData
from .escher import to_escher
from .relayout import subsystem_map
from .split import slug

SETTING = "merge/merging_escher"
SCRIPT = "escher_merge_auto.py"
#: how the subsystems are drawn: {style: description}
STYLES = {"smooth": "Draft map with smooth curves (straight lines bent through the enzyme nodes)",
          "brenda": "Curved as the Escher maps of the BRENDA map (rows, columns, rounded corners; "
                    "ecmapper.brenda_escher)"}
#: the setting of the GUI with the style chosen last
STYLE_SETTING = "merge/style"


def subsystem_maps(m: DraftMap, md: ModelData, scale: float = 2.0, style: str = "smooth") -> dict[str, list]:
    """{subsystem ID: its Escher map} for the subsystems of `m` with reactions (in their positions on `m`),
    drawn in `style` (:data:`STYLES`)."""
    if style not in STYLES:
        raise ValueError(f"unknown style {style!r}, one of {sorted(STYLES)}")
    maps = {}
    for k, s in m.subsystems.items():
        sub = subsystem_map(m, k)
        if not any(n.kind == "enzyme" and n.reactions for n in sub.nodes.values()):
            continue
        if style == "brenda":
            from ..brenda_escher import build, from_draft

            reactions, data = from_draft(sub, md, scale)
            maps[k] = build(reactions, data, s.label or k, description="")
        else:
            maps[k] = to_escher(sub, md, s.label or k, scale=scale, smooth=True)
    return maps


def write_maps(maps: dict[str, list], folder) -> list[Path]:
    """Write the maps to `folder` (its earlier JSON files are removed), one file per map named after it.

    No name starts with "merged" (MergingEscher takes such files for its earlier results and skips them).
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob("*.json"):
        old.unlink()
    paths, used = [], set()
    for k, data in maps.items():
        name = slug(data[0]["map_name"]) or slug(k) or "subsystem"
        if name.startswith("merged"):
            name = "subsystem_" + name
        stem, i = name, 2
        while name in used:
            name, i = f"{stem}_{i}", i + 1
        used.add(name)
        path = folder / f"{name}.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        paths.append(path)
    return paths


def find_tool(location) -> Path | None:
    """MergingEscher's script from `location` (its folder, its src folder or the script); None: not found."""
    if not location:
        return None
    path = Path(location)
    for candidate in (path, path / SCRIPT, path / "src" / SCRIPT):
        if candidate.name == SCRIPT and candidate.is_file():
            return candidate
    return None


def merge_command(script, folder, output, name: str, model=None) -> list[str]:
    """The command that merges the maps in `folder` to `output`.

    The maps are of one model, so their IDs are kept as they are (MergingEscher
    would rewrite g6p_c to g6p[c]). With a COBRA JSON `model`, MergingEscher
    joins groups of maps without common metabolites by model reactions.
    """
    command = [sys.executable, "-u", str(script), str(folder), "--id-style", "none", "--layout", "min-crossings",
               "--name", name, "-o", str(output)]
    if model and Path(model).suffix.lower() == ".json" and Path(model).is_file():
        command += ["--model", str(model)]
    return command
