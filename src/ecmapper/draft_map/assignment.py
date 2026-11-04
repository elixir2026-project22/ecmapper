"""The assignment of subsystems to maps and of model reactions to subsystems as an editable XLSX table (--split).

Sheet "Subsystems": one row per subsystem, one column per map, an "x" marks a map of the subsystem (none: the
subsystem is not drawn; several: it is drawn on each of these maps). Further columns (new maps) may be added; their
order is the order of the maps.
Sheet "Reactions": one row per model reaction (no boundary reactions) with the IDs of the subsystems in which it is
drawn ("Subsystems", separated by ";"); a reaction in several subsystems is drawn in each of them, on their maps.
Sheet "Arrangement" (optional, written by the GUI, :mod:`.gui`): the positions of the subsystem boxes (top left, map
units) on the maps arranged by hand; the other maps are arranged by the pipeline.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .core import Completion, DraftMap, ModelData
from .split import GROUP_NAMES

SUBSYSTEMS = "Subsystems"
REACTIONS = "Reactions"
HELP = "Help"
ARRANGEMENT = "Arrangement"
ARRANGEMENT_COLUMNS = ["Map", "Subsystem ID", "x", "y"]
#: columns of the subsystems sheet before the map columns
SUBSYSTEM_COLUMNS = ["Subsystem ID", "Subsystem", "Colour group", "Set aside (rule)", "Model reactions"]
REACTION_COLUMNS = ["Reaction ID", "Name", "Subsystems", "Maps (when written)", "Equation", "Equation (names)",
                    "Drawn as", "Model subsystem", "EC numbers"]
#: fill of reactions added from the model (not on the BRENDA map), as their enzyme nodes in the SVG (render.py)
ADDED_COLOUR = "CDA0E6"
#: fill of reactions on the BRENDA map that are also added from the model in another subsystem
MIXED_COLOUR = "EEDDF7"
HELP_TEXT = [
    "Assignment of the draft map (ecmapper.draft_map.run --split) to maps; edit it and run again with "
    "--assignment FILE.",
    "",
    "Sheet Subsystems: one row per subsystem of the draft map, one column per map. An x marks the map of the "
    "subsystem.",
    "- Delete the x of a subsystem to leave it out (its reactions are then only drawn where they are in other "
    "subsystems).",
    "- Move the x to put the subsystem on another map. With an x in several columns, the subsystem is drawn on each "
    "of these maps.",
    "- Add a column with a new name (after the existing map columns) for a new map; the order of the columns is the "
    "order of the maps. Maps without subsystems are left out.",
    "- The columns before the maps are for information; changes there are ignored.",
    "",
    "Sheet Reactions: one row per model reaction (no exchange, sink, demand or biomass reactions). Subsystems: the "
    "IDs of the subsystems (first column of sheet Subsystems) in which the reaction is drawn, separated by ';'.",
    "- Remove a subsystem ID to take the reaction out of that subsystem; empty: the reaction is not drawn.",
    "- Add a subsystem ID to draw the reaction there as well (as an added reaction, connected to the metabolites "
    "drawn in that subsystem); with IDs of subsystems on different maps, the reaction is on several maps.",
    "- Rows that are deleted keep their reaction as it is; the other columns are for information.",
    "- Purple rows: reactions added from the model (step 3; not on the BRENDA map, drawn like their purple, dashed "
    "nodes in the SVG). Light purple 'Drawn as': on the BRENDA map in one subsystem, added in another one.",
    "",
    "Sheet Arrangement (written by the GUI, python -m ecmapper.draft_map.gui): the positions (top left of the box, "
    "map units) of the subsystems on the maps arranged by hand. A map with the positions of all its subsystems is "
    "drawn so (boxes closer than --gap are moved apart); the other maps are arranged by the pipeline.",
    "",
    "The table is read at the same point of the pipeline at which it was written (after step 3, the completion), "
    "so run with the same prepared model and the same options of steps 1-3.",
]


@dataclass
class AssignmentTable:
    maps: list[str] = field(default_factory=list)  # in order
    subsystems: dict[str, list[str]] = field(default_factory=dict)  # subsystem ID -> map names (none: not drawn)
    reactions: dict[str, list[str]] = field(default_factory=dict)  # reaction ID -> subsystem IDs
    positions: dict[str, dict[str, tuple[float, float]]] = field(default_factory=dict)  # map -> subsystem -> x, y

    def map_assignment(self, dm: DraftMap) -> tuple[list[str], dict[str, list[str]]]:
        """Map names in order (those with subsystems) and {subsystem ID: map names} of the subsystems of `dm`."""
        assignment = {k: self.subsystems.get(k, []) for k in dm.subsystems}
        used = {m for maps in assignment.values() for m in maps}
        return [m for m in self.maps if m in used], assignment


def reaction_subsystems(dm: DraftMap) -> dict[str, list[str]]:
    """{reaction ID: subsystems in which it is drawn}."""
    result = defaultdict(list)
    for n in dm.nodes.values():
        for r in n.reactions:
            if n.subsystem not in result[r]:
                result[r].append(n.subsystem)
    return result


def origin(statuses) -> str:
    """Where drawn reactions come from by the statuses of their enzyme nodes: "model" (all added from the model,
    step 3), "both" (also added) or "map" (on the BRENDA map or its vicinity); "" if not drawn."""
    statuses = set(statuses)
    if not statuses:
        return ""
    return "model" if statuses == {"added"} else "both" if "added" in statuses else "map"


def equation(md: ModelData, rid: str, names: bool = False) -> str:
    r = md.reactions[rid]

    def side(sign):
        terms = []
        for m, c in r["metabolites"].items():
            if c * sign > 0:
                label = md.metabolites[m]["name"] or m if names else m
                terms.append(label if abs(c) == 1 else f"{abs(c):g} {label}")
        return " + ".join(terms)

    return f"{side(-1)} {'<=>' if r['reversible'] else '-->'} {side(1)}".strip()


def write_assignment(path, dm: DraftMap, md: ModelData, order: list[str], assignment: dict[str, list[str]],
                     positions: dict[str, dict[str, tuple[float, float]]] | None = None) -> None:
    """XLSX table of the assignment of the subsystems of `dm` to the maps and of the model reactions to subsystems.

    `assignment`: {subsystem ID: map names}; `positions`: {map name:
    {subsystem ID: (x, y)}} of the maps arranged by hand (sheet Arrangement).
    """
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    from openpyxl.worksheet.datavalidation import DataValidation

    rank = {m: i for i, m in enumerate(order)}
    in_subsystem = reaction_subsystems(dm)
    counts = Counter(s for subs in in_subsystem.values() for s in subs)
    bold, head_fill = Font(bold=True), PatternFill("solid", fgColor="DDDDDD")
    wb = Workbook()

    ws = wb.active
    ws.title = SUBSYSTEMS
    ws.append(SUBSYSTEM_COLUMNS + order)
    keys = sorted(dm.subsystems, key=lambda k: (min((rank.get(m, len(order)) for m in assignment.get(k, [])),
                                                    default=len(order)), dm.subsystems[k].label.lower()))
    for k in keys:
        s = dm.subsystems[k]
        group = f"{s.color_group} ({GROUP_NAMES[s.color_group]})" if s.color_group in GROUP_NAMES else s.color_group
        maps = assignment.get(k, [])
        ws.append([k, s.label, group, s.aside, counts[k]] + ["x" if m in maps else None for m in order])
    first = len(SUBSYSTEM_COLUMNS) + 1
    for i, width in enumerate([34, 40, 34, 10, 10], 1):
        ws.column_dimensions[get_column_letter(i)].width = width
    for j in range(first, first + len(order) + 10):
        ws.column_dimensions[get_column_letter(j)].width = 14
    for c in ws[1]:
        c.font, c.fill, c.alignment = bold, head_fill, Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[1].height = 45
    ws.freeze_panes = ws.cell(2, first)
    ws.auto_filter.ref = f"A1:{get_column_letter(first + len(order) - 1)}{len(keys) + 1}"
    dv = DataValidation(type="list", formula1='"x"', allow_blank=True)
    ws.add_data_validation(dv)
    dv.add(f"{get_column_letter(first)}2:{get_column_letter(first + len(order) + 9)}{len(keys) + 1}")

    ws = wb.create_sheet(REACTIONS)
    ws.append(REACTION_COLUMNS)
    status = defaultdict(set)
    for n in dm.nodes.values():
        for r in n.reactions:
            status[r].add(n.status)
    label_rank = {k: i for i, k in enumerate(keys)}
    rids = [r for r, v in md.reactions.items() if not v["boundary"]]
    rids.sort(key=lambda r: (min((label_rank[s] for s in in_subsystem.get(r, [])), default=len(keys)), r))
    for rid in rids:
        r, subs = md.reactions[rid], in_subsystem.get(rid, [])
        maps = sorted({m for s in subs for m in assignment.get(s, [])}, key=lambda m: rank.get(m, len(order)))
        ws.append([rid, r["name"], "; ".join(subs), "; ".join(maps), equation(md, rid), equation(md, rid, True),
                   ", ".join(sorted(status[rid])) if subs else ("transport, not drawn" if r["transport"] else
                                                                 "not drawn"),
                   r["subsystem"], ", ".join(r["ecs"])])
    for i, width in enumerate([18, 40, 40, 30, 60, 80, 14, 30, 16], 1):
        ws.column_dimensions[get_column_letter(i)].width = width
    for c in ws[1]:
        c.font, c.fill = bold, head_fill
    edit_fill = PatternFill("solid", fgColor="FFF2CC")
    added_fill, mixed_fill = PatternFill("solid", fgColor=ADDED_COLOUR), PatternFill("solid", fgColor=MIXED_COLOUR)
    i_drawn = REACTION_COLUMNS.index("Drawn as") + 1
    ws.cell(1, 3).fill = PatternFill("solid", fgColor="FFD966")
    for row, rid in enumerate(rids, 2):
        ws.cell(row, 3).fill = edit_fill  # the column to edit keeps its colour
        source = origin(status[rid]) if in_subsystem.get(rid) else ""
        if source == "model":
            for j in range(1, len(REACTION_COLUMNS) + 1):
                if j != 3:
                    ws.cell(row, j).fill = added_fill
        elif source == "both":
            ws.cell(row, i_drawn).fill = mixed_fill
    ws.freeze_panes = "C2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(REACTION_COLUMNS))}{len(rids) + 1}"

    if positions:
        ws = wb.create_sheet(ARRANGEMENT)
        ws.append(ARRANGEMENT_COLUMNS)
        for name in sorted(positions, key=lambda m: rank.get(m, len(order))):
            for k, (x, y) in positions[name].items():
                ws.append([name, k, round(float(x), 1), round(float(y), 1)])
        for i, width in enumerate([30, 40, 10, 10], 1):
            ws.column_dimensions[get_column_letter(i)].width = width
        for c in ws[1]:
            c.font, c.fill = bold, head_fill
        ws.freeze_panes = "A2"

    ws = wb.create_sheet(HELP)
    for line in HELP_TEXT:
        ws.append([line])
    ws.column_dimensions["A"].width = 140
    wb.save(path)


def read_assignment(path) -> AssignmentTable:
    """The (edited) table of :func:`write_assignment`."""
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    table = AssignmentTable()
    rows = wb[SUBSYSTEMS].iter_rows(values_only=True)
    header = [str(c).strip() if c is not None else "" for c in next(rows)]
    columns = [(j, name) for j, name in enumerate(header) if j >= len(SUBSYSTEM_COLUMNS) and name]
    table.maps = list(dict.fromkeys(name for _, name in columns))
    for row in rows:
        if not row or row[0] is None or not str(row[0]).strip():
            continue
        k = str(row[0]).strip()
        marked = [name for j, name in columns if j < len(row) and row[j] is not None and str(row[j]).strip()]
        table.subsystems[k] = list(dict.fromkeys(marked))
    rows = wb[REACTIONS].iter_rows(values_only=True)
    header = [str(c).strip() if c is not None else "" for c in next(rows)]
    i_id, i_subs = header.index("Reaction ID"), header.index("Subsystems")
    for row in rows:
        if not row or row[i_id] is None or not str(row[i_id]).strip():
            continue
        value = row[i_subs] if i_subs < len(row) and row[i_subs] is not None else ""
        subs = [s for s in re.split(r"[;,\s]+", str(value)) if s]
        table.reactions[str(row[i_id]).strip()] = list(dict.fromkeys(subs))
    if ARRANGEMENT in wb.sheetnames:
        rows = wb[ARRANGEMENT].iter_rows(values_only=True)
        header = [str(c).strip() if c is not None else "" for c in next(rows, [])]
        i_map, i_sub, i_x, i_y = (header.index(c) for c in ARRANGEMENT_COLUMNS)
        for row in rows:
            if not row or any(i >= len(row) or row[i] is None or str(row[i]).strip() == ""
                              for i in (i_map, i_sub, i_x, i_y)):
                continue
            name, k = str(row[i_map]).strip(), str(row[i_sub]).strip()
            table.positions.setdefault(name, {})[k] = (float(row[i_x]), float(row[i_y]))
    wb.close()
    return table


def remove_reaction(dm: DraftMap, rid: str, subsystem: str) -> None:
    """Take reaction `rid` out of the enzyme nodes of `subsystem`; nodes left without reactions are removed."""
    gone = set()
    for n in dm.members(subsystem):
        if n.kind == "enzyme" and rid in n.reactions:
            n.reactions.remove(rid)
            if not n.reactions:
                gone.add(n.id)
            elif n.status == "added":
                n.label = n.reactions[0] if len(n.reactions) == 1 else f"{n.reactions[0]} (+{len(n.reactions) - 1})"
    if not gone:
        return
    mets = {e.met for e in dm.edges if e.enzyme in gone}
    dm.remove_nodes(gone)
    used = {e.met for e in dm.edges}
    dm.remove_nodes(mets - used)
    dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]


def drop_subsystems(dm: DraftMap, keys: set[str]) -> None:
    dm.subsystems = {k: v for k, v in dm.subsystems.items() if k not in keys}
    dm.remove_nodes({n.id for n in dm.nodes.values() if n.subsystem in keys})
    dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]


def apply_assignment(dm: DraftMap, md: ModelData, table: AssignmentTable,
                     origin_positions: dict[str, tuple[float, float, str]] | None = None,
                     drop_unassigned: bool = True) -> dict:
    """Apply the reaction sheet and leave out the subsystems without map; after step 3 (completion).

    Reactions added to a subsystem are drawn as in the completion (step 3),
    with the metabolites drawn in the subsystem; reactions taken out of all
    subsystems of an enzyme node remove the node. Subsystems left without
    enzyme nodes are dropped. `drop_unassigned`: False keeps the subsystems
    without map (the GUI shows them). Returns what was changed and warnings.
    """
    warnings = []
    missing = sorted(k for k in dm.subsystems if k not in table.subsystems)
    if missing:
        warnings.append(f"subsystems not in the table, {'left out' if drop_unassigned else 'without map'}: "
                        f"{', '.join(missing)}")
    unknown = sorted(k for k in table.subsystems if k not in dm.subsystems)
    if unknown:
        warnings.append(f"subsystems of the table not in the draft map, ignored: {', '.join(unknown)}")
    no_map = {k for k in dm.subsystems if not table.subsystems.get(k)}
    if drop_unassigned:
        drop_subsystems(dm, no_map)

    current = reaction_subsystems(dm)
    add, remove = [], []
    for rid, wanted in table.reactions.items():
        if rid not in md.reactions:
            warnings.append(f"reaction {rid} not in the model, ignored")
            continue
        for s in wanted:
            if s not in dm.subsystems:
                if s not in table.subsystems:
                    warnings.append(f"reaction {rid}: unknown subsystem {s}, ignored")
            elif s not in current.get(rid, []):
                add.append((rid, s))
        remove.extend((rid, s) for s in current.get(rid, []) if s not in wanted)
    # added first: a subsystem whose reactions are all replaced still has its nodes for placing the new ones
    completion = Completion(dm, md, origin_positions=origin_positions, prefix="user_", anchor_aside=True)
    for rid, s in add:
        mains = completion.main_mets(rid)
        anchors = completion.anchors()
        completion.add(rid, mains, {m: anchors[m] for m in mains if anchors.get(m)}, {}, target=s)
    for rid, s in remove:
        remove_reaction(dm, rid, s)
    empty = {k for k in dm.subsystems if not any(n.kind == "enzyme" for n in dm.members(k))}
    drop_subsystems(dm, empty)
    dm.log.append(
        f"assignment: {len(no_map)} subsystems without map {'left out' if drop_unassigned else 'kept'}; {len(add)} "
        f"reactions added to and {len(remove)} "
        f"taken out of subsystems; {len(empty)} subsystems left without reactions dropped; {len(warnings)} warnings"
    )
    return dict(subsystems_left_out=sorted(no_map) if drop_unassigned else [],
                added=[dict(reaction=r, subsystem=s) for r, s in add],
                removed=[dict(reaction=r, subsystem=s) for r, s in remove], empty_subsystems=sorted(empty),
                warnings=warnings)
