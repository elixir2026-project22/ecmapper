"""The state of the assignment GUI, without Qt: the maps, the subsystems on them and their positions.

A :class:`Session` matches a model to the BRENDA map if it is not prepared
yet (:mod:`..prepare`) and runs the pipeline of :mod:`..run` up to the condensed
subsystems (steps 1-4, with the table of --assignment), so the boxes have
the sizes they have on the final maps. Boxes on a map are always at least
--gap apart: moves are limited with :func:`constrained_delta`, pasted
boxes go to a :func:`free_position`.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import pickle
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from .. import run
from ..ai_layouts import DRAFT, MapLayout, export_draft, load_layouts, make_layout, match_ids, refresh_exit_labels
from ..ai_layouts import apply as apply_layout
from ..ai_layouts import check as check_layout
from ..assignment import equation, origin, reaction_subsystems, read_assignment, write_assignment
from ..core import DraftMap, hub_table, place_free
from ..layout import Ports, anneal_arrangement, arrangement_original, to_absolute, user_arrangement
from ..project import PREP, TABLE, Project, bundled_template, digests, ecmapper_version
from ..relayout import subsystem_map
from ..render import svg_text
from ..split import duplicate_subsystems, slug, split_map

#: the tab of the subsystems on no map
NOT_DRAWN = "Not drawn"
#: margin of the SVG of a subsystem around its box and title (map units)
SVG_MARGIN = 6.0
#: boxes closer than the gap by at most this much (map units, rounding) count as free
TOLERANCE = 1e-3
#: the condensed subsystems of the last start, in OUT_DIR
CACHE = "gui_cache.pkl"
#: the options of :mod:`..run` the drawings of the subsystems depend on (steps 1-4)
DRAWING_OPTIONS = ["keep_subsets", "aside_coverage", "aside_min_matched", "min_anchored", "no_completion", "relayout",
                   "circle", "relayout_margin", "compactness", "pack", "split"]


USAGE = "give MODEL [OUT_DIR] or PREP.pkl BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR (or nothing: a dialog)"


def template() -> tuple[Path, Path]:
    """The BRENDA map of ecmapper (BRENDAmap/ of the project), the default template of the draft maps."""
    return bundled_template()


def build_parser():
    parser = run.build_parser(positional=False)
    parser.prog = "python -m ecmapper.draft_map.gui"
    parser.add_argument("inputs", nargs="*", metavar="FILE",
                        help="MODEL [OUT_DIR]: a model (SBML, or JSON, YAML, MAT by the extension; also compressed: "
                        ".gz, .bz2, .xz, .zip) that is matched to the BRENDA map, or the PREP.pkl of ecmapper.draft_map.prepare (OUT_DIR default: "
                        "MODEL_maps next to the model); or, as ecmapper.draft_map.run, PREP.pkl "
                        "BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR; none: choose in a dialog")
    parser.add_argument("--brenda-json", help="the reactions export of the BRENDA map used as template (default: "
                        "BRENDAmap/BRENDA_map_reactions.json of ecmapper)")
    parser.add_argument("--brenda-svg", help="the SVG of the BRENDA map (default: BRENDAmap/BRENDA_map.svg)")
    parser.add_argument("--subsystems", metavar="JSON", help="with a model: the model in JSON format with the "
                        "subsystems of the reactions (BiGG's SBML files have none)")
    parser.add_argument("--vicinity", default="propagate", choices=["none", "neighbours", "propagate"],
                        help="with a model: search the vicinity of found nodes (see ecmapper.draft_map.prepare)")
    return parser


def set_inputs(args, source=None, out=None, brenda_json=None, brenda_svg=None) -> None:
    """The input files: `source` a model or a PREP.pkl, the BRENDA map (default: :func:`template`), OUT_DIR.

    Sets `model` (None for a PREP.pkl) and `prep` (None for a model until it is prepared) of `args`.
    """
    source = Path(source).resolve() if source else None
    is_prep = source is not None and source.suffix.lower() == ".pkl"
    args.prep, args.model = (source, None) if is_prep else (None, source)
    default_json, default_svg = template()
    args.brenda_json = Path(brenda_json or args.brenda_json or default_json).resolve()
    args.brenda_svg = Path(brenda_svg or args.brenda_svg or default_svg).resolve()
    if out:
        args.out = Path(out).resolve()
    elif source is not None:
        args.out = source.parent / f"{source.name.split('.')[0]}_maps"
    else:
        args.out = None
    for name in ("assignment", "maps", "subsystems"):
        if getattr(args, name):
            setattr(args, name, Path(getattr(args, name)).resolve())


def parse_args(argv=None):
    """The inputs (see :func:`build_parser`), the options of :mod:`..run` (--split implied) and the layouts of step 3b.

    Exits on errors. Without inputs, `model`, `prep` and `out` are None (the GUI asks for them).
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    inputs = args.inputs
    if len(inputs) == 4:
        if args.brenda_json or args.brenda_svg:
            parser.error("the BRENDA map is given twice (positional and --brenda-json/--brenda-svg)")
        set_inputs(args, inputs[0], inputs[3], inputs[1], inputs[2])
    elif len(inputs) <= 2:
        set_inputs(args, *inputs)
    else:
        parser.error(USAGE)
    del args.inputs
    args.split = True
    methods = run.check_args(parser, args)
    args.write_assignment = None
    return args, methods


def _free(moving: np.ndarray, others: np.ndarray, gap: float) -> bool:
    if not len(moving) or not len(others):
        return True
    m, o = moving[:, None, :], others[None, :, :]
    g = gap - TOLERANCE
    hit = ((m[..., 0] < o[..., 2] + g) & (o[..., 0] < m[..., 2] + g)
           & (m[..., 1] < o[..., 3] + g) & (o[..., 1] < m[..., 3] + g))
    return not hit.any()


def boxes_free(moving, others, gap: float) -> bool:
    """Whether the boxes `moving` (rows x0, y0, x1, y1) are at least `gap` apart from the boxes `others`."""
    return _free(np.asarray(moving, float).reshape(-1, 4), np.asarray(others, float).reshape(-1, 4), gap)


def constrained_delta(moving, others, dx: float, dy: float, gap: float) -> tuple[float, float]:
    """The part of the move (dx, dy) of the boxes `moving` that keeps them at least `gap` from `others`.

    The whole move if it ends at a free place; else the boxes move as far as
    they can (up to the gap from an obstacle), then the rest of the move
    along the axis on which they can (they slide along the obstacle).
    `moving` must be free at the start.
    """
    moving = np.asarray(moving, float).reshape(-1, 4)
    others = np.asarray(others, float).reshape(-1, 4)

    def free(ddx, ddy):
        return _free(moving + np.array([ddx, ddy, ddx, ddy]), others, gap)

    def longest(x0, y0, mx, my):
        """The free end of the move (mx, my) from the free place (x0, y0) (bisection)."""
        lo, hi = 0.0, 1.0
        for _ in range(24):
            t = (lo + hi) / 2
            if free(x0 + t * mx, y0 + t * my):
                lo = t
            else:
                hi = t
        return x0 + lo * mx, y0 + lo * my

    if free(dx, dy):
        return dx, dy
    bx, by = longest(0.0, 0.0, dx, dy)
    rx, ry = dx - bx, dy - by
    ends = [longest(bx, by, rx, 0.0), longest(bx, by, 0.0, ry)]
    return max(ends, key=lambda e: abs(e[0] - bx) + abs(e[1] - by))


def free_position(others, w: float, h: float, gap: float, near=None) -> tuple[float, float]:
    """Top left corner of a place for a w x h box at least `gap` from `others` (rows x0, y0, x1, y1).

    Near the centre `near` if given (spiral search, :func:`..core.place_free`), else (or if that finds nothing)
    right of `others`.
    """
    others = np.asarray(others, float).reshape(-1, 4)
    if near is not None:
        cx, cy = place_free(near[0], near[1], w, h, [tuple(b) for b in others], margin=gap,
                            step=max(gap, min(w, h) / 4))
        x, y = cx - w / 2, cy - h / 2
        if boxes_free([(x, y, x + w, y + h)], others, gap):
            return float(x), float(y)
    if not len(others):
        return 0.0, 0.0
    return float(others[:, 2].max() + gap), float(others[:, 1].min())


def shelf_positions(sizes, gap: float, aspect: float = 1.6) -> list[tuple[float, float]]:
    """Top left corners of boxes (w, h) in rows, `gap` apart, about `aspect` times as wide as high."""
    if not sizes:
        return []
    total = sum((w + gap) * (h + gap) for w, h in sizes)
    width = max(max(w for w, _ in sizes), math.sqrt(total * aspect))
    x = y = row = 0.0
    result = []
    for w, h in sizes:
        if x > 0 and x + w > width:
            x, y, row = 0.0, y + row + gap, 0.0
        result.append((x, y))
        x += w + gap
        row = max(row, h)
    return result


class Session:
    """The maps of the draft map of `args` (options of :mod:`..run`), the subsystems on them and their positions.

    `positions[map]` ({subsystem ID: top left (x, y)}) says which subsystems
    are on a map; a subsystem may be on several maps or none. `arranged[map]`:
    the map is arranged by hand (its positions are saved and used by the
    pipeline); the others are annealed by the pipeline. `layouts[map]`: the
    AI layout of a map (:mod:`..ai_layouts`), used by the pipeline while it
    fits the map; `export_svg(escher map) -> SVG or None`: set by the GUI,
    to export layouts whose "→ map" labels changed when saving. `project`:
    the project (:mod:`..project`) whose work folder is OUT_DIR, if any.
    """

    def __init__(self, args, methods: list[str], log=print, mapper=None, project: Project | None = None):
        self.args, self.methods, self.log = args, methods, log
        self.mapper = mapper  # the ECMapper for matching a model (default: all sources)
        self.project = project
        self.path: Path | None = args.assignment
        self.md = self.dm = self.table = None
        self.maps: list[str] = []
        self.positions: dict[str, dict[str, tuple[float, float]]] = {}
        self.arranged: dict[str, bool] = {}
        self.modified = False
        self.counts: dict[str, int] = {}
        self.reactions_of: dict[str, set[str]] = {}
        self.hubs: list[dict] = []
        self.version = 0  # counts the changes of the assignment (the drawings of the subsystems may change)
        self._views = None
        self._svgs: dict = {}
        self.layouts: dict[str, MapLayout] = {}
        self._layout_info: dict = {}
        self.export_svg = None

    @property
    def gap(self) -> float:
        return self.args.gap

    # --- loading and saving ------------------------------------------------------------------

    def cache_key(self) -> str:
        """What the drawings of the subsystems depend on (besides the reaction sheet): inputs, options, code.

        By the contents of the files, so the cache of a project stays valid when it is unpacked elsewhere.
        """
        files = digests([self.args.prep, self.args.brenda_json, self.args.brenda_svg,
                         *sorted(Path(__file__).parents[1].glob("*.py"))])
        options = {o: getattr(self.args, o) for o in DRAWING_OPTIONS}
        return hashlib.sha256(repr((files, options, self.methods)).encode()).hexdigest()

    def prepare_model(self) -> Path:
        """Match the model (`args.model`) to the BRENDA map (:mod:`..prepare`): OUT/MODEL_prep.pkl (in a project:
        OUT/prep.pkl).

        Kept for the next start with the same model, BRENDA map, subsystems
        file, vicinity option and ecmapper version (by contents, OUT/MODEL_prep.key);
        :meth:`rematch` matches again.
        """
        from ...brenda_map import BrendaMap, MapMatcher
        from ...mapper import ECMapper
        from ..prepare import prepare, read_model, read_subsystems

        args = self.args
        model_path = Path(args.model)
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        prep = out / (PREP if self.project is not None else f"{model_path.name.split('.')[0]}_prep.pkl")
        key_file = prep.with_suffix(".key")
        key = hashlib.sha256(repr((digests([model_path, args.brenda_json, args.subsystems]), args.vicinity,
                                   ecmapper_version())).encode()).hexdigest()
        if prep.is_file() and key_file.is_file() and key_file.read_text().strip() == key:
            self.log(f"model matched to the BRENDA map before: {prep}")
            return prep
        self.log(f"reading the model {model_path}")
        model = read_model(model_path)
        self.log(f"matching the {len(model.reactions)} reactions of {model.id} to the BRENDA map (ecmapper; the "
                 "resources are downloaded on first use) …")
        matcher = MapMatcher(BrendaMap.from_json(args.brenda_json), self.mapper or ECMapper(),
                             vicinity=None if args.vicinity == "none" else args.vicinity)
        data = prepare(model, matcher, read_subsystems(args.subsystems) if args.subsystems else None)
        with open(prep, "wb") as f:
            pickle.dump(data, f)
        key_file.write_text(key)
        on_map = sum(m["status"] == "on map" for m in data["matches"].values())
        self.log(f"{len(data['matches'])} reactions, {on_map} on the BRENDA map; written to {prep}")
        return prep

    def rematch(self) -> None:
        """Match the model to the BRENDA map again at the next :meth:`load` (e.g. after the resources were
        updated or the matching changed)."""
        if getattr(self.args, "model", None) and self.args.prep:
            for f in (Path(self.args.prep), Path(self.args.prep).with_suffix(".key")):
                f.unlink(missing_ok=True)

    @staticmethod
    def reusable(cached: dict, table) -> bool:
        """Whether the reaction sheet of `table` would not change the cached drawings (as the table saved here)."""
        if table is None:
            return not cached["changed"]
        dm, md = cached["dm"], cached["md"]
        current = reaction_subsystems(dm)
        for rid, wanted in table.reactions.items():
            if rid not in md.reactions:
                continue
            if any(s in cached["dropped"] for s in wanted):
                return False  # a subsystem left without reactions, which would be drawn again
            if {s for s in wanted if s in dm.subsystems} != set(current.get(rid, [])):
                return False
        return True

    def load(self, use_cache: bool = True) -> None:
        """Steps 1-4 of the pipeline (with the table of --assignment), the maps and positions (from the table).

        The condensed subsystems are kept in OUT/:data:`CACHE` for the next
        start with the same inputs. Maps with the positions of all their
        subsystems in the table are arranged (boxes closer than the gap moved
        apart), the others get the order of the BRENDA map.
        """
        for name in ("brenda_json", "brenda_svg"):
            path = Path(getattr(self.args, name))
            if not path.is_file():
                raise FileNotFoundError(f"the BRENDA map {path} is missing (--{name.replace('_', '-')})")
        if getattr(self.args, "model", None):
            self.args.prep = self.prepare_model()
        table = read_assignment(self.args.assignment) if self.args.assignment else None
        key = self.cache_key()
        cache = Path(self.args.out) / CACHE
        cached = None
        if use_cache and cache.is_file():
            try:
                with open(cache, "rb") as f:
                    cached = pickle.load(f)
                if cached.get("key") != key or not self.reusable(cached, table):
                    cached = None
            except Exception as e:  # noqa: BLE001 (an old or damaged cache is made again)
                self.log(f"cache {cache} not used: {e}")
                cached = None
        if cached is not None:
            self.md, dm = cached["md"], cached["dm"]
            self.log(f"condensed subsystems from {cache}")
        else:
            stats = {}
            self.md, dm, _ = run.draft_until_assignment(self.args, self.log, stats, drop_unassigned=False)
            run.relayout_and_condense(dm, self.md, self.args, self.methods, self.log, stats)
            applied = stats.get("assignment", {})
            try:
                cache.parent.mkdir(parents=True, exist_ok=True)
                with open(cache, "wb") as f:
                    pickle.dump(dict(key=key, md=self.md, dm=dm, changed=bool(applied.get("added") or
                                                                             applied.get("removed")),
                                     dropped=set(applied.get("empty_subsystems", []))), f)
            except OSError as e:
                self.log(f"cache {cache} not written: {e}")
        self.dm, self.table = dm, table
        self.hubs = hub_table(dm)
        self.reactions_of = defaultdict(set)
        statuses = defaultdict(set)
        for n in dm.nodes.values():
            self.reactions_of[n.subsystem].update(n.reactions)
            for r in n.reactions:
                statuses[n.subsystem, r].add(n.status)
        self.counts = {k: len(self.reactions_of[k]) for k in dm.subsystems}
        # reactions added from the model (step 3), not on the BRENDA map
        self.added_counts = {k: sum(origin(statuses[k, r]) == "model" for r in self.reactions_of[k])
                             for k in dm.subsystems}
        order, maps_of = run.map_assignment(dm, self.args, self.table, self.log)
        self.maps = list(dict.fromkeys(self.table.maps)) if self.table is not None else list(order)
        given = self.table.positions if self.table is not None else {}
        self.positions, self.arranged = {}, {}
        for name in self.maps:
            keys = [k for k in dm.subsystems if name in maps_of.get(k, [])]
            pos = given.get(name, {})
            if keys and all(k in pos for k in keys):
                w, h = (np.array(v) for v in zip(*(self.size(k) for k in keys)))
                gx, gy = (np.array(v) for v in zip(*(pos[k] for k in keys)))
                x, y = user_arrangement(w, h, gx, gy, gap=self.gap)
                moved = int(((np.abs(x - gx) > 0.5) | (np.abs(y - gy) > 0.5)).sum())
                if moved:
                    self.log(f"map {name!r}: {moved} boxes of the table closer than the gap moved apart")
                self.positions[name] = {k: (float(a), float(b)) for k, a, b in zip(keys, x, y)}
                self.arranged[name] = True
            else:
                if keys and pos:
                    self.log(f"map {name!r}: the arrangement of the table lacks subsystems, original order")
                self.positions[name] = self.original_order(keys)
                self.arranged[name] = False
        self.layouts = load_layouts(table, Path(self.args.assignment).parent, self.log) if table is not None else {}
        self.modified = False
        self._changed()
        self.log(f"{len(dm.subsystems)} subsystems, {len(self.maps)} maps, {len(self.unassigned())} subsystems "
                 f"without map")
        for name in list(self.layouts):
            info = self.layout_info(name)
            if info is None:
                self.log(f"the AI layout of map {name!r} is left out: there is no such map (with subsystems)")
                del self.layouts[name]
            elif info[2]:
                self.log(f"the AI layout of map {name!r} is outdated: {'; '.join(info[2])}")

    def save(self, path=None) -> Path:
        """The assignment as XLSX table (:func:`..assignment.write_assignment`) with the arranged maps.

        The AI layouts go into folders next to it (TABLE_layouts/MAP); their
        "→ map" labels are updated first (and exported again with
        `export_svg`).
        """
        path = Path(path or self.path)
        assignment = {k: self.maps_of(k) for k in self.dm.subsystems}
        positions = {m: dict(self.positions[m]) for m in self.maps if self.arranged.get(m) and self.positions[m]}
        layouts = {}
        for name, layout in list(self.layouts.items()):
            if name not in self.maps:
                continue
            if name in self.views():
                updated = copy.deepcopy(layout)
                if refresh_exit_labels(updated, self.map_draft(name)):
                    if self.export_svg is not None:
                        updated.svg = self.export_svg(updated.escher)
                    self.layouts[name] = layout = updated
                    self.log(f"the '→ map' labels of the AI layout of {name!r} updated")
            folder = path.with_name(f"{path.stem}_layouts") / slug(name)
            layout.save(folder)
            source = ", ".join(str(layout.source[k]) for k in ("toolkit", "model", "date") if layout.source.get(k))
            layouts[name] = (os.path.relpath(folder, path.parent).replace(os.sep, "/"), source)
        write_assignment(path, self.dm, self.md, self.maps, assignment, positions, layouts)
        self.path = self.args.assignment = path
        self.modified = False
        return path

    def save_project(self, path=None) -> Path:
        """Save the project: the table (:meth:`save`) into its work folder, then the project file `path` (default:
        the one it was opened from or saved as)."""
        table = self.project.work / TABLE
        self.save(table)
        folder = table.with_name(f"{table.stem}_layouts")
        kept = {slug(name) for name in self.layouts if name in self.maps}
        if folder.is_dir():  # the layouts of maps without layout now
            for sub in folder.iterdir():
                if sub.is_dir() and sub.name not in kept:
                    shutil.rmtree(sub, ignore_errors=True)
        return self.project.save(path)

    def build_command(self, out=None, iterations: int | None = None, seeds: int | None = None) -> list[str]:
        """The command line of :mod:`..run` that builds the maps of the saved table (`path`)."""
        args = copy.copy(self.args)
        args.assignment, args.write_assignment, args.central, args.maps = self.path, None, None, None
        args.split = True
        if out is not None:
            args.out = Path(out)
        if iterations is not None:
            args.iterations = iterations
        if seeds is not None:
            args.seeds = seeds
        return [sys.executable, "-m", "ecmapper.draft_map.run", *run.command_line(args)]

    # --- queries -------------------------------------------------------------------------------

    def size(self, k: str) -> tuple[float, float]:
        s = self.dm.subsystems[k]
        return s.w, s.h

    def rect(self, name: str, k: str) -> tuple[float, float, float, float]:
        (x, y), (w, h) = self.positions[name][k], self.size(k)
        return x, y, x + w, y + h

    def rects(self, name: str, exclude=()) -> np.ndarray:
        exclude = set(exclude)
        return np.array([self.rect(name, k) for k in self.positions.get(name, {}) if k not in exclude],
                        float).reshape(-1, 4)

    def maps_of(self, k: str) -> list[str]:
        return [m for m in self.maps if k in self.positions[m]]

    def unassigned(self) -> list[str]:
        """The subsystems on no map, by label."""
        on_map = {k for p in self.positions.values() for k in p}
        return sorted((k for k in self.dm.subsystems if k not in on_map),
                      key=lambda k: self.dm.subsystems[k].label.lower())

    def unassigned_positions(self) -> dict[str, tuple[float, float]]:
        keys = self.unassigned()
        return dict(zip(keys, shelf_positions([self.size(k) for k in keys], self.gap)))

    def find(self, text: str) -> list[str]:
        """Subsystems whose ID or label contains `text`, or with a reaction with this ID (case ignored)."""
        t = text.strip().casefold()
        if not t:
            return []
        return [k for k, s in self.dm.subsystems.items()
                if t in k.casefold() or t in s.label.casefold() or any(r.casefold() == t for r in self.reactions_of[k])]

    def reactions(self, k: str) -> list[tuple[str, str, str, str, str]]:
        """(reaction ID, name, drawn as, equation with names, origin) of the model reactions drawn in subsystem `k`.

        Origin (:func:`..assignment.origin`): "model" for reactions added from the model, "both" for reactions also
        added (in another node of `k`), "map" for those of the BRENDA map.
        """
        status = defaultdict(set)
        for n in self.dm.members(k):
            for r in n.reactions:
                status[r].add(n.status)
        return [(r, self.md.reactions[r]["name"], ", ".join(sorted(st)), equation(self.md, r, names=True), origin(st))
                for r, st in sorted(status.items())]

    def original_order(self, keys) -> dict[str, tuple[float, float]]:
        """Positions of subsystems in the order of the BRENDA map, compacted (:func:`..layout.arrangement_original`)."""
        keys = list(keys)
        if not keys:
            return {}
        sub = DraftMap()
        sub.subsystems = {k: self.dm.subsystems[k] for k in keys}
        sub.orig_centres = {k: self.dm.orig_centres[k] for k in keys}
        x, y = arrangement_original(sub, gap=self.gap)
        return {k: (float(a), float(b)) for k, a, b in zip(keys, x, y)}

    def paste_positions(self, name: str, keys, layout=None, near=None) -> dict[str, tuple[float, float]]:
        """Positions for subsystems added to map `name` as a group at a free place (:func:`free_position`).

        `layout`: their positions relative to each other (e.g. on the map they
        come from), default: in rows. Subsystems already on the map are left
        out.
        """
        keys = [k for k in dict.fromkeys(keys) if k not in self.positions[name]]
        if not keys:
            return {}
        if layout is not None and all(k in layout for k in keys):
            rel = {k: layout[k] for k in keys}
        else:
            rel = dict(zip(keys, shelf_positions([self.size(k) for k in keys], self.gap)))
        x0 = min(rel[k][0] for k in keys)
        y0 = min(rel[k][1] for k in keys)
        w = max(rel[k][0] + self.size(k)[0] for k in keys) - x0
        h = max(rel[k][1] + self.size(k)[1] for k in keys) - y0
        x, y = free_position(self.rects(name), w, h, self.gap, near)
        return {k: (x + rel[k][0] - x0, y + rel[k][1] - y0) for k in keys}

    # --- the maps as the pipeline draws them ----------------------------------------------------

    def views(self) -> dict[str, tuple[DraftMap, dict[str, str]]]:
        """Each map with subsystems as the pipeline draws it: {map: (draft map, {subsystem ID: ID on the map})}.

        Subsystems on several maps are copies (:func:`..split.duplicate_subsystems`), links between maps are
        references (:func:`..split.split_map`), the links inside a map are as with the options
        (:func:`..run.map_links`). Computed again after changes of the assignment.
        """
        if self._views is None:
            maps_of = {k: self.maps_of(k) for k in self.dm.subsystems}
            assigned = {k for k, maps in maps_of.items() if maps}
            dm = DraftMap()
            dm.subsystems = {k: s for k, s in self.dm.subsystems.items() if k in assigned}
            dm.nodes = {i: n for i, n in self.dm.nodes.items() if n.subsystem in assigned}
            dm.edges = [e for e in self.dm.edges if e.enzyme in dm.nodes and e.met in dm.nodes]
            dm.links = [l for l in self.dm.links if l.a in dm.nodes and l.b in dm.nodes]
            dm.orig_centres = dict(self.dm.orig_centres)
            assignment, copies = duplicate_subsystems(dm, {k: maps_of[k] for k in assigned})
            views = {}
            for name, sub in split_map(dm, [m for m in self.maps if self.positions[m]], assignment):
                run.map_links(sub, self.args, self.hubs)
                views[name] = (sub, {k: copies.get((k, name), k) for k in self.positions[name]})
            self._views = views
        return self._views

    def links(self, name: str) -> list[tuple[str, tuple[float, float], str, tuple[float, float]]]:
        """The inter-pathway links of a map: (subsystem, end relative to its box, other subsystem, other end)."""
        if name not in self.views():
            return []
        sub, ids = self.views()[name]
        back = {v: k for k, v in ids.items()}
        result = []
        for l in sub.links:
            a, b = sub.nodes[l.a], sub.nodes[l.b]
            if a.subsystem != b.subsystem:
                result.append((back[a.subsystem], (a.x, a.y), back[b.subsystem], (b.x, b.y)))
        return result

    def subsystem_svg(self, name: str | None, k: str) -> tuple[bytes, tuple[float, float, float, float]]:
        """The SVG of subsystem `k` as drawn on map `name` (with references to other maps) and its view box.

        The view box is in the coordinates of the box (top left at 0, 0). `name`
        None or a map without `k`: the subsystem alone.
        """
        key = (name, k)
        if key not in self._svgs:
            views = self.views()
            if name in views and k in views[name][1]:
                sub, ids = views[name]
                one = subsystem_map(sub, ids[k])
            else:
                one = subsystem_map(self.dm, k)
            text = svg_text(one, md=self.md, margin=SVG_MARGIN)
            self._svgs[key] = (text.encode("utf-8"), one.view_box)
        return self._svgs[key]

    def pack(self, name: str, iterations: int, seed: int = 1, seeds: int = 1, progress=None,
             stop=None) -> dict[str, tuple[float, float]]:
        """Positions of the subsystems of a map arranged by simulated annealing, starting from the current ones.

        As the pipeline arranges the map (:func:`..layout.anneal_arrangement`,
        with the links of the map, without flips); the best of `seeds` runs.
        `progress(fraction)` and `stop()`: see there. Does not change the
        session; call :meth:`views` first when used in another thread.
        """
        keys = list(self.positions[name])
        if len(keys) < 2:
            return dict(self.positions[name])
        sub, ids = self.views()[name]
        back = {v: k for k, v in ids.items()}
        ports = Ports.from_map(sub)
        x = np.array([self.positions[name][back[c]][0] for c in ports.names], float)
        y = np.array([self.positions[name][back[c]][1] for c in ports.names], float)
        best = None
        for i, s in enumerate(range(seed, seed + seeds)):
            report = None if progress is None else (lambda f, _, i=i: progress((i + f) / seeds))
            res = anneal_arrangement(ports, x, y, gap=self.gap, iterations=iterations, seed=s, progress=report,
                                     stop=stop)
            if best is None or res[4]["cost"] < best[4]["cost"]:
                best = res
            if stop is not None and stop():
                break
        return {back[c]: (float(a), float(b)) for c, a, b in zip(ports.names, best[0], best[1])}

    # --- AI layouts (..ai_layouts) --------------------------------------------------------------

    def map_draft(self, name: str) -> DraftMap:
        """Map `name` as the pipeline draws it with the current arrangement (absolute positions)."""
        sub, ids = self.views()[name]
        m = copy.deepcopy(sub)
        back = {v: k for k, v in ids.items()}
        for c, s in m.subsystems.items():
            s.x, s.y = self.positions[name][back[c]]
        to_absolute(m)
        return m

    def export_map(self, name: str, folder, file: str = DRAFT, reuse: dict | None = None) -> tuple[Path, dict]:
        """Write map `name` as Escher map to `folder`/`file` for an AI job; the draft objects of its Escher IDs
        (`reuse`: those of an earlier export, kept)."""
        data, ids = export_draft(self.map_draft(name), self.md, name, scale=self.args.escher_scale, reuse=reuse)
        path = Path(folder) / file
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data), encoding="utf-8")
        return path, ids

    def new_layout(self, name: str, final: list, ids: dict | None, svg: str | None = None, run_log: dict | None = None,
                   source: dict | None = None) -> tuple[MapLayout, list[str]]:
        """The layout of map `name` laid out as the Escher map `final` (exported with `ids`; None: made elsewhere,
        matched by content, :func:`..ai_layouts.match_ids`) and why it does not fit the map (empty: it fits)."""
        m = self.map_draft(name)
        scale = self.args.escher_scale
        if ids is None:
            draft, draft_ids = export_draft(m, self.md, name, scale=scale)
            ids = match_ids(draft, draft_ids, final)
        layout = make_layout(m, self.md, name, final, ids, scale, svg, run_log, source)
        return layout, check_layout(layout, m, self.md)

    def set_layout(self, name: str, layout: MapLayout | None) -> None:
        """Lay out map `name` so (None: draw it automatically again)."""
        if layout is None:
            if self.layouts.pop(name, None) is None:
                return
        else:
            self.layouts[name] = layout
        self.modified = True
        self._layout_info = {}
        self.version += 1

    def layout_info(self, name: str | None):
        """The layout of map `name`, the boxes of its subsystems in it ({subsystem: (x, y, w, h)}) and why it is
        outdated (empty: it fits the map); None without layout."""
        layout = self.layouts.get(name)
        if layout is None or name not in self.views():
            return None
        key = (name, self.version, id(layout))
        if key not in self._layout_info:
            sub, ids = self.views()[name]
            m = copy.deepcopy(sub)
            problems = check_layout(layout, m, self.md)
            apply_layout(m, layout)
            back = {v: k for k, v in ids.items()}
            boxes = {back[c]: (s.x, s.y, s.w, s.h) for c, s in m.subsystems.items() if c in back}
            self._layout_info[key] = (layout, boxes, problems)
        return self._layout_info[key]

    def outdated(self, name: str | None) -> bool:
        info = self.layout_info(name)
        return info is not None and bool(info[2])

    # --- changes -------------------------------------------------------------------------------

    def _changed(self) -> None:
        """The assignment changed: the maps are drawn anew."""
        self._views = None
        self._svgs = {}
        self._layout_info = {}
        self.version += 1

    def state(self):
        return (list(self.maps), {m: dict(p) for m, p in self.positions.items()}, dict(self.arranged),
                dict(self.layouts))

    def restore(self, state) -> None:
        maps, positions, arranged, layouts = state
        members = {m: set(p) for m, p in positions.items()}
        assignment_changed = maps != self.maps or members != {m: set(p) for m, p in self.positions.items()}
        layouts_changed = layouts != self.layouts
        self.maps, self.arranged = list(maps), dict(arranged)
        self.positions = {m: dict(p) for m, p in positions.items()}
        self.layouts = dict(layouts)
        self.modified = True
        if assignment_changed:
            self._changed()
        elif layouts_changed:
            self._layout_info = {}
            self.version += 1

    def add(self, name: str, positions: dict[str, tuple[float, float]]) -> None:
        new = {k: p for k, p in positions.items() if k not in self.positions[name]}
        if new:
            self.positions[name].update(new)
            self.modified = True
            self._changed()

    def remove(self, name: str, keys) -> None:
        keys = [k for k in keys if k in self.positions[name]]
        if keys:
            for k in keys:
                del self.positions[name][k]
            self.modified = True
            self._changed()

    def move(self, name: str, positions: dict[str, tuple[float, float]]) -> None:
        """New positions of subsystems on a map (all of them: an arrangement); the map counts as arranged."""
        self.positions[name].update({k: (float(x), float(y)) for k, (x, y) in positions.items()
                                     if k in self.positions[name]})
        self.arranged[name] = True
        self.modified = True

    def set_arranged(self, name: str, arranged: bool) -> None:
        self.arranged[name] = arranged
        self.modified = True

    def check_name(self, name: str, old: str | None = None) -> str | None:
        """Why `name` can not be the name of a map (None: it can)."""
        name = name.strip()
        if not name:
            return "The name is empty."
        if name == NOT_DRAWN or (name in self.maps and name != old):
            return f"There is a map {name!r} already."
        return None

    def new_map(self, name: str) -> None:
        self.maps.append(name)
        self.positions[name], self.arranged[name] = {}, False
        self.modified = True

    def rename_map(self, old: str, new: str) -> None:
        if old == new:
            return
        self.maps[self.maps.index(old)] = new
        self.positions[new] = self.positions.pop(old)
        self.arranged[new] = self.arranged.pop(old, False)
        if old in self.layouts:
            self.layouts[new] = self.layouts.pop(old)
        self.modified = True
        self._changed()

    def delete_map(self, name: str) -> None:
        self.maps.remove(name)
        del self.positions[name]
        self.arranged.pop(name, None)
        self.layouts.pop(name, None)
        self.modified = True
        self._changed()

    def reorder(self, maps: list[str]) -> None:
        if maps != self.maps:
            assert sorted(maps) == sorted(self.maps)
            self.maps = list(maps)
            self.modified = True
            self._changed()
