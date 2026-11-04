"""The assignment GUI: the geometry, the session (no Qt) and, with qtpy and a Qt binding, the window (offscreen)."""

import gzip
import json
import pickle
import random

import pytest

pytest.importorskip("scipy")
pytest.importorskip("numba")
pytest.importorskip("openpyxl")

from ecmapper.draft_map import run, split
from ecmapper.draft_map.assignment import read_assignment
from ecmapper.draft_map.gui.session import (
    Session,
    boxes_free,
    constrained_delta,
    free_position,
    parse_args,
    shelf_positions,
)
from test_brenda_map import model  # noqa: F401 (fixture)
from test_draft_map import BRENDA_MAP, no_overlaps, prepared


def test_constrained_delta():
    """A dragged box stops at the gap from an obstacle and slides along it."""
    a, b = [(0, 0, 100, 100)], [(200, 0, 300, 100)]
    assert constrained_delta(a, b, 0, 300, 50) == (0, 300)
    dx, dy = constrained_delta(a, b, 200, 0, 50)
    assert dx == pytest.approx(50, abs=0.01) and dy == 0
    dx, dy = constrained_delta(a, b, 200, 30, 50)
    assert dx == pytest.approx(50, abs=0.01) and dy == pytest.approx(30, abs=0.01)
    # a group: the second box hits the obstacle first
    dx, dy = constrained_delta([(0, 0, 50, 50), (0, 150, 120, 200)], [(250, 120, 300, 220)], 200, 0, 50)
    assert dx == pytest.approx(80, abs=0.01)
    # random moves of a box among others: never closer than the gap
    rng = random.Random(1)
    others = [(x, y, x + 80, y + 60) for x in range(0, 1000, 250) for y in range(0, 1000, 250)]
    box = (150, 150, 200, 180)
    for _ in range(300):
        dx, dy = constrained_delta([box], others, rng.uniform(-80, 80), rng.uniform(-80, 80), 40)
        box = (box[0] + dx, box[1] + dy, box[2] + dx, box[3] + dy)
        assert boxes_free([box], others, 40)


def test_free_position():
    others = [(0, 0, 100, 100)]
    x, y = free_position(others, 50, 50, 20, near=(50, 50))
    assert boxes_free([(x, y, x + 50, y + 50)], others, 20)
    assert free_position(others, 50, 50, 20) == (120.0, 0.0)
    assert free_position([], 50, 50, 20) == (0.0, 0.0)
    sizes = [(100, 50), (60, 80), (200, 30), (50, 50)]
    positions = shelf_positions(sizes, 10)
    assert no_overlaps([p[0] for p in positions], [p[1] for p in positions], [s[0] for s in sizes],
                       [s[1] for s in sizes], 10)


@pytest.fixture
def session(sources, model, tmp_path):  # noqa: F811
    _, prep = prepared(sources, model, tmp_path)
    args, methods = parse_args([str(prep), str(BRENDA_MAP / "BRENDA_map_reactions.json"),
                                str(BRENDA_MAP / "BRENDA_map.svg"), str(tmp_path / "out"), "--pack",
                                "--min-map-reactions", "0", "--iterations", "300"])
    s = Session(args, methods, log=lambda message: None)
    s.load()
    return s


def apart(s: Session) -> bool:
    return all(boxes_free([s.rect(m, k)], s.rects(m, exclude=[k]), s.gap) for m in s.maps for k in s.positions[m])


def test_session(session, tmp_path):
    """Copy a subsystem to another map, pack a map, save, build the maps as the GUI does, start again."""
    s = session
    assert len(s.maps) > 1 and apart(s) and not any(s.arranged.values()) and not s.unassigned()
    # reactions added from the model, not on the BRENDA map
    assert sum(s.added_counts.values()) > 0
    for key in s.dm.subsystems:
        rows = s.reactions(key)
        assert s.added_counts[key] == sum(row[4] == "model" for row in rows)
        assert all(row[4] in ("model", "both", "map") and (row[4] == "model") == (row[2] == "added") for row in rows)
    central, other = s.maps[0], s.maps[-1]
    k = next(iter(s.positions[central]))
    positions = s.paste_positions(other, [k], {k: s.positions[central][k]}, near=(0.0, 0.0))
    s.add(other, positions)
    assert s.maps_of(k) == [central, other] and apart(s)
    assert s.paste_positions(other, [k]) == {}
    data, box = s.subsystem_svg(other, k)
    assert data.startswith(b"<svg") and box[2] >= s.size(k)[0]
    assert set(s.views()) == set(s.maps)
    # undo, as the window does it
    state = s.state()
    s.remove(other, [k])
    assert s.maps_of(k) == [central]
    s.restore(state)
    assert s.maps_of(k) == [central, other]
    # pack the central map by annealing
    s.move(central, s.pack(central, 300))
    assert s.arranged[central] and apart(s)
    # a map without subsystems on the table is a map with them later
    s.new_map("Mine")
    assert s.check_name("Mine") and s.check_name(" ") and s.check_name("Other") is None

    path = s.save(tmp_path / "assignment.xlsx")
    table = read_assignment(path)
    assert table.subsystems[k] == [central, other] and "Mine" in table.maps
    assert set(table.positions) == {central}
    for key, (x, y) in s.positions[central].items():
        assert table.positions[central][key] == pytest.approx((x, y), abs=0.1)

    # build the maps with the command line of the GUI (here in this process)
    built = tmp_path / "built"
    command = s.build_command(built)
    assert command[1:3] == ["-m", "ecmapper.draft_map.run"]
    assert "--assignment" in command and "--write-assignment" not in command and "--iterations" in command
    run.main(command[3:])
    index = json.loads((built / "maps.json").read_text())
    assert [m["name"] for m in index["maps"] if k in m["subsystems"]] == [central, other]
    with open(built / f"{split.slug(central)}.pkl", "rb") as f:
        m = pickle.load(f)
    for key, (x, y) in s.positions[central].items():
        assert (m.subsystems[key].x, m.subsystems[key].y) == pytest.approx((x, y), abs=0.1)

    # the next start takes the drawings from the cache, the maps from the table
    messages = []
    again = Session(s.args, s.methods, log=messages.append)
    again.load()
    assert any("gui_cache" in line for line in messages)
    assert again.maps == s.maps and again.maps_of(k) == [central, other] and again.arranged[central]
    assert {key: again.size(key) for key in again.dm.subsystems} == {key: s.size(key) for key in s.dm.subsystems}
    # a reaction taken out of its subsystems in the table: the drawings are made again
    from openpyxl import load_workbook

    wb = load_workbook(path)
    ws = wb["Reactions"]
    rid = next(row[0].value for row in ws.iter_rows(min_row=2) if row[2].value)
    for row in ws.iter_rows(min_row=2):
        if row[0].value == rid:
            row[2].value = None
    wb.save(path)
    messages.clear()
    again = Session(s.args, s.methods, log=messages.append)
    again.load()
    assert not any("gui_cache" in line for line in messages)
    assert all(rid not in again.reactions_of[key] for key in again.dm.subsystems)


def test_parse_args(tmp_path):
    """MODEL [OUT_DIR] with the BRENDA map of ecmapper as template, or the arguments of run; or nothing."""
    template = (BRENDA_MAP / "BRENDA_map_reactions.json").resolve(), (BRENDA_MAP / "BRENDA_map.svg").resolve()
    model = tmp_path / "iML1515.xml"
    args, _ = parse_args([str(model), "--pack"])
    assert args.model == model.resolve() and args.prep is None and args.split and args.pack
    assert (args.brenda_json, args.brenda_svg) == template
    assert args.out == model.resolve().parent / "iML1515_maps"
    args, _ = parse_args([str(model), str(tmp_path / "out"), "--subsystems", "iML1515.json"])
    assert args.out == (tmp_path / "out").resolve() and args.subsystems.name == "iML1515.json"
    args, _ = parse_args(["p.pkl", "m.json", "m.svg", "out"])
    assert args.model is None and args.prep.name == "p.pkl" and args.brenda_svg.name == "m.svg"
    args, _ = parse_args([str(tmp_path / "p.pkl")])
    assert args.prep is not None and args.model is None
    args, _ = parse_args([])
    assert args.model is None and args.prep is None and args.out is None
    with pytest.raises(SystemExit):
        parse_args(["a", "b", "c"])


def test_session_from_model(sources, model, tmp_path):  # noqa: F811
    """A session started with an SBML file matches it to the BRENDA map first, then only once."""
    import cobra

    from ecmapper import ECMapper

    xml = tmp_path / "m.xml"
    cobra.io.write_sbml_model(model, str(xml))
    sbml = tmp_path / "m.xml.gz"  # a compressed model
    sbml.write_bytes(gzip.compress(xml.read_bytes()))
    xml.unlink()
    args, methods = parse_args([str(sbml), str(tmp_path / "out"), "--pack", "--min-map-reactions", "0",
                                "--brenda-json", str(BRENDA_MAP / "BRENDA_map_reactions.json")])
    messages = []
    s = Session(args, methods, log=messages.append, mapper=ECMapper(sources.values()))
    s.load()
    assert args.prep == tmp_path / "out" / "m_prep.pkl" and args.prep.is_file()
    reference, _ = prepared(sources, model, tmp_path)
    assert {r: m["status"] for r, m in s.md.matches.items()} == {r: m["status"] for r, m in reference["matches"].items()}
    assert s.maps and not any("before" in line for line in messages)
    # the build command passes the prepared model to run
    assert s.build_command()[3] == str(args.prep)
    messages.clear()
    again = Session(args, methods, log=messages.append, mapper=ECMapper(sources.values()))
    again.load()
    assert any("before" in line for line in messages) and any("gui_cache" in line for line in messages)


def test_window(session, monkeypatch, tmp_path):
    """The main window (offscreen): tabs, copy and paste with undo, dragging stops at the gap, both zoom levels."""
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("qtpy")
    from qtpy.QtCore import QPointF, Qt
    from qtpy.QtTest import QTest
    from qtpy.QtWidgets import QApplication

    from ecmapper.draft_map.gui.items import BOXES, DRAWINGS
    from ecmapper.draft_map.gui.window import MainWindow

    app = QApplication.instance() or QApplication([])
    s = session
    w = MainWindow(s)
    w.build_tabs()
    w.resize(1200, 800)
    w.show()
    app.processEvents()
    assert w.tabs.count() == len(s.maps) + 1
    central, other = s.maps[0], s.maps[-1]
    k = next(iter(s.positions[central]))
    w.views[central].map_scene().select([k])
    w.copy()
    w.tabs.setCurrentWidget(w.views[other])
    w.paste()
    assert s.maps_of(k) == [central, other] and apart(s)
    assert k in w.views[other].map_scene().boxes
    w.undo.undo()
    assert s.maps_of(k) == [central] and k not in w.views[other].map_scene().boxes
    w.undo.redo()
    assert s.maps_of(k) == [central, other]

    # drag a box onto another one
    name = next(m for m in s.maps if len(s.positions[m]) >= 2)
    w.tabs.setCurrentWidget(w.views[name])
    view = w.views[name]
    view.zoom_to_fit()
    app.processEvents()
    a, b = list(s.positions[name])[:2]
    ra, rb = s.rect(name, a), s.rect(name, b)
    start = view.mapFromScene(QPointF((ra[0] + ra[2]) / 2, (ra[1] + ra[3]) / 2))
    end = view.mapFromScene(QPointF((rb[0] + rb[2]) / 2, (rb[1] + rb[3]) / 2))
    view.map_scene().clearSelection()
    QTest.mousePress(view.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start)
    for i in range(1, 11):
        QTest.mouseMove(view.viewport(), start + (end - start) * i / 10)
    QTest.mouseRelease(view.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, end)
    app.processEvents()
    assert s.rect(name, a) != ra and s.arranged[name] and apart(s)

    for mode in (DRAWINGS, BOXES):
        w.set_mode(mode)
        assert not w.grab().isNull()
    w.change("New map", lambda: s.new_map("Mine"))
    assert w.tabs.count() == len(s.maps) + 1 and "Mine" in w.views
    w.undo.undo()
    assert "Mine" not in w.views

    # the start dialog: a model proposes the output folder, the BRENDA map of ecmapper is the template
    import copy

    from ecmapper.draft_map.gui.window import StartDialog

    args = copy.copy(s.args)
    dialog = StartDialog(w, args)
    assert dialog.brenda_json.path() == args.brenda_json
    sbml = tmp_path / "e_coli.xml"
    sbml.write_text("<sbml/>")
    dialog.model.edit.setText(str(sbml))
    dialog.model_chosen()
    assert dialog.out.path() == tmp_path / "e_coli_maps" and dialog.problem() is None
    dialog.apply()
    assert args.model == sbml.resolve() and args.prep is None and args.out == (tmp_path / "e_coli_maps").resolve()
    dialog.model.edit.setText("")
    assert dialog.problem()
    w.undo.setClean()
    w.close()
