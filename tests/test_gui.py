"""The assignment GUI: the geometry, the session (no Qt) and, with qtpy and a Qt binding, the window (offscreen)."""

import gzip
import json
import pickle
import random
import zipfile
from pathlib import Path

import pytest

pytest.importorskip("scipy")
pytest.importorskip("numba")
pytest.importorskip("openpyxl")

from ecmapper.draft_map import run, split
from ecmapper.draft_map import project as project_module
from ecmapper.draft_map.assignment import read_assignment
from ecmapper.draft_map.project import Project, options_of, work_id
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


def test_project(session, tmp_path, monkeypatch):
    """A project file: saved from a session, opened again (no new match, the drawings from the zipped cache), its
    work folder reused while unchanged, built with run --project; the template; members outside are refused."""
    import copy

    monkeypatch.setenv("ECMAPPER_CACHE", str(tmp_path / "cache"))
    s = session
    prep = s.args.prep
    new = Project.new(prep, None, s.args.brenda_json, s.args.brenda_svg, options_of(s.args))
    assert (new.work / "prep.pkl").is_file() and new.meta["template"]["bundled"] and new.meta["model"] is None
    s.project = new
    methods = new.apply(s.args)
    assert s.args.out == new.work and s.args.pack and s.args.iterations == 300 and methods == s.methods
    assert s.args.brenda_json == (BRENDA_MAP / "BRENDA_map_reactions.json").resolve()
    s.load()
    central = s.maps[0]
    s.move(central, s.pack(central, 300))
    path = s.save_project(tmp_path / "p")
    assert path == (tmp_path / "p.ecmap").resolve()
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
    assert {"project.json", "prep.pkl", "assignment.xlsx", "gui_cache.pkl"} <= names
    assert ".lock" not in names and not any(n.startswith("template/") for n in names)
    new.close()
    assert not new.work.exists()  # the work folder of the unsaved project

    # opened: unpacked, the drawings from the zipped cache, the arrangement from the table
    project, warnings = Project.open(path)
    assert warnings == [] and project.work.name == work_id(path)
    args = copy.copy(s.args)
    messages = []
    again = Session(args, project.apply(args), log=messages.append, project=project)
    again.load()
    assert any("gui_cache" in line for line in messages) and again.arranged[central]
    for key, (x, y) in s.positions[central].items():
        assert again.positions[central][key] == pytest.approx((x, y), abs=0.1)
    # opened again unchanged: the same work folder; after discarded changes: unpacked again
    (project.work / "marker").write_text("x")
    project.close()
    project, _ = Project.open(path)
    assert (project.work / "marker").is_file()
    project.close(discard=True)
    project, _ = Project.open(path)
    assert not (project.work / "marker").exists()
    project.close()

    # built from the project file, the options of the project unless given
    built = tmp_path / "built"
    run.main(["--project", str(path), str(built), "--seeds", "1"])
    with open(built / f"{split.slug(central)}.pkl", "rb") as f:
        m = pickle.load(f)
    for key, (x, y) in s.positions[central].items():
        assert (m.subsystems[key].x, m.subsystems[key].y) == pytest.approx((x, y), abs=0.1)

    # the BRENDA map of ecmapper is referenced, a change of it is reported; another one is copied into the project
    bundled = tmp_path / "bundled"
    bundled.mkdir()
    files = (bundled / "BRENDA_map_reactions.json", bundled / "BRENDA_map.svg")
    for f in files:
        f.write_text("1")
    monkeypatch.setattr(project_module, "bundled_template", lambda: files)
    ref = Project.new(prep, None, None, None, {})
    assert ref.meta["template"]["bundled"] and ref.template() == files and not (ref.work / "template").exists()
    ref_path = ref.save(tmp_path / "ref.ecmap")
    ref.close()
    files[1].write_text("22")
    ref, warnings = Project.open(ref_path)
    assert len(warnings) == 1 and "changed" in warnings[0]
    ref.close()
    other = Project.new(prep, None, s.args.brenda_json, s.args.brenda_svg, {})
    assert not other.meta["template"]["bundled"]
    assert other.template()[1].read_bytes() == s.args.brenda_svg.read_bytes()
    other.close()

    bad = tmp_path / "bad.ecmap"
    with zipfile.ZipFile(bad, "w") as z:
        z.writestr("project.json", "{}")
        z.writestr("../evil.txt", "x")
    with pytest.raises(ValueError):
        Project.open(bad)
    assert not (tmp_path / "cache" / "projects" / "evil.txt").exists()


def test_window(session, monkeypatch, tmp_path):
    """The main window (offscreen): tabs, copy and paste with undo, dragging stops at the gap, both zoom levels; the
    new project dialog and the recent projects."""
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    monkeypatch.setenv("ECMAPPER_CACHE", str(tmp_path / "cache"))
    pytest.importorskip("qtpy")
    from qtpy.QtCore import QPointF, QSettings, Qt
    from qtpy.QtTest import QTest
    from qtpy.QtWidgets import QApplication, QMessageBox

    from ecmapper.draft_map.gui.items import BOXES, DRAWINGS
    from ecmapper.draft_map.gui.window import MainWindow

    app = QApplication.instance() or QApplication(["ecmapper-test"])
    s = session
    w = MainWindow(s, QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    assert w.stack.currentWidget() is w.welcome
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

    assert w.stack.currentWidget() is w.tabs

    # the new project dialog: the BRENDA map of ecmapper is the template, the model is copied into the project
    import copy

    from ecmapper.draft_map.gui.session import template
    from ecmapper.draft_map.gui.window import StartDialog

    args = copy.copy(s.args)
    args.brenda_json = args.brenda_svg = None
    dialog = StartDialog(w, args)
    assert dialog.brenda_json.path() == template()[0]
    sbml = tmp_path / "e_coli.xml"
    sbml.write_text("<sbml/>")
    dialog.model.edit.setText(str(sbml))
    assert dialog.problem() is None
    project = dialog.project()
    assert project.meta["model"] == "model/e_coli.xml" and (project.work / "model" / "e_coli.xml").is_file()
    assert project.meta["template"]["bundled"] and project.meta["options"]["pack"] and project.path is None
    project.close()
    assert not project.work.exists()
    dialog.model.edit.setText("")
    assert dialog.problem()

    # recent projects: the newest first, missing ones dropped
    for name in ("a", "b", "a"):
        w.add_recent(tmp_path / f"{name}.ecmap")
    assert [Path(p).name for p in w.recent_projects()] == ["a.ecmap", "b.ecmap"] and w.welcome.recent.count() == 2
    w.fill_recent_menu()
    assert [a.text() for a in w.recent_menu.actions() if a.text()] == ["&1 a.ecmap", "&2 b.ecmap", "&Clear list"]
    monkeypatch.setattr(QMessageBox, "warning", lambda *args, **kwargs: None)
    assert not w.open_project(tmp_path / "b.ecmap")
    assert [Path(p).name for p in w.recent_projects()] == ["a.ecmap"]
    w.undo.setClean()
    w.close()


def test_ai_dock(session, monkeypatch, tmp_path):
    """The AI layout dock with a scripted AI: it runs the toolkit and submits a layout of the map, which is shown,
    undone and redone, outdated by a paste, saved with the table and discarded."""
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("qtpy")
    import time

    from qtpy.QtCore import QSettings
    from qtpy.QtWidgets import QApplication, QMessageBox

    from ecmapper.draft_map.ai.providers import ToolCall, Turn
    from ecmapper.draft_map.gui import chat as chat_module
    from ecmapper.draft_map.gui.window import MainWindow
    from test_ai import ScriptedChat, fake_toolkit, shifted

    app = QApplication.instance() or QApplication(["ecmapper-test"])
    s = session
    w = MainWindow(s, QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="10" height="10"/></svg>'
    monkeypatch.setattr(w, "export_svg", lambda escher_map: svg)
    s.export_svg = w.export_svg
    w.build_tabs()
    w.resize(1200, 800)
    w.show()
    central, other = s.maps[0], s.maps[-1]
    w.tabs.setCurrentWidget(w.views[central])
    dock = w.ai_dock
    dock.toolkits = [fake_toolkit(tmp_path)]
    dock.fill_toolkits(0)
    dock.approve_all.setChecked(True)

    def submit_shifted():
        folder = dock.job["folder"]
        data = json.loads((folder / "v1" / "draft.escher_final.json").read_text())
        (folder / "v2").mkdir()
        (folder / "v2" / "final.json").write_text(json.dumps(shifted(data, 400, 200)))
        return Turn("Here it is.", [ToolCall("2", "submit_layout", {"path": "v2/final.json"})], "tool_use")

    turns = [Turn("Running the toolkit.", [ToolCall("1", "run_toolkit", {"args": ["run", "draft.escher.json", "v1"]})],
                  "tool_use"),
             submit_shifted, Turn("Done.", [], "end")]

    class Chat(ScriptedChat):
        def step(self, on_text=None, on_thinking=None, stop=None):
            if self.turns and callable(self.turns[0]):
                self.turns[0] = self.turns[0]()
            return super().step(on_text, on_thinking, stop)

    chats = []
    monkeypatch.setattr(chat_module, "make_chat",
                        lambda connection, system, tools: chats.append(Chat(turns)) or chats[-1])
    steps = w.undo.count()
    dock.input.setPlainText("Lay out the map")
    dock.send()
    start = time.time()
    while dock.worker is not None and time.time() - start < 60:
        app.processEvents()
        time.sleep(0.01)
    assert dock.worker is None and chats[0].users == ["Lay out the map"]
    results = [r for batch in chats[0].results for r in batch]
    assert [r.is_error for r in results] == [False, False], [r.text for r in results]
    assert results[1].text.startswith("Accepted")
    assert (dock.job["folder"] / "draft.escher.json").is_file() and (dock.job["folder"] / "transcript.jsonl").is_file()
    text = dock.transcript.toPlainText()
    assert "You: Lay out the map" in text and "AI: Running the toolkit." in text and "toolkit run" in text

    layout = s.layouts[central]
    assert layout.svg == svg and s.layout_info(central)[2] == []
    m = s.map_draft(central)
    for i, (x, y) in layout.positions.items():
        n = m.nodes[i]
        if n.kind == "enzyme" and len(n.reactions) == 1:
            assert (x, y) == pytest.approx((n.x + 200, n.y + 100), abs=0.5)
    scene = w.views[central].map_scene()
    assert scene.escher_item is not None and scene.read_only and w.tabs.tabText(0).startswith("✦")
    assert all(item.outline is not None for item in scene.boxes.values())
    assert not w.grab().isNull()
    assert w.undo.count() == steps + 1
    w.undo.undo()
    assert central not in s.layouts and scene.escher_item is None and not scene.read_only
    w.undo.redo()
    assert s.layouts[central] is layout and scene.escher_item is not None

    # boxes view, paste a subsystem of another map: the layout is outdated (after asking)
    w.toggle_ai_view()
    assert scene.escher_item is None and not scene.read_only
    k = next(k for k in s.positions[other] if k not in s.positions[central])
    w.views[other].map_scene().select([k])
    w.tabs.setCurrentWidget(w.views[other])
    w.copy()
    w.tabs.setCurrentWidget(w.views[central])
    monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: QMessageBox.StandardButton.No)
    w.paste()
    assert k not in s.positions[central] and not s.outdated(central)
    monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: QMessageBox.StandardButton.Yes)
    w.paste()
    assert k in s.positions[central] and s.outdated(central) and w.tabs.tabText(0).startswith("⚠")
    assert dock.update_button.isVisibleTo(dock)
    w.undo.undo()
    assert not s.outdated(central)

    # saved with the table in a project, read again
    monkeypatch.setenv("ECMAPPER_CACHE", str(tmp_path / "cache"))
    s.project = Project.new(s.args.prep, None, s.args.brenda_json, s.args.brenda_svg, options_of(s.args))
    w.dirty = True  # a new project
    assert w.modified()
    assert w.save_to(tmp_path / "a.ecmap") and not w.modified() and w.windowTitle().startswith("a.ecmap")
    folder = s.project.work / "assignment_layouts" / split.slug(central)
    assert (folder / "layout.json").is_file() and (folder / "final.svg").read_text() == svg
    table = s.project.work / "assignment.xlsx"
    assert read_assignment(table).layouts == {central: f"assignment_layouts/{split.slug(central)}"}
    with zipfile.ZipFile(tmp_path / "a.ecmap") as z:
        assert f"assignment_layouts/{split.slug(central)}/final.svg" in z.namelist()
    w.tabs.setCurrentWidget(w.views[central])
    w.discard_layout()
    assert central not in s.layouts

    # a folder for a coding agent; its result is taken by its Escher IDs
    infos = []
    monkeypatch.setattr(QMessageBox, "exec", lambda box: infos.append(box.text()) or 0)
    monkeypatch.setattr(QMessageBox, "information", lambda parent, title, text: infos.append(text))
    w.use_agent_result()
    assert "Prepare folder for a coding agent" in infos[-1]
    w.prepare_agent_folder()
    folder = next((Path(s.args.out) / "ai" / split.slug(central)).glob("*-agent"))
    assert str(folder) in infos[-1] and (folder / "AGENTS.md").is_file() and (folder / "run_toolkit.py").is_file()
    w.use_agent_result()
    assert "no final map" in infos[-1] and central not in s.layouts
    data = json.loads((folder / "draft.escher.json").read_text())
    (folder / "v1").mkdir()
    (folder / "v1" / "draft.escher_final.json").write_text(json.dumps(shifted(data, 400, 200)))
    w.use_agent_result()
    layout = s.layouts[central]
    assert layout.source["connection"] == "coding agent" and layout.source["file"].startswith("v1")
    assert s.layout_info(central)[2] == [] and w.tabs.tabText(0).startswith("✦")

    # a web chat by copy/paste: the user copies the messages into the chat and its replies back
    from ecmapper.draft_map.ai.providers import WEBCHAT, Connection

    dock.connections.append(Connection("Chat", kind=WEBCHAT, model="web chat", transport="clipboard"))
    dock.fill_connections("Chat")
    replies = ['```json\n{"tool": "run_toolkit", "input": {"args": ["run", "draft.escher.json", "v1"]}}\n```\n'
               "@@END 1@@", "It ran.\n@@END 2@@"]
    dock.job = None
    dock.input.setPlainText("Lay out the map")
    dock.send()
    start = time.time()
    while dock.worker is not None and time.time() - start < 60:
        app.processEvents()
        if dock.clipboard_panel.current is not None and replies:
            assert QApplication.clipboard().text().startswith("ecmapper message")
            QApplication.clipboard().setText(replies.pop(0))
            dock.clipboard_panel.paste()
        time.sleep(0.01)
    assert dock.worker is None and not replies and not dock.clipboard_panel.isVisibleTo(dock)
    assert "AI: It ran." in dock.transcript.toPlainText() and dock.status.text() == "Web chat: 2 messages"
    assert (dock.job["folder"] / "v1" / "draft.escher_final.json").is_file()
    w.undo.setClean()
    w.close()


FAKE_CHAT = """<!doctype html><html><head><meta charset="utf-8"></head><body>
<div id="log" style="white-space: pre-wrap"></div>
<textarea id="input" style="width: 600px; height: 80px"></textarea>
<script>
const input = document.getElementById('input'), log = document.getElementById('log');
input.addEventListener('keydown', e => {
  if (e.key !== 'Enter' || e.shiftKey) return;
  e.preventDefault();
  const text = input.value, m = /ecmapper message (\\d+)/.exec(text);
  input.value = '';
  const sent = document.createElement('div');
  sent.textContent = text;
  log.appendChild(sent);
  if (!m) return;
  const answer = document.createElement('div'), words = ['Reply', 'to', 'message', m[1], '\\n@@END ' + m[1] + '@@'];
  log.appendChild(answer);
  let i = 0;
  const timer = setInterval(() => {
    answer.textContent += (i ? ' ' : '') + words[i];
    if (++i === words.length) clearInterval(timer);
  }, 300);
});
</script></body></html>"""


def test_web_chat(tmp_path, monkeypatch):
    """The transports of web chats: the copy/paste panel, and the browser with a fake chat page (skipped without
    QtWebEngine)."""
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("qtpy")
    import threading
    import time

    from ecmapper.draft_map.gui.escher_svg import web_engine_available

    web = web_engine_available()
    from qtpy.QtCore import QTimer, QUrl
    from qtpy.QtWidgets import QApplication, QMessageBox

    from ecmapper.draft_map.ai.providers import ProviderError
    from ecmapper.draft_map.ai.webchat import tail
    from ecmapper.draft_map.gui import webchat as webchat_gui

    app = QApplication.instance() or QApplication(["ecmapper-test"])  # noqa: F841
    infos = []
    monkeypatch.setattr(QMessageBox, "information", lambda parent, title, text: infos.append(text))
    panel = webchat_gui.ClipboardPanel(lambda fn: fn())
    message = f"ecmapper message 1\n\nHello\n\n{tail(1)}"

    def user():
        assert QApplication.clipboard().text() == message and panel.isVisible()
        panel.paste()  # the message itself: not a reply
        QApplication.clipboard().setText("Hi.\n@@END 1@@")
        panel.paste()

    QTimer.singleShot(200, user)
    assert panel.exchange(message, 1, []) == "Hi." and "does not hold a reply" in infos[0]
    assert panel.current is None and not panel.isVisible()
    assert panel.exchange(message, 1, [], stop=lambda: True) is None
    if not web:
        pytest.skip("no QtWebEngine")
    page = tmp_path / "chat.html"
    page.write_text(FAKE_CHAT, encoding="utf-8")
    browser = webchat_gui.BrowserWindow(lambda fn: fn())
    loaded = []
    browser.view.loadFinished.connect(loaded.append)
    browser.open(QUrl.fromLocalFile(str(page)).toString())
    start = time.time()
    while not loaded and time.time() - start < 30:
        app.processEvents()
        time.sleep(0.01)
    assert loaded == [True]
    # a job with the same site keeps the open page; another site (a changed connection) is loaded
    calls = []
    browser.view.load = calls.append
    browser.open(QUrl.fromLocalFile(str(page)).toString())
    browser.open(" chat.example.org ")
    assert [u.host() for u in calls] == ["chat.example.org"]
    del browser.view.load
    # message 1 waits for Start (logged in, a new conversation open); then the window pastes, sends and reads
    assert browser.exchange("ecmapper message 1\n\nHi", 1, [], stop=lambda: time.time() - start > 3) is None
    assert "Stopped" in browser.status.text() and not browser.started
    QTimer.singleShot(300, browser.start)
    QApplication.clipboard().clear()  # an empty clipboard stays empty
    start = time.time()
    reply = browser.exchange(f"ecmapper message 1\n\nLay out\n\n{tail(1)}", 1, [],
                             stop=lambda: time.time() - start > 60)
    assert reply == "Reply to message 1", browser.status.text()
    assert QApplication.clipboard().text() == ""
    QApplication.clipboard().setText("the user's text")  # is restored after the paste
    start = time.time()
    assert browser.exchange(f"ecmapper message 2\n\nResults\n\n{tail(2)}", 2, [],
                            stop=lambda: time.time() - start > 60) == "Reply to message 2"
    assert QApplication.clipboard().text() == "the user's text"
    browser.current = dict(number=4, done=threading.Event())
    browser.read("text")  # an error in a callback ends the exchange, not the program
    assert "the browser failed" in browser.status.text() and browser.current is None
    monkeypatch.setattr(webchat_gui, "blocked", lambda url: "not here")
    with pytest.raises(ProviderError, match="not here"):
        browser.exchange(f"ecmapper message 3\n\nMore\n\n{tail(3)}", 3, [])
    browser.close()


def test_escher_svg_export():
    """Escher exports a map as SVG in a hidden web view (skipped without QtWebEngine or escher.min.js)."""
    pytest.importorskip("qtpy")
    from ecmapper.draft_map.gui.escher_svg import EscherSvg, find_script, web_engine_available

    if not web_engine_available():
        pytest.skip("no QtWebEngine")
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(["ecmapper-test"])  # noqa: F841
    import os
    from pathlib import Path

    script = Path(os.environ["ECMAPPER_ESCHER_JS"]) if os.environ.get("ECMAPPER_ESCHER_JS") else find_script()
    if script is None:
        pytest.skip("escher.min.js not found (ECMAPPER_ESCHER_JS, setting ai/escher_js, escher package or download)")
    body = dict(reactions={"1": dict(name="R", bigg_id="R", reversibility=False, label_x=60, label_y=40,
                                     gene_reaction_rule="", genes=[],
                                     metabolites=[dict(bigg_id="a", coefficient=-1), dict(bigg_id="b", coefficient=1)],
                                     segments={"5": dict(from_node_id="2", to_node_id="3", b1=None, b2=None),
                                               "6": dict(from_node_id="3", to_node_id="4", b1=None, b2=None)})},
                nodes={"2": dict(node_type="metabolite", x=0, y=0, bigg_id="a", name="A", label_x=5, label_y=-10,
                                 node_is_primary=True),
                       "3": dict(node_type="midmarker", x=100, y=0),
                       "4": dict(node_type="metabolite", x=200, y=0, bigg_id="b", name="B", label_x=205,
                                 label_y=-10, node_is_primary=True)},
                text_labels={}, canvas=dict(x=-100, y=-100, width=400, height=200))
    try:
        exporter = EscherSvg(script)
    except RuntimeError as e:
        pytest.skip(f"Escher could not be run here: {e}")
    svg = exporter.export([dict(map_name="test", map_id="test", map_description="", homepage="", schema=""), body])
    exporter.close()
    assert "<svg" in svg and 'width="400"' in svg and 'height="200"' in svg
