"""The main window of the assignment GUI: a tab per map, editing with undo, packing, building the maps."""

from __future__ import annotations

import html
import json
import shlex
import traceback
from pathlib import Path

from qtpy.QtCore import QProcess, QProcessEnvironment, QRectF, QThread, QTimer, QUrl, Qt, Signal
from qtpy.QtGui import QAction, QActionGroup, QColor, QCursor, QDesktopServices, QKeySequence, QUndoCommand, QUndoStack
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMenu,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..assignment import ADDED_COLOUR, MIXED_COLOUR
from ..split import GROUP_NAMES
from .items import AUTO, BOXES, DETAIL_SCALE, DRAWINGS, MapScene, MapView
from .session import NOT_DRAWN, Session, set_inputs, template

TITLE = "ecmapper draft map assignment"


def subsystems(n: int) -> str:
    return f"{n} subsystem{'' if n == 1 else 's'}"


class Snapshot(QUndoCommand):
    """An undoable change: the states of the session before and after it (pushed after the change was made)."""

    def __init__(self, window: MainWindow, before, after, text: str):
        super().__init__(text)
        self.window, self.before, self.after = window, before, after
        self.first = True

    def redo(self) -> None:
        if self.first:  # the change is made already
            self.first = False
            return
        self.window.session.restore(self.after)
        self.window.refresh()

    def undo(self) -> None:
        self.window.session.restore(self.before)
        self.window.refresh()


class LoadWorker(QThread):
    """Steps 1-4 of the pipeline (:meth:`.session.Session.load`) in the background."""

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session, self.error = session, None

    def run(self) -> None:
        try:
            self.session.load()
        except Exception:  # noqa: BLE001 (shown to the user)
            self.error = traceback.format_exc()


class PackWorker(QThread):
    """The simulated annealing of a map (:meth:`.session.Session.pack`) in the background."""

    progress = Signal(float)

    def __init__(self, session: Session, name: str, iterations: int, seed: int, seeds: int, parent=None):
        super().__init__(parent)
        self.session, self.name = session, name
        self.iterations, self.seed, self.seeds = iterations, seed, seeds
        self.stopped = False
        self.result = self.error = None

    def stop(self) -> None:
        self.stopped = True

    def run(self) -> None:
        try:
            self.result = self.session.pack(self.name, self.iterations, self.seed, self.seeds,
                                            progress=self.progress.emit, stop=lambda: self.stopped)
        except Exception:  # noqa: BLE001 (shown to the user)
            self.error = traceback.format_exc()


def spin_box(value: int, maximum: int, step: int = 1) -> QSpinBox:
    box = QSpinBox()
    box.setRange(1, maximum)
    box.setSingleStep(step)
    box.setValue(value)
    box.setGroupSeparatorShown(True)
    return box


class PackDialog(QDialog):
    def __init__(self, parent, name: str, args):
        super().__init__(parent)
        self.setWindowTitle("Pack map")
        layout = QFormLayout(self)
        layout.addRow(QLabel(f"Arrange the subsystems of <b>{html.escape(name)}</b> by simulated annealing, as the "
                             "pipeline does, starting from their current places."))
        self.iterations = spin_box(args.iterations, 100_000_000, 10000)
        self.seeds = spin_box(args.seeds, 64)
        self.seed = spin_box(args.seed, 1_000_000)
        layout.addRow("Iterations", self.iterations)
        layout.addRow("Runs (seeds), the best is kept", self.seeds)
        layout.addRow("First seed", self.seed)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)


class BuildDialog(QDialog):
    def __init__(self, parent, session: Session):
        super().__init__(parent)
        self.setWindowTitle("Build maps")
        args = session.args
        layout = QFormLayout(self)
        by_pipeline = [m for m in session.maps if session.positions[m] and not session.arranged.get(m)]
        note = (f"Builds the maps of <b>{html.escape(session.path.name)}</b> with "
                "<code>ecmapper.draft_map.run --assignment</code>. Maps arranged here are drawn as arranged")
        note += (f"; the pipeline arranges {html.escape(', '.join(by_pipeline))}." if by_pipeline else ".")
        label = QLabel(note)
        label.setWordWrap(True)
        layout.addRow(label)
        row = QHBoxLayout()
        self.out = QLineEdit(str(args.out))
        browse = QPushButton("…")
        browse.clicked.connect(self.browse)
        row.addWidget(self.out)
        row.addWidget(browse)
        layout.addRow("Output folder", row)
        self.iterations = spin_box(args.iterations, 100_000_000, 10000)
        self.seeds = spin_box(args.seeds, 64)
        layout.addRow("Iterations (maps arranged by the pipeline)", self.iterations)
        layout.addRow("Runs (seeds) per map", self.seeds)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Build")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)
        self.resize(640, self.sizeHint().height())

    def browse(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Output folder", self.out.text())
        if path:
            self.out.setText(path)


class FileField(QWidget):
    """A path with a button that opens a file dialog (`kind`: "open", "save" or "folder")."""

    def __init__(self, path, kind: str = "open", filter: str = "", parent=None):
        super().__init__(parent)
        self.kind, self.filter = kind, filter
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.edit = QLineEdit(str(path) if path else "")
        browse = QPushButton("…")
        browse.clicked.connect(self.browse)
        layout.addWidget(self.edit)
        layout.addWidget(browse)

    def path(self) -> Path | None:
        text = self.edit.text().strip()
        return Path(text) if text else None

    def browse(self) -> None:
        start = self.edit.text()
        if self.kind == "folder":
            path = QFileDialog.getExistingDirectory(self, "Folder", start)
        else:
            path, _ = QFileDialog.getOpenFileName(self, "File", start, self.filter)
        if path:
            self.edit.setText(path)
            self.edited()

    def edited(self) -> None:
        """Hook for dialogs (a new model proposes an output folder)."""


class StartDialog(QDialog):
    """The input files of a session: the model (or a PREP.pkl), the BRENDA map as template, the output folder."""

    def __init__(self, parent, args):
        super().__init__(parent)
        self.setWindowTitle("New draft maps from a model")
        self.args = args
        layout = QFormLayout(self)
        intro = QLabel("The model is matched to the BRENDA map (ecmapper; the resources are downloaded on first use), "
                       "which is the template of the draft maps: its subsystems with model reactions, completed with "
                       "further model reactions, split into maps.")
        intro.setWordWrap(True)
        layout.addRow(intro)
        default_json, default_svg = template()
        self.model = FileField(args.model or args.prep, "open",
                               "Models (*.xml *.sbml *.json *.yml *.yaml *.mat *.gz *.bz2 *.xz *.zip *.pkl);;"
                               "All files (*)")
        self.model.edited = self.model_chosen
        self.subsystems = FileField(getattr(args, "subsystems", None), "open",
                                    "JSON models (*.json *.json.gz *.json.bz2 *.json.xz *.zip)")
        self.assignment = FileField(args.assignment, "open", "Excel tables (*.xlsx)")
        self.out = FileField(args.out, "folder")
        self.brenda_json = FileField(args.brenda_json or default_json, "open", "JSON (*.json)")
        self.brenda_svg = FileField(args.brenda_svg or default_svg, "open", "SVG (*.svg)")
        self.pack = QCheckBox("Pack the maps without inter-pathway links (--pack)")
        self.pack.setChecked(bool(args.pack) or args.model is None and args.prep is None)
        layout.addRow("Model (SBML, JSON, MAT, also zipped; or PREP.pkl)", self.model)
        layout.addRow("Subsystems (JSON model, optional)", self.subsystems)
        layout.addRow("Assignment table (XLSX, optional)", self.assignment)
        layout.addRow("Output folder", self.out)
        layout.addRow("BRENDA map: reactions (JSON)", self.brenda_json)
        layout.addRow("BRENDA map: SVG", self.brenda_svg)
        layout.addRow(self.pack)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Start")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)
        self.resize(860, self.sizeHint().height())

    def model_chosen(self) -> None:
        model = self.model.path()
        if model is not None:
            self.out.edit.setText(str(model.parent / f"{model.name.split('.')[0]}_maps"))

    def problem(self) -> str | None:
        model, out = self.model.path(), self.out.path()
        if model is None or not model.is_file():
            return "Choose the model file."
        if out is None:
            return "Choose the output folder."
        for field in (self.brenda_json, self.brenda_svg, self.subsystems, self.assignment):
            if field.path() is not None and not field.path().is_file():
                return f"{field.path()} does not exist."
        if self.pack.isChecked() and (self.args.hub_subsystems or self.args.flips != "none"):
            return "--pack does not go with --hub-subsystems and --flips."
        return None

    def accept(self) -> None:
        problem = self.problem()
        if problem:
            QMessageBox.warning(self, self.windowTitle(), problem)
            return
        super().accept()

    def apply(self) -> None:
        """Set the inputs of the dialog in the options."""
        args = self.args
        args.subsystems, args.assignment = self.subsystems.path(), self.assignment.path()
        args.pack = self.pack.isChecked()
        set_inputs(args, self.model.path(), self.out.path(), self.brenda_json.path(), self.brenda_svg.path())


class InfoPanel(QWidget):
    """The selected subsystem: what it is, its maps and its model reactions (read only)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Reaction", "Name", "Drawn as", "Equation"])
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.table.horizontalHeader().setStretchLastSection(True)
        legend = QLabel(f'<span style="background-color:#{ADDED_COLOUR}; color:black">&nbsp;Reaction&nbsp;</span> '
                        "added from the model (step 3; not on the BRENDA map) &nbsp; "
                        f'<span style="background-color:#{MIXED_COLOUR}; color:black">&nbsp;Reaction&nbsp;</span> '
                        "on the map and added")
        legend.setWordWrap(True)
        layout.addWidget(self.text)
        layout.addWidget(self.table, 1)
        layout.addWidget(legend)

    def show_subsystems(self, session: Session, keys: list[str]) -> None:
        self.table.setRowCount(0)
        if not keys:
            self.text.setText("<i>No subsystem selected.</i>")
            return
        if len(keys) > 1:
            labels = ", ".join(html.escape(session.dm.subsystems[k].label) for k in keys[:12])
            self.text.setText(f"<b>{subsystems(len(keys))} selected</b>: {labels}{' ...' if len(keys) > 12 else ''}")
            return
        k = keys[0]
        s = session.dm.subsystems[k]
        group = f"{s.color_group} ({GROUP_NAMES[s.color_group]})" if s.color_group in GROUP_NAMES else s.color_group
        maps = session.maps_of(k)
        rows = [f"<b>{html.escape(s.label)}</b>", f"ID: {html.escape(k)}", f"Colour group: {html.escape(group)}",
                f"Maps: {html.escape(', '.join(maps)) if maps else '<i>not drawn</i>'}",
                f"Box: {s.w:.0f} × {s.h:.0f} map units, {session.counts.get(k, 0)} model reactions"
                + (f", {session.added_counts[k]} of them added from the model" if session.added_counts.get(k) else "")]
        if s.aside:
            rows.append(f"Set aside (rule {html.escape(s.aside)})")
        self.text.setText("<br>".join(rows))
        reactions = session.reactions(k)
        self.table.setRowCount(len(reactions))
        fills = {"model": QColor(f"#{ADDED_COLOUR}"), "both": QColor(f"#{MIXED_COLOUR}")}
        for i, (*row, source) in enumerate(reactions):
            for j, value in enumerate(row):
                item = QTableWidgetItem(value)
                if source in fills:
                    item.setBackground(fills[source])
                    item.setForeground(QColor("black"))
                    item.setToolTip("added from the model (step 3; not on the BRENDA map)" if source == "model" else
                                    "on the BRENDA map, and added from the model in another node")
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()


class MainWindow(QMainWindow):
    """A tab per map (and one for the subsystems on no map); see the module :mod:`ecmapper.draft_map.gui`."""

    log_message = Signal(str)

    def __init__(self, session: Session):
        super().__init__()
        self.session = session
        session.log = self.log_message.emit
        self.undo = QUndoStack(self)
        self.undo.cleanChanged.connect(lambda clean: self.setWindowModified(not clean))
        self.clipboard = None  # (subsystem IDs, their positions where they were copied)
        self.views: dict[str | None, MapView] = {}
        self.mode = AUTO
        self.found: list[str] = []
        self.process = self.build_progress = None
        self._reordering = False
        self.tabs = QTabWidget()
        self.tabs.setMovable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabBar().tabMoved.connect(lambda *_: QTimer.singleShot(0, self.tab_moved))
        self.tabs.currentChanged.connect(self.current_changed)
        self.tabs.tabBarDoubleClicked.connect(self.rename_map)
        self.setCentralWidget(self.tabs)
        self.status_label = QLabel()
        self.statusBar().addPermanentWidget(self.status_label)
        self.make_docks()
        self.make_actions()
        self.log_message.connect(self.show_log)
        self.update_title()
        self.update_actions()
        self.resize(1600, 1000)

    # --- setup ---------------------------------------------------------------------------------

    def make_docks(self) -> None:
        self.info = InfoPanel()
        self.info.show_subsystems(self.session, [])
        self.find_edit = QLineEdit()
        self.find_edit.setPlaceholderText("Subsystem name or ID, reaction ID")
        self.find_edit.setClearButtonEnabled(True)
        self.find_edit.textChanged.connect(self.find)
        self.find_list = QListWidget()
        self.find_list.itemActivated.connect(lambda item: self.go_to(item.data(Qt.ItemDataRole.UserRole)))
        find = QWidget()
        layout = QVBoxLayout(find)
        layout.addWidget(self.find_edit)
        layout.addWidget(self.find_list)
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.docks = []
        for title, widget, area in (("Find", find, Qt.DockWidgetArea.RightDockWidgetArea),
                                    ("Subsystem", self.info, Qt.DockWidgetArea.RightDockWidgetArea),
                                    ("Log", self.log_view, Qt.DockWidgetArea.BottomDockWidgetArea)):
            dock = QDockWidget(title, self)
            dock.setObjectName(title)
            dock.setWidget(widget)
            self.addDockWidget(area, dock)
            self.docks.append(dock)
        self.resizeDocks(self.docks[:2], [260, 600], Qt.Orientation.Vertical)
        self.resizeDocks(self.docks[:1], [480], Qt.Orientation.Horizontal)
        self.resizeDocks(self.docks[2:], [150], Qt.Orientation.Vertical)

    def action(self, text: str, slot, shortcut=None, tip: str = "") -> QAction:
        a = QAction(text, self)
        a.triggered.connect(lambda *_: slot())
        if shortcut is not None:
            a.setShortcut(QKeySequence(shortcut))
        if tip:
            a.setStatusTip(tip)
            a.setToolTip(tip)
        return a

    def make_actions(self) -> None:
        SK = QKeySequence.StandardKey
        self.a_new = self.action("&New from model…", self.new_from_model, SK.New,
                                 "Match a model (SBML, JSON, MAT) to the BRENDA map and draw its maps")
        self.a_open = self.action("&Open table…", self.open_table, SK.Open, "Open another assignment table (XLSX)")
        self.a_save = self.action("&Save", self.save, SK.Save, "Save the assignment table (XLSX)")
        self.a_save_as = self.action("Save &as…", self.save_as, SK.SaveAs)
        self.a_build = self.action("&Build maps…", self.build_maps, "Ctrl+B",
                                   "Save the table and build the final maps with ecmapper.draft_map.run")
        self.a_quit = self.action("&Quit", self.close, SK.Quit)
        self.a_undo = self.undo.createUndoAction(self, "&Undo")
        self.a_undo.setShortcut(QKeySequence(SK.Undo))
        self.a_redo = self.undo.createRedoAction(self, "&Redo")
        self.a_redo.setShortcut(QKeySequence(SK.Redo))
        self.a_cut = self.action("Cu&t", self.cut, SK.Cut, "Take the selected subsystems off this map, to paste them")
        self.a_copy = self.action("&Copy", self.copy, SK.Copy,
                                  "Copy the selected subsystems (pasted, they are drawn on both maps)")
        self.a_paste = self.action("&Paste", self.paste, SK.Paste,
                                   "Put the copied subsystems on this map (at the mouse pointer or at a free place)")
        self.a_remove = self.action("&Remove from map", self.remove, SK.Delete)
        self.a_select_all = self.action("Select &all", self.select_all, SK.SelectAll)
        self.a_fit = self.action("Zoom to &fit", lambda: self.current_view() and self.current_view().zoom_to_fit(),
                                 "Ctrl+0")
        self.a_zoom_in = self.action("Zoom &in", lambda: self.current_view() and self.current_view().zoom(1.25),
                                     SK.ZoomIn)
        self.a_zoom_out = self.action("Zoom &out", lambda: self.current_view() and self.current_view().zoom(0.8),
                                      SK.ZoomOut)
        self.a_zoom_sel = self.action("Zoom to &selection", self.zoom_to_selection, "Ctrl+E")
        group = QActionGroup(self)
        self.mode_actions = {}
        for mode, text, key, tip in ((AUTO, "&Auto (by zoom)", "F2", "Boxes when zoomed out, drawings when zoomed in"),
                                     (BOXES, "&Boxes", "F3", "Subsystems as boxes with their names"),
                                     (DRAWINGS, "&Drawings", "F4", "Subsystems as drawn on the map (SVG)")):
            a = self.action(text, lambda m=mode: self.set_mode(m), key, tip)
            a.setCheckable(True)
            a.setChecked(mode == self.mode)
            group.addAction(a)
            self.mode_actions[mode] = a
        self.a_pack = self.action("&Pack map…", self.pack_map, "Ctrl+Shift+P",
                                  "Arrange the subsystems of this map by simulated annealing")
        self.a_original = self.action("&Original order", self.original_order, None,
                                      "Arrange the subsystems of this map as on the BRENDA map, compacted")
        self.a_pipeline = self.action("&Let the pipeline arrange this map", self.pipeline_arranges, None,
                                      "Do not save the arrangement of this map; the pipeline anneals it")
        self.a_new_map = self.action("&New map…", self.new_map, "Ctrl+Shift+N")
        self.a_rename_map = self.action("&Rename map…", lambda: self.rename_map(self.tabs.currentIndex()), "F6")
        self.a_delete_map = self.action("&Delete map", self.delete_map)

        bar = self.menuBar()
        m = bar.addMenu("&File")
        for a in (self.a_new, self.a_open, self.a_save, self.a_save_as, None, self.a_build, None, self.a_quit):
            m.addSeparator() if a is None else m.addAction(a)
        m = bar.addMenu("&Edit")
        for a in (self.a_undo, self.a_redo, None, self.a_cut, self.a_copy, self.a_paste, self.a_remove, None,
                  self.a_select_all):
            m.addSeparator() if a is None else m.addAction(a)
        self.copy_menu = m.addMenu("Copy &to map")
        self.move_menu = m.addMenu("&Move to map")
        self.copy_menu.aboutToShow.connect(lambda: self.fill_map_menu(self.copy_menu, False))
        self.move_menu.aboutToShow.connect(lambda: self.fill_map_menu(self.move_menu, True))
        m = bar.addMenu("&View")
        for a in (self.a_fit, self.a_zoom_in, self.a_zoom_out, self.a_zoom_sel, None, *self.mode_actions.values(),
                  None):
            m.addSeparator() if a is None else m.addAction(a)
        for dock in self.docks:
            m.addAction(dock.toggleViewAction())
        m = bar.addMenu("&Arrange")
        for a in (self.a_pack, self.a_original, None, self.a_pipeline):
            m.addSeparator() if a is None else m.addAction(a)
        m = bar.addMenu("&Maps")
        for a in (self.a_new_map, self.a_rename_map, self.a_delete_map):
            m.addAction(a)
        tools = self.addToolBar("Tools")
        tools.setObjectName("Tools")
        for a in (self.a_save, self.a_undo, self.a_redo, None, self.a_cut, self.a_copy, self.a_paste, None, self.a_fit,
                  *self.mode_actions.values(), None, self.a_pack, None, self.a_build):
            tools.addSeparator() if a is None else tools.addAction(a)

    # --- loading -------------------------------------------------------------------------------

    def start_loading(self) -> None:
        """Run steps 1-4 of the pipeline in the background, then show the maps."""
        self.setEnabled(False)
        text = "Drawing the subsystems (steps 1-4 of the pipeline) …"
        if getattr(self.session.args, "model", None):
            text = "Matching the model to the BRENDA map, drawing the subsystems …"
        dialog = QProgressDialog(text, "", 0, 0, self)
        dialog.setWindowTitle(TITLE)
        dialog.setCancelButton(None)
        dialog.setMinimumDuration(0)
        dialog.setMinimumWidth(700)
        dialog.setWindowModality(Qt.WindowModality.WindowModal)
        last_line = lambda text: dialog.setLabelText(text[:120])  # noqa: E731
        self.log_message.connect(last_line)
        worker = LoadWorker(self.session, self)

        def done():
            self.log_message.disconnect(last_line)
            dialog.close()
            self.setEnabled(True)
            if worker.error:
                self.show_log(worker.error)
                QMessageBox.critical(self, TITLE, f"The draft map could not be made:\n\n{worker.error[-1500:]}")
                return
            self.build_tabs()
            self.undo.clear()
            self.undo.setClean()
            self.update_title()

        worker.finished.connect(done)
        dialog.show()
        worker.start()
        self.load_worker = worker

    def build_tabs(self) -> None:
        self.tabs.blockSignals(True)
        while self.tabs.count():
            self.tabs.widget(0).deleteLater()
            self.tabs.removeTab(0)
        self.views = {}
        for name in [*self.session.maps, None]:
            self.add_view(name)
        self.tabs.blockSignals(False)
        self.tabs.setCurrentIndex(0)
        self.refresh_tab_titles()
        self.current_changed()
        self.find(self.find_edit.text())

    def add_view(self, name: str | None) -> MapView:
        scene = MapScene(self.session, name, self)
        scene.set_mode(self.mode)
        scene.sync()
        scene.selectionChanged.connect(self.selection_changed)
        view = MapView(scene)
        view.moved.connect(lambda positions, v=view: self.move(v, positions))
        view.menu.connect(self.context_menu)
        view.zoom_changed.connect(self.update_status)
        self.tabs.addTab(view, name or NOT_DRAWN)
        self.views[name] = view
        return view

    # --- state of the window -----------------------------------------------------------------

    def current_view(self) -> MapView | None:
        return self.tabs.currentWidget()

    def current_scene(self) -> MapScene | None:
        view = self.current_view()
        return view.map_scene() if view is not None else None

    def current_map(self) -> str | None:
        scene = self.current_scene()
        return scene.name if scene is not None else None

    def update_title(self) -> None:
        name = self.session.path.name if self.session.path else "new assignment"
        self.setWindowTitle(f"{name}[*] – {TITLE}")
        self.setWindowModified(not self.undo.isClean())

    def refresh(self) -> None:
        """Tabs, items and texts as in the session (after a change, undo or redo)."""
        names = self.session.maps
        for name in [n for n in self.views if n is not None and n not in names]:
            view = self.views.pop(name)
            self.tabs.removeTab(self.tabs.indexOf(view))
            view.deleteLater()
        for name in names:
            if name not in self.views:
                self.add_view(name)
        self._reordering = True
        for i, name in enumerate([*names, None]):
            j = self.tabs.indexOf(self.views[name])
            if j != i:
                self.tabs.tabBar().moveTab(j, i)
        self._reordering = False
        for view in self.views.values():
            view.map_scene().sync()
        if self.found:
            self.find(self.find_edit.text())
        self.refresh_tab_titles()
        self.selection_changed()
        self.update_status()

    def refresh_tab_titles(self) -> None:
        for i in range(self.tabs.count()):
            scene = self.tabs.widget(i).map_scene()
            n = len(scene.positions())
            found = sum(1 for k in self.found if k in scene.boxes)
            if scene.name is None:
                text, tip = f"{NOT_DRAWN} ({n})", "Subsystems on no map; copy them onto a map"
            else:
                arranged = self.session.arranged.get(scene.name)
                text = f"{'✎ ' if arranged else ''}{scene.name} ({n})"
                tip = "arranged here (the pipeline keeps it)" if arranged else "arranged by the pipeline when building"
            self.tabs.setTabText(i, text + (f" [{found} found]" if found else ""))
            self.tabs.setTabToolTip(i, tip)

    def current_changed(self, *_) -> None:
        self.selection_changed()
        self.update_status()

    def selection_changed(self) -> None:
        scene = self.current_scene()
        self.info.show_subsystems(self.session, scene.selected_keys() if scene is not None else [])
        self.update_actions()

    def update_actions(self) -> None:
        scene = self.current_scene()
        loaded = scene is not None
        on_map = loaded and not scene.read_only
        selected = bool(loaded and scene.selected_keys())
        for a in (self.a_copy, self.a_zoom_sel):
            a.setEnabled(selected)
        self.a_cut.setEnabled(selected)
        self.a_remove.setEnabled(selected and on_map)
        self.a_paste.setEnabled(on_map and self.clipboard is not None)
        self.copy_menu.setEnabled(selected)
        self.move_menu.setEnabled(selected)
        for a in (self.a_pack, self.a_original, self.a_rename_map, self.a_delete_map):
            a.setEnabled(on_map)
        self.a_pipeline.setEnabled(on_map and bool(self.session.arranged.get(scene.name)))
        self.a_open.setEnabled(bool(self.session.args.prep or getattr(self.session.args, "model", None)))
        for a in (self.a_save, self.a_save_as, self.a_build, self.a_new_map, self.a_select_all, self.a_fit,
                  self.a_zoom_in, self.a_zoom_out):
            a.setEnabled(loaded)

    def update_status(self) -> None:
        view = self.current_view()
        if view is None:
            self.status_label.setText("")
            return
        scene, scale = view.map_scene(), view.scale_factor()
        drawn = "drawings" if scene.mode == DRAWINGS or (scene.mode == AUTO and scale >= DETAIL_SCALE) else "boxes"
        if scene.name is None:
            text = f"{subsystems(len(scene.positions()))} on no map"
        else:
            how = "arranged here" if self.session.arranged.get(scene.name) else "arranged by the pipeline"
            text = f"{scene.name}: {subsystems(len(scene.positions()))}, {how}"
        self.status_label.setText(f"{text}  |  zoom {scale:.2f} ({drawn})")

    def show_log(self, text: str) -> None:
        self.log_view.appendPlainText(text.rstrip())
        self.statusBar().showMessage(text.strip().splitlines()[-1][:200] if text.strip() else "", 8000)

    def status(self, text: str) -> None:
        self.statusBar().showMessage(text, 6000)

    # --- changes -------------------------------------------------------------------------------

    def change(self, text: str, apply, select=None, scene: MapScene | None = None) -> None:
        """Make a change of the session (`apply()`) undoable and show it; select `select` on `scene`."""
        before = self.session.state()
        apply()
        after = self.session.state()
        if after == before:
            return
        self.undo.push(Snapshot(self, before, after, text))
        self.refresh()
        if select:
            (scene or self.current_scene()).select(select)

    def move(self, view: MapView, positions: dict) -> None:
        name = view.map_scene().name
        self.change(f"Move {subsystems(len(positions))}", lambda: self.session.move(name, positions))

    def copy(self) -> bool:
        scene = self.current_scene()
        keys = scene.selected_keys() if scene is not None else []
        if not keys:
            self.status("No subsystem selected")
            return False
        positions = scene.positions()
        self.clipboard = (keys, {k: positions[k] for k in keys})
        self.status(f"{subsystems(len(keys))} copied")
        self.update_actions()
        return True

    def cut(self) -> None:
        if not self.copy():
            return
        name, keys = self.current_map(), self.clipboard[0]
        if name is not None:
            self.change(f"Cut {subsystems(len(keys))}", lambda: self.session.remove(name, keys))
            self.status(f"{subsystems(len(keys))} cut (subsystems on no map are in the tab {NOT_DRAWN!r})")

    def paste(self) -> None:
        name, view = self.current_map(), self.current_view()
        if name is None or self.clipboard is None:
            return
        keys, layout = self.clipboard
        positions = self.session.paste_positions(name, keys, layout, view.cursor_scene_pos(QCursor.pos()))
        if not positions:
            self.status("The copied subsystems are on this map already")
            return
        self.change(f"Paste {subsystems(len(positions))}", lambda: self.session.add(name, positions),
                    select=list(positions))
        self.show_positions(view, positions)

    def show_positions(self, view: MapView, positions: dict) -> None:
        """Zoom out if the boxes at `positions` are outside the visible part of the map."""
        rect = QRectF()
        for k, (x, y) in positions.items():
            w, h = self.session.size(k)
            rect = rect.united(QRectF(x, y, w, h))
        if not view.mapToScene(view.viewport().rect()).boundingRect().contains(rect):
            view.zoom_to_fit()

    def remove(self) -> None:
        name, scene = self.current_map(), self.current_scene()
        keys = scene.selected_keys() if scene is not None else []
        if name is not None and keys:
            self.change(f"Remove {subsystems(len(keys))} from {name}", lambda: self.session.remove(name, keys))

    def put_on(self, target: str, move: bool) -> None:
        """Copy (or move) the selected subsystems to map `target`, at a free place."""
        source, scene = self.current_map(), self.current_scene()
        keys = scene.selected_keys()
        positions = self.session.paste_positions(target, keys, {k: scene.positions()[k] for k in keys})

        def apply():
            self.session.add(target, positions)
            if move and source is not None:
                self.session.remove(source, keys)

        verb = "Move" if move and source is not None else "Copy"
        self.change(f"{verb} {subsystems(len(keys))} to {target}", apply)
        self.views[target].map_scene().select(positions)
        self.status(f"{subsystems(len(positions))} put on {target}"
                    + (f", {len(keys) - len(positions)} were there already" if len(positions) < len(keys) else ""))

    def fill_map_menu(self, menu: QMenu, move: bool) -> None:
        menu.clear()
        for name in self.session.maps:
            if name != self.current_map():
                menu.addAction(name, lambda n=name: self.put_on(n, move))

    def select_all(self) -> None:
        scene = self.current_scene()
        if scene is not None:
            scene.select(scene.boxes)

    def context_menu(self, global_pos) -> None:
        scene = self.current_scene()
        menu = QMenu(self)
        for a in (self.a_cut, self.a_copy, self.a_paste, self.a_remove):
            menu.addAction(a)
        if scene.selected_keys():
            menu.addSeparator()
            copy_menu = menu.addMenu("Copy to map" if not scene.read_only else "Put on map")
            self.fill_map_menu(copy_menu, False)
            if not scene.read_only:
                self.fill_map_menu(menu.addMenu("Move to map"), True)
        menu.addSeparator()
        menu.addAction(self.a_zoom_sel)
        menu.addAction(self.a_fit)
        menu.exec(global_pos)

    def zoom_to_selection(self) -> None:
        scene, view = self.current_scene(), self.current_view()
        rect = QRectF()
        for k in scene.selected_keys():
            item = scene.boxes[k]
            rect = rect.united(item.mapRectToScene(item.boundingRect()))
        view.zoom_to(rect)

    def set_mode(self, mode: str) -> None:
        self.mode = mode
        for view in self.views.values():
            view.map_scene().set_mode(mode)
        self.mode_actions[mode].setChecked(True)
        self.update_status()

    # --- arranging -----------------------------------------------------------------------------

    def pack_map(self) -> None:
        name = self.current_map()
        if name is None:
            return
        if len(self.session.positions[name]) < 2:
            self.status("Nothing to pack")
            return
        dialog = PackDialog(self, name, self.session.args)
        if not dialog.exec():
            return
        self.session.views()  # here: the worker only reads them
        worker = PackWorker(self.session, name, dialog.iterations.value(), dialog.seed.value(), dialog.seeds.value(),
                            self)
        progress = QProgressDialog(f"Packing {name} …", "Stop", 0, 1000, self)
        progress.setWindowTitle("Pack map")
        progress.setToolTip("Stop: keep the best arrangement so far")
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.setMinimumDuration(0)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        worker.progress.connect(lambda f: progress.setValue(int(1000 * f)))
        progress.canceled.connect(worker.stop)

        def done():
            progress.canceled.disconnect(worker.stop)  # closing the dialog cancels it
            progress.close()
            if worker.error:
                self.show_log(worker.error)
                QMessageBox.critical(self, "Pack map", f"Packing failed:\n\n{worker.error[-1500:]}")
            elif worker.result:
                self.change(f"Pack {name}", lambda: self.session.move(name, worker.result))
                self.views[name].zoom_to_fit()
                self.status(f"{name} packed" + (" (stopped)" if worker.stopped else ""))

        worker.finished.connect(done)
        progress.show()
        worker.start()
        self.pack_worker = worker

    def original_order(self) -> None:
        name = self.current_map()
        if name is not None and self.session.positions[name]:
            positions = self.session.original_order(self.session.positions[name])
            self.change(f"Original order of {name}", lambda: self.session.move(name, positions))
            self.views[name].zoom_to_fit()

    def pipeline_arranges(self) -> None:
        name = self.current_map()
        if name is not None:
            self.change(f"Let the pipeline arrange {name}", lambda: self.session.set_arranged(name, False))

    # --- maps ----------------------------------------------------------------------------------

    def ask_name(self, title: str, label: str, text: str = "", old: str | None = None) -> str | None:
        while True:
            name, ok = QInputDialog.getText(self, title, label, QLineEdit.EchoMode.Normal, text)
            if not ok:
                return None
            problem = self.session.check_name(name, old)
            if problem is None:
                return name.strip()
            QMessageBox.warning(self, title, problem)
            text = name

    def new_map(self) -> None:
        name = self.ask_name("New map", "Name of the new map:")
        if name is not None:
            self.change(f"New map {name}", lambda: self.session.new_map(name))
            self.tabs.setCurrentWidget(self.views[name])

    def rename_map(self, index: int) -> None:
        if index < 0:
            return
        scene = self.tabs.widget(index).map_scene()
        old = scene.name
        if old is None:
            return
        new = self.ask_name("Rename map", "New name of the map:", old, old)
        if new is not None and new != old:
            self.views[new] = self.views.pop(old)
            scene.name = new
            self.change(f"Rename {old} to {new}", lambda: self.session.rename_map(old, new))

    def delete_map(self) -> None:
        name = self.current_map()
        if name is None:
            return
        keys = list(self.session.positions[name])
        if keys:
            only_here = sum(1 for k in keys if self.session.maps_of(k) == [name])
            answer = QMessageBox.question(self, "Delete map",
                                          f"Delete the map {name!r} with {subsystems(len(keys))}? {only_here} of them "
                                          f"are on no other map and go to {NOT_DRAWN!r}.")
            if answer != QMessageBox.StandardButton.Yes:
                return
        self.change(f"Delete map {name}", lambda: self.session.delete_map(name))

    def tab_moved(self) -> None:
        """A tab was dragged: the order of the maps (the tab of the subsystems on no map stays last)."""
        if self._reordering or None not in self.views:
            return
        self._reordering = True
        last = self.tabs.indexOf(self.views[None])
        if last != self.tabs.count() - 1:
            self.tabs.tabBar().moveTab(last, self.tabs.count() - 1)
        self._reordering = False
        order = [self.tabs.widget(i).map_scene().name for i in range(self.tabs.count())]
        order = [n for n in order if n is not None]
        if order != self.session.maps:
            self.change("Reorder maps", lambda: self.session.reorder(order))

    # --- finding -------------------------------------------------------------------------------

    def find(self, text: str) -> None:
        if self.session.dm is None:
            return
        keys = self.session.find(text) if len(text.strip()) >= 2 else []
        keys.sort(key=lambda k: self.session.dm.subsystems[k].label.lower())
        self.found = keys
        for view in self.views.values():
            view.map_scene().set_found(keys)
        self.find_list.clear()
        for k in keys:
            maps = self.session.maps_of(k)
            item = QListWidgetItem(f"{self.session.dm.subsystems[k].label}  —  {', '.join(maps) or NOT_DRAWN}")
            item.setData(Qt.ItemDataRole.UserRole, k)
            item.setToolTip(k)
            self.find_list.addItem(item)
        self.refresh_tab_titles()

    def go_to(self, key: str) -> None:
        """Show subsystem `key`: on the next of its maps (or on the tab of the subsystems on no map)."""
        maps = self.session.maps_of(key) or [None]
        current = self.current_map()
        name = maps[(maps.index(current) + 1) % len(maps)] if current in maps else maps[0]
        view = self.views[name]
        self.tabs.setCurrentWidget(view)
        scene = view.map_scene()
        scene.select([key])
        item = scene.boxes[key]
        view.zoom_to(item.mapRectToScene(item.boundingRect()))

    # --- files ---------------------------------------------------------------------------------

    def save(self) -> bool:
        if self.session.dm is None:
            return False
        if self.session.path is None:
            return self.save_as()
        try:
            path = self.session.save()
        except PermissionError:
            QMessageBox.warning(self, TITLE, f"{self.session.path} can not be written. Is it open in Excel?")
            return False
        self.undo.setClean()
        self.update_title()
        self.status(f"Saved {path}")
        return True

    def save_as(self) -> bool:
        start = str(self.session.path or Path(self.session.args.out) / "assignment.xlsx")
        path, _ = QFileDialog.getSaveFileName(self, "Save assignment table", start, "Excel tables (*.xlsx)")
        if not path:
            return False
        path = Path(path)
        if path.suffix.lower() != ".xlsx":
            path = path.with_suffix(".xlsx")
        self.session.path = path
        return self.save()

    def maybe_save(self) -> bool:
        """True if the changes are saved or may be discarded."""
        if self.undo.isClean():
            return True
        answer = QMessageBox.question(
            self, TITLE, "Save the changes of the assignment?",
            QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
        if answer == QMessageBox.StandardButton.Save:
            return self.save()
        return answer == QMessageBox.StandardButton.Discard

    def new_from_model(self) -> bool:
        """Choose the model, the template and the output folder (:class:`StartDialog`), then start."""
        if self.process is not None:
            self.status("The maps are being built")
            return False
        if not self.maybe_save():
            return False
        dialog = StartDialog(self, self.session.args)
        if not dialog.exec():
            return False
        dialog.apply()
        self.session.path = self.session.args.assignment
        self.update_title()
        self.start_loading()
        return True

    def open_table(self) -> None:
        if not self.maybe_save():
            return
        start = str(self.session.path or self.session.args.out)
        path, _ = QFileDialog.getOpenFileName(self, "Open assignment table", start, "Excel tables (*.xlsx)")
        if path:
            self.session.args.assignment = self.session.path = Path(path)
            self.start_loading()

    def closeEvent(self, event) -> None:
        if self.maybe_save():
            if self.process is not None:
                self.process.kill()
                self.process.waitForFinished(3000)
            event.accept()
        else:
            event.ignore()

    # --- building the maps ---------------------------------------------------------------------

    def build_maps(self) -> None:
        if self.process is not None:
            self.status("The maps are being built")
            return
        if (not self.undo.isClean() or self.session.path is None or not self.session.path.exists()) and not self.save():
            return
        dialog = BuildDialog(self, self.session)
        if not dialog.exec():
            return
        out = Path(dialog.out.text())
        command = self.session.build_command(out, dialog.iterations.value(), dialog.seeds.value())
        self.show_log("$ " + " ".join(shlex.quote(c) for c in command))
        process = QProcess(self)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        process.setProcessEnvironment(env)
        process.readyReadStandardOutput.connect(self.read_build_output)
        process.finished.connect(lambda code, status: self.build_finished(code, status, out))
        self.build_progress = QProgressDialog("Building the maps …", "Cancel", 0, 0, self)
        self.build_progress.setWindowTitle("Build maps")
        self.build_progress.setMinimumWidth(700)
        self.build_progress.setMinimumDuration(0)
        self.build_progress.canceled.connect(process.kill)
        self.build_progress.show()
        self.process = process
        process.start(command[0], command[1:])

    def read_build_output(self) -> None:
        text = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        lines = [line for line in text.splitlines() if line.strip()]
        if lines:
            self.log_view.appendPlainText("\n".join(lines))
            if self.build_progress is not None:
                self.build_progress.setLabelText(lines[-1][:120])

    def build_finished(self, code: int, status, out: Path) -> None:
        self.read_build_output()
        self.process = None
        if self.build_progress is not None:
            self.build_progress.canceled.disconnect()  # closing the dialog cancels it
            self.build_progress.close()
            self.build_progress = None
        if status != QProcess.ExitStatus.NormalExit or code != 0:
            crashed = status != QProcess.ExitStatus.NormalExit
            self.show_log("cancelled" if crashed else f"ecmapper.draft_map.run stopped with exit code {code}")
            if not crashed:
                QMessageBox.warning(self, "Build maps", f"Building the maps failed (exit code {code}); see the log.")
            return
        try:
            stats = json.loads((out / "stats.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            stats = {}
        split, cov = stats.get("split", {}), stats.get("coverage", {})
        rows = [f"<b>{split.get('maps', '?')} maps</b> with {split.get('reactions', '?')} reactions "
                f"({split.get('reactions_on_several_maps', 0)} on several maps; "
                f"{subsystems(split.get('subsystem_copies', 0))} drawn on further maps)"]
        if cov:
            rows.append(f"Drawn: {cov.get('drawn', '?')} of {cov.get('reactions', '?')} model reactions; "
                        f"{cov.get('coverage_other', 0):.0%} of those that are neither boundary nor transport "
                        "reactions")
        for name, s in stats.get("maps", {}).items():
            how = "as arranged here" if s.get("arranged_by") == "table" else "arranged by the pipeline"
            rows.append(f"{html.escape(name)}: {how}")
        box = QMessageBox(self)
        box.setWindowTitle("Maps built")
        box.setText(f"The maps are in {html.escape(str(out))}.<br><br>" + "<br>".join(rows))
        overview = box.addButton("Open overview", QMessageBox.ButtonRole.ActionRole)
        folder = box.addButton("Open folder", QMessageBox.ButtonRole.ActionRole)
        box.addButton(QMessageBox.StandardButton.Close)
        box.exec()
        if box.clickedButton() is overview:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(out / "overview.svg")))
        elif box.clickedButton() is folder:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(out)))
