"""The AI layout dock: a chat with an AI of the user's choice that lays out the current map with a toolkit.

The AI (a connection, :mod:`..ai.providers`; or a web chat,
:mod:`..ai.webchat`, carried by :mod:`.webchat`) works in a folder
OUT/ai/MAP/DATE-TIME with the map exported as Escher map and runs the
toolkit (:mod:`..ai.toolkits`; each run approved by the user unless allowed
for the job); what it submits becomes the map's AI layout (undoable, see
:meth:`.window.MainWindow.apply_ai_layout`). The model runs in a thread
(:class:`AgentWorker`); approvals and submissions are done in the GUI
thread (:class:`GuiCall`).
"""

from __future__ import annotations

import html
import json
import shlex
import time
import traceback
from pathlib import Path

from qtpy.QtCore import QObject, QSettings, QThread, QUrl, Qt, Signal
from qtpy.QtGui import QColor, QFont, QKeySequence, QShortcut, QTextCharFormat, QTextCursor
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QTextBrowser,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..ai import coding_agent
from ..ai.agent import Agent, Events, ToolError
from ..ai.providers import (
    ANTHROPIC,
    ANTHROPIC_MODELS,
    EFFORTS,
    OPENAI,
    PRESETS,
    WEBCHAT,
    Connection,
    ProviderError,
    list_models,
    make_chat,
)
from ..ai.toolkits import LayoutJob, Toolkit
from ..ai.webchat import BROWSER, CLIPBOARD, WebChat, check_connection
from ..ai_layouts import DRAFT, map_stubs
from ..split import slug
from .escher_svg import web_engine_available
from .webchat import BrowserWindow, ClipboardPanel

KEYRING = "ecmapper-draft-map"
FIRST_REQUEST = "Lay out this map."
GREY, RED, BLUE = QColor(110, 110, 110), QColor(180, 30, 30), QColor(30, 80, 170)


def keyring_module():
    try:
        import keyring
    except ImportError:
        return None
    return keyring


class GuiCall(QObject):
    """Run a function in the thread of this object (the GUI thread) and wait for it, from any thread."""

    request = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.request.connect(self._run, Qt.ConnectionType.BlockingQueuedConnection)

    def _run(self, box) -> None:
        fn, out = box
        try:
            out["value"] = fn()
        except BaseException as e:  # noqa: BLE001 (raised again in the calling thread)
            out["error"] = e

    def __call__(self, fn):
        if QThread.currentThread() is self.thread():
            return fn()
        box = (fn, {})
        self.request.emit(box)
        if "error" in box[1]:
            raise box[1]["error"]
        return box[1].get("value")


class AgentWorker(QThread):
    """Runs the agent with a message; reports what happens with signals. `reason`: why it stopped."""

    text = Signal(str)
    thinking = Signal(str)
    turn = Signal(object)
    tool_started = Signal(object)
    tool_finished = Signal(object, object)
    failed = Signal(str)

    def __init__(self, agent: Agent, message: str | None, approve, parent=None):
        super().__init__(parent)
        self.agent, self.message, self.approve = agent, message, approve
        self.stopped = False
        self.reason = None

    def stop(self) -> None:
        self.stopped = True

    def run(self) -> None:
        worker = self

        class WorkerEvents(Events):
            def text(self, delta):
                worker.text.emit(delta)

            def thinking(self, delta):
                worker.thinking.emit(delta)

            def turn(self, turn):
                worker.turn.emit(turn)

            def tool_started(self, call):
                worker.tool_started.emit(call)

            def tool_finished(self, call, result):
                worker.tool_finished.emit(call, result)

            def approve(self, call):
                return worker.approve(call)

            def cancelled(self):
                return worker.stopped

        try:
            self.reason = self.agent.run(self.message, WorkerEvents())
        except ProviderError as e:
            self.reason = "error"
            self.failed.emit(str(e))
        except Exception:  # noqa: BLE001 (shown to the user)
            self.reason = "error"
            self.failed.emit(traceback.format_exc())


class ConnectionDialog(QDialog):
    """Edit a connection: kind and URL (presets), model, API key, effort, images; a web chat: its URL and
    transport."""

    def __init__(self, parent, connection: Connection, key: str = ""):
        super().__init__(parent)
        self.setWindowTitle("AI connection")
        self.setMinimumWidth(620)
        form = QFormLayout(self)
        self.name = QLineEdit(connection.name)
        self.preset = QComboBox()
        self.preset.addItems(["(choose to fill in)", *PRESETS])
        self.preset.activated.connect(self.preset_chosen)
        self.kind = QComboBox()
        self.kind.addItem("Anthropic API", ANTHROPIC)
        self.kind.addItem("OpenAI-compatible API", OPENAI)
        self.kind.addItem("Web chat (browser or copy/paste)", WEBCHAT)
        self.kind.setCurrentIndex(max(self.kind.findData(connection.kind), 0))
        self.kind.currentIndexChanged.connect(self.kind_changed)
        self.base_url = QLineEdit(connection.base_url)
        self.url_label = QLabel("Base URL:")
        self.transport = QComboBox()
        self.transport.addItem("Browser in ecmapper (sends and reads the messages)", BROWSER)
        self.transport.addItem("Copy/paste by hand (any chat)", CLIPBOARD)
        self.transport.setCurrentIndex(max(self.transport.findData(connection.transport), 0))
        self.model = QComboBox()
        self.model.setEditable(True)
        self.models_button = QPushButton("List models")
        self.models_button.clicked.connect(self.fetch_models)
        model_box = QWidget()
        model_row = QHBoxLayout(model_box)
        model_row.setContentsMargins(0, 0, 0, 0)
        model_row.addWidget(self.model, 1)
        model_row.addWidget(self.models_button)
        self.key_env = QLineEdit(connection.key_env)
        self.key_env.setPlaceholderText("e.g. ANTHROPIC_API_KEY")
        self.key = QLineEdit(key)
        self.key.setEchoMode(QLineEdit.EchoMode.Password)
        stored = ("stored with keyring" if keyring_module() else
                  "kept for this session only (install keyring to store it)")
        self.key.setPlaceholderText(f"or the key itself ({stored})")
        self.effort = QComboBox()
        self.effort.addItems(EFFORTS)
        self.effort.setCurrentText(connection.effort)
        self.images = QCheckBox("Send images (renders of the toolkit) to the model")
        self.images.setChecked(connection.images)
        form.addRow("Name:", self.name)
        form.addRow("Preset:", self.preset)
        form.addRow("Kind:", self.kind)
        form.addRow(self.url_label, self.base_url)
        self.api_rows = [("Model:", model_box), ("Key from variable:", self.key_env), ("API key:", self.key),
                         ("Effort (Anthropic):", self.effort)]
        self.web_rows = [("Transport:", self.transport)]
        for label, widget in self.web_rows + self.api_rows:
            form.addRow(label, widget)
        form.addRow("", self.images)
        self.api_note = QLabel("Anthropic: without a key, the SDK's own credentials are used (ANTHROPIC_API_KEY, ant "
                               "auth). OpenAI-compatible: the model must support tool calls; local servers need no "
                               "key.")
        self.web_note = QLabel("A web chat gets the tools as text and calls them in its replies. Use the browser only "
                               "for chats whose terms allow automated use (e.g. your own or your institute's Open "
                               "WebUI); claude.ai, ChatGPT and Gemini do not, so ecmapper does not send messages "
                               "there: copy/paste by hand is fine everywhere. If the chat has an API key "
                               "(OpenAI-compatible), that connection is more reliable.")
        for note in (self.api_note, self.web_note):
            note.setWordWrap(True)
            form.addRow(note)
        self.form = form
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)
        self.kind_changed()
        self.model.setCurrentText(connection.model)

    def kind_changed(self) -> None:
        kind = self.kind.currentData()
        anthropic, web = kind == ANTHROPIC, kind == WEBCHAT
        current = self.model.currentText()
        self.model.clear()
        if anthropic:
            self.model.addItems(ANTHROPIC_MODELS)
        self.model.setCurrentText(current or (ANTHROPIC_MODELS[0] if anthropic else ""))
        self.effort.setEnabled(anthropic)
        for rows, shown in ((self.api_rows, not web), (self.web_rows, web)):
            for _, widget in rows:
                widget.setVisible(shown)
                self.form.labelForField(widget).setVisible(shown)
        self.api_note.setVisible(not web)
        self.web_note.setVisible(web)
        self.url_label.setText("Chat URL:" if web else "Base URL:")
        self.base_url.setPlaceholderText("the page of the chat (for the browser)" if web else
                                         "default of the API (e.g. http://localhost:11434/v1 for Ollama)")

    def preset_chosen(self, index: int) -> None:
        if index <= 0:
            return
        name = self.preset.currentText()
        kind, url, env, model = PRESETS[name]
        self.kind.setCurrentIndex(self.kind.findData(kind))
        self.base_url.setText(url)
        self.key_env.setText(env)
        self.model.setCurrentText(model)
        if not self.name.text().strip() or self.name.text() in PRESETS:
            self.name.setText(name)

    def connection(self) -> Connection:
        kind, url = self.kind.currentData(), self.base_url.text().strip()
        model = self.model.currentText().strip()
        if kind == WEBCHAT:
            host = QUrl.fromUserInput(url).host() if url else ""
            model = f"web chat {host}".strip()
        return Connection(name=self.name.text().strip() or "AI", kind=kind, model=model, base_url=url,
                          key_env=self.key_env.text().strip(), key=self.key.text().strip(),
                          effort=self.effort.currentText(), images=self.images.isChecked(),
                          transport=self.transport.currentData())

    def fetch_models(self) -> None:
        try:
            models = list_models(self.connection())
        except ProviderError as e:
            QMessageBox.warning(self, "List models", str(e))
            return
        current = self.model.currentText()
        self.model.clear()
        self.model.addItems(models)
        self.model.setCurrentText(current if current in models else (models[0] if models else ""))

    def accept(self) -> None:
        connection = self.connection()
        if connection.kind == WEBCHAT:
            try:
                check_connection(connection)
            except ProviderError as e:
                QMessageBox.warning(self, "AI connection", str(e))
                return
        elif not connection.model:
            QMessageBox.warning(self, "AI connection", "Give a model.")
            return
        super().accept()


class AIDock(QDockWidget):
    """The chat panel; `main`: the :class:`.window.MainWindow`."""

    def __init__(self, main, settings: QSettings):
        super().__init__("AI layout", main)
        self.setObjectName("AI layout")
        self.main, self.settings = main, settings
        self.keys: dict[str, str] = {}  # connection name -> key (this session)
        self.connections = self.load_connections()
        self.toolkits: list[Toolkit] = []
        for folder in json.loads(settings.value("ai/toolkits", "[]") or "[]"):
            try:
                self.toolkits.append(Toolkit.detect(folder))
            except (OSError, ValueError):
                pass
        self.job = None  # the current job: dict(map, folder, job, agent, ids, ...)
        self.worker = None
        self.browser = None  # the window of web chats (made when first used)
        self.gui_call = GuiCall(self)
        self.usage = dict(input_tokens=0, output_tokens=0, cache_read_input_tokens=0)
        self.in_answer = False
        self.make_widgets()

    # --- settings ------------------------------------------------------------------------------

    def load_connections(self) -> list[Connection]:
        try:
            saved = json.loads(self.settings.value("ai/connections", "") or "[]")
            connections = [Connection(**{k: v for k, v in c.items() if k in Connection.__dataclass_fields__})
                           for c in saved]
        except (ValueError, TypeError):
            connections = []
        return connections or [Connection("Anthropic")]

    def save_settings(self) -> None:
        self.settings.setValue("ai/connections", json.dumps([c.settings() for c in self.connections]))
        self.settings.setValue("ai/connection", self.connection_box.currentText())
        self.settings.setValue("ai/toolkits", json.dumps([str(t.folder) for t in self.toolkits]))
        self.settings.setValue("ai/toolkit", self.toolkit_box.currentIndex())
        self.settings.setValue("ai/approve_all", self.approve_all.isChecked())

    def key_of(self, connection: Connection) -> str:
        if connection.name in self.keys:
            return self.keys[connection.name]
        kr = keyring_module()
        if kr is not None:
            try:
                return kr.get_password(KEYRING, connection.name) or ""
            except Exception:  # noqa: BLE001 (no keyring backend)
                return ""
        return ""

    def store_key(self, name: str, key: str) -> None:
        self.keys[name] = key
        kr = keyring_module()
        if kr is not None and key:
            try:
                kr.set_password(KEYRING, name, key)
            except Exception:  # noqa: BLE001 (no keyring backend: kept for the session)
                pass

    # --- widgets -------------------------------------------------------------------------------

    def make_widgets(self) -> None:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        row = QHBoxLayout()
        row.addWidget(QLabel("AI:"))
        self.connection_box = QComboBox()
        self.connection_box.setToolTip("The AI the chat goes to")
        row.addWidget(self.connection_box, 1)
        edit = QToolButton()
        edit.setText("⚙")
        edit.setToolTip("Edit this connection (model, key, …)")
        edit.clicked.connect(lambda: self.edit_connection())
        add = QToolButton()
        add.setText("+")
        add.setToolTip("New connection")
        add.clicked.connect(lambda: self.new_connection())
        row.addWidget(edit)
        row.addWidget(add)
        layout.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel("Toolkit:"))
        self.toolkit_box = QComboBox()
        self.toolkit_box.setToolTip("The layout toolkit the AI works with (a folder, e.g. escher_handlayout_toolkit)")
        row.addWidget(self.toolkit_box, 1)
        add_toolkit = QToolButton()
        add_toolkit.setText("…")
        add_toolkit.setToolTip("Add a toolkit folder")
        add_toolkit.clicked.connect(self.add_toolkit)
        row.addWidget(add_toolkit)
        layout.addLayout(row)
        self.target = QLabel()
        self.target.setWordWrap(True)
        layout.addWidget(self.target)
        self.transcript = QTextBrowser()
        self.transcript.setOpenExternalLinks(True)
        layout.addWidget(self.transcript, 1)
        self.clipboard_panel = ClipboardPanel(self.gui_call)  # web chats by copy/paste
        layout.addWidget(self.clipboard_panel)
        self.input = QPlainTextEdit()
        self.input.setPlaceholderText("Layout instructions for this map, e.g. “glycolysis as the central column, the "
                                      "citric acid cycle below it as a ring” (Ctrl+Enter sends)")
        self.input.setMaximumHeight(110)
        layout.addWidget(self.input)
        QShortcut(QKeySequence("Ctrl+Return"), self.input, activated=self.send)
        row = QHBoxLayout()
        self.send_button = QPushButton("Send")
        self.send_button.clicked.connect(self.send)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop)
        self.update_button = QPushButton("Update with AI")
        self.update_button.setToolTip("The map changed since it was laid out: let the AI update its layout")
        self.update_button.clicked.connect(self.update_layout)
        self.restart_button = QPushButton("New conversation")
        self.restart_button.setToolTip("Continue the job in its work folder with a fresh context")
        self.restart_button.clicked.connect(self.new_conversation)
        for b in (self.send_button, self.stop_button, self.update_button, self.restart_button):
            row.addWidget(b)
        row.addStretch(1)
        layout.addLayout(row)
        self.approve_all = QCheckBox("Run toolkit commands without asking")
        self.approve_all.setToolTip("The toolkit runs Python code the AI wrote (its spec), with your rights")
        self.approve_all.setChecked(str(self.settings.value("ai/approve_all", "false")).lower() == "true")
        layout.addWidget(self.approve_all)
        self.status = QLabel()
        self.status.setStyleSheet("color: gray")
        layout.addWidget(self.status)
        self.setWidget(panel)
        self.fill_connections(self.settings.value("ai/connection", ""))
        self.fill_toolkits(int(self.settings.value("ai/toolkit", 0) or 0))
        self.connection_box.currentIndexChanged.connect(lambda *_: self.save_settings())
        self.toolkit_box.currentIndexChanged.connect(lambda *_: self.save_settings())
        self.approve_all.toggled.connect(lambda *_: self.save_settings())
        self.update_state()

    def fill_connections(self, current: str = "") -> None:
        self.connection_box.blockSignals(True)
        self.connection_box.clear()
        for c in self.connections:
            self.connection_box.addItem(c.name)
            self.connection_box.setItemData(self.connection_box.count() - 1, f"{c.kind}: {c.model}",
                                            Qt.ItemDataRole.ToolTipRole)
        names = [c.name for c in self.connections]
        self.connection_box.setCurrentIndex(names.index(current) if current in names else 0)
        self.connection_box.blockSignals(False)

    def fill_toolkits(self, current: int = 0) -> None:
        self.toolkit_box.blockSignals(True)
        self.toolkit_box.clear()
        for t in self.toolkits:
            self.toolkit_box.addItem(f"{t.name} ({t.folder.name})")
            self.toolkit_box.setItemData(self.toolkit_box.count() - 1, str(t.folder), Qt.ItemDataRole.ToolTipRole)
        if not self.toolkits:
            self.toolkit_box.addItem("(no toolkit: add a folder with …)")
        self.toolkit_box.setCurrentIndex(min(max(current, 0), self.toolkit_box.count() - 1))
        self.toolkit_box.blockSignals(False)

    def current_connection(self) -> Connection:
        c = self.connections[max(self.connection_box.currentIndex(), 0)]
        c.key = self.key_of(c)
        return c

    def current_toolkit(self) -> Toolkit | None:
        i = self.toolkit_box.currentIndex()
        return self.toolkits[i] if 0 <= i < len(self.toolkits) else None

    def edit_connection(self, connection: Connection | None = None) -> None:
        index = self.connection_box.currentIndex()
        old = connection or self.connections[index]
        dialog = ConnectionDialog(self, old, self.key_of(old))
        if not dialog.exec():
            return
        new = dialog.connection()
        if new.key:
            self.store_key(new.name, new.key)
        new.key = ""
        if connection is None:
            self.connections[index] = new
        else:
            self.connections.append(new)
        self.fill_connections(new.name)
        self.save_settings()

    def new_connection(self) -> None:
        self.edit_connection(Connection("New connection", model=""))

    def add_toolkit(self) -> None:
        start = str(self.toolkits[-1].folder.parent) if self.toolkits else ""
        folder = QFileDialog.getExistingDirectory(self, "Toolkit folder (e.g. escher_handlayout_toolkit)", start)
        if not folder:
            return
        try:
            toolkit = Toolkit.detect(folder)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Toolkit", str(e))
            return
        self.toolkits = [t for t in self.toolkits if t.folder != toolkit.folder] + [toolkit]
        self.fill_toolkits(len(self.toolkits) - 1)
        self.save_settings()
        commands = ", ".join(toolkit.commands) or "any"
        self.note(f"Toolkit {toolkit.name} in {toolkit.folder}: runs python -m {toolkit.package} ({commands}); "
                  f"instructions: {', '.join(toolkit.instructions)}")

    # --- state ---------------------------------------------------------------------------------

    def update_state(self) -> None:
        """The target (the current map) and the buttons."""
        name = self.main.current_map() if self.main.session.dm is not None else None
        running = self.worker is not None
        session = self.main.session
        if name is None:
            self.target.setText("<i>Choose a map tab to lay it out.</i>")
        else:
            info = session.layout_info(name)
            state = ("drawn automatically" if info is None else
                     "<b>AI layout, outdated</b> (Build draws it automatically)" if info[2] else "AI layout")
            job = " — this chat's map" if self.job is not None and self.job["map"] == name else ""
            self.target.setText(f"Map: <b>{html.escape(name)}</b> ({state}){job}")
        self.send_button.setEnabled(not running and name is not None)
        self.stop_button.setEnabled(running)
        self.restart_button.setEnabled(not running and self.job is not None)
        self.update_button.setVisible(name is not None and session.outdated(name))
        self.update_button.setEnabled(not running)
        for w in (self.connection_box, self.toolkit_box):
            w.setEnabled(not running)

    # --- transcript ----------------------------------------------------------------------------

    def _insert(self, text: str, colour: QColor | None = None, bold=False, italic=False, new_block=False) -> None:
        cursor = self.transcript.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        if new_block and not cursor.atStart():
            cursor.insertBlock()
        fmt = QTextCharFormat()
        if colour is not None:
            fmt.setForeground(colour)
        if bold:
            fmt.setFontWeight(QFont.Weight.Bold)
        fmt.setFontItalic(italic)
        cursor.insertText(text, fmt)
        self.transcript.setTextCursor(cursor)
        self.transcript.ensureCursorVisible()

    def note(self, text: str, colour: QColor = GREY) -> None:
        self.in_answer = False
        self._insert(text, colour, italic=True, new_block=True)

    def show_user(self, text: str) -> None:
        self.in_answer = False
        self._insert("You: ", BLUE, bold=True, new_block=True)
        self._insert(text)

    def on_text(self, delta: str) -> None:
        if not self.in_answer:
            self._insert("AI: ", bold=True, new_block=True)
            self.in_answer = True
        self._insert(delta)

    def on_thinking(self, delta: str) -> None:
        if self.in_answer:
            self.in_answer = False
            self._insert("", new_block=True)
        self._insert(delta, GREY, italic=True)

    def on_tool_started(self, call) -> None:
        self.in_answer = False
        self._insert(f"▸ {describe(call)}", GREY, new_block=True)

    def on_tool_finished(self, call, result) -> None:
        if result.is_error:
            first = (result.text.strip().splitlines() or [""])[0][:160]
            if not self.transcript.toPlainText().rstrip().endswith(describe(call)):
                self._insert(f"▸ {describe(call)}", GREY, new_block=True)
            self._insert(f"  ✗ {first}", RED)
        else:
            self._insert("  ✓", GREY)

    def on_turn(self, turn) -> None:
        self.in_answer = False
        for k in self.usage:
            self.usage[k] += turn.usage.get(k, 0) or 0
        u = self.usage
        chat = self.job["agent"].chat if self.job is not None else None
        if isinstance(chat, WebChat):
            self.status.setText(f"Web chat: {chat.number} messages")
            return
        cached = f", {u['cache_read_input_tokens'] / 1000:.1f}k from the cache" if u["cache_read_input_tokens"] else ""
        self.status.setText(f"Tokens of this chat: {u['input_tokens'] / 1000:.1f}k in{cached}, "
                            f"{u['output_tokens'] / 1000:.1f}k out")

    # --- jobs ----------------------------------------------------------------------------------

    def start_job(self, name: str) -> dict:
        """A new job for map `name`: export it into a new work folder; raises ValueError if it can not start."""
        toolkit = self.current_toolkit()
        if toolkit is None:
            raise ValueError("Add a toolkit folder first (…), e.g. escher_handlayout_toolkit.")
        session = self.main.session
        if name not in session.views():
            raise ValueError(f"The map {name!r} has no subsystems to lay out.")
        folder = Path(session.args.out) / "ai" / slug(name) / time.strftime("%Y%m%d-%H%M%S")
        _, ids = session.export_map(name, folder)
        subsystems, stubs, facts = self.map_facts(name)
        job = dict(map=name, folder=folder, toolkit=toolkit, ids=ids, drafts=1, key=object(), approve_all=False)
        job["job"] = LayoutJob(toolkit, folder, name, facts, subsystems, stubs, self.submit)
        self.job = job
        try:
            self.new_chat()
        except ProviderError:
            self.job = None
            raise
        self.note(f"New job: {name} with {toolkit.name}; work folder {folder}")
        return job

    def map_facts(self, name: str) -> tuple[list[tuple[str, str, int]], dict[str, list[str]], str]:
        """The subsystems of map `name` (ID, name, reactions), its metabolites that continue on other maps and a
        summary, for the AI's task."""
        session = self.main.session
        sub, on_map = session.views()[name]
        subsystems = [(k, session.dm.subsystems[k].label, session.counts.get(k, 0)) for k in on_map]
        stubs = {}
        for node, maps in map_stubs(sub).items():
            stubs.setdefault(sub.nodes[node].label, maps)
        reactions = sum(n for *_, n in subsystems)
        return subsystems, stubs, f"{len(subsystems)} subsystems, {reactions} reactions"

    def prepare_agent_folder(self, name: str) -> Path:
        """A work folder for a coding agent on map `name` (:func:`..ai.coding_agent.prepare`); ValueError if it can
        not be made."""
        toolkit = self.current_toolkit()
        if toolkit is None:
            raise ValueError("Add a toolkit folder first (…), e.g. escher_handlayout_toolkit.")
        session = self.main.session
        if name not in session.views():
            raise ValueError(f"The map {name!r} has no subsystems to lay out.")
        folder = Path(session.args.out) / "ai" / slug(name) / (time.strftime("%Y%m%d-%H%M%S") + "-agent")
        _, ids = session.export_map(name, folder)
        subsystems, stubs, facts = self.map_facts(name)
        return coding_agent.prepare(folder, toolkit, name, facts, subsystems, stubs, ids)

    def new_chat(self) -> None:
        """A fresh conversation (connection as chosen now) for the current job."""
        job = self.job
        connection = self.current_connection()
        tools = job["job"].tools()
        if connection.kind == WEBCHAT:
            chat = WebChat(connection, job["job"].system_prompt(), tools, self.web_transport(connection))
        else:
            chat = make_chat(connection, job["job"].system_prompt(), tools)
        job["agent"] = Agent(chat, tools, transcript=job["folder"] / "transcript.jsonl")
        job["connection"] = connection
        self.usage = dict.fromkeys(self.usage, 0)

    def web_transport(self, connection: Connection):
        """The transport of a web chat: the browser window (opened at the chat's page) or the copy/paste panel."""
        check_connection(connection)
        if connection.transport != BROWSER:
            return self.clipboard_panel
        if self.browser is None:
            if not web_engine_available():
                raise ProviderError("The browser needs QtWebEngine (e.g. pip install PyQt6-WebEngine); choose "
                                    "copy/paste in the connection (⚙) instead.")
            self.browser = BrowserWindow(self.gui_call, self)
        self.browser.open(connection.base_url)
        self.note("Web chat: log in in the browser window, open a new conversation there and click Start.")
        return self.browser

    def send(self) -> None:
        if self.worker is not None:
            return
        name = self.main.current_map()
        if name is None:
            return
        text = self.input.toPlainText().strip()
        try:
            if self.job is None or self.job["map"] != name or self.job["toolkit"] is not self.current_toolkit():
                self.start_job(name)
                message = text or FIRST_REQUEST
            elif self.job["connection"].name != self.connection_box.currentText():
                self.new_chat()
                message = resume(text)
            else:
                message = text or "Go on."
        except (ValueError, ProviderError, OSError) as e:
            QMessageBox.warning(self, "AI layout", str(e))
            return
        self.input.clear()
        self.show_user(text or message)
        self.run(message)

    def new_conversation(self) -> None:
        if self.job is None or self.worker is not None:
            return
        try:
            self.new_chat()
        except ProviderError as e:
            QMessageBox.warning(self, "AI layout", str(e))
            return
        text = self.input.toPlainText().strip()
        self.input.clear()
        self.note("New conversation (the files of the work folder are kept)")
        self.show_user(text or "(continue)")
        self.run(resume(text))

    def update_layout(self) -> None:
        """The map changed since it was laid out: give the AI the map as it is now."""
        name = self.main.current_map()
        session = self.main.session
        info = session.layout_info(name)
        if self.worker is not None or info is None:
            return
        problems = "; ".join(info[2])
        try:
            if self.job is not None and self.job["map"] == name:
                self.job["drafts"] += 1
                file = f"draft_v{self.job['drafts']}.escher.json"
                _, self.job["ids"] = session.export_map(name, self.job["folder"], file, reuse=self.job["ids"])
                message = (f"The map changed since you laid it out: {problems}. It is now in {file}: the objects that "
                           f"were in {DRAFT} have the same IDs there. Update your layout for it (keep what still fits) "
                           "and submit it again.")
            else:
                self.start_job(name)
                message = (f"This map was laid out before, but it has changed: {problems}. Lay it out (the user may "
                           "tell you how).")
        except (ValueError, ProviderError, OSError) as e:
            QMessageBox.warning(self, "AI layout", str(e))
            return
        self.show_user(message)
        self.run(message)

    def run(self, message: str | None) -> None:
        worker = AgentWorker(self.job["agent"], message, lambda call: self.gui_call(lambda: self.approve(call)), self)
        worker.text.connect(self.on_text)
        worker.thinking.connect(self.on_thinking)
        worker.turn.connect(self.on_turn)
        worker.tool_started.connect(self.on_tool_started)
        worker.tool_finished.connect(self.on_tool_finished)
        worker.failed.connect(lambda text: self.note(text, RED))
        worker.finished.connect(self.finished)
        self.worker = worker
        self.update_state()
        worker.start()

    def stop(self) -> None:
        if self.worker is not None:
            self.worker.stop()
            self.note("Stopping …")

    def finished(self) -> None:
        reason, self.worker = self.worker.reason, None
        texts = {"stopped": "Stopped.", "max_steps": "Stopped after the maximum number of steps; send a message to "
                 "go on.", "max_tokens": "The answer reached the output limit."}
        if reason and reason.startswith("refusal"):
            self.note(f"The model declined ({reason}).", RED)
        elif reason in texts:
            self.note(texts[reason])
        self.update_state()

    # --- called by the agent (GUI thread) -------------------------------------------------------

    def approve(self, call) -> bool:
        job = self.job
        if self.approve_all.isChecked() or job["approve_all"]:
            return True
        args = " ".join(shlex.quote(str(a)) for a in call.input.get("args", []))
        box = QMessageBox(self)
        box.setWindowTitle("Run the toolkit?")
        box.setIcon(QMessageBox.Icon.Question)
        box.setText(f"The AI wants to run\n\npython -m {job['toolkit'].package} {args}\n\nin {job['folder']}.")
        box.setInformativeText("The toolkit runs Python code the AI wrote (its spec files), with your rights.")
        run = box.addButton("Run", QMessageBox.ButtonRole.AcceptRole)
        run_all = box.addButton("Run all of this job", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Decline", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        if box.clickedButton() is run_all:
            job["approve_all"] = True
        return box.clickedButton() in (run, run_all)

    def submit(self, path: Path) -> str:
        """submit_layout of the job (called in the agent's thread)."""
        return self.gui_call(lambda: self._submit(path))

    def _submit(self, path: Path) -> str:
        job = self.job
        try:
            final = json.loads(Path(path).read_text(encoding="utf-8"))
            final[1]["nodes"], final[1]["reactions"]  # noqa: B018 (an Escher map)
        except (OSError, ValueError, KeyError, IndexError, TypeError) as e:
            raise ToolError(f"{path.name}: not an Escher map ({e})") from None
        run_log = None
        if (path.parent / "run_log.json").is_file():
            try:
                run_log = json.loads((path.parent / "run_log.json").read_text(encoding="utf-8"))
            except ValueError:
                run_log = None
        connection = job["connection"]
        source = dict(toolkit=job["toolkit"].name, connection=connection.name, model=connection.model,
                      date=time.strftime("%Y-%m-%d %H:%M"), folder=str(job["folder"]), file=path.name)
        session = self.main.session
        if job["map"] not in session.views():
            raise ToolError(f"the map {job['map']!r} is gone (the user changed the maps)")
        layout, problems = session.new_layout(job["map"], final, job["ids"], None, run_log, source)
        if problems:
            raise ToolError("ecmapper did not accept the map: " + "; ".join(problems) + ". The reactions and their "
                            "stoichiometry must be those of the draft; fix it and submit again.")
        layout.svg = self.main.export_svg(layout.escher)  # with the "→ map" labels of the map now
        self.main.apply_ai_layout(job["map"], layout, job["key"])
        self.update_state()
        x, y, w, h = layout.extent()
        drawn = sum(layout.reactions.values())
        missing = sum(1 for i in session.views()[job["map"]][0].nodes if i not in layout.positions)
        note = "" if layout.svg is not None else " (no SVG: Escher could not be run here; the map is shown as a note)"
        self.note(f"Layout of {job['map']} from {path.name} shown{note}.")
        return (f"Accepted and shown to the user: {drawn} reactions, {w:.0f} × {h:.0f} map units"
                + (f"; {missing} draft nodes not found in the map (placed next to their neighbours)" if missing else "")
                + note + ".")


def describe(call) -> str:
    """A short description of a tool call."""
    if call.name == "run_toolkit":
        return "toolkit " + " ".join(str(a) for a in call.input.get("args", []))
    arg = call.input.get("path")
    return f"{call.name} {arg}" if arg else call.name


def resume(text: str) -> str:
    message = ("Your conversation was reset; the files of the work folder are kept. Look at them (notes, spec files, "
               "the newest run) and go on with the job.")
    return message + (f"\n\nThe user: {text}" if text else "")
