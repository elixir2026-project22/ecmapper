"""The transports of web chats (:mod:`..ai.webchat`): a browser window that sends ecmapper's messages to a chat and
reads the replies (:class:`BrowserWindow`, needs QtWebEngine), and a panel of the dock for copying them by hand
(:class:`ClipboardPanel`).

The browser keeps its logins (the profile :data:`PROFILE` in the application's
data folder). It pastes a message into the chat's input field (the lowest
visible text field of the page), as the user would, presses Enter (or clicks
a send button if that did not send it) and reads the text of the page every
second until the reply is complete. It sends nothing to chats whose terms
forbid automated use (:func:`..ai.webchat.blocked`).
"""

from __future__ import annotations

import functools
import threading
import time

from qtpy.QtCore import QEvent, QEventLoop, QMimeData, QThread, QTimer, QUrl, Qt
from qtpy.QtGui import QImage, QKeyEvent
from qtpy.QtWidgets import (
    QApplication,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..ai.providers import ProviderError
from ..ai.webchat import blocked, end_line, reply_of, tail

PROFILE = "ecmapper-webchat"
#: focus the chat's input field (the lowest visible text field); false if there is none
FOCUS_JS = """(function () {
  const shown = e => { const r = e.getBoundingClientRect(), s = getComputedStyle(e);
                       return r.width > 20 && r.height > 8 && s.visibility !== 'hidden' && s.display !== 'none'; };
  const fields = [...document.querySelectorAll('textarea, [contenteditable="true"], [contenteditable=""]')]
    .filter(e => shown(e) && !e.disabled && !e.readOnly);
  if (!fields.length) return false;
  fields.sort((a, b) => b.getBoundingClientRect().bottom - a.getBoundingClientRect().bottom);
  fields[0].focus();
  return true;
})()"""
#: click the send button nearest to the focused field; false if there is none
SEND_JS = """(function () {
  const field = document.activeElement, at = field ? field.getBoundingClientRect() : null;
  const name = b => [b.getAttribute('aria-label'), b.title, b.type === 'submit' ? 'submit' : '',
                     (b.textContent || '').trim().slice(0, 20)].join(' ');
  const buttons = [...document.querySelectorAll('button, [role="button"]')]
    .filter(b => !b.disabled && /send|submit|senden/i.test(name(b)));
  if (!buttons.length) return false;
  const away = b => { if (!at) return 0; const r = b.getBoundingClientRect();
                      return Math.hypot(r.left - at.right, r.top - at.top); };
  buttons.sort((a, b) => away(a) - away(b));
  buttons[0].click();
  return true;
})()"""
TEXT_JS = "document.body ? document.body.innerText : ''"


def copy_mime(data: QMimeData | None) -> QMimeData | None:
    """A copy of the clipboard's data (it is restored after a paste); None: the clipboard is empty."""
    if data is None:
        return None
    copy = QMimeData()
    for f in data.formats():
        copy.setData(f, data.data(f))
    return copy


def guarded(method):
    """A slot or callback of the browser window: an error ends the exchange (not the program, as PyQt would)."""

    @functools.wraps(method)
    def run(self, *args):
        try:
            return method(self, *args)
        except Exception as e:  # noqa: BLE001 (reported to the agent's thread as the error of the exchange)
            self.finish(error=f"the browser failed: {type(e).__name__}: {e}")

    return run


class Transport:
    """What the transports share: :meth:`exchange` (from any thread) runs `begin` in the GUI thread (by
    `gui_call`) and waits until it ends with :meth:`finish`, or calls `cancel` when `stop()` is true."""

    images = False
    reply = error = current = None

    def exchange(self, message: str, number: int, images: list[bytes], stop=None) -> str | None:
        done = threading.Event()
        self.gui_call(lambda: self.begin(message, number, list(images), done))
        in_gui = QThread.currentThread() is QApplication.instance().thread()
        while not done.is_set():
            if stop is not None and stop():
                self.gui_call(self.cancel)
                return None
            if in_gui:
                QApplication.processEvents(QEventLoop.ProcessEventsFlag.AllEvents, 50)
                time.sleep(0.01)
            else:
                done.wait(0.2)
        if self.error:
            raise ProviderError(self.error)
        return self.reply

    def finish(self, reply: str | None = None, error: str | None = None) -> None:
        current, self.current = self.current, None
        self.reply, self.error = reply, error
        self.ended(current)
        if current is not None:
            current["done"].set()

    def ended(self, current: dict | None) -> None:
        """The exchange `current` ended (or was stopped: None)."""


class BrowserWindow(QWidget, Transport):
    """The web chat in a browser window (see the module); `gui_call` runs a function in the GUI thread."""

    images = True

    def __init__(self, gui_call, parent=None):
        super().__init__(parent, Qt.WindowType.Window)
        from qtpy.QtWebEngineCore import QWebEnginePage, QWebEngineProfile
        from qtpy.QtWebEngineWidgets import QWebEngineView

        self.gui_call = gui_call
        self.paste_action = QWebEnginePage.WebAction.Paste
        self.setWindowTitle("Web chat (ecmapper)")
        self.resize(1000, 850)
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        back, reload = QToolButton(), QToolButton()
        back.setText("◀")
        back.setToolTip("Back")
        reload.setText("⟳")
        reload.setToolTip("Reload")
        self.address = QLineEdit()
        self.address.returnPressed.connect(lambda: self.open(self.address.text(), again=True))
        for w in (back, reload, self.address):
            row.addWidget(w)
        layout.addLayout(row)
        self.view = QWebEngineView(self)
        self.profile = QWebEngineProfile(PROFILE, self)  # after the view: deleted after its page
        self.profile.setPersistentCookiesPolicy(QWebEngineProfile.PersistentCookiesPolicy.ForcePersistentCookies)
        self.view.setPage(make_page(self.profile, self.view, self))
        back.clicked.connect(self.view.back)
        reload.clicked.connect(self.view.reload)
        self.view.urlChanged.connect(self.url_changed)
        layout.addWidget(self.view, 1)
        row = QHBoxLayout()
        self.status = QLabel("Log in to the chat and open a new conversation.")
        self.status.setWordWrap(True)
        self.start_button = QPushButton("Start")
        self.start_button.setToolTip("You are logged in and a new conversation is open: send ecmapper's message")
        self.start_button.clicked.connect(self.start)
        self.take_button = QPushButton("Use the reply as it is")
        self.take_button.setToolTip("The reply is complete, but the model did not end it with the end line")
        self.take_button.clicked.connect(lambda: self.take())
        row.addWidget(self.status, 1)
        row.addWidget(self.start_button)
        row.addWidget(self.take_button)
        layout.addLayout(row)
        self.timer = QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.poll)
        self.started = False  # the user clicked Start (logged in, a conversation open); again for each message 1
        self.saved = None  # [the clipboard's data before a paste (None: empty)]
        self.update_buttons()

    def open(self, url: str, again: bool = False) -> None:
        """Show the window with the page `url` unless a page of its site is open (`again`: load it in any case)."""
        target = QUrl.fromUserInput(url.strip())
        if again or self.view.url().isEmpty() or self.view.url().host().lower() != target.host().lower():
            self.view.load(target)
        self.show()
        self.raise_()
        self.activateWindow()

    def say(self, text: str) -> None:
        self.status.setText(text)

    def url_changed(self, url) -> None:
        self.address.setText(url.toString())
        reason = blocked(url.toString())
        if reason:
            self.say(f"ecmapper sends nothing here: {reason}")

    def update_buttons(self) -> None:
        cur = self.current
        self.start_button.setEnabled(cur is not None and cur["sent"] is None and not cur["pasting"])
        self.take_button.setEnabled(cur is not None and cur["sent"] is not None)

    # --- an exchange (GUI thread) --------------------------------------------------------------------------------

    def begin(self, message: str, number: int, images: list[bytes], done) -> None:
        self.reply = self.error = None
        self.current = dict(message=message, number=number, images=images, done=done, sent=None, pasting=False,
                            text=None, changed=time.time(), clicked=False)
        if number == 1:
            self.started = False
        self.show()
        self.raise_()
        if self.started:
            self.send()
        else:
            self.say(f"Log in to the chat, open a new (empty) conversation and click Start: ecmapper then sends "
                     f"message {number} and reads the replies." if number == 1 else
                     f"Open the conversation of this job and click Start to send message {number}.")
        self.update_buttons()

    def start(self) -> None:
        self.started = True
        if self.current is not None and self.current["sent"] is None:
            self.send()

    @guarded
    def send(self) -> None:
        reason = blocked(self.view.url().toString())
        if reason:
            self.finish(error=reason)
            return
        self.current["pasting"] = True
        self.update_buttons()
        self.view.setFocus()
        self.view.page().runJavaScript(FOCUS_JS, self.focused)

    @guarded
    def focused(self, found) -> None:
        cur = self.current
        if cur is None:
            return
        if not found:
            cur["pasting"] = self.started = False
            self.say("No input field found on this page: open the chat (logged in), then click Start.")
            self.update_buttons()
            return
        self.saved = [copy_mime(QApplication.clipboard().mimeData())]
        self.paste_next()

    @guarded
    def paste_next(self) -> None:
        """Paste the images of the message, then its text, then send it."""
        cur = self.current
        if cur is None:
            return
        clipboard = QApplication.clipboard()
        if cur["images"]:
            clipboard.setImage(QImage.fromData(cur["images"].pop(0)))
            self.view.page().triggerAction(self.paste_action)
            QTimer.singleShot(1500, self.paste_next)  # the chat uploads it
            return
        clipboard.setText(cur["message"])
        self.view.page().triggerAction(self.paste_action)
        QTimer.singleShot(700, self.press_enter)

    @guarded
    def press_enter(self) -> None:
        cur = self.current
        if self.saved is not None:
            (data,), self.saved = self.saved, None
            if data is None:
                QApplication.clipboard().clear()
            else:
                QApplication.clipboard().setMimeData(data)
        if cur is None:
            return
        target = self.view.focusProxy() or self.view
        for kind in (QEvent.Type.KeyPress, QEvent.Type.KeyRelease):
            QApplication.sendEvent(target, QKeyEvent(kind, Qt.Key.Key_Return, Qt.KeyboardModifier.NoModifier, "\r"))
        cur["sent"], cur["pasting"] = time.time(), False
        self.say(f"Message {cur['number']} sent; waiting for the reply …")
        self.update_buttons()
        self.timer.start()

    @guarded
    def poll(self) -> None:
        if self.current is not None and self.current["sent"] is not None:
            self.view.page().runJavaScript(TEXT_JS, self.read)

    @guarded
    def read(self, text) -> None:
        cur = self.current
        if cur is None or not isinstance(text, str):
            return
        now, n = time.time(), cur["number"]
        if text != cur["text"]:
            cur["text"], cur["changed"] = text, now
        reply, complete = reply_of(text, n)
        if complete and now - cur["changed"] >= 2:
            self.finish(reply=reply)
        elif tail(n) not in text and not cur["clicked"] and now - cur["sent"] > 10:
            cur["clicked"] = True  # Enter did not send it (or the chat shows it as an attachment): its send button
            self.view.page().runJavaScript(SEND_JS)
            self.say(f"Message {n}: clicked the send button; if the message is still in the input field, send it "
                     "yourself.")
        elif not complete and reply and now - cur["changed"] > 60:
            self.say(f"The reply to message {n} has not ended with {end_line(n)} yet: if it is complete, click "
                     "“Use the reply as it is”.")

    @guarded
    def take(self) -> None:
        cur = self.current
        if cur is not None and cur["text"]:
            reply, _ = reply_of(cur["text"], cur["number"])
            if reply:
                self.finish(reply=reply)

    def cancel(self) -> None:
        self.current = None
        self.ended(None)

    def ended(self, current: dict | None) -> None:
        self.timer.stop()
        if self.error:
            self.say(f"Stopped: {self.error}")
        elif current is None:
            self.say("Stopped.")
        else:
            self.say(f"Reply to message {current['number']} read; ecmapper works on it.")
        self.update_buttons()

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt)
        self.started = False
        if self.current is not None:
            self.finish(error="the browser window was closed")
        super().closeEvent(event)


_PAGE = None


def make_page(profile, parent, owner):
    """A page of `profile` that opens the chat's pop-ups (e.g. for a login) in windows of `owner`."""
    global _PAGE
    if _PAGE is None:
        from qtpy.QtWebEngineCore import QWebEnginePage
        from qtpy.QtWebEngineWidgets import QWebEngineView

        class Page(QWebEnginePage):
            def __init__(self, profile, parent, owner):
                super().__init__(profile, parent)
                self.owner = owner

            def createWindow(self, kind):  # noqa: N802 (Qt)
                view = QWebEngineView(self.owner)
                view.setWindowFlags(Qt.WindowType.Window)
                view.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
                page = Page(self.profile(), view, self.owner)
                page.windowCloseRequested.connect(view.close)
                view.setPage(page)
                view.resize(700, 800)
                view.show()
                return page

        _PAGE = Page
    return _PAGE(profile, parent, owner)


class ClipboardPanel(QFrame, Transport):
    """A panel of the dock for a web chat by hand: the user copies ecmapper's messages into the chat and the
    model's replies back; `gui_call` runs a function in the GUI thread."""

    def __init__(self, gui_call, parent=None):
        super().__init__(parent)
        self.gui_call = gui_call
        self.setFrameShape(QFrame.Shape.StyledPanel)
        layout = QVBoxLayout(self)
        self.label = QLabel()
        self.label.setWordWrap(True)
        layout.addWidget(self.label)
        row = QHBoxLayout()
        self.copy_button = QPushButton("Copy message")
        self.copy_button.clicked.connect(self.copy)
        self.paste_button = QPushButton("Paste reply")
        self.paste_button.setToolTip("Take the reply of the model from the clipboard")
        self.paste_button.clicked.connect(self.paste)
        row.addWidget(self.copy_button)
        row.addWidget(self.paste_button)
        row.addStretch(1)
        layout.addLayout(row)
        self.hide()

    def begin(self, message: str, number: int, images: list[bytes], done) -> None:
        self.reply = self.error = None
        self.current = dict(message=message, number=number, done=done)
        self.copy_button.setText(f"Copy message {number}")
        new = " in a new conversation" if number == 1 else ""
        self.label.setText(f"<b>Web chat:</b> message {number} is on the clipboard. Paste it into the chat{new} and "
                           "send it; when the model is done, copy its whole reply and click <i>Paste reply</i>.")
        self.copy()
        self.show()

    def copy(self) -> None:
        if self.current is not None:
            QApplication.clipboard().setText(self.current["message"])

    def paste(self) -> None:
        cur = self.current
        if cur is None:
            return
        text = QApplication.clipboard().text()
        reply, _ = reply_of(text, cur["number"])
        if not text.strip() or text.strip() == cur["message"].strip() or not reply:
            QMessageBox.information(self, "Web chat", "The clipboard does not hold a reply of the model: copy its "
                                    "whole reply (e.g. with the chat's copy button), then click Paste reply.")
            return
        self.finish(reply=reply)

    def cancel(self) -> None:
        self.current = None
        self.ended(None)

    def ended(self, current: dict | None) -> None:
        self.hide()
