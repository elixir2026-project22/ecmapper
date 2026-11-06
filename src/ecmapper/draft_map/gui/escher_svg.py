"""SVG export of Escher maps by Escher itself (escher.min.js in a hidden QWebEngineView), for the AI layouts.

escher.min.js is taken from (in this order) the setting `ai/escher_js`, the
escher Python package, a copy in the application's data folder, or
downloaded once (:data:`ESCHER_URL`, after asking). Without QtWebEngine
(e.g. PyQt6-WebEngine not installed) there is no export.
"""

from __future__ import annotations

import json
import os
import tempfile
import time
import urllib.request
from pathlib import Path

from qtpy.QtCore import QEventLoop, QSettings, QStandardPaths, QTimer, QUrl, Qt

ESCHER_VERSION = "1.7.3"
ESCHER_URL = f"https://unpkg.com/escher@{ESCHER_VERSION}/dist/escher.min.js"
SETTING = "ai/escher_js"

PAGE = """<!doctype html><html><head><meta charset="utf-8"><script src="{script}"></script></head>
<body style="margin:0"><div id="map" style="width:1200px;height:800px"></div>
<script>
window.__svg = null; window.__error = null;
const createObjectURL = URL.createObjectURL.bind(URL);
URL.createObjectURL = function (blob) {{  // Escher's save_svg hands the SVG to FileSaver: keep its text
  if (blob && blob.type === 'image/svg+xml') {{ blob.text().then(t => {{ window.__svg = t; }}); }}
  return createObjectURL(blob);
}};
window.exportMap = function (data) {{
  window.__svg = null; window.__error = null;
  document.getElementById('map').innerHTML = '';
  try {{
    escher.Builder(data, null, null, document.getElementById('map'), {{
      menu: 'zoom', enable_editing: false, enable_keys: false, enable_search: false, never_ask_before_quit: true,
      scroll_behavior: 'none', full_screen_button: false,
      first_load_callback: function (builder) {{ builder.map.save_svg(); }}
    }});
  }} catch (e) {{ window.__error = String(e && e.stack || e); }}
}};
</script></body></html>"""


def web_engine_available() -> bool:
    """Whether QtWebEngine can be used (imported here: call it before the QApplication is made)."""
    os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu")  # the export needs no GPU
    try:
        from qtpy.QtWebEngineWidgets import QWebEngineView  # noqa: F401
    except Exception:  # noqa: BLE001 (ImportError, or a binding without QtWebEngine)
        return False
    return True


def cached_script() -> Path:
    location = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    folder = Path(location or tempfile.gettempdir())
    return folder / f"escher-{ESCHER_VERSION}.min.js"


def find_script(settings: QSettings | None = None) -> Path | None:
    """escher.min.js from the setting, the escher package or the data folder (None: not found)."""
    setting = (settings or QSettings()).value(SETTING, "")
    if setting and Path(setting).is_file():
        return Path(setting)
    try:
        import escher

        path = Path(escher.__file__).parent / "static" / "escher.min.js"
        if path.is_file():
            return path
    except ImportError:
        pass
    path = cached_script()
    return path if path.is_file() else None


def download_script() -> Path:
    """Download escher.min.js (:data:`ESCHER_URL`) into the data folder."""
    path = cached_script()
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(ESCHER_URL, timeout=60) as response:
        data = response.read()
    if b"escher" not in data[:2000] and b"Builder" not in data:
        raise OSError(f"{ESCHER_URL} did not give escher.min.js")
    path.write_bytes(data)
    return path


class EscherSvg:
    """Exports Escher maps as SVG; made in the GUI thread (needs QtWebEngine, see :func:`web_engine_available`)."""

    def __init__(self, script: Path):
        from qtpy.QtWebEngineWidgets import QWebEngineView
        from qtpy.QtWidgets import QApplication

        app = QApplication.instance()
        if app is None or not app.arguments():  # QtWebEngine crashes without a program name
            raise RuntimeError("QtWebEngine needs a QApplication made with a program name (argv[0])")
        self.folder = Path(tempfile.mkdtemp(prefix="ecmapper_escher_"))
        page = self.folder / "export.html"
        page.write_text(PAGE.format(script=QUrl.fromLocalFile(str(script)).toString()), encoding="utf-8")
        self.view = QWebEngineView()
        self.view.setAttribute(Qt.WidgetAttribute.WA_DontShowOnScreen)  # shown (so scripts and timers run), not seen
        self.view.resize(1200, 800)
        self.view.show()
        loop = QEventLoop()
        self.loaded = None
        self.view.loadFinished.connect(lambda ok: (setattr(self, "loaded", ok), loop.quit()))
        self.view.load(QUrl.fromLocalFile(str(page)))
        QTimer.singleShot(30000, loop.quit)
        loop.exec()
        if not self.loaded:
            raise RuntimeError(f"the page with {script} could not be loaded")

    def _js(self, code: str, timeout: float = 5.0):
        result = {}
        loop = QEventLoop()
        self.view.page().runJavaScript(code, lambda value: (result.__setitem__("value", value), loop.quit()))
        QTimer.singleShot(int(timeout * 1000), loop.quit)
        loop.exec()
        return result.get("value")

    def export(self, escher_map: list, timeout: float = 120.0) -> str:
        """The SVG of an Escher map ([head, body]) as Escher exports it (the canvas)."""
        self.view.page().runJavaScript(f"window.exportMap({json.dumps(escher_map)});")
        start = time.time()
        while time.time() - start < timeout:
            value = self._js("window.__svg || (window.__error ? 'ERROR ' + window.__error : null)")
            if value:
                if value.startswith("ERROR "):
                    raise RuntimeError(f"Escher could not draw the map: {value[6:500]}")
                return value
            loop = QEventLoop()
            QTimer.singleShot(100, loop.quit)
            loop.exec()
        raise RuntimeError(f"Escher did not export the map within {timeout:.0f} s")

    def close(self) -> None:
        self.view.close()
        self.view.deleteLater()
