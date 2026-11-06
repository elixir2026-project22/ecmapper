"""Projects of the draft map GUI: all files of a layout task in one zip file (``.ecmap``).

A project holds the model (and the JSON model with its subsystems), the
result of matching it to the BRENDA map (``prep.pkl``), the assignment table
with the AI layouts, the AI work folders and the condensed subsystems of
the GUI (``gui_cache.pkl``, so a project opens quickly). The BRENDA map of
ecmapper (``BRENDAmap/``) is only referenced by its SHA-256; another
template is copied into the project. The options of :mod:`.run` are kept in
``project.json``.

An open project is unpacked into a work folder in the ecmapper cache
(``projects/``), where the GUI and the AI jobs write their files; saving
zips the work folder again. A work folder is reused when the project is
opened again unchanged.
"""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
import shutil
import sys
import time
import uuid
import zipfile
from pathlib import Path

#: the extension of project files
EXT = ".ecmap"
#: the version of the layout of project files
FORMAT = 1
META, PREP, PREP_KEY, TABLE = "project.json", "prep.pkl", "prep.key", "assignment.xlsx"
MODEL_DIR, TEMPLATE_DIR = "model", "template"
TEMPLATE_FILES = ("BRENDA_map_reactions.json", "BRENDA_map.svg")
#: files of the work folder that are not saved: the lock and the stamp of the unpacked project file
LOCK, SOURCE = ".lock", ".source"
EXCLUDE = {LOCK, SOURCE}
#: options of :mod:`.run` that are not kept (the GUI always splits, the table replaces --central and --maps)
NOT_KEPT = {"split", "assignment", "write_assignment", "central", "maps", "project"}
#: options of the GUI besides those of :mod:`.run`, with their defaults
GUI_OPTIONS = {"vicinity": "propagate"}

_digests: dict = {}


def bundled_template() -> tuple[Path, Path]:
    """The BRENDA map of ecmapper (BRENDAmap/ of the project), the default template of the draft maps."""
    folder = Path(__file__).resolve().parents[3] / "BRENDAmap"
    return folder / TEMPLATE_FILES[0], folder / TEMPLATE_FILES[1]


def work_root() -> Path:
    """The folder of the work folders of open projects (``projects/`` in the ecmapper cache)."""
    from ..sources.base import default_cache_dir

    return default_cache_dir() / "projects"


def file_digest(path) -> str:
    """SHA-256 of the contents of a file (kept per path, size and modification time in this process)."""
    path = Path(path)
    st = path.stat()
    key = (str(path.resolve()), st.st_size, st.st_mtime_ns)
    if key not in _digests:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        _digests[key] = h.hexdigest()
    return _digests[key]


def digests(files) -> list:
    """(name, SHA-256) of the files that are given (None: left out), for keys of caches."""
    return [(Path(f).name, file_digest(f)) for f in files if f]


def ecmapper_version() -> str:
    try:
        from importlib.metadata import version

        return version("ecmapper")
    except Exception:  # noqa: BLE001 (not installed)
        return "unknown"


def option_defaults() -> dict:
    """The kept options of :mod:`.run` and of the GUI with their defaults."""
    from .run import build_parser

    defaults = {a.dest: a.default for a in build_parser(positional=False)._actions
                if a.option_strings and a.dest != "help" and a.dest not in NOT_KEPT}
    return {**defaults, **GUI_OPTIONS}


def options_of(args) -> dict:
    """The kept options of `args` (JSON values)."""
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in
            ((k, getattr(args, k, d)) for k, d in option_defaults().items())}


def relayout_methods(args) -> list[str]:
    """The layouts of step 3b of `args` (as :func:`.run.check_args`)."""
    return [] if args.relayout == "none" else [m for m in args.relayout.split(",") if m]


def work_id(path) -> str:
    """The name of the work folder of the project file `path`."""
    return hashlib.sha1(str(Path(path).resolve()).lower().encode()).hexdigest()[:16]


def _stamp(path: Path) -> dict:
    st = path.stat()
    return dict(path=str(path.resolve()), size=st.st_size, mtime_ns=st.st_mtime_ns)


def _alive(pid: int) -> bool:
    if pid == os.getpid():
        return True
    if sys.platform == "win32":
        handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        code = ctypes.c_ulong()
        ctypes.windll.kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
        ctypes.windll.kernel32.CloseHandle(handle)
        return code.value == 259  # STILL_ACTIVE
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def extract(path, folder) -> None:
    """Unpack the project file `path` into `folder` (members outside it are refused)."""
    folder = Path(folder).resolve()
    with zipfile.ZipFile(path) as z:
        for info in z.infolist():
            target = (folder / info.filename).resolve()
            if info.filename.startswith(("/", "\\")) or ":" in info.filename or not target.is_relative_to(folder):
                raise ValueError(f"{path}: the member {info.filename!r} lies outside the project")
        z.extractall(folder)


class Project:
    """A project file (`path`, None until it is saved) and its work folder (`work`); `meta`: ``project.json``."""

    def __init__(self, work: Path, path: Path | None, meta: dict):
        self.work, self.path, self.meta = Path(work), Path(path) if path else None, meta

    # --- making and opening ---------------------------------------------------------------------

    @classmethod
    def new(cls, model=None, subsystems=None, brenda_json=None, brenda_svg=None, options=None,
            assignment=None) -> Project:
        """A new project of `model` (a model file or a PREP.pkl) in a new work folder.

        The files are copied into it; the template (default:
        :func:`bundled_template`) is referenced if it is the bundled one, else
        copied. `assignment`: a table (XLSX) to start from; it is read when the
        project is loaded and saved into the project then.
        """
        work = work_root() / f"new-{uuid.uuid4().hex[:12]}"
        work.mkdir(parents=True)
        meta = dict(format=FORMAT, ecmapper=ecmapper_version(), created=time.strftime("%Y-%m-%d %H:%M:%S"),
                    model=None, subsystems=None, options=dict(options or {}), build_out=None,
                    start_assignment=str(Path(assignment).resolve()) if assignment else None)
        if model is not None:
            model = Path(model)
            meta["start_model"] = str(model.resolve())
            if model.suffix.lower() == ".pkl":
                shutil.copyfile(model, work / PREP)
            else:
                (work / MODEL_DIR).mkdir()
                shutil.copyfile(model, work / MODEL_DIR / model.name)
                meta["model"] = f"{MODEL_DIR}/{model.name}"
        if subsystems is not None:
            subsystems = Path(subsystems)
            name = subsystems.name
            if f"{MODEL_DIR}/{name}" == meta["model"]:
                name = f"subsystems_{name}"
            (work / MODEL_DIR).mkdir(exist_ok=True)
            shutil.copyfile(subsystems, work / MODEL_DIR / name)
            meta["subsystems"] = f"{MODEL_DIR}/{name}"
        bundled = bundled_template()
        given = (Path(brenda_json or bundled[0]), Path(brenda_svg or bundled[1]))
        if all(g.resolve() == b.resolve() for g, b in zip(given, bundled)):
            meta["template"] = dict(bundled=True, sha256={f.name: file_digest(f) for f in bundled})
        else:
            (work / TEMPLATE_DIR).mkdir()
            for g, name in zip(given, TEMPLATE_FILES):
                shutil.copyfile(g, work / TEMPLATE_DIR / name)
            meta["template"] = dict(bundled=False, source=[str(g.resolve()) for g in given])
        project = cls(work, None, meta)
        project.write_meta()
        project.lock()
        return project

    @classmethod
    def open(cls, path) -> tuple[Project, list[str]]:
        """The project file `path` in its work folder (reused if the file did not change since it was unpacked
        or saved there) and warnings."""
        path = Path(path).resolve()
        work = work_root() / work_id(path)
        warnings = []
        lock = work / LOCK
        if lock.is_file():
            try:
                pid = int(lock.read_text().strip())
            except ValueError:
                pid = None
            if pid is not None and pid != os.getpid() and _alive(pid):
                warnings.append(f"{path.name} is open in another window too (process {pid}); changes saved there "
                                "and here overwrite each other")
                work = work.with_name(f"{work.name}-{os.getpid()}")  # its own work folder, unpacked again
        try:
            reuse = json.loads((work / SOURCE).read_text(encoding="utf-8")) == _stamp(path)
        except (OSError, ValueError):
            reuse = False
        if not reuse:
            if work.exists():
                shutil.rmtree(work)
            work.mkdir(parents=True)
            extract(path, work)
            (work / SOURCE).write_text(json.dumps(_stamp(path)), encoding="utf-8")
        try:
            meta = json.loads((work / META).read_text(encoding="utf-8"))
        except (OSError, ValueError) as e:
            raise ValueError(f"{path} is no ecmapper project: {e}") from e
        if meta.get("format", 0) > FORMAT:
            raise ValueError(f"{path} was saved by a newer ecmapper (project format {meta['format']})")
        project = cls(work, path, meta)
        template = meta.get("template", {})
        if template.get("bundled"):
            for f in bundled_template():
                if not f.is_file():
                    raise FileNotFoundError(f"the BRENDA map {f} of ecmapper is missing")
                if template.get("sha256", {}).get(f.name) not in (None, file_digest(f)):
                    warnings.append(f"the BRENDA map of ecmapper ({f.name}) changed since the project was made; "
                                    "the subsystems are drawn from the current one")
        project.lock()
        return project, warnings

    # --- files -----------------------------------------------------------------------------------

    @property
    def name(self) -> str:
        return self.path.name if self.path else "new project"

    def file(self, key: str) -> Path | None:
        """The file of the entry `key` ("model", "subsystems") of the project, None if it has none."""
        return self.work / self.meta[key] if self.meta.get(key) else None

    def template(self) -> tuple[Path, Path]:
        if self.meta.get("template", {}).get("bundled", True):
            return bundled_template()
        return self.work / TEMPLATE_DIR / TEMPLATE_FILES[0], self.work / TEMPLATE_DIR / TEMPLATE_FILES[1]

    @property
    def build_out(self) -> Path | None:
        """The folder of the maps built from the project: the last one used, else PROJECT_maps next to it."""
        if self.meta.get("build_out"):
            return Path(self.meta["build_out"])
        if self.path is not None:
            return self.path.with_name(f"{self.path.stem}_maps")
        model = self.meta.get("start_model")
        return Path(model).parent / f"{Path(model).name.split('.')[0]}_maps" if model else None

    def apply_files(self, args) -> None:
        """Set the input files of `args` to those of the project; its work folder is OUT_DIR."""
        args.model, args.subsystems = self.file("model"), self.file("subsystems")
        args.prep = self.work / PREP if (self.work / PREP).is_file() or args.model is None else None
        args.brenda_json, args.brenda_svg = self.template()
        start = self.meta.get("start_assignment")
        args.assignment = (self.work / TABLE if (self.work / TABLE).is_file() else Path(start) if start else None)
        args.out = self.work
        args.split = True
        args.central = args.maps = args.write_assignment = None

    def apply(self, args) -> list[str]:
        """Set the options and the input files of `args` to those of the project; returns the layouts of step 3b."""
        for k, v in {**option_defaults(), **self.meta.get("options", {})}.items():
            setattr(args, k, v)
        self.apply_files(args)
        return relayout_methods(args)

    def write_meta(self) -> None:
        (self.work / META).write_text(json.dumps(self.meta, indent=1), encoding="utf-8")

    def save(self, path=None) -> Path:
        """Zip the work folder into the project file `path` (default: `path`); written to a temporary file first,
        then put in place."""
        path = Path(path or self.path).resolve()
        if path.suffix.lower() != EXT:
            path = path.with_name(path.name + EXT)
        self.meta.update(format=FORMAT, ecmapper=ecmapper_version(), saved=time.strftime("%Y-%m-%d %H:%M:%S"))
        self.meta.pop("start_assignment", None)  # the table is in the project now (if it was loaded)
        self.write_meta()
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        try:
            with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
                for f in sorted(self.work.rglob("*")):
                    rel = f.relative_to(self.work)
                    if f.is_file() and rel.parts[0] not in EXCLUDE:
                        z.write(f, rel.as_posix())
            os.replace(tmp, path)
        finally:
            if tmp.exists():
                tmp.unlink()
        self.path = path
        if self.work.name == work_id(path):
            (self.work / SOURCE).write_text(json.dumps(_stamp(path)), encoding="utf-8")
        else:  # saved as another file: this folder no longer shows the file it was unpacked from
            (self.work / SOURCE).unlink(missing_ok=True)
        return path

    # --- the work folder -------------------------------------------------------------------------

    def lock(self) -> None:
        (self.work / LOCK).write_text(str(os.getpid()))

    def close(self, discard: bool = False) -> None:
        """The project is not used any more: release the work folder; remove it if it is not the work folder of a
        saved project file (or `discard`: changes since it was saved are dropped, it is unpacked again next time)."""
        try:
            (self.work / LOCK).unlink(missing_ok=True)
            if self.path is None or self.work.name != work_id(self.path):
                shutil.rmtree(self.work, ignore_errors=True)
            elif discard:
                (self.work / SOURCE).unlink(missing_ok=True)
        except OSError:
            pass
