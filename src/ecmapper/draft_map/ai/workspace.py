"""The work folder of an AI job and the tools on it: list, read, write and edit files there (the toolkit's files are
readable as toolkit/...), view images, run commands (with a timeout, the output shortened)."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .agent import Tool, ToolError

#: name under which the toolkit's folder is readable
TOOLKIT = "toolkit"
SKIP = {"__pycache__", ".git", ".pytest_cache"}


def clip(text: str, limit: int) -> str:
    """`text` with the middle left out if it is longer than `limit` characters."""
    if len(text) <= limit:
        return text
    head = int(limit * 0.6)
    tail = limit - head
    return f"{text[:head]}\n… [{len(text) - limit} characters left out] …\n{text[-tail:]}"


class Workspace:
    """Files of an AI job: everything in `folder`; `toolkit` (a folder) readable as toolkit/..."""

    def __init__(self, folder, toolkit=None, max_output: int = 20000):
        self.folder = Path(folder).resolve()
        self.folder.mkdir(parents=True, exist_ok=True)
        self.toolkit = Path(toolkit).resolve() if toolkit else None
        self.max_output = max_output

    def path(self, name: str, write: bool = False) -> Path:
        """The file `name` (relative to the work folder, or toolkit/... for reading)."""
        name = str(name).replace("\\", "/").strip()
        parts = Path(name).parts
        if self.toolkit is not None and parts and parts[0] == TOOLKIT:
            if write:
                raise ToolError(f"{name}: the toolkit's files are read-only; write in the work folder")
            base, rest = self.toolkit, Path(*parts[1:]) if len(parts) > 1 else Path(".")
        else:
            base, rest = self.folder, Path(name)
        path = (base / rest).resolve()
        if path != base and base not in path.parents:
            raise ToolError(f"{name}: outside the work folder")
        return path

    def show(self, path: Path) -> str:
        """The name of `path` for the model."""
        if self.toolkit is not None and (path == self.toolkit or self.toolkit in path.parents):
            return (Path(TOOLKIT) / path.relative_to(self.toolkit)).as_posix()
        return path.relative_to(self.folder).as_posix() or "."

    def list_files(self, path: str = ".", depth: int = 2) -> str:
        root = self.path(path)
        if not root.is_dir():
            raise ToolError(f"{path}: not a folder")
        lines = []

        def walk(folder: Path, level: int):
            for p in sorted(folder.iterdir(), key=lambda p: (p.is_file(), p.name.lower())):
                if p.name in SKIP or p.name.startswith("."):
                    continue
                if len(lines) >= 300:
                    return
                indent = "  " * level
                if p.is_dir():
                    lines.append(f"{indent}{p.name}/")
                    if level + 1 < depth:
                        walk(p, level + 1)
                else:
                    lines.append(f"{indent}{p.name} ({p.stat().st_size} bytes)")

        walk(root, 0)
        if len(lines) >= 300:
            lines.append("… (more files)")
        head = f"{self.show(root)}/" + ("" if root != self.folder or self.toolkit is None else
                                         f"  (the toolkit's files: {TOOLKIT}/)")
        return "\n".join([head, *lines]) if lines else f"{head} (empty)"

    def read_file(self, path: str, offset: int = 0, limit: int = 400) -> str:
        p = self.path(path)
        if not p.is_file():
            raise ToolError(f"{path}: no such file")
        data = p.read_bytes()
        if b"\0" in data[:4096]:
            raise ToolError(f"{path}: not a text file")
        lines = data.decode("utf-8", errors="replace").splitlines()
        offset, limit = max(0, offset), max(1, min(limit, 2000))
        part = [f"{i + 1:5d}  {line if len(line) <= 2000 else line[:2000] + ' …'}"
                for i, line in enumerate(lines[offset:offset + limit], start=offset)]
        more = len(lines) - offset - len(part)
        text = "\n".join(part)
        text = clip(text, self.max_output)
        note = f"\n… {more} more lines (offset={offset + len(part)})" if more > 0 else ""
        return f"{self.show(p)}: {len(lines)} lines\n{text}{note}"

    def write_file(self, path: str, content: str) -> str:
        p = self.path(path, write=True)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding="utf-8")
        return f"wrote {self.show(p)} ({len(content.splitlines())} lines)"

    def edit_file(self, path: str, old: str, new: str) -> str:
        p = self.path(path, write=True)
        if not p.is_file():
            raise ToolError(f"{path}: no such file")
        text = p.read_text(encoding="utf-8")
        count = text.count(old)
        if count != 1:
            raise ToolError(f"{path}: the old text occurs {count} times; it must occur exactly once (give more "
                            "context)")
        p.write_text(text.replace(old, new), encoding="utf-8")
        return f"edited {self.show(p)}"

    def view_image(self, path: str):
        p = self.path(path)
        if not p.is_file() or p.suffix.lower() != ".png":
            raise ToolError(f"{path}: no PNG file")
        if p.stat().st_size > 5_000_000:
            raise ToolError(f"{path}: larger than 5 MB; render a region or a smaller scale")
        return f"{self.show(p)}", [p.read_bytes()]

    def run(self, command: list[str], env: dict | None = None, timeout: float = 300.0) -> str:
        """Run `command` in the work folder; its exit code and output (shortened)."""
        try:
            done = subprocess.run(command, cwd=self.folder, env=env, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise ToolError(f"timed out after {timeout:.0f} s") from None
        except OSError as e:
            raise ToolError(f"could not run {command[0]}: {e}") from None
        out = done.stdout.decode("utf-8", errors="replace")
        err = done.stderr.decode("utf-8", errors="replace")
        text = out + (("\n[stderr]\n" + err) if err.strip() else "")
        return f"exit code {done.returncode}\n{clip(text.strip(), self.max_output)}"

    def tools(self) -> list[Tool]:
        """The file tools."""
        where = "the work folder" + (f" (the toolkit's files are readable as {TOOLKIT}/...)" if self.toolkit else "")
        return [
            Tool("list_files", f"List the files of a folder in {where}, two levels deep.",
                 {"type": "object", "properties": {"path": {"type": "string", "description": "folder, default ."}},
                  "additionalProperties": False}, self.list_files),
            Tool("read_file", f"Read a text file of {where} with line numbers, `limit` lines from `offset` (0-based).",
                 {"type": "object", "properties": {"path": {"type": "string"}, "offset": {"type": "integer"},
                                                   "limit": {"type": "integer"}},
                  "required": ["path"], "additionalProperties": False}, self.read_file),
            Tool("write_file", "Write (create or replace) a text file in the work folder.",
                 {"type": "object", "properties": {"path": {"type": "string"}, "content": {"type": "string"}},
                  "required": ["path", "content"], "additionalProperties": False}, self.write_file),
            Tool("edit_file", "Replace the text `old`, which must occur exactly once, by `new` in a file of the work "
                 "folder.",
                 {"type": "object", "properties": {"path": {"type": "string"}, "old": {"type": "string"},
                                                   "new": {"type": "string"}},
                  "required": ["path", "old", "new"], "additionalProperties": False}, self.edit_file),
            Tool("view_image", "Look at a PNG image (e.g. a render of the toolkit) of the work folder.",
                 {"type": "object", "properties": {"path": {"type": "string"}}, "required": ["path"],
                  "additionalProperties": False}, self.view_image),
        ]
