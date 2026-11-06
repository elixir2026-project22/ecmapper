"""Layout toolkits for AI jobs and the jobs themselves.

A toolkit is a folder with a Python package that an AI runs to lay out an
Escher map (e.g. escher_handlayout) and instructions for the AI. It is
described by `ecmapper_toolkit.json` in the folder (:data:`MANIFEST`)::

    {"name": "my_toolkit", "package": "my_toolkit", "instructions": ["PROMPT.md"],
     "commands": ["inspect", "run"], "result": "**/*_final.json", "notes": "..."}

`package` is run as ``python -m package ARGS`` in the work folder (with the
toolkit folder on PYTHONPATH), `commands` are the allowed first arguments
(empty: any), `instructions` are put into the AI's system prompt, `result`
finds the final map in the work folder. escher_handlayout and
escher_handlayout_kit are recognised without a manifest
(:data:`ESCHER_HANDLAYOUT`, :data:`ESCHER_HANDLAYOUT_KIT`).
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

from ..ai_layouts import DRAFT
from .agent import Tool, ToolError
from .workspace import TOOLKIT, Workspace

MANIFEST = "ecmapper_toolkit.json"
ESCHER_HANDLAYOUT = dict(
    name="escher_handlayout", package="escher_handlayout", instructions=["PROMPT.md"],
    commands=["inspect", "run", "issues", "query", "patterns", "validate", "render"],
    notes=f"Parameters of PROMPT.md: INPUT_MAP = {DRAFT} (in the work folder), OUTPUT_DIR = . (the work folder; use a "
          "new subfolder per run, e.g. v1, v2, as PROMPT.md says). The commands written `python -m escher_handlayout "
          "...` or `eh ...` there are run with the tool run_toolkit (its arguments after the package name). The other "
          f"documents (SPEC_REFERENCE.md, PATTERNS.md, examples/) are readable as {TOOLKIT}/... The map's subsystems "
          "(the title labels) are the zone proposals.",
)
#: escher_handlayout_kit (run.py with inspect, check-spec, build, validate; JSON layout specs)
ESCHER_HANDLAYOUT_KIT = dict(
    name="escher_handlayout_kit", package="run", instructions=["PROMPT.md"],
    commands=["inspect", "check-spec", "build", "validate"], result="**/*_handlayout.json",
    notes=f"Parameters of PROMPT.md: KIT_DIR = {TOOLKIT} (read-only: write specs into the work folder, e.g. "
          f"specs/MAP.layout.json, starting from a copy of {TOOLKIT}/specs/layout_template.json), INPUT_MAP = {DRAFT}, "
          "OUT_DIR = . (the work folder), MAP_NAME = draft. The commands written `python run.py ...` there are run "
          "with the tool run_toolkit (its arguments after run.py, e.g. [\"inspect\", \"--input\", "
          f"\"{DRAFT}\", \"--out\", \"analysis\"]). Skip the setup step (the requirements are installed; there is no "
          "selftest command). Only JSON specs: never a Python spec or --unsafe-python-spec. The final map is the "
          "best build's *_handlayout.json; submit that file. The map's subsystems (the title labels) are the zone "
          f"proposals. The other documents (docs/GUIDE.md, schema/, specs/, examples/) are readable as {TOOLKIT}/....",
)


@dataclass
class Toolkit:
    name: str
    folder: Path
    package: str
    instructions: list[str] = field(default_factory=lambda: ["PROMPT.md"])
    commands: list[str] = field(default_factory=list)
    result: str = "**/*_final.json"
    notes: str = ""

    @classmethod
    def detect(cls, folder) -> Toolkit:
        """The toolkit in `folder` (its manifest, or escher_handlayout); ValueError if there is none."""
        folder = Path(folder)
        manifest = folder / MANIFEST
        if manifest.is_file():
            data = json.loads(manifest.read_text(encoding="utf-8"))
            if "package" not in data:
                raise ValueError(f"{manifest}: 'package' is missing")
            data.setdefault("name", folder.name)
            return cls(folder=folder, **{k: v for k, v in data.items() if k in cls.__dataclass_fields__})
        if (folder / "escher_handlayout" / "__main__.py").is_file():
            return cls(folder=folder, **ESCHER_HANDLAYOUT)
        if (folder / "run.py").is_file() and (folder / "escher_handlayout" / "declspec.py").is_file():
            return cls(folder=folder, **ESCHER_HANDLAYOUT_KIT)
        raise ValueError(f"no toolkit in {folder}: it needs {MANIFEST} (see DRAFT_MAP.md) or the escher_handlayout "
                         "package")

    def settings(self) -> dict:
        return dict(folder=str(self.folder))

    def instructions_text(self) -> str:
        parts = []
        for name in self.instructions:
            p = self.folder / name
            if p.is_file():
                parts.append(f"=== {name} (toolkit/{name}) ===\n{p.read_text(encoding='utf-8', errors='replace')}")
        return "\n\n".join(parts)

    def command(self, args: list[str]) -> list[str]:
        if self.commands and (not args or args[0] not in self.commands):
            raise ToolError(f"allowed commands: {', '.join(self.commands)}")
        return [sys.executable, "-m", self.package, *args]

    def env(self) -> dict:
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(p for p in (str(self.folder), env.get("PYTHONPATH", "")) if p)
        env["PYTHONIOENCODING"] = "utf-8"
        env["PYTHONUTF8"] = "1"  # files written as UTF-8 also on Windows (the maps have "→" labels)
        env["MPLBACKEND"] = "Agg"
        return env

    def newest_result(self, folder) -> Path | None:
        """The newest final map in `folder` (none of *_PARTIAL.json)."""
        found = [p for p in Path(folder).glob(self.result) if "_PARTIAL" not in p.name]
        return max(found, key=lambda p: p.stat().st_mtime) if found else None


PREAMBLE = """You lay out a metabolic map for ecmapper, a program that draws maps of metabolic models from the BRENDA \
map. The user works in ecmapper's GUI and talks to you in a small chat panel: keep your messages short.

The map "{map}" ({facts}) is in your work folder as {draft}, an Escher map as ecmapper draws it now: the \
subsystems are boxes, titled by text labels, arranged by the user (take that arrangement as a rough plan unless the \
user says otherwise). Lay out the whole map with the toolkit {toolkit} as its instructions below say; you may \
interconnect the subsystems (e.g. merge metabolites they share). The user's requests in the chat take precedence over \
the toolkit's defaults. Never change the biochemistry: the same reactions, stoichiometry and reversibility.

Tools: list_files, read_file, write_file and edit_file work in the work folder (the toolkit's files are readable as \
{toolkit_dir}/...); run_toolkit runs `python -m {package} ARGS` in the work folder (the user may have to approve each \
run); view_image shows a PNG. When a final map is written, call submit_layout with its path: ecmapper checks it \
against the draft and shows it in the GUI; the user may then ask for changes, and you may submit again. Work on your \
own until the map is submitted (or you are stuck), then report briefly.

Subsystems of the map (ID: name, reactions):
{subsystems}
{stubs}
=== Notes for this toolkit ===
{notes}

{instructions}"""


def map_facts(subsystems: list[tuple[str, str, int]], stubs: dict[str, list[str]]) -> tuple[str, str]:
    """The lists of the subsystems and of the metabolites that continue on other maps (empty if none) for a task."""
    text = "\n".join(f"- {k}: {label}, {n} reactions" for k, label, n in subsystems)
    stub_text = ""
    if stubs:
        stub_text = ("\nMetabolites that continue on other maps (ecmapper keeps their '→ map' labels up to date):\n"
                     + "\n".join(f"- {m}: {', '.join(maps)}" for m, maps in sorted(stubs.items())) + "\n")
    return text, stub_text


class LayoutJob:
    """An AI job: lay out map `map_name` (exported as :data:`DRAFT` into `folder`) with `toolkit`.

    `submit(path)` is called by the tool submit_layout: it checks the map
    and shows it (raises :class:`ToolError` with the problems otherwise) and
    returns a summary. `approve_all`: run the toolkit without asking.
    """

    def __init__(self, toolkit: Toolkit, folder, map_name: str, facts: str, subsystems: list[tuple[str, str, int]],
                 stubs: dict[str, list[str]], submit):
        self.toolkit, self.map_name, self.facts = toolkit, map_name, facts
        self.subsystems, self.stubs, self.submit = subsystems, stubs, submit
        self.workspace = Workspace(folder, toolkit.folder)
        self.folder = self.workspace.folder
        self.submitted: list[Path] = []

    def system_prompt(self) -> str:
        subsystems, stubs = map_facts(self.subsystems, self.stubs)
        return PREAMBLE.format(map=self.map_name, facts=self.facts, draft=DRAFT, toolkit=self.toolkit.name,
                               toolkit_dir=TOOLKIT, package=self.toolkit.package, subsystems=subsystems, stubs=stubs,
                               notes=self.toolkit.notes or "(none)", instructions=self.toolkit.instructions_text())

    def run_toolkit(self, args: list[str]) -> str:
        return self.workspace.run(self.toolkit.command([str(a) for a in args]), env=self.toolkit.env())

    def submit_layout(self, path: str) -> str:
        p = self.workspace.path(path)
        if not p.is_file():
            newest = self.toolkit.newest_result(self.folder)
            hint = f" (the newest result is {self.workspace.show(newest)})" if newest else ""
            raise ToolError(f"{path}: no such file{hint}")
        summary = self.submit(p)
        self.submitted.append(p)
        return summary

    def tools(self) -> list[Tool]:
        commands = f" Allowed commands: {', '.join(self.toolkit.commands)}." if self.toolkit.commands else ""
        return [
            *self.workspace.tools(),
            Tool("run_toolkit", f"Run `python -m {self.toolkit.package} ARGS` in the work folder; returns the exit "
                 f"code and the output (shortened).{commands}",
                 {"type": "object", "properties": {"args": {"type": "array", "items": {"type": "string"},
                                                            "description": "the arguments after the package name"}},
                  "required": ["args"], "additionalProperties": False}, self.run_toolkit, needs_approval=True),
            Tool("submit_layout", "Submit a final Escher map (path in the work folder) as the layout of the map: "
                 "ecmapper checks it against the draft and shows it to the user. Returns a summary or the problems.",
                 {"type": "object", "properties": {"path": {"type": "string"}}, "required": ["path"],
                  "additionalProperties": False}, self.submit_layout),
        ]
