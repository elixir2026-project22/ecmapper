"""Work folders for coding agents (Claude Code, pi, …): the AI job without ecmapper's chat.

:func:`prepare` writes, next to the exported map (:data:`..ai_layouts.DRAFT`),
the task as AGENTS.md (read by pi, Codex and others when they are started in
the folder; CLAUDE.md imports it for Claude Code), a script that runs the toolkit (:data:`RUNNER`) and the job's record
(:data:`JOB`: the map and the draft objects of its Escher IDs). The GUI then
takes the agent's final map as the map's layout (:func:`find_job`,
:func:`newest_job`), like a map submitted in the chat.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from ..ai_layouts import DRAFT
from .toolkits import Toolkit, map_facts

TASK = "AGENTS.md"
CLAUDE = "CLAUDE.md"
RUNNER = "run_toolkit.py"
JOB = "ecmapper_job.json"
FIRST_PROMPT = "Lay out this map."
AGENTS = {  # name: its start command (TOOLKIT: the toolkit's folder); Claude Code asks before it runs commands, pi not
    "Claude Code": ["claude", "--add-dir", "TOOLKIT", FIRST_PROMPT],
    "pi": ["pi", FIRST_PROMPT],
}

TASK_TEXT = """# Lay out the map "{map}" for ecmapper

ecmapper draws maps of metabolic models from the BRENDA map. The map "{map}" ({facts}) is in this folder as \
`{draft}`, an Escher map as ecmapper draws it now: the subsystems are boxes, titled by text labels, arranged by the \
user (take that arrangement as a rough plan unless the user says otherwise). Lay out the whole map with the toolkit \
{toolkit} as its instructions say; you may interconnect the subsystems (e.g. merge metabolites they share). The \
user's requests take precedence over the toolkit's defaults. Never change the biochemistry: the same reactions, \
stoichiometry and reversibility.

## The toolkit

- Its folder is `{toolkit_dir}`; read its instructions first: {instructions}.
- Run it only with `{python} {runner} ARGS`, which runs `python -m {package} ARGS` in this folder with the toolkit \
on the path and UTF-8 output. Where the instructions write `python -m {package} ...` (or a short alias of it), use \
this command.{commands}
- Notes of ecmapper for this toolkit (written for ecmapper's own chat: "run_toolkit" there is `{runner}` here, and \
`toolkit/...` is the toolkit's folder): {notes}

## The map

Subsystems (ID: name, reactions):
{subsystems}
{stubs}
## When you are done

Write the final map into this folder (the toolkit's final `*_final.json`, not a partial one) and tell the user its \
path. The user takes it into ecmapper with *AI ▸ Use coding agent result*: ecmapper checks it against `{draft}` \
(the same reactions and stoichiometry) and uses the newest final map in this folder. Do not edit `{draft}` or \
`{job}`, `{claude}` or this file.
"""

RUNNER_TEXT = '''"""Runs the toolkit {name} for this work folder (written by ecmapper): python {runner} ARGS."""

import os
import subprocess
import sys

TOOLKIT = {folder!r}
PACKAGE = {package!r}
PYTHON = {python!r}

env = dict(os.environ)
env["PYTHONPATH"] = os.pathsep.join(p for p in (TOOLKIT, env.get("PYTHONPATH", "")) if p)
env.update(PYTHONIOENCODING="utf-8", PYTHONUTF8="1", MPLBACKEND="Agg")
python = PYTHON if os.path.isfile(PYTHON) else sys.executable
folder = os.path.dirname(os.path.abspath(__file__))
sys.exit(subprocess.call([python, "-m", PACKAGE, *sys.argv[1:]], env=env, cwd=folder))
'''


def prepare(folder, toolkit: Toolkit, map_name: str, facts: str, subsystems: list[tuple[str, str, int]],
            stubs: dict[str, list[str]], ids: dict, python: str = sys.executable) -> Path:
    """Write the task, the toolkit's runner and the job record into `folder` (with the exported map in it)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    instructions = ", ".join(f"`{toolkit.folder / name}`" for name in toolkit.instructions
                             if (toolkit.folder / name).is_file()) or "its README"
    commands = f" Its commands: {', '.join(toolkit.commands)}." if toolkit.commands else ""
    subsystem_text, stub_text = map_facts(subsystems, stubs)
    (folder / TASK).write_text(TASK_TEXT.format(
        map=map_name, facts=facts, draft=DRAFT, toolkit=toolkit.name, toolkit_dir=toolkit.folder,
        instructions=instructions, python=f'"{python}"',
        runner=RUNNER, package=toolkit.package, commands=commands, notes=toolkit.notes or "(none)",
        subsystems=subsystem_text, stubs=stub_text, job=JOB, claude=CLAUDE), encoding="utf-8")
    (folder / CLAUDE).write_text(f"@{TASK}\n", encoding="utf-8")  # Claude Code reads CLAUDE.md, imports AGENTS.md
    (folder / RUNNER).write_text(RUNNER_TEXT.format(name=toolkit.name, runner=RUNNER, folder=str(toolkit.folder),
                                                    package=toolkit.package, python=python), encoding="utf-8")
    job = dict(map=map_name, toolkit=str(toolkit.folder), toolkit_name=toolkit.name, ids=ids)
    (folder / JOB).write_text(json.dumps(job, indent=1), encoding="utf-8")
    return folder


def agent_args(agent: str, toolkit: Toolkit) -> list[str]:
    """The command (program and arguments) that starts coding agent `agent` (:data:`AGENTS`) on a job."""
    return [str(toolkit.folder) if a == "TOOLKIT" else a for a in AGENTS[agent]]


def command(folder, toolkit: Toolkit, agent: str) -> str:
    """The commands that start coding agent `agent` on the job in `folder` (two lines, for a terminal)."""
    args = " ".join(f'"{a}"' if " " in a or "\\" in a or a == FIRST_PROMPT else a for a in agent_args(agent, toolkit))
    return f'cd "{folder}"\n{args}'


def start_terminal(folder, toolkit: Toolkit, agent: str) -> None:
    """Start coding agent `agent` on the job in a new terminal window (Windows: cmd; elsewhere OSError)."""
    if sys.platform != "win32":
        raise OSError("starting a terminal is only implemented for Windows; run the command yourself")
    subprocess.Popen(["cmd", "/k", *agent_args(agent, toolkit)], cwd=str(folder),
                     creationflags=subprocess.CREATE_NEW_CONSOLE)


def read_job(folder) -> dict | None:
    """The job record in `folder` (None: no job there)."""
    path = Path(folder) / JOB
    try:
        job = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return job if isinstance(job, dict) and "map" in job and "ids" in job else None


def find_job(path) -> tuple[Path, dict] | None:
    """The job folder (and its record) that file `path` lies in (in a subfolder of it, too)."""
    for folder in Path(path).resolve().parents:
        job = read_job(folder)
        if job is not None:
            return folder, job
    return None


def newest_job(base, map_name: str) -> tuple[Path, dict, Path | None] | None:
    """The newest prepared job folder of map `map_name` below `base` (OUT/ai/MAP), its record and the newest final
    map in it (None if it has none yet); None without job."""
    jobs = []
    for path in Path(base).glob(f"*/{JOB}"):
        job = read_job(path.parent)
        if job is not None and job["map"] == map_name:
            jobs.append((path.stat().st_mtime, path.parent, job))
    if not jobs:
        return None
    _, folder, job = max(jobs, key=lambda j: j[0])
    try:
        result = Toolkit.detect(job["toolkit"]).newest_result(folder)
    except (ValueError, OSError):
        found = [p for p in folder.glob("**/*_final.json") if "_PARTIAL" not in p.name]
        result = max(found, key=lambda p: p.stat().st_mtime) if found else None
    return folder, job, result
