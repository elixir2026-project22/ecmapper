"""AI layouts of whole maps: export, positions, outdated layouts, the pipeline (ai_layouts.py), the harness (ai/)."""

import copy
import json
import pickle
import subprocess
import sys
from collections import Counter

import pytest

pytest.importorskip("scipy")
pytest.importorskip("numba")
pytest.importorskip("openpyxl")

from ecmapper.draft_map import ai_layouts, split
from ecmapper.draft_map.ai_layouts import MapLayout, check, export_draft, make_layout, map_positions, match_ids
from ecmapper.draft_map.core import ModelData
from test_brenda_map import model  # noqa: F401 (fixture)
from test_draft_map import prepared, run_draft_map


@pytest.fixture
def built(sources, model, tmp_path):  # noqa: F811
    """The maps of the test model (--split --pack) built from an assignment table: (prep, table, out, maps.json)."""
    _, prep = prepared(sources, model, tmp_path)
    xlsx = tmp_path / "assignment.xlsx"
    run_draft_map(prep, tmp_path / "first", "--split", "--pack", "--min-map-reactions", "0", "--write-assignment",
                  str(xlsx))
    out = tmp_path / "reference"
    run_draft_map(prep, out, "--split", "--pack", "--assignment", str(xlsx))
    return prep, xlsx, out, json.loads((out / "maps.json").read_text())


def final_map(out, name):
    with open(out / f"{split.slug(name)}.pkl", "rb") as f:
        return pickle.load(f)


def shifted(data, dx, dy):
    """A copy of an Escher map with all positions moved by (dx, dy), as a toolkit might lay it out."""
    data = copy.deepcopy(data)
    body = data[1]
    for n in body["nodes"].values():
        n["x"] += dx
        n["y"] += dy
        if "label_x" in n:
            n["label_x"] += dx
            n["label_y"] += dy
    for r in body["reactions"].values():
        r["label_x"] += dx
        r["label_y"] += dy
        for s in r["segments"].values():
            for b in ("b1", "b2"):
                if s.get(b):
                    s[b] = dict(x=s[b]["x"] + dx, y=s[b]["y"] + dy)
    for t in body["text_labels"].values():
        t["x"] += dx
        t["y"] += dy
    body["canvas"] = dict(body["canvas"], x=body["canvas"]["x"] + dx, y=body["canvas"]["y"] + dy)
    return data


def merge(data, keep, other):
    """Merge Escher node `other` into `keep` (as escher_handlayout does)."""
    body = data[1]
    for r in body["reactions"].values():
        for s in r["segments"].values():
            for end in ("from_node_id", "to_node_id"):
                if s[end] == other:
                    s[end] = keep
    del body["nodes"][other]


def test_export_and_positions(built):
    prep, _, out, index = built
    md = ModelData(str(prep))
    name = index["maps"][0]["name"]
    m = final_map(out, name)
    draft, ids = export_draft(m, md, name)
    kinds = Counter(key[0] for key in ids.values())
    assert set(ids) <= set(draft[1]["nodes"]) | set(draft[1]["reactions"]) | set(draft[1]["text_labels"])
    enzymes = [n for n in m.nodes.values() if n.kind == "enzyme" and n.reactions]
    assert kinds["reaction"] == kinds["mid"] == sum(len(n.reactions) for n in enzymes) == len(draft[1]["reactions"])
    assert kinds["title"] == len(m.subsystems)
    drawn = {key[1] for key in ids.values() if key[0] == "node"}
    assert drawn and all(m.nodes[i].kind != "enzyme" for i in drawn)

    # an export with the IDs of the first one: the same IDs for the same objects (markers and segments renumbered);
    # without a subsystem, the IDs of the others are kept
    again, ids2 = export_draft(m, md, name, reuse=ids)
    assert ids2 == ids
    for i, key in ids.items():
        if key[0] in ("node", "mid", "free"):
            assert again[1]["nodes"][i] == draft[1]["nodes"][i]
    smaller = copy.deepcopy(m)
    k = next(iter(smaller.subsystems))
    gone = {i for i, n in smaller.nodes.items() if n.subsystem == k}
    smaller.subsystems.pop(k)
    smaller.nodes = {i: n for i, n in smaller.nodes.items() if i not in gone}
    smaller.edges = [e for e in smaller.edges if e.enzyme not in gone and e.met not in gone]
    _, ids3 = export_draft(smaller, md, name, reuse=ids)
    assert {i: key for i, key in ids3.items() if i in ids} == {i: key for i, key in ids.items() if i in ids3}
    assert {tuple(key) for key in ids3.values()} == {tuple(key) for key in ids.values() if key[1] not in gone | {k}}

    # laid out elsewhere: everything moved; a metabolite node merged into another one of the same metabolite
    at_draft = map_positions(draft, ids, m, 2.0)
    assert set(at_draft) == drawn | {n.id for n in enzymes}
    for n in enzymes:
        assert at_draft[n.id] == pytest.approx((n.x, n.y), abs=80)
    final = shifted(draft, 1000, 500)
    by_metabolite = {}
    pair = None
    for i, key in ids.items():
        if key[0] == "node":
            if key[2] in by_metabolite and by_metabolite[key[2]][1] != key[1]:
                pair = (by_metabolite[key[2]], (i, key[1]))
                break
            by_metabolite[key[2]] = (i, key[1])
    run_log = None
    assert pair is not None  # the test map has a metabolite drawn twice
    if pair is not None:
        (keep, keep_node), (other, other_node) = pair
        merge(final, keep, other)
        run_log = dict(merges=[[ids[keep][2], keep, other]])
    layout = make_layout(m, md, name, final, ids, 2.0, svg="<svg/>", run_log=run_log)
    for node, (x, y) in at_draft.items():
        if pair is not None and node == other_node:
            continue
        assert layout.positions[node] == pytest.approx((x + 500, y + 250), abs=1e-6), node
    if pair is not None:
        assert layout.positions[other_node] == pytest.approx(layout.positions[keep_node])

    # made elsewhere (other Escher IDs): matched by content
    renumbered = copy.deepcopy(final)
    body = renumbered[1]
    new = {i: str(int(i) + 100000) for i in [*body["nodes"], *body["reactions"], *body["text_labels"]]}
    body["nodes"] = {new[i]: n for i, n in body["nodes"].items()}
    body["text_labels"] = {new[i]: t for i, t in body["text_labels"].items()}
    body["reactions"] = {new[i]: r for i, r in body["reactions"].items()}
    for r in body["reactions"].values():
        for s in r["segments"].values():
            s["from_node_id"], s["to_node_id"] = new[s["from_node_id"]], new[s["to_node_id"]]
    matched = match_ids(draft, ids, renumbered)
    positions = map_positions(renumbered, matched, m, 2.0)
    for node, p in layout.positions.items():
        if pair is None or node != other_node:
            assert positions[node] == pytest.approx(p, abs=1e-6), node
    assert {tuple(k) for k in matched.values() if k[0] == "stub"} == {tuple(k) for k in ids.values() if k[0] == "stub"}

    # saved and read again
    layout.save(out / "layout")
    again = MapLayout.load(out / "layout")
    assert again.positions == layout.positions and again.escher == layout.escher and again.svg == "<svg/>"
    assert again.ids == layout.ids and again.reactions == layout.reactions


def test_check_and_exit_labels(built):
    prep, _, out, index = built
    md = ModelData(str(prep))
    name = index["maps"][0]["name"]
    m = final_map(out, name)
    draft, ids = export_draft(m, md, name)
    layout = make_layout(m, md, name, shifted(draft, 10, 10), ids, 2.0)
    assert check(layout, m, md) == []
    # a reaction less in the layout, a coefficient changed
    broken = copy.deepcopy(layout)
    body = broken.escher[1]
    rid, r = next(iter(body["reactions"].items()))
    del body["reactions"][rid]
    other = next(iter(body["reactions"].values()))
    other["metabolites"][0]["coefficient"] *= 2
    problems = check(broken, m, md)
    assert any(r["bigg_id"] in p and "not in the layout" in p for p in problems)
    assert any(other["bigg_id"] in p and "stoichiometry" in p for p in problems)
    # a subsystem more on the map
    m2 = copy.deepcopy(m)
    k = next(iter(m2.subsystems))
    m2.subsystems["extra"] = copy.deepcopy(m2.subsystems[k])
    assert any("subsystems added" in p for p in check(layout, m2, md))

    # "→ map" labels follow the stubs of the map: changed, removed, added
    metabolites = [i for i, n in m.nodes.items() if n.kind != "enzyme" and any(key[0] == "node" and key[1] == i
                                                                                for key in ids.values())]
    m3 = copy.deepcopy(m)
    with_stub = [i for i in metabolites if getattr(m3.nodes[i], "stubs", None)]
    assert with_stub
    for i in with_stub[:1]:
        m3.nodes[i].stubs = []
    without = [i for i in metabolites if not getattr(m3.nodes[i], "stubs", None) and i not in with_stub][:2]
    for i in without:
        m3.nodes[i].stubs = ["Some other map"]
    assert ai_layouts.refresh_exit_labels(layout, m3)
    texts = {key[1]: layout.escher[1]["text_labels"][i]["text"] for i, key in layout.ids.items()
             if key[0] == "stub" and i in layout.escher[1]["text_labels"]}
    for i in with_stub[:1]:
        assert i not in texts
    for i in without:
        assert texts[i] == "→ Some other map"
    assert not ai_layouts.refresh_exit_labels(layout, m3)


def test_pipeline_with_ai_layout(built, tmp_path):
    """A map with an AI layout is written as the layout (Escher map, SVG, node positions); the other maps as before;
    an outdated layout is not used."""
    from openpyxl import load_workbook

    prep, xlsx, out, index = built
    md = ModelData(str(prep))
    name = index["maps"][0]["name"]
    m = final_map(out, name)
    draft, ids = export_draft(m, md, name)
    layout = make_layout(m, md, name, shifted(draft, 300, -200), ids, 2.0, svg="<svg>AI</svg>",
                         source=dict(toolkit="test"))
    folder = tmp_path / "assignment_layouts" / split.slug(name)
    layout.save(folder)
    wb = load_workbook(xlsx)
    ws = wb.create_sheet("AI layouts")
    ws.append(["Map", "Folder", "Source"])
    ws.append([name, f"assignment_layouts/{split.slug(name)}", "test"])
    wb.save(xlsx)

    edited = tmp_path / "edited"
    run_draft_map(prep, edited, "--split", "--pack", "--assignment", str(xlsx))
    prefix = split.slug(name)
    stats = json.loads((edited / "stats.json").read_text())
    assert stats["maps"][name]["arranged_by"] == "AI layout"
    assert (edited / f"{prefix}.svg").read_text() == "<svg>AI</svg>"
    written = json.loads((edited / f"{prefix}.escher.json").read_text())
    assert written[1] == layout.escher[1] and written[0]["map_name"].endswith("(BRENDA map draft)")
    laid_out = final_map(edited, name)
    for i, (x, y) in layout.positions.items():
        assert (laid_out.nodes[i].x, laid_out.nodes[i].y) == pytest.approx((x, y))
    for s in laid_out.subsystems.values():
        members = laid_out.members(s.id)
        assert all(s.x <= n.x <= s.x + s.w and s.y <= n.y <= s.y + s.h for n in members)
    for other in index["maps"][1:]:
        p = split.slug(other["name"])
        assert (edited / f"{p}.escher.json").read_bytes() == (out / f"{p}.escher.json").read_bytes()

    # outdated: a reaction less in the layout
    body = layout.escher[1]
    del body["reactions"][next(iter(body["reactions"]))]
    layout.save(folder)
    run_draft_map(prep, edited, "--split", "--pack", "--assignment", str(xlsx))
    stats = json.loads((edited / "stats.json").read_text())
    assert stats["maps"][name]["arranged_by"] != "AI layout" and name in stats["ai_layouts_outdated"]
    assert (edited / f"{prefix}.escher.json").read_bytes() == (out / f"{prefix}.escher.json").read_bytes()


# --- the harness (no Qt, no SDKs: fake clients) ---------------------------------------------------------------------

from types import SimpleNamespace as NS  # noqa: E402

from ecmapper.draft_map.ai import coding_agent, webchat  # noqa: E402
from ecmapper.draft_map.ai.agent import Agent, Events, Tool, ToolError, schema_problem  # noqa: E402
from ecmapper.draft_map.ai.providers import (  # noqa: E402
    FALLBACK_BETA,
    WEBCHAT,
    AnthropicChat,
    Chat,
    Connection,
    OpenAIChat,
    ProviderError,
    ToolCall,
    ToolResult,
    Turn,
)
from ecmapper.draft_map.ai.toolkits import LayoutJob, Toolkit  # noqa: E402
from ecmapper.draft_map.ai.webchat import WebChat  # noqa: E402
from ecmapper.draft_map.ai.workspace import Workspace, clip  # noqa: E402

FAKE_TOOLKIT = """import shutil, sys
from pathlib import Path
command, args = sys.argv[1], sys.argv[2:]
if command == "run":
    out = Path(args[1]); out.mkdir(parents=True, exist_ok=True)
    shutil.copy(args[0], out / "draft.escher_final.json")
print("fake toolkit", command, " ".join(args))
"""


def fake_toolkit(tmp_path):
    folder = tmp_path / "toolkit"
    (folder / "fake_tk").mkdir(parents=True)
    (folder / "fake_tk" / "__init__.py").write_text("")
    (folder / "fake_tk" / "__main__.py").write_text(FAKE_TOOLKIT)
    (folder / "PROMPT.md").write_text("Lay out the map nicely.")
    (folder / "ecmapper_toolkit.json").write_text(json.dumps(dict(name="fake", package="fake_tk",
                                                                  commands=["run", "inspect"])))
    return Toolkit.detect(folder)


def test_workspace(tmp_path):
    toolkit = fake_toolkit(tmp_path)
    w = Workspace(tmp_path / "work", toolkit.folder, max_output=500)
    assert "wrote" in w.write_file("a/b.txt", "one\ntwo\nthree\n")
    assert "2  two" in w.read_file("a/b.txt") and "3 lines" in w.read_file("a/b.txt")
    assert "1 more lines" in w.read_file("a/b.txt", 0, 2)
    w.edit_file("a/b.txt", "two", "2")
    assert (tmp_path / "work" / "a" / "b.txt").read_text() == "one\n2\nthree\n"
    with pytest.raises(ToolError, match="exactly once"):
        w.edit_file("a/b.txt", "e", "E")
    assert "Lay out" in w.read_file("toolkit/PROMPT.md")
    assert "fake_tk/" in w.list_files("toolkit")
    for bad in ("../x.txt", str(tmp_path / "x.txt"), "a/../../x.txt"):
        with pytest.raises(ToolError, match="outside"):
            w.write_file(bad, "x")
    with pytest.raises(ToolError, match="read-only"):
        w.write_file("toolkit/PROMPT.md", "x")
    with pytest.raises(ToolError, match="outside"):
        w.read_file("toolkit/../../secret.txt")
    assert len(clip("x" * 2000, 500)) < 600 and "left out" in clip("x" * 2000, 500)
    assert schema_problem({"path": 1}, {"type": "object", "properties": {"path": {"type": "string"}}}) == \
        "input.path must be of type string"
    assert "missing" in schema_problem({}, {"type": "object", "required": ["path"]})


def test_layout_job(tmp_path):
    toolkit = fake_toolkit(tmp_path)
    assert Toolkit.detect(toolkit.folder).package == "fake_tk"
    with pytest.raises(ValueError, match="no toolkit"):
        Toolkit.detect(tmp_path)
    kit = tmp_path / "escher_handlayout_kit"  # recognised by run.py and escher_handlayout/declspec.py
    (kit / "escher_handlayout").mkdir(parents=True)
    (kit / "run.py").write_text("")
    (kit / "escher_handlayout" / "declspec.py").write_text("")
    found = Toolkit.detect(kit)
    assert (found.package, found.result) == ("run", "**/*_handlayout.json") and "build" in found.commands
    submitted = []

    def submit(path):
        submitted.append(path)
        return "accepted"

    job = LayoutJob(toolkit, tmp_path / "work", "Central", "2 subsystems, 5 reactions",
                    [("glycolysis", "Glycolysis", 3)], {"pyr": ["Amino acids"]}, submit)
    (job.folder / "draft.escher.json").write_text("[{}, {}]")
    system = job.system_prompt()
    assert "Central" in system and "Lay out the map nicely." in system and "pyr: Amino acids" in system
    tools = {t.name: t for t in job.tools()}
    assert tools["run_toolkit"].needs_approval and not tools["submit_layout"].needs_approval
    with pytest.raises(ToolError, match="allowed commands"):
        job.run_toolkit(["render", "x"])
    out = job.run_toolkit(["run", "draft.escher.json", "v1"])
    assert out.startswith("exit code 0") and "fake toolkit run" in out
    final = toolkit.newest_result(job.folder)
    assert final == job.folder / "v1" / "draft.escher_final.json"
    with pytest.raises(ToolError, match="newest result is v1/draft.escher_final.json"):
        job.submit_layout("v2/x.json")
    assert job.submit_layout("v1/draft.escher_final.json") == "accepted" and submitted == [final]


def test_agent_folder(tmp_path):
    toolkit = fake_toolkit(tmp_path)
    base = tmp_path / "ai" / "central"
    ids = {"1": ["node", "pyr_n", "pyr_c"], "2": ["reaction", "e1", "PYK"]}
    folder = base / "20261104-120000-agent"
    folder.mkdir(parents=True)
    (folder / "draft.escher.json").write_text("[{}, {}]")
    assert coding_agent.newest_job(base, "Central") is None
    coding_agent.prepare(folder, toolkit, "Central", "1 subsystems, 3 reactions", [("glycolysis", "Glycolysis", 3)],
                        {"pyr": ["Amino acids"]}, ids)
    task = (folder / coding_agent.TASK).read_text(encoding="utf-8")
    assert "Central" in task and "PROMPT.md" in task and "glycolysis: Glycolysis, 3 reactions" in task
    assert "pyr: Amino acids" in task and "run_toolkit.py" in task and "Use coding agent result" in task
    assert coding_agent.newest_job(base, "Central") == (folder, coding_agent.read_job(folder), None)
    assert coding_agent.newest_job(base, "Other") is None
    # the runner runs the toolkit in the folder, from anywhere
    done = subprocess.run([sys.executable, str(folder / coding_agent.RUNNER), "run", "draft.escher.json", "v1"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0 and "fake toolkit run" in done.stdout
    final = folder / "v1" / "draft.escher_final.json"
    assert final.is_file()
    found_folder, job, result = coding_agent.newest_job(base, "Central")
    assert (found_folder, result, job["ids"], job["map"]) == (folder, final, ids, "Central")
    assert coding_agent.find_job(final) == (folder.resolve(), job)
    assert coding_agent.find_job(tmp_path / "x.json") is None
    assert (folder / coding_agent.CLAUDE).read_text(encoding="utf-8").strip() == "@AGENTS.md"
    claude, pi = (coding_agent.command(folder, toolkit, agent) for agent in ("Claude Code", "pi"))
    assert claude.startswith(f'cd "{folder}"') and f'--add-dir "{toolkit.folder}"' in claude
    assert pi.splitlines()[1] == 'pi "Lay out this map."'


class ScriptedChat(Chat):
    """Returns the given turns; records what it is sent."""

    def __init__(self, turns):
        super().__init__(Connection("test"), "system", [])
        self.turns, self.users, self.results = list(turns), [], []

    def add_user(self, text):
        self.users.append(text)

    def add_results(self, results):
        self.results.append(results)

    def step(self, on_text=None, on_thinking=None, stop=None):
        turn = self.turns.pop(0) if self.turns else Turn("again", [ToolCall("x", "echo", {"text": "x"})], "tool_use")
        if on_text:
            on_text(turn.text)
        return turn


class Recorder(Events):
    def __init__(self, approve=True, cancel_after=None):
        self.texts, self.finished, self.approve_answer, self.cancel_after = [], [], approve, cancel_after

    def text(self, delta):
        self.texts.append(delta)

    def tool_finished(self, call, result):
        self.finished.append((call.name, result.is_error))

    def approve(self, call):
        return self.approve_answer

    def cancelled(self):
        return self.cancel_after is not None and len(self.finished) >= self.cancel_after


def test_agent(tmp_path):
    ran = []
    tools = [Tool("echo", "", {"type": "object", "properties": {"text": {"type": "string"}}, "required": ["text"]},
                  lambda text: (ran.append(text), f"echo {text}")[1]),
             Tool("danger", "", {"type": "object"}, lambda: ran.append("danger") or "done", needs_approval=True),
             Tool("broken", "", {"type": "object"}, lambda: 1 / 0)]
    calls = [ToolCall("1", "echo", {"text": "a"}), ToolCall("2", "danger", {}), ToolCall("3", "echo", {"text": 5}),
             ToolCall("4", "nope", {}), ToolCall("5", "broken", {}), ToolCall("6", "echo", {}, error="cut off")]
    chat = ScriptedChat([Turn("working", calls, "tool_use"), Turn("done", [], "end")])
    events = Recorder(approve=False)
    agent = Agent(chat, tools, transcript=tmp_path / "t.jsonl")
    assert agent.run("lay it out", events) == "done"
    assert chat.users == ["lay it out"] and ran == ["a"] and events.texts == ["working", "done"]
    (results,) = chat.results  # all results of a turn in one message
    assert [r.id for r in results] == ["1", "2", "3", "4", "5", "6"]
    assert [r.is_error for r in results] == [False, True, True, True, True, True]
    assert "declined" in results[1].text and "type string" in results[2].text and "Unknown tool" in results[3].text
    assert "ZeroDivisionError" in results[4].text and "cut off" in results[5].text
    assert len((tmp_path / "t.jsonl").read_text().splitlines()) == 4
    # approved; stopped by the user after the first tool: the other calls get results, the run stops
    chat = ScriptedChat([Turn("", [ToolCall("1", "danger", {}), ToolCall("2", "echo", {"text": "b"})], "tool_use")])
    assert Agent(chat, tools).run("go", Recorder(cancel_after=1)) == "stopped"
    assert [r.is_error for r in chat.results[0]] == [False, True] and "Stopped" in chat.results[0][1].text
    assert Agent(ScriptedChat([]), tools, max_steps=3).run("forever") == "max_steps"
    assert Agent(ScriptedChat([Turn("no", [], "refusal", detail="cyber")]), tools).run("x") == "refusal: cyber"


class FakeStream:
    def __init__(self, events, message):
        self.events, self.message = events, message

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        return iter(self.events)

    def get_final_message(self):
        return self.message


def test_anthropic_chat():
    sent = []
    content = [NS(type="thinking", thinking="plan"), NS(type="text", text="Hi"),
               NS(type="tool_use", id="t1", name="read_file", input={"path": "a"})]
    message = NS(stop_reason="tool_use", content=content,
                 usage=NS(input_tokens=10, output_tokens=5, cache_read_input_tokens=3, cache_creation_input_tokens=0))
    events = [NS(type="content_block_delta", delta=NS(type="thinking_delta", thinking="plan")),
              NS(type="content_block_delta", delta=NS(type="text_delta", text="Hi"))]

    def stream(**kwargs):
        sent.append(dict(kwargs, messages=list(kwargs["messages"])))
        return FakeStream(events, message)

    client = NS(beta=NS(messages=NS(stream=stream)))
    tool = Tool("read_file", "Read", {"type": "object"}, lambda path: path)
    chat = AnthropicChat(Connection("c", model="claude-opus-5-5", effort="high"), "SYSTEM", [tool], client)
    chat.add_user("lay out")
    texts, thoughts = [], []
    turn = chat.step(texts.append, thoughts.append)
    assert turn.stop == "tool_use" and turn.text == "Hi" and turn.thinking == "plan" and texts == ["Hi"]
    assert turn.calls == [ToolCall("t1", "read_file", {"path": "a"})] and turn.usage["cache_read_input_tokens"] == 3
    request = sent[0]
    assert request["model"] == "claude-opus-5-5" and request["system"] == "SYSTEM"
    assert request["thinking"] == {"type": "adaptive", "display": "summarized"}
    assert request["output_config"] == {"effort": "high"} and request["cache_control"] == {"type": "ephemeral"}
    assert request["betas"] == [FALLBACK_BETA] and request["fallbacks"] == "default"
    assert request["tools"][0]["eager_input_streaming"] and request["tools"][0]["input_schema"] == {"type": "object"}
    assert chat.messages[1] == {"role": "assistant", "content": message.content}  # appended unchanged
    chat.add_results([ToolResult("t1", "text of a", [b"PNG"]), ToolResult("t2", "", is_error=True)])
    results = chat.messages[-1]
    assert results["role"] == "user" and [b["tool_use_id"] for b in results["content"]] == ["t1", "t2"]
    assert results["content"][0]["content"][1]["type"] == "image" and results["content"][1]["is_error"]
    # no fallbacks for other models; refusals and cut-off calls
    chat = AnthropicChat(Connection("c", model="claude-haiku-5-5"), "S", [tool], client)
    assert "fallbacks" not in chat.request() and "betas" not in chat.request()
    message.stop_reason, message.stop_details = "refusal", NS(category="cyber")
    turn = chat.step()
    assert turn.stop == "refusal" and turn.detail == "cyber" and chat.messages == []
    message.stop_reason = "max_tokens"
    assert chat.step().calls[0].error


def test_openai_chat():
    sent = []

    def chunk(content=None, tool_calls=None, finish=None, usage=None, choices=True):
        delta = NS(content=content, tool_calls=tool_calls)
        return NS(choices=[NS(delta=delta, finish_reason=finish)] if choices else [], usage=usage)

    def call(index, id=None, name=None, arguments=None):
        return NS(index=index, id=id, function=NS(name=name, arguments=arguments))

    chunks = [chunk("Look"), chunk(tool_calls=[call(0, "c1", "read_file", '{"pa')]),
              chunk(tool_calls=[call(0, None, None, 'th": "a"}')]),
              chunk(tool_calls=[call(1, "c2", "view_image", "{bad")], finish="tool_calls"),
              chunk(choices=False, usage=NS(prompt_tokens=7, completion_tokens=3))]

    def create(**kwargs):
        sent.append(kwargs)
        return iter(chunks)

    client = NS(chat=NS(completions=NS(create=create)))
    tool = Tool("read_file", "Read", {"type": "object"}, lambda path: path)
    chat = OpenAIChat(Connection("o", kind="openai", model="local-model"), "SYSTEM", [tool], client)
    chat.add_user("lay out")
    turn = chat.step()
    assert turn.stop == "tool_use" and turn.text == "Look" and turn.usage == dict(input_tokens=7, output_tokens=3)
    assert turn.calls[0] == ToolCall("c1", "read_file", {"path": "a"}) and turn.calls[1].error
    assert sent[0]["tools"][0]["function"]["name"] == "read_file" and sent[0]["stream_options"]
    assert chat.messages[0] == {"role": "system", "content": "SYSTEM"}
    assistant = chat.messages[-1]
    assert assistant["tool_calls"][0]["function"]["arguments"] == '{"path": "a"}'
    chat.add_results([ToolResult("c1", "text", [b"PNG"]), ToolResult("c2", "bad", is_error=True)])
    assert chat.messages[-3] == {"role": "tool", "tool_call_id": "c1", "content": "text"}
    assert chat.messages[-2]["content"].startswith("Error:")
    assert chat.messages[-1]["role"] == "user" and chat.messages[-1]["content"][1]["type"] == "image_url"
    chat.connection.images = False
    chat.add_results([ToolResult("c3", "t", [b"PNG"])])
    assert chat.messages[-1]["role"] == "tool"


def test_webchat_protocol():
    reply = ('I look first.\n```json\n{"tool": "list_files", "input": {"path": "."}}\n```\njson\nCopy code\n'
             '{"tool": "read_file", "input": {"path": "a.txt", "limit": 5}} and {"x": 1}\n'
             '```json\n{"tool": "run_toolkit", "input": {"args": ["run", }\n```\n\n'
             '{"tool": "edit_file", "input": "nope"}')
    text, calls = webchat.parse(reply, 3)
    assert [c.id for c in calls] == ["3.1", "3.2", "3.3", "3.4"]
    assert calls[0] == ToolCall("3.1", "list_files", {"path": "."})
    assert calls[1].input == {"path": "a.txt", "limit": 5} and calls[1].error is None
    assert calls[2].name == "run_toolkit" and "not valid JSON" in calls[2].error
    assert calls[3].name == "edit_file" and "must be a JSON object" in calls[3].error
    assert text.startswith("I look first.") and '"tool"' not in text and '{"x": 1}' in text and "```" not in text
    assert webchat.parse("All done.", 4) == ("All done.", [])
    # the reply in the text of the page: after the message's last line, up to the end line
    page = (f"ecmapper message 1\n{webchat.PROTOCOL[:200]}\n{webchat.tail(1)}\nCalling\njson\nCopy\n{{...}}\n"
            f"@@END 1@@\necmapper message 2\nresults\n{webchat.tail(2)}\nDone so far")
    assert webchat.reply_of(page, 1) == ("Calling\njson\nCopy\n{...}", True)
    assert webchat.reply_of(page, 2) == ("Done so far", False)
    assert webchat.reply_of(page + " @@END 2@@\nCopy", 2) == ("Done so far", True)
    assert webchat.reply_of("Hi @@END 1@@ pasted.txt\nHello\n@@END 2@@", 2) == ("pasted.txt\nHello", True)
    assert webchat.reply_of("Hello\n@@END 12@@", 1) == ("Hello\n@@END 12@@", False)
    # chats whose terms forbid automated use
    assert "Anthropic" in webchat.blocked("https://claude.ai/new") and webchat.blocked("chatgpt.com")
    assert webchat.blocked("https://x.gemini.google.com/app") and webchat.blocked("https://notclaude.ai") is None
    assert webchat.blocked("https://chat.example.org/c/1") is None
    with pytest.raises(ProviderError, match="Anthropic"):
        webchat.check_connection(Connection("w", kind=WEBCHAT, base_url="https://claude.ai", transport="browser"))
    with pytest.raises(ProviderError, match="URL"):
        webchat.check_connection(Connection("w", kind=WEBCHAT, transport="browser"))
    webchat.check_connection(Connection("w", kind=WEBCHAT, base_url="https://claude.ai", transport="clipboard"))


class FakeTransport:
    """Gives the replies in turn (N@@ in them: the number of the message); records the messages."""

    images = False

    def __init__(self, replies):
        self.replies, self.messages = list(replies), []

    def exchange(self, message, number, images, stop=None):
        self.messages.append(message)
        return self.replies.pop(0).replace("N@@", f"{number}@@") if self.replies else None


def test_webchat(tmp_path):
    toolkit = fake_toolkit(tmp_path)
    submitted = []
    job = LayoutJob(toolkit, tmp_path / "work", "Central", "1 subsystems, 3 reactions",
                    [("glycolysis", "Glycolysis", 3)], {}, lambda path: submitted.append(path) or "accepted")
    (job.folder / "draft.escher.json").write_text("[{}, {}]")
    transport = FakeTransport([
        'I run it.\n```json\n{"tool": "run_toolkit", "input": {"args": ["run", "draft.escher.json", "v1"]}}\n```\n'
        "@@END N@@",
        '```json\n{"tool": "submit_layout", "input": {"path": "v1/draft.escher_final.json"}}\n```',
        "Done: it is submitted.\n@@END N@@"])
    chat = WebChat(Connection("w", kind=WEBCHAT, transport="clipboard"), job.system_prompt(), job.tools(), transport)
    events = Recorder()
    assert Agent(chat, job.tools()).run("Lay out this map.", events) == "done"
    first, second, third = transport.messages
    assert first.startswith("ecmapper message 1\n") and first.endswith(webchat.tail(1))
    assert "Lay out the map nicely." in first and '"tool": "list_files"' in first and "- submit_layout:" in first
    assert "=== The user ===\nLay out this map." in first and "How to call the tools" not in second
    assert "[call 1.1: run_toolkit]\nexit code 0" in second and second.endswith(webchat.tail(2))
    assert "[call 2.1: submit_layout]\naccepted" in third
    assert submitted == [job.folder / "v1" / "draft.escher_final.json"]
    assert events.texts == ["I run it.", "Done: it is submitted."] and chat.history[2][1] == "Done: it is submitted."
    chat.add_results([ToolResult("4.1", "a render", [b"PNG"])])
    assert "1 image(s)" in chat.pending[-1] and chat.images == []
    assert Agent(chat, job.tools()).run("Go on.") == "stopped"  # the transport stopped (no reply)
