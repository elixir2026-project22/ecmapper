import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ecmapper.draft_map import merge_escher
from ecmapper.draft_map.core import DraftMap, Edge, ModelData, Node, Subsystem


def two_pathways():
    """Subsystems S1 (a -R1-> b) and S2 (b -R2-> c), side by side with b in both, and S3 without reactions."""
    md = ModelData(dict(
        metabolites={m: dict(name=m.upper(), compartment="c", keys={m.upper()}, currency=False)
                     for m in ("a_c", "b_c", "c_c")},
        reactions={"R1": dict(name="R1", reversible=False, metabolites={"a_c": -1, "b_c": 1}),
                   "R2": dict(name="R2", reversible=False, metabolites={"b_c": -1, "c_c": 1})},
        node_met_keys={}, matches={}))
    dm = DraftMap()
    for k, label, x in (("S1", "Upper part", 0), ("S2", "Merged pathway", 400), ("S3", "Empty", 800)):
        dm.subsystems[k] = Subsystem(k, label, "", "#ffffff", x, 0, 300, 100)

    def node(i, kind, k, x, mets=(), reactions=()):
        dm.nodes[i] = Node(i, kind, k, x, 50, 20, 20, i, reactions=list(reactions), model_mets=set(mets))

    for k, x, (s, p), rid in (("S1", 0, ("a", "b"), "R1"), ("S2", 400, ("b", "c"), "R2")):
        node(f"{k}s", "main", k, x + 40, {s})
        node(f"{k}e", "enzyme", k, x + 150, reactions=[rid])
        node(f"{k}p", "main", k, x + 260, {p})
        dm.edges += [Edge(f"{k}e", f"{k}s", "main", False, False), Edge(f"{k}e", f"{k}p", "main", False, True)]
    node("S3m", "main", "S3", 850, {"a"})
    return dm, md


def test_subsystem_maps(tmp_path):
    dm, md = two_pathways()
    maps = merge_escher.subsystem_maps(dm, md)
    assert list(maps) == ["S1", "S2"]  # S3 has no reactions
    for k, data in maps.items():
        head, body = data
        assert head["map_name"] == dm.subsystems[k].label
        assert len(body["reactions"]) == 1
        # smooth: the main lines are curves
        for r in body["reactions"].values():
            for s in r["segments"].values():
                ends = {body["nodes"][s[e]]["node_type"] for e in ("from_node_id", "to_node_id")}
                if "metabolite" in ends:
                    assert s["b1"] and s["b2"]
    paths = merge_escher.write_maps(maps, tmp_path / "maps")
    assert [p.name for p in paths] == ["upper_part.json", "subsystem_merged_pathway.json"]
    assert json.loads(paths[0].read_text(encoding="utf-8")) == maps["S1"]
    (tmp_path / "maps" / "old.json").write_text("[]")
    merge_escher.write_maps({"S1": maps["S1"], "S9": maps["S1"]}, tmp_path / "maps")
    assert sorted(p.name for p in (tmp_path / "maps").iterdir()) == ["upper_part.json", "upper_part_2.json"]



def test_subsystem_maps_brenda_style():
    dm, md = two_pathways()
    # a second model reaction at the enzyme node of S1 (drawn beside R1), with a proton
    md.metabolites["h_c"] = dict(name="H+", compartment="c", keys={"H"}, currency=True)
    md.reactions["R1b"] = dict(name="R1 (proton)", reversible=True, metabolites={"a_c": -1, "h_c": -1, "b_c": 1},
                               genes=[("g1", "G1")], gpr="g1")
    md.base["h_c"] = "h"
    dm.nodes["S1e"].reactions.append("R1b")
    with pytest.raises(ValueError):
        merge_escher.subsystem_maps(dm, md, style="other")
    maps = merge_escher.subsystem_maps(dm, md, style="brenda")
    assert list(maps) == ["S1", "S2"]
    head, body = maps["S1"]
    assert head["map_name"] == "Upper part"
    reactions = {r["bigg_id"]: r for r in body["reactions"].values()}
    assert sorted(reactions) == ["R1", "R1b"]
    r = reactions["R1b"]
    assert r["reversibility"] and r["gene_reaction_rule"] == "g1" and r["genes"] == [dict(bigg_id="g1", name="G1")]
    assert {m["bigg_id"]: m["coefficient"] for m in r["metabolites"]} == {"a_c": -1, "h_c": -1, "b_c": 1}
    nodes = body["nodes"]
    drawn = {n["bigg_id"] for n in nodes.values() if n["node_type"] == "metabolite"}
    assert drawn == {"a_c", "b_c"}  # one primary node each, shared by both reactions; no proton
    assert all(n["node_is_primary"] for n in nodes.values() if n["node_type"] == "metabolite")
    mids = [(n["x"], n["y"]) for n in nodes.values() if n["node_type"] == "midmarker"]
    assert len(mids) == 2 and mids[0][0] == mids[1][0] and abs(mids[0][1] - mids[1][1]) == 48  # side by side
    assert [t["text"] for t in body["text_labels"].values()] == ["Upper part"]


def test_find_tool_and_command(tmp_path):
    assert merge_escher.find_tool("") is None and merge_escher.find_tool(tmp_path) is None
    script = tmp_path / "src" / merge_escher.SCRIPT
    script.parent.mkdir()
    script.write_text("")
    for location in (tmp_path, tmp_path / "src", script):
        assert merge_escher.find_tool(str(location)) == script
    command = merge_escher.merge_command(script, tmp_path / "maps", tmp_path / "out.json", "Map")
    assert command[:3] == [sys.executable, "-u", str(script)] and "--model" not in command
    assert command[command.index("--id-style") + 1] == "none" and command[command.index("-o") + 1].endswith("out.json")
    model = tmp_path / "model.json"
    model.write_text("{}")
    assert "--model" in merge_escher.merge_command(script, "maps", "out.json", "Map", model)
    assert "--model" not in merge_escher.merge_command(script, "maps", "out.json", "Map", tmp_path / "model.xml")


def test_merge_with_merging_escher(tmp_path):
    """The real tool (its folder in ECMAPPER_MERGING_ESCHER): the two pathways joined at b."""
    script = merge_escher.find_tool(os.environ.get("ECMAPPER_MERGING_ESCHER"))
    if script is None:
        pytest.skip("MergingEscher not found (set ECMAPPER_MERGING_ESCHER to its folder)")
    dm, md = two_pathways()
    folder, output = tmp_path / "maps", tmp_path / "merged" / "map_merged.json"
    merge_escher.write_maps(merge_escher.subsystem_maps(dm, md), folder)
    done = subprocess.run(merge_escher.merge_command(script, folder, output, "Two pathways"), capture_output=True,
                          text=True, timeout=300)
    assert done.returncode == 0, done.stdout + done.stderr
    head, body = json.loads(Path(output).read_text())
    assert head["map_name"] == "Two pathways"
    assert sorted(r["bigg_id"] for r in body["reactions"].values()) == ["R1", "R2"]
    primaries = [n["bigg_id"] for n in body["nodes"].values() if n.get("node_is_primary")]
    assert primaries.count("b_c") == 1 and {"a_c", "c_c"} <= set(primaries)
