import json
import re
import shutil

import pytest

from ecmapper.pwl import PWLError, detect_motifs, parse
from ecmapper.pwl.__main__ import main as pwl_main
from ecmapper.pwl.compounds import CompoundRegistry, parse_equation, slug, split_location
from ecmapper.pwl.motifs import check_motifs
from ecmapper.sources import BKMSSource

EXAMPLE = """\
pwl 2.0;
subsystem upper-glycolysis 'upper glycolysis' {
  flow down;
  compartment in 'cytosol';
  compartment out;
  compound Glc 'D-glucose' brenda:35 kegg:C00031;
  compound G6P 'D-glucose 6-phosphate';
  compound F6P;
  compound FBP;
  compound DHAP;
  compound GAP;
  compound ATP;
  compound ADP;
  compound H2O;
  compound Pi;
  compound H 'H+';
  port G6P side R to pentose-phosphate-pathway;

  Glc@out -[t1: 2.7.1.199 | -H@out, +H@in]-> Glc#2@in;  // a transporter
  Glc#2 -[r1: 2.7.1.1 | ATP>ADP {br=BR1, ev=cof+eq}]=> G6P
    <-R[r2: 5.3.1.9]-> F6P
    -[r3: -> 2.7.1.11 | ATP>ADP || r4: <- 3.1.3.11 | -H2O, +Pi ^L]- FBP;
  FBP -[r5: 4.1.2.13]-> (DHAP + GAP);
  ring cw {
    DHAP -[r6: 5.3.1.1]- GAP;
    chord GAP -[r7: ~? | -2*H2O]- DHAP;
    react [r8: 1.1.1.1] {GAP, DHAP, F6P};
  }
  star F6P arc 90 order cw {
    -[r9: 1.1.1.2]- _;
  }
  motif substrate_cycle { F6P FBP } forward { r3 } backward { r4 };
  motif transport { Glc@out Glc@in } via { t1 };
  motif pathway name='upper part' { r1 r2 r3 };
  effector ATP inhibits r3, r1;
}
"""


def test_parse_example():
    doc = parse(EXAMPLE)
    assert doc.id == "upper-glycolysis" and doc.title == "upper glycolysis" and doc.flow == "down"
    assert doc.compounds["Glc"].xref("kegg") == ["C00031"] and doc.compounds["Glc"].label == "D-glucose"
    rs = doc.reaction_index()
    assert len(rs) == 10
    assert rs["t1"].substrates[0].species == "Glc@out" and rs["t1"].products[0].node == "Glc#2"
    assert rs["r1"].reversible is False and [p.token for p in rs["r1"].substrates] == ["Glc"]
    assert rs["r1"].attrs == {"br": "BR1", "ev": "cof+eq"}
    assert rs["r2"].reversible is True
    assert [p.token for p in rs["r3"].substrates] == ["F6P"] and [p.token for p in rs["r4"].substrates] == ["FBP"]
    assert rs["r4"].cofactor_side == "L" and [p.token for p in rs["r4"].consumed] == ["H2O"]
    assert [p.token for p in rs["r5"].products] == ["DHAP", "GAP"]
    assert rs["r6"].direction is None and rs["r6"].substrates is None and rs["r6"].sided
    assert rs["r7"].spontaneous and rs["r7"].consumed[0].coef == "2" and rs["r7"].context == "chord"
    assert not rs["r8"].sided and [p.token for p in rs["r8"].unsided] == ["GAP", "DHAP", "F6P"]
    assert rs["r9"].right == [] and rs["r9"].context == "star"
    assert rs["t1"].comment == "a transporter"
    assert doc.effectors[0].effect == "inhibits" and doc.effectors[0].reactions == ["r3", "r1"]
    assert doc.motifs[2].params == {"name": "upper part"}
    assert check_motifs(doc) == []


def test_detect_motifs_in_hand_written_document():
    kinds = {m.kind: m for m in detect_motifs(parse(EXAMPLE))}
    assert kinds["substrate_cycle"].roles == {"forward": ["r3"], "backward": ["r4"]}
    assert kinds["transport"].members == ["Glc@in", "Glc@out", "H@in", "H@out"]  # a symport


@pytest.mark.parametrize(
    "change, message",
    [
        (("pwl 2.0;\n", ""), "starts with 'pwl 2.0;'"),
        (("pwl 2.0", "pwl 3.0"), "unsupported PWL version"),
        (("-[r6: 5.3.1.1]-", "<=[r6: 5.3.1.1]->"), "undefined link arrows"),
        (("-[r6: 5.3.1.1]-", "-[r6: 5.3.1.1"), "malformed link"),
        (("]=> G6P", "]=> G6P ]"), "unexpected character ']'"),
        (("-[r3: -> 2.7.1.11 | ATP>ADP || r4: <- 3.1.3.11 | -H2O, +Pi ^L]- FBP",
          "-[r3: -> 2.7.1.11 | ATP>ADP]-> FBP"), "plain link"),
        (("compound F6P;", "compound ring;"), "keyword"),
        (("compound F6P;", ""), "'F6P' is not declared"),
        (("compound F6P;", "compound F6P;\n  compound F6P;"), "declared twice"),
        (("backward { r4 }", "backward { r42 }"), "unknown reaction r42"),
        (("react [r8: 1.1.1.1]", "react [r8: -> 1.1.1.1]"), "no direction"),
        (("[r9: 1.1.1.2]", "[r9: 1.1.1.2 ^L]"), "without cofactors"),
        (("ring cw {\n", "ring cw {\n    flow up;\n"), "not allowed inside a ring"),
        (("r6: 5.3.1.1", "r1: 5.3.1.1"), "used twice"),
        (("-2*H2O", "-x*H2O"), "bad coefficient"),
        (("Glc#2 -[r1", "Glc#1 -[r1"), "copy numbers start at #2"),
        (("(DHAP + GAP)", "(DHAP)"), "at least two members"),
        (("H@in]", "H@nowhere]"), "compartment 'nowhere' is not declared"),
        (("\n}\n", "\n}\nflow up;\n"), "after the end"),
        (("ev=cof+eq}", "ev cof}"), "bad attribute"),
    ],
)
def test_parse_errors(change, message):
    old, new = change
    assert old in EXAMPLE
    with pytest.raises(PWLError, match=re.escape(message)):
        parse(EXAMPLE.replace(old, new, 1))


def test_error_position():
    with pytest.raises(PWLError) as e:
        parse(EXAMPLE.replace("-[r6: 5.3.1.1]-", "-[r6: 5.3.1.1]- $"))
    assert e.value.line == EXAMPLE.split("\n").index("    DHAP -[r6: 5.3.1.1]- GAP;") + 1


def test_motif_check_reports_unsupported_motif():
    extra = "  motif cycle { F6P FBP DHAP };\n  motif repeat period=2 count=2 { r1 r2 r3 };\n"
    text = EXAMPLE.replace("  effector", extra + "  effector")
    problems = check_motifs(parse(text))
    assert any("no reaction converts DHAP into F6P" in p for p in problems)
    assert any("expected period*count = 4" in p for p in problems)


# ---------------------------------------------------------------------------------- compounds


def test_names_and_equations():
    assert slug("NAD+") == "NAD" and slug("D-Glc 6-P") == "D_Glc_6_P" and slug("ring") == "ring_c"
    assert split_location("H+[side 2]") == ("H+", "side2") and split_location("Na+/in") == ("Na+", "in")
    lhs, rhs = parse_equation("2 ubiquinol + O2 + n H+[side 1] = 2 ubiquinone + 2 H2O + n H+[side 2]")
    assert [(t.coef, t.name, t.location) for t in lhs] == [
        ("2", "ubiquinol", None), (None, "O2", None), ("n", "H+", "side1")]
    assert rhs[2].location == "side2"
    assert parse_equation("A <=> B") is None


def _met(node, display, brenda=None, lig=None):
    return {"node_id": node, "display_name": display, "brenda_name": brenda or display, "brenda_ligand_id": lig,
            "node_class": "nodeHM"}


def test_registry_merges_ligand_ids_and_keeps_names_apart():
    reactions = []
    for k in range(5):
        reactions.append({"main_metabolites": [_met(f"m{k}", "H+", lig=12), _met(f"n{k}", "NADH", lig=8)],
                          "cofactors": {"consumed": [_met(f"o{k}", "O2", lig=9 if k else 1180)], "produced": []}})
    # one drawn NADH carries the ligand ID of H+
    reactions.append({"main_metabolites": [_met("x", "NADH", lig=12)], "cofactors": {"consumed": [], "produced": []}})
    reg = CompoundRegistry({"reactions": reactions})
    nadh = reg.by_token("NADH")
    assert nadh.ligand_ids == {8} and nadh.reassigned == {(12, "NADH"): 1}
    assert reg.by_token("H").ligand_ids == {12} and reg.by_token("H").label == "H+"
    assert reg.by_token("O2").ligand_ids == {9, 1180}


def test_bkms_equations(cache_dir):
    bkms = BKMSSource(cache_dir=cache_dir)
    assert bkms.equations("BR19236")[0] == ("3.5.1.14", "N2-acetyl-L-ornithine + H2O = acetate + L-ornithine")
    assert "N2-acetyl-L-ornithine + H2O = acetate + L-ornithine" in bkms.ec_equations("3.5.1.16")
    assert bkms.equations("nothing") == []


# ---------------------------------------------------------------------------------- generator


class FakeBKMS:
    def __init__(self, equations):
        self._eq = equations  # {BRENDA reaction ID: [(EC, equation)]}

    def equations(self, br):
        return list(self._eq.get(br, []))

    def ec_equations(self, ec):
        return [eq for items in self._eq.values() for e, eq in items if e == ec]


def _reaction(node, ec, mains, consumed=(), produced=(), brs=(), x=0, y=0, name=""):
    return {
        "node_id": node, "subsystem": "mini", "ec_label": ec, "ec_numbers": [ec], "enzyme_name": name or ec,
        "brenda_reaction_ids": list(brs), "layout": {"x": x, "y": y, "width": 60, "height": 30},
        "main_metabolites": [dict(_met(n, d, b, lig), position={"x": px, "y": py}) for n, d, b, lig, px, py in mains],
        "cofactors": {"consumed": [_met(node + "c" + d, d, lig=lig) for d, lig in consumed],
                      "produced": [_met(node + "p" + d, d, lig=lig) for d, lig in produced], "undirected": []},
    }


def mini_map():
    # (node, display name, BRENDA name, ligand ID, x, y); the equations call DHAP "glycerone phosphate"
    F6P = ("m1", "Fru 6-P", "D-fructose 6-phosphate", 87, 0, 0)
    FBP = ("m2", "FBP", "fructose-1,6-bisphosphate", 6523, 0, 200)
    DHAP = ("m3", "DHAP", "dihydroxyacetone phosphate", 278, -100, 400)
    GAP = ("m4", "GAP", "D-glyceraldehyde 3-phosphate", 147, 100, 400)
    A, B, C = ("m5", "A", "A", 501, 400, 0), ("m6", "B", "B", 502, 600, 200), ("m7", "C", "C", 503, 400, 400)
    X, Y = ("m8", "X", "X", 504, 800, 0), ("m9", "Y", "Y", 505, 800, 200)
    atp, adp, nad, nadh, h2o, pi = ("ATP", 4), ("ADP", 13), ("NAD+", 7), ("NADH", 8), ("H2O", 1), ("Pi", 16)
    reactions = [
        _reaction("node10", "2.7.1.11", [F6P, FBP], [atp], [adp], ["BR1"], -30, 100, "6-phosphofructokinase"),
        _reaction("node11", "3.1.3.11", [FBP, F6P], [h2o], [pi], ["BR2"], 30, 100),
        _reaction("node12", "4.1.2.13", [FBP, DHAP, GAP], brs=["BR3"], x=0, y=300),
        _reaction("node13", "5.3.1.1", [GAP, DHAP], brs=["BR4"], x=0, y=400),
        _reaction("node14", "1.1.1.1", [X, Y], x=800, y=100),
        _reaction("node15", "1.1.1.10", [A, B], [nad], [nadh], ["BR6"], 500, 100),
        _reaction("node16", "1.1.1.11", [B, C], [nadh], [nad], ["BR7"], 500, 300),
        _reaction("node17", "4.2.1.1", [C, A], [], [h2o], ["BR8"], 400, 200),
    ]
    nodes = {}
    for r in reactions:
        for m in r["main_metabolites"]:
            nodes[m["node_id"]] = dict(m, subsystem="mini", x=m["position"]["x"], y=m["position"]["y"], width=40,
                                       height=30)
    return {
        "source": "test",
        "subsystems": [{"id": "mini", "label_text": "mini map", "bbox": {"x": -200, "y": -100, "width": 1200,
                                                                          "height": 600}},
                       {"id": "other", "label_text": "other", "bbox": {"x": 0, "y": -900, "width": 10, "height": 10}}],
        "reactions": reactions,
        "metabolite_nodes": list(nodes.values()),
        "metabolite_links": [{"start": [0, 0], "end": [0, -500], "start_node_id": "m1", "end_node_id": "zz",
                              "start_subsystem": "mini", "end_subsystem": "other"}],
        "unassigned_links": [],
    }


MINI_BKMS = FakeBKMS({
    "BR1": [("2.7.1.11", "ATP + D-fructose 6-phosphate = ADP + D-fructose 1,6-bisphosphate")],
    "BR2": [("3.1.3.11", "D-fructose 1,6-bisphosphate + H2O = D-fructose 6-phosphate + phosphate")],
    "BR3": [("4.1.2.13", "D-fructose 1,6-bisphosphate = glycerone phosphate + D-glyceraldehyde 3-phosphate")],
    "BR4": [("5.3.1.1", "D-glyceraldehyde 3-phosphate = glycerone phosphate")],
    "BR4r": [("5.3.1.1", "glycerone phosphate = D-glyceraldehyde 3-phosphate")],
    "BR6": [("1.1.1.10", "A + NAD+ = B + NADH + H+")],
    "BR7": [("1.1.1.11", "B + NADH = C + NAD+")],
    "BR8": [("4.2.1.1", "C = A + H2O")],
})


@pytest.fixture
def generated(tmp_path):
    pytest.importorskip("networkx")
    pytest.importorskip("numpy")
    from ecmapper.pwl.brenda import generate

    data = mini_map()
    (tmp_path / "map.json").write_text(json.dumps(data))
    out = tmp_path / "pwl"
    summary = generate(data, out, MINI_BKMS)
    return data, out, summary


def test_generate_mini_map(generated):
    data, out, summary = generated
    text = (out / "subsystems" / "mini.pwl").read_text()
    doc = parse(text)
    rs = doc.reaction_index()
    # DHAP is "glycerone phosphate" in the equations: the only drawn metabolite and term left, so they match
    assert [p.token for p in rs["node10"].substrates] == ["Fru_6_P"] and rs["node10"].attrs["ev"] == "cof+eq"
    assert [p.token for p in rs["node11"].substrates] == ["FBP"]
    assert "FBP -[node12: 4.1.2.13 {br=BR3, ev=eq, match=" in text and "]-> (DHAP + GAP);" in text
    assert rs["node13"].reversible is True and rs["node13"].attrs["rev"] == "eq"
    assert rs["node14"].direction is None and "br" not in rs["node14"].attrs
    assert doc.compounds["FBP"].xref("brenda") == ["6523"] and doc.compounds["NAD"].label == "NAD+"
    assert "port Fru_6_P side U to other;" in text
    kinds = {}
    for m in doc.motifs:
        kinds.setdefault(m.kind, []).append(m)
    assert kinds["substrate_cycle"][0].roles == {"forward": ["node11"], "backward": ["node10"]}
    assert kinds["cycle"][0].members == ["A", "B", "C"]
    assert kinds["carrier"][0].members == ["NAD", "NADH"]
    assert summary["totals"]["oriented"] == 7 and summary["totals"]["reversible"] == 1
    index = json.loads((out / "node_index.json").read_text())
    assert index["mini"]["FBP"] == "m2"


def _check(data, out):
    from ecmapper.pwl.check import check

    return check(data, out, MINI_BKMS)


def test_check_mini_map(generated):
    data, out, _ = generated
    report = _check(data, out)
    assert report["parse_errors"] == [] and report["mismatches"] == [] and report["direction_mismatches"] == []
    assert report["reactions_checked"] == 8 and report["direction_checked"] == 7 and report["motif_problems"] == []


@pytest.mark.parametrize(
    "old, new, key",
    [
        ("2.7.1.11 | ATP>ADP", "2.7.1.12 | ATP>ADP", "mismatches"),  # EC number
        (" | ATP>ADP", " | ATP>AMP", "mismatches"),  # cofactors
        ("compound FBP brenda:6523", "compound FBP brenda:6524", "mismatches"),  # ligand ID
        ("]-> (DHAP + GAP)", "]- (DHAP + GAP)", None),  # no direction is not an error ...
        ("FBP -[node12", "(DHAP + GAP) -[node12", "mismatches"),  # ... but the drawn nodes must match
        ("compound X", "compound ring", "parse_errors"),
    ],
)
def test_check_detects_corruptions(generated, old, new, key):
    data, out, _ = generated
    path = out / "subsystems" / "mini.pwl"
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new, 1))
    report = _check(data, out)
    if key is None:
        assert report["mismatches"] == [] and report["parse_errors"] == []
    else:
        assert report[key]


def test_check_detects_wrong_direction(generated):
    data, out, _ = generated
    path = out / "subsystems" / "mini.pwl"
    text = path.read_text()
    # turn the aldolase around: (DHAP + GAP) -> FBP contradicts BR3 with ev=eq
    text = text.replace("FBP -[node12: 4.1.2.13", "FBP <-[node12: 4.1.2.13")
    text = text.replace("]-> (DHAP + GAP)", "]- (DHAP + GAP)")
    path.write_text(text)
    assert any("node12" in m for m in _check(data, out)["direction_mismatches"])


def test_check_detects_missing_reaction(generated):
    data, out, _ = generated
    path = out / "subsystems" / "mini.pwl"
    lines = [line for line in path.read_text().split("\n") if "node14" not in line]
    path.write_text("\n".join(lines))
    report = _check(data, out)
    assert any("missing" in m for m in report["mismatches"])


def test_cli(generated, tmp_path, capsys, cache_dir):
    data, out, _ = generated
    assert pwl_main(["parse", str(out / "subsystems" / "mini.pwl")]) == 0
    assert "8 reactions, 7 with direction" in capsys.readouterr().out
    assert pwl_main(["motifs", str(out / "subsystems" / "mini.pwl")]) == 0
    assert "motif cycle { A B C };" in capsys.readouterr().out
    bad = tmp_path / "bad.pwl"
    bad.write_text("subsystem x 'x' {}")
    assert pwl_main(["parse", str(bad)]) == 1
    # generate without BKMS and KEGG: no directions, the files are still complete and valid
    out2 = tmp_path / "plain"
    args = ["generate", str(tmp_path / "map.json"), str(out2), "--no-bkms", "--no-kegg", "--cache-dir", str(cache_dir)]
    assert pwl_main(args) == 0
    assert json.loads((out2 / "validation_report.json").read_text())["mismatches"] == 0
    assert "-[node12: 4.1.2.13]-" not in (out2 / "subsystems" / "mini.pwl").read_text()
    assert "react [node12: 4.1.2.13] {FBP, DHAP, GAP};" in (out2 / "subsystems" / "mini.pwl").read_text()
    shutil.rmtree(out2)
