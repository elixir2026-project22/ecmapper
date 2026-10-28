import cobra
import pandas as pd
from cobra import Metabolite, Model, Reaction

from ecmapper.cli import main


def test_cli(tmp_path, cache_dir, capsys):
    model = Model("cli_test")
    a, b = Metabolite("a_c", compartment="c"), Metabolite("b_c", compartment="c")
    for rid, annotation in [("PGI", {"bigg.reaction": "PGI"}), ("PYK", {"kegg.reaction": "R00200"})]:
        r = Reaction(rid)
        r.add_metabolites({a: -1, b: 1})
        r.annotation = annotation
        model.add_reactions([r])
    sbml = tmp_path / "model.xml"
    cobra.io.write_sbml_model(model, str(sbml))

    out = tmp_path / "mapping.tsv"
    annotated = tmp_path / "annotated.xml"
    assert main([str(sbml), "-o", str(out), "--cache-dir", str(cache_dir), "--write-model", str(annotated)]) == 0
    printed = capsys.readouterr().out
    assert "resources   2" in printed

    table = pd.read_csv(out, sep="\t", keep_default_na=False).set_index("reaction_id")
    assert table.loc["PGI", "ec_numbers"] == "5.3.1.9"
    assert table.loc["PYK", "ec_numbers"] == "2.7.1.40"
    assert cobra.io.read_sbml_model(str(annotated)).reactions.PGI.annotation["ec-code"] == "5.3.1.9"

    assert main([str(sbml), "-o", str(out), "--long", "--sources", "kegg", "--cache-dir", str(cache_dir)]) == 0
    table = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert list(table.columns) == ["reaction_id", "origin", "ec", "complete", "source", "identifier", "via"]
    assert set(table.source) == {"", "kegg"}
