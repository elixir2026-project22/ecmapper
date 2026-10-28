import pytest
from cobra import Metabolite, Reaction

from ecmapper.metabolites import MetaboliteComparison
from ecmapper.sources import BiGGrMetabolites, EnzymeNomenclature, KEGGReactionCompounds, KEGGSource, check_updates


@pytest.fixture
def enzyme(cache_dir):
    return EnzymeNomenclature(cache_dir=cache_dir)


@pytest.fixture
def comparison(cache_dir, enzyme):
    return MetaboliteComparison(
        KEGGSource(cache_dir=cache_dir),
        KEGGReactionCompounds(cache_dir=cache_dir),
        BiGGrMetabolites(cache_dir=cache_dir),
        enzyme,
    )


def test_enzyme(enzyme):
    assert enzyme.current("4.2.1.52") == {"4.3.3.7"}
    # transfer to several EC numbers, written over two lines
    assert enzyme.current("1.7.99.4") == {"1.7.1.1", "1.7.1.2", "1.7.1.3", "1.7.5.1", "1.7.7.2", "1.9.6.1"}
    assert enzyme.current("1.1.1.74") == set()
    assert enzyme.current("4.3.3.7") == {"4.3.3.7"}
    # not in ENZYME: unchanged
    assert enzyme.current("9.9.9.9") == {"9.9.9.9"}
    assert enzyme.current("2.7.1.-") == {"2.7.1.-"}
    assert enzyme.entry("5.3.1.9").name == "glucose-6-phosphate isomerase"
    assert enzyme.entry("1.1.1.74").deleted
    assert "ECOLI" in enzyme.organisms("2.7.1.40")
    assert "ECOLI" not in enzyme.organisms("3.5.1.14")
    assert enzyme.organisms("9.9.9.9") == frozenset()


def test_kegg_compounds(cache_dir):
    compounds = KEGGReactionCompounds(cache_dir=cache_dir)
    assert compounds.compounds("R00200") == {"C00002", "C00008", "C00022", "C00074"}
    assert compounds.compounds("R99999") == frozenset()


def test_biggr_metabolites(cache_dir):
    metabolites = BiGGrMetabolites(cache_dir=cache_dir)
    assert "C00074" in metabolites.xrefs("pep")["kegg.compound"]
    assert metabolites.xrefs("pep")["chebi"]
    # legacy IDs with single underscores
    assert metabolites.xrefs("glyc_R") == metabolites.xrefs("glyc__R") != {}
    assert metabolites.xrefs("unknown") == {}


def pyk():
    """Pyruvate kinase with BiGG metabolites, one with a KEGG annotation."""
    mets = {m: Metabolite(f"{m}_c", compartment="c") for m in ("adp", "h", "pep", "atp", "pyr")}
    mets["pyr"].annotation = {"kegg.compound": "C00022"}
    r = Reaction("PYK")
    r.add_metabolites({mets["adp"]: -1, mets["h"]: -1, mets["pep"]: -1, mets["atp"]: 1, mets["pyr"]: 1})
    return r


def test_metabolite_compounds(comparison):
    m = Metabolite("pyr_c", compartment="c")
    assert comparison.metabolite_compounds(m) == {"C00022"}
    m = Metabolite("x_c", compartment="c")
    m.annotation = {"kegg.compound": ["C00022", "kegg.compound:C00074"], "bigg.metabolite": "atp"}
    assert comparison.metabolite_compounds(m) == {"C00022", "C00074", "C00002"}
    # protons are ignored
    assert comparison.metabolite_compounds(Metabolite("h_c", compartment="c")) == set()
    # compartment IDs that differ from the suffix of the metabolite ID
    assert comparison.metabolite_compounds(Metabolite("pyr_c", compartment="Cytosol")) == {"C00022"}


def test_scores(comparison):
    scores = comparison.scores(pyk(), ["2.7.1.40", "2.7.9.2", "9.9.9.9"])
    assert scores["2.7.1.40"] == 1.0
    # R00199 (ATP + pyruvate + H2O = AMP + PEP + phosphate): ADP missing, H2O, AMP and phosphate extra
    assert scores["2.7.9.2"] == pytest.approx(2.5 / 4.5)
    assert scores["9.9.9.9"] is None
    assert comparison.reference_reactions("2.7.9.2") == {"R00199"}
    # no metabolite with KEGG compound
    r = Reaction("R")
    r.add_metabolites({Metabolite("x_c", compartment="c"): -1})
    assert comparison.scores(r, ["2.7.1.40"]) == {}


def test_check_updates_includes_references(cache_dir, server):
    statuses = check_updates(["enzyme", "kegg-compounds", "biggr-metabolites"], cache_dir=cache_dir)
    # the fixture files have no metadata and the fake server does not know the URLs
    assert [s.file for s in statuses] == ["enzyme.dat", "kegg_link_cpd_rn.tsv", "biggr_metabolites.json"]
    with pytest.raises(ValueError, match="unknown resource"):
        check_updates(["chebi"])
