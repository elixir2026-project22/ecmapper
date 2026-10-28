import json

import cobra
import pytest
from cobra import Metabolite, Model, Reaction

from ecmapper import ECMapper, check_updates, update_sources
from ecmapper.cli import main
from ecmapper.sources import KEGGSource, ModelSEEDSource

KEGG_URL = KEGGSource.files["kegg_link_ec_rn.tsv"]
SEED_FILE = "Unique_ModelSEED_Reaction_ECs.txt"
SEED_URL = ModelSEEDSource.files[SEED_FILE]
SEED_V1 = b"ModelSEED ID\tExternal ID\tSource\nrxn00001\t3.6.1.1\tEnzyme Class\n"
SEED_V2 = SEED_V1 + b"rxn00002\t3.5.1.54\tEnzyme Class\n"


def statuses(results):
    return {s.file: s.status for s in results}


def test_etag(tmp_path, server):
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"', "Last-Modified": "Mon, 01 Sep 2025 10:00:00 GMT"})
    seed = ModelSEEDSource(cache_dir=tmp_path)
    assert statuses(seed.check_updates()) == {SEED_FILE: "missing"}
    assert seed.lookup("rxn00001") == {"3.6.1.1"}
    meta = json.loads((tmp_path / f"{SEED_FILE}.meta.json").read_text())
    assert meta["etag"] == '"v1"'

    server.requests.clear()
    (status,) = seed.check_updates()
    assert status.status == "current"
    assert status.cached == meta["downloaded"]
    assert server.requests[0][1]["If-None-Match"] == '"v1"'

    # a new version on the server is detected without downloading it
    server.files[SEED_URL] = (SEED_V2, {"ETag": '"v2"', "Last-Modified": "Mon, 01 Dec 2025 10:00:00 GMT"})
    (status,) = seed.check_updates()
    assert status.status == "outdated" and status.update_available
    assert status.remote_modified == "Mon, 01 Dec 2025 10:00:00 GMT"
    assert (tmp_path / SEED_FILE).read_bytes() == SEED_V1
    assert seed.lookup("rxn00002") == set()

    # update downloads it and the source is reloaded
    assert statuses(seed.update()) == {SEED_FILE: "updated"}
    assert (tmp_path / SEED_FILE).read_bytes() == SEED_V2
    assert seed.lookup("rxn00002") == {"3.5.1.54"}
    assert statuses(seed.check_updates()) == {SEED_FILE: "current"}
    assert statuses(seed.update()) == {SEED_FILE: "current"}


def test_without_validators(tmp_path, server):
    """KEGG sends neither ETag nor Last-Modified, the content is compared."""
    server.files[KEGG_URL] = (b"rn:R00200\tec:2.7.1.40\n", {})
    kegg = KEGGSource(cache_dir=tmp_path)
    assert kegg.lookup("R00200") == {"2.7.1.40"}
    assert statuses(kegg.check_updates()) == {"kegg_link_ec_rn.tsv": "current"}

    server.files[KEGG_URL] = (b"rn:R00200\tec:2.7.1.40\nrn:R00004\tec:3.6.1.1\n", {})
    assert statuses(kegg.check_updates()) == {"kegg_link_ec_rn.tsv": "outdated"}
    assert kegg.lookup("R00004") == set()
    assert statuses(kegg.update()) == {"kegg_link_ec_rn.tsv": "updated"}
    assert kegg.lookup("R00004") == {"3.6.1.1"}
    # no temporary files are left behind
    assert sorted(p.name for p in tmp_path.iterdir()) == ["kegg_link_ec_rn.tsv", "kegg_link_ec_rn.tsv.meta.json"]


def test_cache_without_metadata(tmp_path, server):
    """Files cached without metadata are compared by content, then the validators are stored."""
    (tmp_path / SEED_FILE).write_bytes(SEED_V1)
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"'})
    assert statuses(check_updates([ModelSEEDSource(cache_dir=tmp_path)])) == {SEED_FILE: "current"}
    assert json.loads((tmp_path / f"{SEED_FILE}.meta.json").read_text())["etag"] == '"v1"'
    server.requests.clear()
    assert statuses(check_updates([ModelSEEDSource(cache_dir=tmp_path)])) == {SEED_FILE: "current"}
    assert server.requests[0][1]["If-None-Match"] == '"v1"'


def test_modified_cache_file(tmp_path, server):
    """A cached file that was changed locally is compared by content, not by its stale validators."""
    server.files[SEED_URL] = (SEED_V2, {"ETag": '"v2"'})
    seed = ModelSEEDSource(cache_dir=tmp_path)
    seed.update()
    (tmp_path / SEED_FILE).write_bytes(SEED_V1)
    server.requests.clear()
    assert statuses(seed.check_updates()) == {SEED_FILE: "outdated"}
    assert "If-None-Match" not in server.requests[0][1]
    assert statuses(seed.update()) == {SEED_FILE: "updated"}
    assert (tmp_path / SEED_FILE).read_bytes() == SEED_V2


def test_force(tmp_path, server):
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"'})
    seed = ModelSEEDSource(cache_dir=tmp_path)
    seed.update()
    server.requests.clear()
    assert statuses(seed.update(force=True)) == {SEED_FILE: "current"}
    assert "If-None-Match" not in server.requests[0][1]


def test_missing_local_and_error(tmp_path, cache_dir, server):
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"'})
    results = update_sources(["modelseed"], cache_dir=tmp_path)
    assert statuses(results) == {SEED_FILE: "downloaded"}
    assert results[0].cached

    local = KEGGSource(local_files={"kegg_link_ec_rn.tsv": cache_dir / "kegg_link_ec_rn.tsv"})
    assert statuses(local.check_updates()) == {"kegg_link_ec_rn.tsv": "local"}

    (tmp_path / "kegg_link_ec_rn.tsv").write_text("rn:R00200\tec:2.7.1.40\n")
    (status,) = KEGGSource(cache_dir=tmp_path).check_updates()  # the server does not know the URL
    assert status.status == "error"
    assert status.detail


@pytest.fixture
def model_file(tmp_path):
    model = Model("m")
    a, b = Metabolite("a_c", compartment="c"), Metabolite("b_c", compartment="c")
    r = Reaction("R1")
    r.add_metabolites({a: -1, b: 1})
    r.annotation = {"seed.reaction": "rxn00002"}
    model.add_reactions([r])
    path = tmp_path / "model.xml"
    cobra.io.write_sbml_model(model, str(path))
    return path


def test_mapper_update(tmp_path, server, model_file):
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"'})
    mapper = ECMapper(["modelseed"], cache_dir=tmp_path / "cache", current_ec=False, compare_metabolites=False)
    model = cobra.io.read_sbml_model(str(model_file))
    assert mapper.map_model(model)["R1"].origin == "none"
    server.files[SEED_URL] = (SEED_V2, {"ETag": '"v2"'})
    assert statuses(mapper.check_updates()) == {SEED_FILE: "outdated"}
    assert statuses(mapper.update()) == {SEED_FILE: "updated"}
    assert mapper.map_model(model)["R1"].ec_numbers == ["3.5.1.54"]


def test_cli(tmp_path, server, model_file, capsys):
    cache = tmp_path / "cache"
    server.files[SEED_URL] = (SEED_V1, {"ETag": '"v1"', "Last-Modified": "Mon, 01 Sep 2025 10:00:00 GMT"})
    assert main(["--check-updates", "--sources", "modelseed", "--no-current-ec", "--no-metabolites", "--cache-dir", str(cache)]) == 0
    out = capsys.readouterr().out
    assert "missing" in out and "1 file(s) can be updated" in out

    assert main(["--update", "--sources", "modelseed", "--no-current-ec", "--no-metabolites", "--cache-dir", str(cache)]) == 0
    assert "downloaded" in capsys.readouterr().out

    assert main(["--check-updates", "--sources", "modelseed", "--no-current-ec", "--no-metabolites", "--cache-dir", str(cache)]) == 0
    assert "all files are up to date" in capsys.readouterr().out

    # update before mapping
    server.files[SEED_URL] = (SEED_V2, {"ETag": '"v2"'})
    out_file = tmp_path / "out.tsv"
    args = [str(model_file), "--update", "--sources", "modelseed", "--no-current-ec", "--no-metabolites", "--cache-dir", str(cache), "-o", str(out_file)]
    assert main(args) == 0
    assert "updated" in capsys.readouterr().out
    assert "3.5.1.54" in out_file.read_text()

    # a failed check is reported with exit code 1
    del server.files[SEED_URL]
    assert main(["--check-updates", "--sources", "modelseed", "--no-current-ec", "--no-metabolites", "--cache-dir", str(cache)]) == 1
