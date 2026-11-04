"""Match a model to the BRENDA map with ecmapper and store what the draft map needs.

usage: python -m ecmapper.draft_map.prepare MODEL.xml BRENDA_map_reactions.json OUT.pkl [--vicinity propagate]
       [--subsystems MODEL.json]

The model and the subsystems file may be compressed (.gz, .bz2, .xz, or a .zip with one model file).
"""

from __future__ import annotations

import argparse
import bz2
import gzip
import json
import lzma
import pickle
import tempfile
import zipfile
from pathlib import Path

import cobra

from ..brenda_map import CURRENCY_KEYS, BrendaMap, MapMatcher
from ..mapper import ECMapper


def prepare(model: cobra.Model, matcher: MapMatcher, subsystems: dict[str, str] | None = None) -> dict:
    """Matches of the model's reactions to the map, with its reactions and metabolites (see :class:`.core.ModelData`).

    ``subsystems``: {reaction ID: subsystem} for models without subsystems and groups.
    """
    matches = matcher.match_model(model)
    met_keys = {m.id: matcher.metabolite_keys(m) for m in model.metabolites}
    node_met_keys = {}  # (enzyme node, metabolite node) -> keys, for all drawn metabolites
    inferred = matcher._index["inferred"]
    for node in matcher.map.nodes.values():
        for m in node.metabolites:
            keys = matcher.ligand_keys(m.names)
            if keys and m.ligand_id in inferred:
                keys = keys | {inferred[m.ligand_id]}
            node_met_keys[(node.node_id, m.node_id)] = keys
    boundary = {r.id for r in model.boundary}
    group_of = {}
    for group in model.groups:
        for member in group.members:
            group_of.setdefault(member.id, group.name)
    group_of.update(subsystems or {})
    reactions = {
        r.id: dict(
            name=r.name,
            metabolites={m.id: c for m, c in r.metabolites.items()},
            subsystem=r.subsystem or group_of.get(r.id, ""),
            reversible=r.reversibility,
            boundary=r.id in boundary,
            transport=(r.id not in boundary) and matcher.is_transport(r),
            ecs=sorted(matcher.reaction_ecs(r)) if r.id not in boundary else [],
            gpr=r.gene_reaction_rule,
            genes=sorted((g.id, g.name or g.id) for g in r.genes),
        )
        for r in model.reactions
    }
    metabolites = {
        m.id: dict(name=m.name, compartment=m.compartment, formula=m.formula,
                   currency=bool(met_keys[m.id] & CURRENCY_KEYS), keys=met_keys[m.id])
        for m in model.metabolites
    }
    return dict(
        model_id=model.id,
        reactions=reactions,
        metabolites=metabolites,
        node_met_keys=node_met_keys,
        matches={
            rid: dict(
                status=m.status,
                candidates=[
                    dict(node_id=c.node.node_id, score=c.score, via=c.via, accepted=c.accepted,
                         vicinity_of=c.vicinity_of)
                    for c in m.candidates
                ],
            )
            for rid, m in matches.items()
        },
    )


#: extensions of the model formats (:func:`read_model`)
MODEL_EXTENSIONS = (".xml", ".sbml", ".json", ".yml", ".yaml", ".mat")
#: single compressed files by extension; a .zip holds one model file
DECOMPRESS = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}
COMPRESSED = (*DECOMPRESS, ".zip")


def uncompressed(path) -> tuple[bytes, str] | None:
    """The content and name of a compressed model file (None if `path` is not compressed).

    A .zip archive must hold exactly one file with a model extension (:data:`MODEL_EXTENSIONS`).
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in DECOMPRESS:
        with DECOMPRESS[suffix](path, "rb") as f:
            return f.read(), path.stem
    if suffix != ".zip":
        return None
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if not n.endswith("/") and not n.startswith("__MACOSX/")
                 and not Path(n).name.startswith("._") and Path(n).suffix.lower() in MODEL_EXTENSIONS]
        if len(names) != 1:
            raise ValueError(f"{path}: expected one model file in the archive, found {len(names)} "
                             f"({', '.join(names) or 'extensions ' + ', '.join(MODEL_EXTENSIONS)})")
        return z.read(names[0]), Path(names[0]).name


def _read_model(path: str) -> cobra.Model:
    name = path.lower()
    if name.endswith(".json"):
        return cobra.io.load_json_model(path)
    if name.endswith((".yml", ".yaml")):
        return cobra.io.load_yaml_model(path)
    if name.endswith(".mat"):
        return cobra.io.load_matlab_model(path)
    return cobra.io.read_sbml_model(path)


def read_model(path) -> cobra.Model:
    """A model in SBML (default), JSON, YAML or MATLAB format (by the file extension), possibly compressed
    (:data:`COMPRESSED`, the format by the extension inside, e.g. MODEL.xml.gz)."""
    content = uncompressed(path)
    if content is None:
        return _read_model(str(path))
    data, name = content
    with tempfile.TemporaryDirectory() as folder:
        file = Path(folder) / name
        file.write_bytes(data)
        return _read_model(str(file))


def read_subsystems(path) -> dict[str, str]:
    """{reaction ID: subsystem} of a model in JSON format (for SBML files without subsystems, e.g. BiGG's), possibly
    compressed."""
    content = uncompressed(path)
    data = Path(path).read_bytes() if content is None else content[0]
    return {r["id"]: r["subsystem"] for r in json.loads(data.decode("utf-8"))["reactions"] if r.get("subsystem")}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("brenda_map")
    parser.add_argument("out")
    parser.add_argument("--vicinity", default="propagate", choices=["none", "neighbours", "propagate"])
    parser.add_argument("--subsystems", help="model in JSON format with the subsystems of the reactions "
                        "(BiGG's SBML files have none)")
    args = parser.parse_args(argv)

    model = read_model(args.model)
    brenda_map = BrendaMap.from_json(args.brenda_map)
    matcher = MapMatcher(brenda_map, ECMapper(), vicinity=None if args.vicinity == "none" else args.vicinity)
    subsystems = read_subsystems(args.subsystems) if args.subsystems else None
    result = prepare(model, matcher, subsystems)
    with open(args.out, "wb") as f:
        pickle.dump(result, f)
    on_map = sum(m["status"] == "on map" for m in result["matches"].values())
    print(f"{len(result['matches'])} reactions, {on_map} on map; written to {args.out}")


if __name__ == "__main__":
    main()
