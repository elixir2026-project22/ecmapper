"""Match a model to the BRENDA map with ecmapper and store what the draft map needs.

usage: python -m ecmapper.draft_map.prepare MODEL.xml BRENDA_map_reactions.json OUT.pkl [--vicinity propagate]
       [--subsystems MODEL.json]
"""

from __future__ import annotations

import argparse
import json
import pickle

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


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("brenda_map")
    parser.add_argument("out")
    parser.add_argument("--vicinity", default="propagate", choices=["none", "neighbours", "propagate"])
    parser.add_argument("--subsystems", help="model in JSON format with the subsystems of the reactions "
                        "(BiGG's SBML files have none)")
    args = parser.parse_args(argv)

    model = cobra.io.read_sbml_model(args.model)
    brenda_map = BrendaMap.from_json(args.brenda_map)
    matcher = MapMatcher(brenda_map, ECMapper(), vicinity=None if args.vicinity == "none" else args.vicinity)
    subsystems = None
    if args.subsystems:
        with open(args.subsystems, encoding="utf-8") as f:
            subsystems = {r["id"]: r["subsystem"] for r in json.load(f)["reactions"] if r.get("subsystem")}
    result = prepare(model, matcher, subsystems)
    with open(args.out, "wb") as f:
        pickle.dump(result, f)
    on_map = sum(m["status"] == "on map" for m in result["matches"].values())
    print(f"{len(result['matches'])} reactions, {on_map} on map; written to {args.out}")


if __name__ == "__main__":
    main()
