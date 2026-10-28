"""Evaluate the matching of model reactions to a BRENDA map on models from BiGG.

usage: python evaluate_map_models.py BRENDA_map_reactions.json OUT_DIR [--cache-dir DIR] [MODEL ...]

The SBML files of the models are downloaded from BiGG (bigg.ucsd.edu, which
has more annotations than bigg.bio) into OUT_DIR/models.
For each model, OUT_DIR gets the matches (MODEL_matches.tsv) and the matches
found in the vicinity of others (MODEL_vicinity.tsv, both modes), and the
summary tables coverage.tsv and holdout.tsv:

coverage
    reactions per status of the regular matching and the reactions and nodes
    only found by each vicinity mode (full annotation);
holdout
    the annotation and notes of 30 % of the reactions are removed (3
    repetitions); for the removed reactions that are on the map with full
    annotation: how many are found again on one of their nodes (`recovered`)
    by the regular matching and by both vicinity modes, how many nodes the
    vicinity search accepts that are not among their nodes (`other_node`), and
    how many vicinity matches the removed reactions get that were not on the
    map with full annotation (`new`).
"""

import argparse
import dataclasses
import gzip
import random
from pathlib import Path

import cobra
import pandas as pd
import requests

from ecmapper import ECMapper
from ecmapper.brenda_map import VICINITY_MODES, BrendaMap, MapMatcher, map_matches_to_dataframe, map_summary

#: manually curated prokaryote models with cytosol, periplasm and extracellular space only (see BRENDA_MAP_EVALUATION.md)
MODELS = (
    "iYO844", "iYS854", "iCN900", "iJN1463", "iCN718", "STM_v1_0", "iYL1228", "iAF987", "iIT341", "iEK1008", "iAF692",
)
URL = "https://bigg.ucsd.edu/static/models/{}.xml.gz"
HOLD_OUT, REPETITIONS = 0.3, 3


def download(model_id: str, directory: Path) -> Path:
    path = directory / f"{model_id}.xml"
    if not path.exists():
        response = requests.get(URL.format(model_id), timeout=120)
        response.raise_for_status()
        directory.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.decompress(response.content))
    return path


def copy(matches):
    return {
        rid: dataclasses.replace(m, candidates=[dataclasses.replace(c) for c in m.candidates])
        for rid, m in matches.items()
    }


def on_map(matches):
    return {rid: {c.node.node_id for c in m.nodes} for rid, m in matches.items() if m.nodes}


def evaluate(matcher: MapMatcher, model, model_id: str, out_dir: Path) -> tuple[dict, list[dict]]:
    brenda_map = matcher.map
    full = matcher.match_model(model)
    stats = map_summary(brenda_map, full)
    coverage = {
        "model": model_id,
        "reactions": stats["reactions"] - stats["status"]["boundary"],
        **{s.replace(" ", "_"): stats["status"][s] for s in ("on map", "related", "not on map")},
        "transporters": sum(matcher.is_transport(r) for r in model.reactions if full[r.id].status == "not on map"),
        "nodes": stats["nodes_with_model_reactions"],
        "subsystems": stats["subsystems_with_model_reactions"],
    }
    finds = []
    for mode in VICINITY_MODES:
        matches = copy(full)
        new = matcher.search_vicinity(model.reactions, matches, mode)
        s = map_summary(brenda_map, matches)
        coverage.update({f"{mode}_pairs": len(new), f"{mode}_reactions": s["vicinity_reactions"],
                         f"{mode}_nodes": s["vicinity_nodes"]})
        table = map_matches_to_dataframe(matches)
        finds.append(table[table.via.str.contains("vicinity")].assign(mode=mode))
    pd.concat(finds).to_csv(out_dir / f"{model_id}_vicinity.tsv", sep="\t", index=False)
    map_matches_to_dataframe(full).to_csv(out_dir / f"{model_id}_matches.tsv", sep="\t", index=False)

    reference = on_map(full)
    candidates = [r.id for r in model.reactions if full[r.id].status != "boundary"]
    holdout = []
    for repetition in range(REPETITIONS):
        held = set(random.Random(repetition).sample(candidates, round(HOLD_OUT * len(candidates))))
        degraded = model.copy()
        for rid in held:
            degraded.reactions.get_by_id(rid).annotation = {}
            degraded.reactions.get_by_id(rid).notes = {}
        base = matcher.match_model(degraded)
        target = {rid for rid in held if rid in reference}
        row = {"model": model_id, "repetition": repetition, "held_out": len(held), "on_map_reference": len(target)}
        for label in ("regular", *VICINITY_MODES):
            matches = copy(base)
            if label != "regular":
                matcher.search_vicinity(degraded.reactions, matches, label)
            found = on_map(matches)
            row[f"{label}_recovered"] = sum(bool(found.get(rid, set()) & reference[rid]) for rid in target)
            row[f"{label}_other_node"] = sum(len(found.get(rid, set()) - reference[rid]) for rid in target)
            row[f"{label}_new"] = sum(1 for rid in held - target for c in matches[rid].nodes if "vicinity" in c.via)
        holdout.append(row)
    return coverage, holdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("brenda_map", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("models", nargs="*", default=MODELS)
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    matcher = MapMatcher(BrendaMap.from_json(args.brenda_map), ECMapper(cache_dir=args.cache_dir),
                         cache_dir=args.cache_dir)
    coverage, holdout = [], []
    for model_id in args.models:
        model = cobra.io.read_sbml_model(str(download(model_id, args.out_dir / "models")))
        c, h = evaluate(matcher, model, model_id, args.out_dir)
        coverage.append(c)
        holdout.extend(h)
        print(f"{model_id}: {c['on_map']} of {c['reactions']} reactions on the map", flush=True)
    pd.DataFrame(coverage).to_csv(args.out_dir / "coverage.tsv", sep="\t", index=False)
    pd.DataFrame(holdout).to_csv(args.out_dir / "holdout.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
