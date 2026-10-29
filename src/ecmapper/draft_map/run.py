"""Build a draft map of a model from the BRENDA map.

usage: python -m ecmapper.draft_map.run PREP.pkl BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR [options]

PREP.pkl is written by :mod:`ecmapper.draft_map.prepare`. See DRAFT_MAP.md for the steps.
"""

from __future__ import annotations

import argparse
import copy
import json
import pickle
import time
from collections import Counter
from pathlib import Path

import numpy as np

from .core import (
    Completion,
    DraftMap,
    ModelData,
    anchor_metabolites,
    disconnect_hubs,
    drop_empty_subsystems,
    drop_redundant_subsystems,
    hub_table,
    localise_edges,
    prune_unassociated,
)
from .escher import write_escher
from .layout import (
    Ports,
    SequencePair,
    anneal_arrangement,
    apply_arrangement,
    arrangement_original,
    condense_subsystems,
    evaluate_arrangement,
    fit_boxes_original,
    internal_metrics,
    to_absolute,
    to_relative,
)
from .relayout import METHODS as RELAYOUT_METHODS
from .relayout import available_methods, relayout_subsystems
from .render import write_cnapy_maps, write_reaction_positions, write_svg
from .routes import route_edges
from .split import (
    ASIDE,
    DEFAULT_CENTRAL,
    assign_maps,
    cross_map_metabolites,
    overview,
    slug,
    split_map,
    write_overview_svg,
)


def coverage(md, dm) -> dict:
    """How many model reactions the draft map shows."""
    drawn = {r for n in dm.nodes.values() for r in n.reactions}
    added = {r for n in dm.nodes.values() if n.status == "added" for r in n.reactions}
    total = len(md.reactions)
    boundary = sum(r["boundary"] for r in md.reactions.values())
    transport = sum(r["transport"] for r in md.reactions.values())
    other = total - boundary - transport
    drawn_other = sum(1 for r in drawn if not md.reactions[r]["boundary"] and not md.reactions[r]["transport"])
    return dict(reactions=total, boundary=boundary, transporters=transport, other=other, drawn=len(drawn),
                drawn_from_map=len(drawn - added), drawn_added=len(added), drawn_other=drawn_other,
                coverage_other=drawn_other / other)


def arrange(dm, args, final, log):
    """Step 6: the original order (compacted, and as packed sequence pair) and the annealed arrangement.

    Returns the ports, the compacted original order (x, y), {name: (x, y,
    flip_x, flip_y, evaluation)} and the evaluations of the annealing seeds.
    """
    ports = Ports.from_map(dm)
    results = {}
    x, y = arrangement_original(dm, gap=args.gap)
    results["original order (compacted)"] = (x, y, None, None, evaluate_arrangement(ports, x, y))
    log(f"original order: {results['original order (compacted)'][4]}")
    sp_ = SequencePair(ports.w, ports.h, args.gap)
    gp, gn = SequencePair.from_positions(x + ports.w / 2, y + ports.h / 2)
    xs, ys = sp_.pack(gp, gn)
    results["original order (sequence pair packing)"] = (xs, ys, None, None, evaluate_arrangement(ports, xs, ys))
    log(f"original order, packed: {results['original order (sequence pair packing)'][4]}")
    if len(ports.names) < 2:
        results[final] = results["original order (sequence pair packing)"]
        return ports, x, y, results, []
    runs = []
    for seed in range(args.seed, args.seed + args.seeds):
        res = anneal_arrangement(ports, x, y, gap=args.gap, iterations=args.iterations, seed=seed,
                                 flips=False if args.flips == "none" else args.flips,
                                 log=log if args.seeds == 1 else None)
        runs.append(res)
        log(f"{final} (seed {seed}): {res[4]}")
    res = min(runs, key=lambda r: r[4]["cost"])
    results[final] = res
    log(f"{final}: {res[4]}")
    return ports, x, y, results, [r[4] for r in runs]


def write_outputs(m, md, out: Path, prefix: str, title: str, args, stats: dict) -> None:
    """Files of a final map (absolute positions): SVG, CNApy box positions, Escher map, pickle, reaction positions.

    Main edges that pass through other nodes are routed around them first (:mod:`.routes`).
    """
    stats["routed_edges"] = route_edges(m) if available_methods(["fdp"]) else 0  # needs pygraphviz
    write_svg(m, str(out / f"{prefix}.svg"), md=md)
    write_cnapy_maps(m, str(out / f"{prefix}.maps"))
    stats["escher"] = write_escher(m, md, str(out / f"{prefix}.escher.json"), title, scale=args.escher_scale,
                                   subsystem_labels=not args.escher_no_labels,
                                   description="Draft map from the BRENDA map, made with ecmapper.draft_map")
    with open(out / f"{prefix}.pkl", "wb") as f:
        pickle.dump(m, f)
    positions = "reaction_positions.json" if prefix == "draft_map" else f"{prefix}.reaction_positions.json"
    write_reaction_positions(m, str(out / positions))
    stats["final"] = m.summary()
    stats["coverage"] = coverage(md, m)


def split_maps(dm, md, args, final, log, stats) -> None:
    """--split: a central map and auxiliary maps, each arranged on its own, and an overview of them."""
    explicit = None
    if args.maps:
        with open(args.maps, encoding="utf-8") as f:
            explicit = json.load(f)
    central = DEFAULT_CENTRAL if args.central is None else [s for s in args.central.split(",") if s]
    n_log = len(dm.log)
    order, assignment = assign_maps(dm, central, explicit, args.min_map_reactions)
    for line in dm.log[n_log:]:
        log(line)
    cross = cross_map_metabolites(dm, assignment)
    maps = split_map(dm, order, assignment)
    log(dm.log[-1])
    hubs = hub_table(dm)
    if args.hub_subsystems:
        # the links between maps (of hubs, too) are references now; inside the maps, hub links are removed
        stats["hubs_disconnected"] = [r for r in hubs if r["subsystems"] >= args.hub_subsystems]
    stats["log"] = dm.log
    stats["maps"] = {}
    arranged = []
    index = dict(maps=[], cross_map_metabolites=cross)
    for name, sub in maps:
        prefix = slug(name)
        log(f"map {name!r}: {len(sub.subsystems)} subsystems, {len(sub.links)} inter-pathway links")
        if args.pack:
            sub.links = []
        elif args.hub_subsystems:
            disconnect_hubs(sub, args.hub_subsystems, hubs)
            log(sub.log[-1])
        ports, _, _, results, seeds = arrange(sub, args, final, log)
        ax, ay, fx, fy, _ = results[final]
        m = copy.deepcopy(sub)
        apply_arrangement(m, ports.names, ax, ay, fx, fy)
        to_absolute(m)
        s = dict(arrangements={k: v[4] for k, v in results.items()}, anneal_seeds=seeds)
        write_outputs(m, md, args.out, prefix, f"{md.model_id} {name} (BRENDA map draft)", args, s)
        del s["coverage"]
        s["references"] = sum(1 for n in m.nodes.values() if getattr(n, "stubs", None))
        stats["maps"][name] = s
        arranged.append((name, m))
        index["maps"].append(dict(name=name, file=f"{prefix}.escher.json", subsystems=list(m.subsystems),
                                  reactions=sorted({r for n in m.nodes.values() for r in n.reactions})))
    layout = overview(arranged, dm, assignment, gap=args.gap, seed=args.seed,
                      hubs={r["key"] for r in stats.get("hubs_disconnected", [])})
    write_overview_svg(layout, arranged, str(args.out / "overview.svg"),
                       {name: f"{slug(name)}.svg" for name, _ in arranged})
    with open(args.out / "maps.json", "w", encoding="utf-8") as f:
        json.dump(index, f, indent=1)
    on_maps = Counter(r for entry in index["maps"] for r in entry["reactions"])
    stats["split"] = dict(maps=len(arranged), reactions=len(on_maps),
                          reactions_on_several_maps=sum(1 for c in on_maps.values() if c > 1),
                          cross_map_metabolites=len(cross))
    stats["coverage"] = coverage(md, dm)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("prep")
    parser.add_argument("brenda_json")
    parser.add_argument("brenda_svg")
    parser.add_argument("out", type=Path)
    parser.add_argument("--keep-subsets", action="store_true",
                        help="keep subsystems whose model reactions are a proper subset of another subsystem's (rule A)")
    parser.add_argument("--aside-coverage", type=float, default=0.9,
                        help="set aside subsystems of which at least this fraction of the model reactions are in other "
                        "subsystems (rules B and C; 1: only fully covered ones, above 1: off)")
    parser.add_argument("--aside-min-matched", type=float, default=0.0,
                        help="set aside subsystems in which less than this fraction of the enzyme nodes have model "
                        "reactions (rule D)")
    parser.add_argument("--min-anchored", type=int, default=2)
    parser.add_argument("--no-completion", action="store_true")
    parser.add_argument("--relayout", default=",".join(m for m in RELAYOUT_METHODS if m != "current"),
                        help="step 3b: layouts tried for subsystems with added reactions, comma-separated from "
                        "spring, fdp, dot, neato (default: all); 'none': keep the positions of step 3")
    parser.add_argument("--circle", default="", metavar="IDS",
                        help="comma-separated BRENDA subsystem IDs whose cycle of main metabolites is drawn on a "
                        "circle in step 3b (e.g. citric-acid-cycle)")
    parser.add_argument("--relayout-margin", type=float, default=0.25,
                        help="step 3b: take dot or neato (which redraw the BRENDA part) only with at least this "
                        "fraction fewer crossings")
    parser.add_argument("--iterations", type=int, default=20000)
    parser.add_argument("--gap", type=float, default=80.0)
    parser.add_argument("--compactness", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--seeds", type=int, default=1, help="anneal with this many seeds, keep the best")
    parser.add_argument("--hub-subsystems", type=int, default=0,
                        help="remove the inter-pathway links of metabolites drawn in at least this many subsystems")
    parser.add_argument("--flips", choices=["none", "x", "xy"], default="none",
                        help="allow mirroring subsystems in the final arrangement (x: left-right, xy: also upside down)")
    parser.add_argument("--escher-scale", type=float, default=2.0, help="scale of the Escher map")
    parser.add_argument("--escher-no-labels", action="store_true", help="no subsystem names in the Escher map")
    parser.add_argument("--compare-flips", action="store_true", help="also anneal with mirroring (x and xy)")
    parser.add_argument("--extra", nargs="*", default=[], help="further arrangements: graphviz elk user-order")
    parser.add_argument("--pack", action="store_true",
                        help="after step 4, drop the inter-pathway links and pack the subsystems into as little space "
                        "as possible (instead of steps 5 and 6)")
    parser.add_argument("--split", action="store_true",
                        help="split into a central map and auxiliary maps (BRENDA colour groups), each arranged on "
                        "its own, with an overview")
    parser.add_argument("--central", help="with --split: comma-separated BRENDA subsystem IDs of the central map "
                        "(default: glycolysis, gluconeogenesis, pentose phosphate and Entner-Doudoroff pathways, "
                        "pyruvate dehydrogenase, citric acid cycle, oxidative phosphorylation, fermentations)")
    parser.add_argument("--maps", help="with --split: JSON file {map name: [BRENDA subsystem IDs]} of further maps")
    parser.add_argument("--min-map-reactions", type=int, default=25,
                        help="with --split: merge colour group maps with fewer model reactions into another map")
    args = parser.parse_args(argv)
    methods = [] if args.relayout == "none" else [m for m in args.relayout.split(",") if m]
    if set(methods) - set(RELAYOUT_METHODS) - {"current"}:
        parser.error(f"--relayout: unknown layouts {sorted(set(methods) - set(RELAYOUT_METHODS))}")
    if not args.split and (args.central is not None or args.maps):
        parser.error("--central and --maps need --split")
    if args.split and (args.compare_flips or args.extra):
        parser.error("--split: --compare-flips and --extra do not apply")
    if args.pack and (args.hub_subsystems or args.flips != "none" or args.compare_flips or args.extra):
        parser.error("--pack arranges the subsystems without links: --hub-subsystems, --flips, --compare-flips "
                     "and --extra do not apply")
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    def log(msg):
        print(f"[{time.time() - t0:6.1f}s] {msg}", flush=True)

    md = ModelData(args.prep)
    dm = DraftMap.from_brenda(args.brenda_json, args.brenda_svg)
    stats = {"original": dm.summary()}
    origin_positions = {n.id: (n.x, n.y, n.subsystem) for n in dm.nodes.values()}

    # step 1: subsystems with model reactions
    drop_empty_subsystems(dm, md)
    log(dm.log[-1])
    stats["redundant_subsystems"] = drop_redundant_subsystems(dm, not args.keep_subsets, args.aside_coverage,
                                                              args.aside_min_matched)
    log(dm.log[-1])
    anchor_metabolites(dm, md)
    stats["step1_subsystems"] = dm.summary()
    write_svg(dm, str(args.out / "step1_subsystems.svg"), md=md)

    # step 2: enzyme nodes without model reactions (the arrangement follows after condensation, see DRAFT_MAP.md)
    prune_unassociated(dm)
    stats["step2_pruned"] = dm.summary()
    write_svg(dm, str(args.out / "step2_pruned.svg"), md=md)
    log(dm.log[-1])

    # step 3: model reactions that connect drawn metabolites
    added = []
    if not args.no_completion:
        completion = Completion(dm, md, min_anchored=args.min_anchored, origin_positions=origin_positions)
        added = completion.run()
        log(dm.log[-1])
    stats["step3_completed"] = dm.summary()
    write_svg(dm, str(args.out / "step3_completed.svg"), md=md)
    with open(args.out / "added_reactions.json", "w") as f:
        json.dump(added, f, indent=1)

    aside = {k for k, s in dm.subsystems.items() if s.aside}
    if args.pack or args.split:
        # no edges between subsystems, which would cross the packed boxes
        localise_edges(dm)
        log(dm.log[-1])
    elif aside:
        # no edges to set-aside subsystems, which are drawn on a separate map
        localise_edges(dm, aside)
        log(dm.log[-1])

    # step 3b: lay out the subsystems with added reactions again
    if methods:
        missing = set(methods) - set(available_methods(methods))
        if missing:
            log(f"step 3b: layouts {sorted(missing)} skipped (networkx or pygraphviz not installed)")
        stats["relayout"] = relayout_subsystems(dm, methods, margin=args.relayout_margin,
                                                circle=[k for k in args.circle.split(",") if k],
                                                compactness=args.compactness)
        log(dm.log[-1])
        write_svg(dm, str(args.out / "step3b_relayout.svg"), md=md)

    # original subsystem centres (for the arrangement that keeps the original order)
    orig_centres = {k: (s.x + s.w / 2, s.y + s.h / 2) for k, s in dm.subsystems.items()}

    # step 4: condense each subsystem
    before_condense = copy.deepcopy(dm)
    stats["internal_before_condensing"] = internal_metrics(dm)
    condense_subsystems(dm, compactness=args.compactness)
    dm.orig_centres = orig_centres
    stats["internal_after_condensing"] = internal_metrics(dm)
    log(f"step 4: subsystems condensed {stats['internal_before_condensing']} -> {stats['internal_after_condensing']}")

    # step 5: hub metabolites
    stats["hubs_all"] = hub_table(dm)[:40]
    if args.hub_subsystems and not args.split:
        stats["hubs_disconnected"] = disconnect_hubs(dm, args.hub_subsystems)
        log(dm.log[-1])
    final = "packed" if args.pack else "annealed"
    if args.split:
        # the links between the maps become references; with --pack, the links inside the maps are dropped there
        split_maps(dm, md, args, final, log, stats)
        with open(args.out / "stats.json", "w") as f:
            json.dump(stats, f, indent=1, default=float)
        log("done")
        return
    aside_map = None
    if aside:
        # the set-aside subsystems on a separate map; links to them become references
        main_name = f"{md.model_id} draft map"
        whole = dm
        (_, dm), (_, aside_map) = split_map(whole, [main_name, ASIDE],
                                            {k: ASIDE if k in aside else main_name for k in whole.subsystems})
        dm.log.append(whole.log[-1])
        log(dm.log[-1])
    if args.pack:
        # without links, the annealing below only minimises the area (and the deviation from the aspect ratio)
        dm.log.append(f"pack: removed all {len(dm.links)} inter-pathway links")
        dm.links = []
        log(dm.log[-1])

    # step 6: arrangement
    ports, x, y, results, stats["anneal_seeds"] = arrange(dm, args, final, log)
    sp_ = SequencePair(ports.w, ports.h, args.gap)
    if args.compare_flips:
        for mode in ("x", "xy"):
            res = anneal_arrangement(ports, x, y, gap=args.gap, iterations=args.iterations, seed=args.seed, flips=mode)
            results[f"annealed with flips {mode}"] = res
            log(f"annealed with flips {mode}: {res[4]}")
    if "user-order" in args.extra:
        # the order of the user's scheme: arrange the uncondensed (pruned and completed) subsystems, then condense
        um = copy.deepcopy(before_condense)
        if args.hub_subsystems:
            disconnect_hubs(um, args.hub_subsystems)
        fit_boxes_original(um)
        to_relative(um)
        uports = Ports.from_map(um)
        ux = np.array([um.subsystems[k].x for k in uports.names])
        uy = np.array([um.subsystems[k].y for k in uports.names])
        ures = anneal_arrangement(uports, ux, uy, gap=args.gap, iterations=args.iterations, seed=args.seed, log=log)
        log(f"user order, before condensing: {ures[4]}")
        # same sequence pair, condensed boxes
        gpu, gnu = SequencePair.from_positions(ures[0] + uports.w / 2, ures[1] + uports.h / 2)
        xs, ys = sp_.pack(gpu, gnu)
        results["user order (arrange, then condense)"] = (xs, ys, None, None, evaluate_arrangement(ports, xs, ys))
        log(f"user order, condensed: {results['user order (arrange, then condense)'][4]}")
    for name in args.extra:
        if name == "graphviz":
            from .external import graphviz_arrangement

            for prog in ("neato", "fdp"):
                gx, gy = graphviz_arrangement(dm, ports, prog=prog, gap=args.gap)
                results[f"graphviz {prog}"] = (gx, gy, None, None, evaluate_arrangement(ports, gx, gy))
                log(f"graphviz {prog}: {results[f'graphviz {prog}'][4]}")
        elif name == "elk":
            from .external import elk_arrangement

            for algorithm in ("layered", "stress"):
                ex, ey = elk_arrangement(dm, ports, algorithm=algorithm, gap=args.gap)
                results[f"elk {algorithm}"] = (ex, ey, None, None, evaluate_arrangement(ports, ex, ey))
                log(f"elk {algorithm}: {results[f'elk {algorithm}'][4]}")

    stats["arrangements"] = {k: v[4] for k, v in results.items()}
    stats["log"] = dm.log
    for name, (ax, ay, fx, fy, _) in results.items():
        m = copy.deepcopy(dm)
        apply_arrangement(m, ports.names, ax, ay, fx, fy)
        to_absolute(m)
        fname = name.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")
        write_svg(m, str(args.out / f"map_{fname}.svg"), md=md)
        if name == final:
            write_outputs(m, md, args.out, "draft_map", f"{md.model_id} draft map (BRENDA map)", args, stats)
    if aside_map is not None:
        log(f"set-aside map: {len(aside_map.subsystems)} subsystems")
        if args.pack:
            aside_map.links = []
        ports, _, _, results, _ = arrange(aside_map, args, final, log)
        m = copy.deepcopy(aside_map)
        apply_arrangement(m, ports.names, *results[final][:4])
        to_absolute(m)
        s = dict(arrangement=results[final][4])
        write_outputs(m, md, args.out, "set_aside", f"{md.model_id} set-aside subsystems (BRENDA map draft)", args, s)
        del s["coverage"]
        stats["set_aside_map"] = s
    with open(args.out / "stats.json", "w") as f:
        json.dump(stats, f, indent=1, default=float)
    log("done")


if __name__ == "__main__":
    main()
