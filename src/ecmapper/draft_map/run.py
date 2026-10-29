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
from pathlib import Path

import numpy as np

from .core import (
    Completion,
    DraftMap,
    ModelData,
    anchor_metabolites,
    disconnect_hubs,
    drop_empty_subsystems,
    hub_table,
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
from .render import write_cnapy_maps, write_reaction_positions, write_svg


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


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("prep")
    parser.add_argument("brenda_json")
    parser.add_argument("brenda_svg")
    parser.add_argument("out", type=Path)
    parser.add_argument("--min-anchored", type=int, default=2)
    parser.add_argument("--no-completion", action="store_true")
    parser.add_argument("--iterations", type=int, default=20000)
    parser.add_argument("--gap", type=float, default=80.0)
    parser.add_argument("--compactness", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--seeds", type=int, default=1, help="anneal with this many seeds, keep the best")
    parser.add_argument("--hub-subsystems", type=int, default=0,
                        help="remove the inter-pathway links of metabolites drawn in at least this many subsystems")
    parser.add_argument("--flips", choices=["none", "x", "xy"], default="none",
                        help="allow mirroring subsystems in the final arrangement (x: left-right, xy: also upside down)")
    parser.add_argument("--escher-scale", type=float, default=2.5, help="scale of the Escher map")
    parser.add_argument("--escher-no-labels", action="store_true", help="no subsystem names in the Escher map")
    parser.add_argument("--compare-flips", action="store_true", help="also anneal with mirroring (x and xy)")
    parser.add_argument("--extra", nargs="*", default=[], help="further arrangements: graphviz elk user-order")
    args = parser.parse_args(argv)
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
    anchor_metabolites(dm, md)
    stats["step1_subsystems"] = dm.summary()
    write_svg(dm, str(args.out / "step1_subsystems.svg"), md=md)
    log(dm.log[-1])

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
    if args.hub_subsystems:
        stats["hubs_disconnected"] = disconnect_hubs(dm, args.hub_subsystems)
        log(dm.log[-1])

    # step 6: arrangement
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
    runs = []
    for seed in range(args.seed, args.seed + args.seeds):
        res = anneal_arrangement(ports, x, y, gap=args.gap, iterations=args.iterations, seed=seed,
                                 flips=False if args.flips == "none" else args.flips,
                                 log=log if args.seeds == 1 else None)
        runs.append(res)
        log(f"annealed (seed {seed}): {res[4]}")
    stats["anneal_seeds"] = [r[4] for r in runs]
    res = min(runs, key=lambda r: r[4]["cost"])
    results["annealed"] = res
    log(f"annealed: {res[4]}")
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
        slug = name.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")
        write_svg(m, str(args.out / f"map_{slug}.svg"), md=md)
        if name == "annealed":
            write_svg(m, str(args.out / "draft_map.svg"), md=md)
            write_cnapy_maps(m, str(args.out / "draft_map.maps"))
            stats["escher"] = write_escher(m, md, str(args.out / "draft_map.escher.json"),
                                           f"{md.model_id} draft map (BRENDA map)", scale=args.escher_scale,
                                           subsystem_labels=not args.escher_no_labels,
                                           description="Draft map from the BRENDA map, made with ecmapper.draft_map")
            with open(args.out / "draft_map.pkl", "wb") as f:
                pickle.dump(m, f)
            write_reaction_positions(m, str(args.out / "reaction_positions.json"))
            stats["final"] = m.summary()
            stats["coverage"] = coverage(md, m)
    with open(args.out / "stats.json", "w") as f:
        json.dump(stats, f, indent=1, default=float)
    log("done")


if __name__ == "__main__":
    main()
