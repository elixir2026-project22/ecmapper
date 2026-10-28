"""Command line interface: ``python -m ecmapper model.xml -o ec_mapping.tsv``."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .brenda_map import BrendaMap, MapMatcher, highlight_svg, map_matches_to_dataframe, map_summary, node_coverage
from .evaluation import evaluate, evaluation_summary
from .mapper import ECMapper, annotate_model, summary, to_dataframe
from .sources import SOURCES, FileStatus, default_cache_dir


def load_model(path: Path):
    """Load a model with cobrapy, the format is chosen by the file extension."""
    import cobra

    suffixes = "".join(path.suffixes[-2:]).lower()
    if suffixes.endswith((".json", ".json.gz")):
        return cobra.io.load_json_model(path)
    if suffixes.endswith((".yml", ".yaml", ".yml.gz", ".yaml.gz")):
        return cobra.io.load_yaml_model(path)
    if suffixes.endswith(".mat"):
        return cobra.io.load_matlab_model(path)
    return cobra.io.read_sbml_model(str(path))


def _print_summary(stats: dict, mapper: ECMapper) -> None:
    print(f"reactions: {stats['reactions']}")
    for origin, n in stats["origin"].items():
        print(f"  {origin:<11} {n}")
    print(f"with EC number: {stats['with_ec']}")
    print("selection of the EC number:")
    for selection, n in stats["selection"].items():
        if selection not in ("none", "boundary"):
            print(f"  {selection:<11} {n}")
    print("reactions with EC numbers per source:")
    for source, n in stats["reactions_per_source"].items():
        print(f"  {source:<11} {n}")
    unavailable = [s.name for s in mapper.sources.values() if s.loaded and not s.available]
    if unavailable:
        print(f"unavailable sources (skipped): {', '.join(unavailable)}")


def _print_file_statuses(statuses: list[FileStatus]) -> None:
    print(f"{'resource':<17} {'file':<34} {'status':<10} {'cached':<26} remote last modified")
    for s in statuses:
        line = f"{s.source:<17} {s.file:<34} {s.status:<10} {s.cached:<26} {s.remote_modified}"
        print(line + (f" ({s.detail})" if s.detail else ""))


def _print_map_summary(stats: dict) -> None:
    print(f"reactions: {stats['reactions']}")
    for status, n in stats["status"].items():
        print(f"  {status:<11} {n}")
    print(
        f"map nodes with model reactions: {stats['nodes_with_model_reactions']} of {stats['nodes']}, "
        f"in {stats['subsystems_with_model_reactions']} subsystems"
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="ecmapper",
        description="Map the reactions of a constraint-based model to EC numbers.",
    )
    parser.add_argument("model", type=Path, nargs="?", help="model file (SBML, JSON, YAML or MAT)")
    parser.add_argument("-o", "--output", type=Path, help="write the mapping to this TSV file")
    parser.add_argument("--long", action="store_true", help="one row per reaction, EC number and evidence")
    parser.add_argument(
        "--sources",
        default=",".join(SOURCES),
        help=f"comma separated resources to query (default: {','.join(SOURCES)})",
    )
    parser.add_argument("--always-query", action="store_true", help="also query reactions that have EC numbers")
    parser.add_argument("--no-notes", action="store_true", help="ignore EC numbers in reaction notes")
    parser.add_argument("--no-xrefs", action="store_true", help="do not follow cross-references between resources")
    parser.add_argument("--include-boundary", action="store_true", help="also map boundary reactions")
    parser.add_argument(
        "--no-bkms-consensus", action="store_true", help="do not select the EC number on which BKMS-react agrees"
    )
    parser.add_argument(
        "--no-current-ec", action="store_true", help="do not replace transferred EC numbers (ExPASy ENZYME)"
    )
    parser.add_argument(
        "--no-metabolites", action="store_true", help="do not select EC numbers by comparing metabolites"
    )
    parser.add_argument("--min-score", type=float, default=0.6, help="minimum score of the comparison (0.6)")
    parser.add_argument("--min-margin", type=float, default=0.1, help="minimum margin to the next EC number (0.1)")
    parser.add_argument(
        "--tie-breakers", default="", help="comma separated tie-breakers: sources, organism (default: none)"
    )
    parser.add_argument("--organism", help="Swiss-Prot organism mnemonic for the organism tie-breaker, e.g. ECOLI")
    parser.add_argument("--cache-dir", type=Path, default=None, help=f"download cache (default: {default_cache_dir()})")
    parser.add_argument("--refresh", action="store_true", help="download the resource files again")
    parser.add_argument(
        "--check-updates", action="store_true", help="check whether newer resource files are available"
    )
    parser.add_argument(
        "--update", action="store_true", help="download newer resource files into the cache (before mapping)"
    )
    parser.add_argument("--write-model", type=Path, help="write the model with the found EC numbers as SBML")
    parser.add_argument(
        "--evaluate",
        action="store_true",
        help="remove the model's EC numbers, map it and compare (-o writes one row per reaction)",
    )
    parser.add_argument(
        "--bigg-only", action="store_true", help="with --evaluate: keep only the BiGG annotations of the model"
    )
    group = parser.add_argument_group("BRENDA map", "find the reactions of the model on a BRENDA pathway map")
    group.add_argument("--brenda-map", type=Path, help="reactions of the map (JSON), -o writes the matches")
    group.add_argument("--node-table", type=Path, help="write one row per map node with its model reactions")
    group.add_argument("--svg", type=Path, help="the map as SVG, with --svg-out")
    group.add_argument("--svg-out", type=Path, help="write the map SVG with the model's reactions highlighted")
    group.add_argument("--all-candidates", action="store_true", help="list all candidate nodes of each reaction")
    group.add_argument("--min-map-score", type=float, default=0.6, help="minimum score of a node match (0.6)")
    group.add_argument(
        "--min-map-score-unlinked",
        type=float,
        default=0.8,
        help="minimum score of nodes found only by metabolites (0.8)",
    )
    group.add_argument(
        "--max-map-gap", type=float, default=0.1, help="maximum distance to the best score of a reaction (0.1)"
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.model is None and not (args.check_updates or args.update):
        parser.error("a model file, --check-updates or --update is required")
    if bool(args.svg) != bool(args.svg_out):
        parser.error("--svg and --svg-out are used together")

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s: %(message)s")
    mapper = ECMapper(
        [s.strip() for s in args.sources.split(",") if s.strip()],
        use_notes=not args.no_notes,
        always_query=args.always_query,
        skip_boundary=not args.include_boundary,
        follow_xrefs=not args.no_xrefs,
        bkms_consensus=not args.no_bkms_consensus,
        current_ec=not args.no_current_ec,
        compare_metabolites=not args.no_metabolites,
        min_score=args.min_score,
        min_margin=args.min_margin,
        tie_breakers=[t.strip() for t in args.tie_breakers.split(",") if t.strip()],
        organism=args.organism,
        cache_dir=args.cache_dir,
        refresh=args.refresh,
    )
    matcher = None
    if args.brenda_map:
        matcher = MapMatcher(
            BrendaMap.from_json(args.brenda_map),
            mapper,
            min_score=args.min_map_score,
            min_score_unlinked=args.min_map_score_unlinked,
            max_gap=args.max_map_gap,
            cache_dir=args.cache_dir,
            refresh=args.refresh,
        )
    if args.update:
        if matcher is not None:
            statuses = [s for r in matcher.resources for s in r.update(force=args.refresh)]
            matcher.reset()
        else:
            statuses = mapper.update(force=args.refresh)
        _print_file_statuses(statuses)
        if any(s.status == "error" for s in statuses):
            return 1
    elif args.check_updates:
        resources = matcher.resources if matcher is not None else mapper.resources
        statuses = [s for r in resources for s in r.check_updates()]
        _print_file_statuses(statuses)
        n = sum(s.update_available for s in statuses)
        print(f"{n} file(s) can be updated with --update" if n else "all files are up to date")
        if any(s.status == "error" for s in statuses):
            return 1
    if args.model is None:
        return 0

    model = load_model(args.model)
    if matcher is not None:
        matches = matcher.match_model(model, skip_boundary=not args.include_boundary)
        _print_map_summary(map_summary(matcher.map, matches))
        if args.output:
            map_matches_to_dataframe(matches, all_candidates=args.all_candidates).to_csv(
                args.output, sep="\t", index=False
            )
            print(f"matches written to {args.output}")
        if args.node_table:
            node_coverage(matcher.map, matches).to_csv(args.node_table, sep="\t", index=False)
            print(f"map nodes written to {args.node_table}")
        if args.svg:
            n = highlight_svg(args.svg, args.svg_out, matches)
            print(f"{n} map nodes highlighted in {args.svg_out}")
        return 0
    if args.evaluate:
        table = evaluate(model, mapper, bigg_only=args.bigg_only)
        print(evaluation_summary(table).to_string())
        if args.output:
            table.to_csv(args.output, sep="\t", index=False)
            print(f"evaluation written to {args.output}")
        return 0
    mappings = mapper.map_model(model)
    _print_summary(summary(mappings), mapper)

    if args.output:
        to_dataframe(mappings, long=args.long).to_csv(args.output, sep="\t", index=False)
        print(f"mapping written to {args.output}")
    if args.write_model:
        import cobra

        changed = annotate_model(model, mappings)
        cobra.io.write_sbml_model(model, str(args.write_model))
        print(f"EC numbers added to {changed} reactions, model written to {args.write_model}")
    return 0
