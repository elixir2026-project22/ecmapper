"""Command line interface: ``python -m ecmapper model.xml -o ec_mapping.tsv``."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

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
    print(f"with EC number: {stats['with_ec']}, of these with more than one: {stats['ambiguous']}")
    print("reactions with EC numbers per source:")
    for source, n in stats["reactions_per_source"].items():
        print(f"  {source:<11} {n}")
    unavailable = [s.name for s in mapper.sources.values() if s.loaded and not s.available]
    if unavailable:
        print(f"unavailable sources (skipped): {', '.join(unavailable)}")


def _print_file_statuses(statuses: list[FileStatus]) -> None:
    print(f"{'source':<10} {'file':<34} {'status':<10} {'cached':<26} remote last modified")
    for s in statuses:
        line = f"{s.source:<10} {s.file:<34} {s.status:<10} {s.cached:<26} {s.remote_modified}"
        print(line + (f" ({s.detail})" if s.detail else ""))


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
    parser.add_argument("--cache-dir", type=Path, default=None, help=f"download cache (default: {default_cache_dir()})")
    parser.add_argument("--refresh", action="store_true", help="download the resource files again")
    parser.add_argument(
        "--check-updates", action="store_true", help="check whether newer resource files are available"
    )
    parser.add_argument(
        "--update", action="store_true", help="download newer resource files into the cache (before mapping)"
    )
    parser.add_argument("--write-model", type=Path, help="write the model with the found EC numbers as SBML")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.model is None and not (args.check_updates or args.update):
        parser.error("a model file, --check-updates or --update is required")

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s: %(message)s")
    mapper = ECMapper(
        [s.strip() for s in args.sources.split(",") if s.strip()],
        use_notes=not args.no_notes,
        always_query=args.always_query,
        skip_boundary=not args.include_boundary,
        follow_xrefs=not args.no_xrefs,
        cache_dir=args.cache_dir,
        refresh=args.refresh,
    )
    if args.update:
        statuses = mapper.update(force=args.refresh)
        _print_file_statuses(statuses)
        if any(s.status == "error" for s in statuses):
            return 1
    elif args.check_updates:
        statuses = mapper.check_updates()
        _print_file_statuses(statuses)
        n = sum(s.update_available for s in statuses)
        print(f"{n} file(s) can be updated with --update" if n else "all files are up to date")
        if any(s.status == "error" for s in statuses):
            return 1
    if args.model is None:
        return 0

    model = load_model(args.model)
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
