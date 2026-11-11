"""Command line: ``python -m ecmapper.pwl {generate,check,parse,motifs} ...``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _resources(args):
    from ..sources import BKMSSource, KEGGCompoundNames

    bkms = None if args.no_bkms else BKMSSource(cache_dir=args.cache_dir)
    if bkms is not None:
        bkms.data  # noqa: B018, loads (downloads) the file
        if bkms.error:
            print(f"BKMS-react not available ({bkms.error}); reactions get no direction", file=sys.stderr)
            bkms = None
    kegg = None
    if not getattr(args, "no_kegg", True):
        names = KEGGCompoundNames(cache_dir=args.cache_dir)
        kegg = names.data or None
        if names.error:
            print(f"KEGG compound names not available ({names.error})", file=sys.stderr)
    return bkms, kegg


def _check(args, data, bkms) -> int:
    from .check import check, summary

    report = check(data, args.outdir, bkms)
    s = summary(report)
    with open(Path(args.outdir) / "validation_report.json", "w", encoding="utf-8") as f:
        json.dump(s, f, indent=1)
        f.write("\n")
    print(json.dumps({k: v for k, v in s.items() if not k.endswith("_examples")}, indent=1))
    for k in ("parse_errors", "mismatches", "direction_mismatches", "motif_problems"):
        for x in report[k][:10]:
            print(f"{k}: {x}")
    return 1 if report["parse_errors"] or report["mismatches"] else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m ecmapper.pwl", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    g = sub.add_parser("generate", help="write PWL files for the subsystems of a BRENDA map export and check them")
    g.add_argument("map_json", type=Path, help="BRENDA_map_reactions.json")
    g.add_argument("outdir", type=Path)
    g.add_argument("--subsystem", action="append", help="only these subsystems (repeatable)")
    g.add_argument("--no-bkms", action="store_true", help="do not use BKMS-react (no sides and directions)")
    g.add_argument("--no-kegg", action="store_true", help="do not use KEGG compound names (synonyms, KEGG IDs)")
    g.add_argument("--cache-dir", type=Path)
    c = sub.add_parser("check", help="check PWL files against the BRENDA map export")
    c.add_argument("map_json", type=Path)
    c.add_argument("outdir", type=Path)
    c.add_argument("--no-bkms", action="store_true")
    c.add_argument("--cache-dir", type=Path)
    p = sub.add_parser("parse", help="parse and validate PWL files, print a summary")
    p.add_argument("files", nargs="+", type=Path)
    m = sub.add_parser("motifs", help="detect motifs in PWL files")
    m.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args(argv)

    if args.command in ("generate", "check"):
        data = json.loads(args.map_json.read_text(encoding="utf-8"))
        bkms, kegg = _resources(args)
        if args.command == "generate":
            from .brenda import generate

            summary = generate(data, args.outdir, bkms, kegg, args.subsystem)
            print(json.dumps(summary["totals"], indent=1))
            if args.subsystem:
                return 0
        return _check(args, data, bkms)

    from .motifs import check_motifs, detect_motifs, format_motif
    from .syntax import PWLError, parse_file

    status = 0
    for path in args.files:
        try:
            doc = parse_file(path)
        except PWLError as e:
            print(f"{path}: {e}")
            status = 1
            continue
        if args.command == "parse":
            rs = doc.reactions()
            print(f"{path}: {doc.id}: {len(rs)} reactions, {sum(r.substrates is not None for r in rs)} with direction, "
                  f"{sum(not r.sided for r in rs)} without sides, {len(doc.compounds)} compounds, "
                  f"{len(doc.motifs)} motifs")
            for problem in check_motifs(doc):
                print(f"  {problem}")
                status = 1
        else:
            for motif in detect_motifs(doc):
                print(f"{path}: {format_motif(motif)}")
    return status


if __name__ == "__main__":
    sys.exit(main())
