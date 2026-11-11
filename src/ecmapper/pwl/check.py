"""Check generated PWL files against the BRENDA map export they were made from.

The check parses every file (syntax and :func:`~ecmapper.pwl.syntax.validate`)
and compares it with ``BRENDA_map_reactions.json``:

- every reaction appears exactly once, under its node ID, in the file of its subsystem;
- EC label, spontaneous flag, main-metabolite nodes and consumed/produced cofactors match;
- every compound is declared with the BRENDA ligand IDs of its drawn metabolites
  (IDs that were not merged into the compound are listed in ``compounds.json``);
- with BKMS-react: substrates and products of every reaction with a direction are
  on the sides of its BRENDA equation (``br=``) that the evidence (``ev=``) says;
- every motif is supported by the reactions (:func:`~ecmapper.pwl.motifs.check_motif`).
"""

from __future__ import annotations

import collections
import json
import re
from pathlib import Path

from .compounds import CompoundRegistry, loose_key, parse_equation, split_location
from .motifs import check_motifs
from .syntax import Document, PWLError, parse_file


def _cofactor_tokens(registry: CompoundRegistry, entries) -> collections.Counter:
    return collections.Counter(registry.token(e) for e in entries)


def _equation_side(names: set[str], equation) -> set[str]:
    keys = {loose_key(n) for n in names}
    return {side for side, terms in zip("LR", equation) for t in terms if loose_key(t.name) in keys}


def check(data: dict, pwl_dir, bkms=None) -> dict:
    """Check the files in `pwl_dir` (``subsystems/``, ``node_index.json``, ``compounds.json``); returns a report."""
    pwl_dir = Path(pwl_dir)
    node_index = json.loads((pwl_dir / "node_index.json").read_text(encoding="utf-8"))
    compounds_json = json.loads((pwl_dir / "compounds.json").read_text(encoding="utf-8"))
    registry = CompoundRegistry(data)
    rx = {r["node_id"]: r for r in data["reactions"]}
    by_sub = collections.defaultdict(list)
    for r in data["reactions"]:
        by_sub[r["subsystem"]].append(r)
    report = dict(files=0, reactions_checked=0, parse_errors=[], mismatches=[], direction_checked=0,
                  direction_mismatches=[], motif_problems=[], statements=collections.Counter(),
                  reactions_in=collections.Counter(), motifs=collections.Counter())
    seen = collections.Counter()
    subsystem_ids = {s["id"] for s in data.get("subsystems", [])} | set(by_sub)

    def bad(sid, text):
        report["mismatches"].append(f"{sid}: {text}")

    # compound declarations against the source (ligand IDs)
    for info in registry:
        entry = compounds_json.get(info.token)
        if entry is None:
            bad("compounds.json", f"compound {info.token} missing")
            continue
        if sorted(info.ligand_ids) != entry["brenda_ligand_ids"]:
            bad("compounds.json", f"{info.token}: ligand IDs {entry['brenda_ligand_ids']} != {sorted(info.ligand_ids)}")
    explained = {(x["ligand_id"], x["display_name"]) for v in compounds_json.values() for x in v["reassigned"]}

    for sid in sorted(by_sub):
        path = pwl_dir / "subsystems" / f"{sid}.pwl"
        report["files"] += 1
        try:
            doc: Document = parse_file(path)
        except (PWLError, OSError) as e:
            report["parse_errors"].append(f"{sid}: {e}")
            continue
        if doc.id != sid:
            bad(sid, f"subsystem ID {doc.id}")
        tokmap = node_index.get(sid, {})
        reactions = doc.reactions()
        for r in reactions:
            report["reactions_in"][r.context] += 1
            src = rx.get(r.rid)
            if src is None or src["subsystem"] != sid:
                bad(sid, f"unknown reaction {r.rid}")
                continue
            seen[r.rid] += 1
            want = sorted({m["node_id"] for m in src["main_metabolites"]})
            got = sorted({tokmap.get(p.node, "?" + p.node) for p in r.participants})
            if want != got:
                bad(sid, f"{r.rid}: main metabolites {got} != {want}")
            ec = src["ec_label"].strip()
            spont = ec.lower() == "spontaneous"
            expected = "?" if (spont or ec.startswith("?")) else ec
            if r.enzyme_label != expected or r.spontaneous != spont:
                bad(sid, f"{r.rid}: enzyme {r.enzyme_label} != {expected}")
            if collections.Counter(p.token for p in r.consumed) != _cofactor_tokens(
                    registry, src["cofactors"]["consumed"] + src["cofactors"].get("undirected", [])):
                bad(sid, f"{r.rid}: consumed cofactors differ")
            produced = _cofactor_tokens(registry, src["cofactors"]["produced"])
            if collections.Counter(p.token for p in r.produced) != produced:
                bad(sid, f"{r.rid}: produced cofactors differ")
            # declared ligand IDs cover the drawn metabolites
            for m in src["main_metabolites"] + src["cofactors"]["consumed"] + src["cofactors"]["produced"]:
                token = registry.token(m)
                compound = doc.compounds.get(token)
                lig = m.get("brenda_ligand_id")
                if compound is None:
                    bad(sid, f"{r.rid}: compound {token} not declared")
                elif lig is not None and str(lig) not in compound.xref("brenda") and (
                        lig, split_location(m["display_name"])[0]) not in explained:
                    bad(sid, f"{r.rid}: ligand ID {lig} of {m['display_name']} not declared")
            if bkms is not None and r.substrates is not None and "br" in r.attrs:
                _check_direction(r, doc, compounds_json, bkms, report, sid)
        for s in doc.statements:
            report["statements"][type(s).__name__.lower()] += 1
        for m in doc.motifs:
            report["motifs"][m.kind] += 1
        report["motif_problems"] += [f"{sid}: {p}" for p in check_motifs(doc)]
        drawn = {tokmap.get(n) for n in doc.nodes()}
        if set(tokmap.values()) - drawn:
            bad(sid, "drawn nodes not used: " + ", ".join(sorted(set(tokmap.values()) - drawn)[:5]))
        for port in doc.ports:
            if port.target and port.target not in subsystem_ids:
                bad(sid, f"port to unknown subsystem {port.target}")
        report["reactions_checked"] += len(by_sub[sid])
    missing = [rid for rid in rx if seen[rid] == 0]
    twice = [rid for rid, n in seen.items() if n > 1]
    if missing:
        bad("all", f"{len(missing)} reactions missing, e.g. {missing[:5]}")
    if twice:
        bad("all", f"{len(twice)} reactions written more than once, e.g. {twice[:5]}")
    for key in ("statements", "reactions_in", "motifs"):
        report[key] = dict(sorted(report[key].items()))
    return report


def _check_direction(r, doc, compounds_json, bkms, report, sid):
    """Substrates and products on the sides of the BRENDA equation that the evidence says."""
    ev = r.attrs.get("ev")
    equations = bkms.equations(r.attrs["br"]) or bkms.equations(re.sub(r"-[a-z]$", "", r.attrs["br"]))
    parsed = [e for e in (parse_equation(eq) for _ec, eq in equations) if e]
    if not parsed:
        report["direction_mismatches"].append(f"{sid}: {r.rid}: no equation for {r.attrs['br']}")
        return

    def names(p):
        c = compounds_json.get(p.token, {})
        return set(c.get("display_names", {})) | set(c.get("brenda_names", {}))

    expected_sub = "R" if ev == "cof-eq" else "L"
    report["direction_checked"] += 1
    if r.reversible:
        # <->: either direction, only the sides must be opposite
        for equation in parsed:
            left = set().union(set(), *(_equation_side(names(p), equation) for p in r.left))
            right = set().union(set(), *(_equation_side(names(p), equation) for p in r.right))
            if not (len(left) == 1 and left == right):
                return
        report["direction_mismatches"].append(f"{sid}: {r.rid}: both sides on one side of {r.attrs['br']}")
        return
    for equation in parsed:
        ok = True
        for p, want in [(p, expected_sub) for p in r.substrates] + [
                (p, "L" if expected_sub == "R" else "R") for p in r.products]:
            sides = _equation_side(names(p), equation)
            if sides and want not in sides:
                ok = False
        for p, consumed in [(p, True) for p in r.consumed] + [(p, False) for p in r.produced]:
            sides = _equation_side(names(p), equation)
            want = "L" if consumed == (ev != "cof-eq") else "R"
            if sides == {"L", "R"} or not sides or ev == "eq":
                continue
            if want not in sides:
                ok = False
        if ok:
            return
    report["direction_mismatches"].append(f"{sid}: {r.rid}: sides contradict {r.attrs['br']} ({ev})")


def summary(report: dict) -> dict:
    """The report with the lists replaced by their lengths (and the first entries)."""
    out = {}
    for k, v in report.items():
        if isinstance(v, list):
            out[k] = len(v)
            if v:
                out[k + "_examples"] = v[:10]
        else:
            out[k] = v
    return out
