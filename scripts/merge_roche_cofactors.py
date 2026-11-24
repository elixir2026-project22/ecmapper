"""Integrate the Roche-chart cofactor combinations into the BRENDA-map table
cofactor_combinations_without_H2O_H+.csv.

Rules (agreed with the user):
* Roche count = number of grid squares listed in example_grid_squares. A Roche row
  with alternatives counts its grid squares for every variant; several Roche rows
  falling on one combination are summed. n_total = BRENDA reactions + Roche squares.
* Alternatives are expanded: "X or Y" (paired by position across the row, so
  "ADP or GDP" -> "ATP or GTP" gives ADP->ATP and GDP->GTP) and explicit lists in
  brackets ("NTPs (ATP, CTP, GTP, UTP)"). Generic classes (UDP-sugar, acyl-CoA,
  NAD(P)+, acceptor ...) only match BRENDA's own generic entries.
* Exact matches only are merged (if synonymous BRENDA labels give several exact
  matches, the Roche count goes to the one with most reactions and the others are
  listed in also_exact_matches_brenda); a combination and its reverse are one entry.
  For Roche variants without an exact match, BRENDA combinations that contain the
  Roche one or are contained in it (in either direction) are listed in near_matches.
* As on the BRENDA side, H2O, H+, O2, CO2, NH3, NH4+ and H2O2 are removed.

Synonymous BRENDA cofactor labels are unified with BRENDA_SYNONYMS (brenda_synonyms.py)
already in brenda_map_cofactors.py; this script applies the mapping again as a no-op
safeguard (merged_brenda_rows would list any rows merged here).

Roche compound names are translated to the BRENDA cofactor labels with ROCHE_NAMES;
a Roche name may stand for several BRENDA labels (synonyms with different BRENDA
ligand IDs, e.g. reduced ferredoxin). Names not listed are used as written.
"""
import csv
import itertools
import re
from collections import defaultdict
from pathlib import Path

BASE = Path("BRENDAmap/cofactor_analysis")
EXCLUDED = {"H2O", "H+", "O2", "CO2", "NH3", "NH4+", "H2O2"}
ROCHE_NAMES = {
    "acetyl-CoA": ["Ac-CoA"], "CoA-SH": ["CoA"], "glycine": ["Gly"],
    "ferredoxin (red)": ["Fd_red", "reduced ferredoxin [Fe-S] cluster"],
    "ferredoxin (ox)": ["Fd_ox", "oxidized ferredoxin [Fe-S] cluster"],
    "propionyl-CoA": ["propanoyl-CoA"], "propionate": ["propanoate"], "acetoacetate": ["AcAc"],
    "FMN·E": ["FMN"], "FMNH2·E": ["FMNH2"],
    "Acceptor": ["acceptor", "an acceptor", "oxidized acceptor", "electron acceptor",
                 "oxidized electron acceptor", "an oxidized electron acceptor", "A"],
    "Acceptor·H2": ["reduced acceptor", "a reduced acceptor", "AH2", "reduced electron acceptor",
                    "a reduced electron acceptor"],
    "ubiquinone": ["UQ"], "ubiquinol": ["UQH2"],
    "thioredoxin-(SH)2": ["red. thioredoxin", "thioredoxin"],
    "thioredoxin-S2": ["ox. thioredoxin", "thioredoxin disulfide"],
    "PRPP": ["5-P-alpha-D-Rib 1-PP"], "ribose 1-P": ["alpha-D-Rib 1-P"], "HCO3-": ["HCO3"],
    "L-aspartate": ["Asp"], "L-glutamate": ["Glu"], "L-glutamine": ["Gln"], "L-alanine": ["Ala"],
    "2-oxoglutarate": ["2OG"], "fumarate": ["FUM"],
    "SAM": ["S-adenosyl-L-methionine"], "SAH": ["S-adenosyl-L-homocysteine"],
    "THF": ["THF", "(6S)-THFA", "5,6,7,8-tetrahydrofolate"],
    "PAPS": ["3'-phosphoadenylyl sulfate", "3'-phosphoadenylylsulfate"],
    "PAP": ["adenosine 3',5'-bisphosphate"], "H2SO4": ["sulfate"],
    "polyprenyl-PP": ["an all-trans-polyprenyl PP"], "phytyl-PP": ["phytyl-PP"],
}


def fmt(names):
    names = list(names)
    return "∅" if not names else names[0] if len(names) == 1 else "(" + ", ".join(names) + ")"


def split_items(field):
    """'a; b or c; d' -> [['a'], ['b', 'c'], ['d']]; bracket lists become alternatives."""
    items = []
    for it in filter(None, (x.strip() for x in field.split(";"))):
        m = re.match(r"^[^()]*\(([^()]*,[^()]*)\)$", it)  # "NTPs (ATP, CTP, GTP, UTP)"
        alts = [a.strip() for a in m.group(1).split(",")] if m else [a.strip() for a in it.split(" or ")]
        items.append(alts)
    return items


def variants(consumed, released):
    c, p = split_items(consumed), split_items(released)
    n = max([len(a) for a in c + p] + [1])
    out = []
    for k in range(n):
        pick = lambda items: tuple(sorted({a[k] if len(a) == n else a[0] for a in items} - EXCLUDED, key=str.lower))
        out.append((pick(c), pick(p)))
    return list(dict.fromkeys(out))


def labels(name):
    return ROCHE_NAMES.get(name, [name])


def side_matches(roche_side, brenda_side):
    """Exact match allowing synonyms: a one-to-one assignment of Roche names to BRENDA labels."""
    if len(roche_side) != len(brenda_side):
        return False
    return any(all(b in labels(r) for r, b in zip(roche_side, perm))
               for perm in itertools.permutations(brenda_side))


def covered(small, big, small_is_roche):
    """Every compound of `small` occurs in `big`."""
    if small_is_roche:
        return all(any(b in labels(r) for b in big) for r in small)
    return all(any(b in labels(r) for r in big) for b in small)


def shared(rs, bs):
    return sum(any(b in labels(r) for b in bs) for r in rs)


# Synonymous BRENDA cofactor labels -> one label (see brenda_synonyms.py). The BRENDA
# table is already built with them; applying them again here is a no-op safeguard.
from brenda_synonyms import BRENDA_SYNONYMS  # noqa: E402


def relabel(side):
    return tuple(sorted({BRENDA_SYNONYMS.get(x, x) for x in side}, key=str.lower))


def merge_brenda_synonyms(raw):
    groups = {}
    for b in raw:
        A = relabel(filter(None, b["side_a"].split(" + ")))
        B = relabel(filter(None, b["side_b"].split(" + ")))
        g = groups.setdefault(tuple(sorted((A, B))), dict(A=A, B=B, fwd=0, rev=0, src=[], ec=set(),
                                                           subs=set(), var=0))
        same = (A, B) == (g["A"], g["B"])
        g["fwd"] += int(b["n_a_to_b"] if same else b["n_b_to_a"])
        g["rev"] += int(b["n_b_to_a"] if same else b["n_a_to_b"])
        g["src"].append(b)
        g["ec"] |= set(filter(None, b["ec_numbers"].split("; ")))
        g["subs"] |= set(filter(None, b["subsystems"].split("; ")))
        g["var"] += int(b["n_variants_with_H2O_H+"])
    out = []
    for g in groups.values():
        A, B, fwd, rev = g["A"], g["B"], g["fwd"], g["rev"]
        if rev > fwd:
            A, B, fwd, rev = B, A, rev, fwd
        out.append(dict(
            A=A, B=B, side_a=" + ".join(A), side_b=" + ".join(B), pair=f"({fmt(A)}, {fmt(B)})",
            n_reactions=str(fwd + rev), n_a_to_b=str(fwd), n_b_to_a=str(rev),
            n_subsystems=str(len(g["subs"])), ec_numbers="; ".join(sorted(g["ec"])),
            subsystems="; ".join(sorted(g["subs"])), **{"n_variants_with_H2O_H+": str(g["var"])},
            merged_brenda_rows="; ".join(f"{x['pair']} [{x['n_reactions']}]" for x in g["src"])
            if len(g["src"]) > 1 or g["src"][0]["pair"] != f"({fmt(A)}, {fmt(B)})" else "",
            roche=[]))
    out.sort(key=lambda x: (-int(x["n_reactions"]), x["pair"].lower()))
    return out


brenda_raw = list(csv.DictReader(open(BASE / "cofactor_combinations_without_H2O_H+.csv")))
brenda = merge_brenda_synonyms(brenda_raw)
roche = list(csv.DictReader(open(BASE / "roche_cofactor_combinations.csv")))

roche_only = {}  # unordered key -> entry
for r in roche:
    squares = [s.strip() for s in r["example_grid_squares"].split(";") if s.strip()]
    for C, P in variants(r["consumed"], r["released"]):
        hits = []
        for b in brenda:  # sorted by reaction count, so the first hit is the largest
            if side_matches(C, b["A"]) and side_matches(P, b["B"]):
                hits.append((b, False))
            elif side_matches(C, b["B"]) and side_matches(P, b["A"]):
                hits.append((b, True))
        hit, rev = hits[0] if hits else (None, None)
        # Synonymous BRENDA labels can give several exact matches: the Roche count goes to
        # the largest one, the others are listed.
        rec = dict(row=r, C=C, P=P, n=len(squares), squares=squares,
                   also=[f"{b['pair']} [{b['n_reactions']}]" for b, _ in hits[1:]])
        if hit:
            rec["reverse"] = rev
            hit["roche"].append(rec)
            continue
        key = tuple(sorted((C, P)))
        e = roche_only.setdefault(key, dict(A=C, B=P, roche=[]))
        rec["reverse"] = (C, P) != (e["A"], e["B"])
        e["roche"].append(rec)


def near_matches(C, P):
    out = []
    size = len(C) + len(P)
    for b in brenda:
        if not b["A"] and not b["B"]:
            continue
        for X, Y in ((b["A"], b["B"]), (b["B"], b["A"])):
            sub = covered(X, C, False) and covered(Y, P, False)   # BRENDA within Roche
            sup = covered(C, X, True) and covered(P, Y, True)     # Roche within BRENDA
            if (sub or sup) and shared(C, X) + shared(P, Y) >= min(2, size):
                out.append((int(b["n_reactions"]), b["pair"], "BRENDA ⊂ Roche" if sub else "Roche ⊂ BRENDA"))
                break
    out.sort(key=lambda t: -t[0])
    return "; ".join(f"{p} [{n}, {kind}]" for n, p, kind in out[:5])


def roche_cols(recs):
    def j(f):
        return " | ".join(dict.fromkeys(f(x) for x in recs if f(x)))
    return dict(
        n_grid_squares_roche=sum(x["n"] for x in recs),
        roche_as_drawn=j(lambda x: f"({fmt(x['C'])} → {fmt(x['P'])})"),
        roche_rows=j(lambda x: f"{x['row']['consumed']} → {x['row']['released']}"),
        roche_group=j(lambda x: x["row"]["group"]),
        roche_grid_squares=j(lambda x: "; ".join(x["squares"])),
        roche_example_reactions=j(lambda x: x["row"]["example_reactions"]),
        roche_drawn_reversible=j(lambda x: x["row"]["drawn_reversible"]),
        roche_uncertain=j(lambda x: x["row"]["uncertain"]),
        roche_note=j(lambda x: x["row"]["note"]),
        also_exact_matches_brenda=j(lambda x: "; ".join(x.get("also", []))),
    )


rows = []
for b in brenda:
    rc = roche_cols(b["roche"]) if b["roche"] else {}
    rows.append(dict(
        pair=b["pair"], side_a=b["side_a"], side_b=b["side_b"],
        source="both" if b["roche"] else "BRENDA only",
        n_reactions_brenda=int(b["n_reactions"]), n_a_to_b=b["n_a_to_b"], n_b_to_a=b["n_b_to_a"],
        n_grid_squares_roche=rc.get("n_grid_squares_roche", 0), **{k: v for k, v in rc.items() if k != "n_grid_squares_roche"},
        near_matches="", n_subsystems_brenda=b["n_subsystems"], ec_numbers=b["ec_numbers"],
        subsystems=b["subsystems"], n_variants_with_H2O_H=b["n_variants_with_H2O_H+"],
        merged_brenda_rows=b["merged_brenda_rows"]))
for e in roche_only.values():
    rc = roche_cols(e["roche"])
    rows.append(dict(
        pair=f"({fmt(e['A'])}, {fmt(e['B'])})", side_a=" + ".join(e["A"]), side_b=" + ".join(e["B"]),
        source="Roche only", n_reactions_brenda=0, n_a_to_b="", n_b_to_a="", **rc,
        near_matches=near_matches(e["A"], e["B"])))
for x in rows:
    x["n_total"] = x["n_reactions_brenda"] + x["n_grid_squares_roche"]
rows.sort(key=lambda x: (-x["n_total"], -x["n_grid_squares_roche"], x["pair"].lower()))

cols = ["pair", "side_a", "side_b", "source", "n_reactions_brenda", "n_grid_squares_roche", "n_total",
        "n_a_to_b", "n_b_to_a", "roche_as_drawn", "roche_rows", "roche_group", "roche_grid_squares",
        "roche_example_reactions", "roche_drawn_reversible", "roche_uncertain", "roche_note",
        "also_exact_matches_brenda", "near_matches", "merged_brenda_rows", "n_subsystems_brenda", "ec_numbers", "subsystems", "n_variants_with_H2O_H"]
out = BASE / "cofactor_combinations_without_H2O_H+_with_roche.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows({c: x.get(c, "") for c in cols} for x in rows)

from collections import Counter
print(Counter(x["source"] for x in rows), "BRENDA rows", len(brenda_raw), "->", len(brenda))
for b in brenda:
    if b["merged_brenda_rows"].count("; ") >= 1 or (b["merged_brenda_rows"] and False):
        print("MERGED", b["pair"], b["n_reactions"], "<=", b["merged_brenda_rows"])
n_var = sum(len(variants(r["consumed"], r["released"])) for r in roche)
print("Roche rows", len(roche), "variants", n_var)
for x in rows:
    if x["source"] != "BRENDA only":
        print(f"{x['n_total']:>4} {x['n_reactions_brenda']:>4} {x['n_grid_squares_roche']:>3} {x['source']:<10} {x['pair'][:70]:<70} <- {x['roche_rows'][:60]}")
        if x["near_matches"]:
            print("        near:", x["near_matches"][:230])
