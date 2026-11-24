"""List cofactor (secondary-metabolite) combinations and cross-subsystem
metabolite connections on the BRENDA map.

Primary metabolites = main nodes (nodeHM / nodeSM); secondary metabolites =
cofactor nodes (nodeCM), split into consumed/produced by arrow direction.
Compounds are identified by BRENDA ligand ID; each ID is shown under its most
frequent cofactor label on the map.
"""
import collections
import csv
import json
import sys
from pathlib import Path

src = Path(sys.argv[1] if len(sys.argv) > 1 else "BRENDAmap/BRENDA_map_reactions.json")
out = Path(sys.argv[2] if len(sys.argv) > 2 else "BRENDAmap/cofactor_analysis")
out.mkdir(parents=True, exist_ok=True)
d = json.loads(src.read_text())
sub_label = {s["id"]: s["label_text"] for s in d["subsystems"]}

# Labels whose ligand ID belongs to a different compound (label kept, ID ignored).
MISLINKED = {(12, "NADH"), (17, "dimethylallyl-PP"), (157, "PPi"), (5, "NAD(P)H"),
             (3252, "reduced acceptor"), (811, "malonyl-CoA")}


def key(m):
    lid = m.get("brenda_ligand_id")
    return ("label", m["display_name"]) if (lid, m["display_name"]) in MISLINKED or lid is None else lid


# canonical label per ligand ID
labels = collections.defaultdict(collections.Counter)
for r in d["reactions"]:
    for k in ("consumed", "produced", "undirected"):
        for m in r["cofactors"][k]:
            labels[key(m)][m["display_name"]] += 1
canon = {lid: c.most_common(1)[0][0] for lid, c in labels.items()}
canon[1], canon[12] = "H2O", "H+"
aliases = {lid: sorted(c) for lid, c in labels.items() if len(c) > 1}
WATER_PROTON = {canon.get(1), canon.get(12)}  # H2O, H+


def fmt(names):
    names = list(names)
    if not names:
        return "∅"
    return names[0] if len(names) == 1 else "(" + ", ".join(names) + ")"


# Small inorganic molecules left out of the combinations (matched on the canonical label;
# O2 has two ligand IDs and H2O2 two, NH4+ is the ammonium form of NH3).
EXCLUDED = {"O2", "CO2", "NH3", "NH4+", "H2O2"}


def side(ms):
    return tuple(sorted({canon[key(m)] for m in ms} - EXCLUDED, key=str.lower))


combos = collections.defaultdict(list)
no_cof = only_excluded = 0
for r in d["reactions"]:
    C, P = side(r["cofactors"]["consumed"]), side(r["cofactors"]["produced"])
    if not C and not P:
        if any(r["cofactors"][k] for k in ("consumed", "produced")):
            only_excluded += 1
        else:
            no_cof += 1
        continue
    combos[(C, P)].append(r)

def merge_reverse(directed):
    """Merge (A, B) and (B, A) into one entry, shown in its more frequent direction.

    directed: {(consumed, produced): [reactions]}. Returns rows with side_a/side_b,
    n_fwd (side_a consumed -> side_b produced), n_rev (reverse) and the reactions.
    """
    groups = collections.defaultdict(dict)
    for (C, P), rs in directed.items():
        groups[tuple(sorted((C, P)))][(C, P)] = rs
    out_rows = []
    for dirs in groups.values():
        # orientation: most reactions first, then alphabetical
        (A, B) = min(dirs, key=lambda cp: (-len(dirs[cp]), fmt(cp[0]).lower(), fmt(cp[1]).lower()))
        fwd = dirs.get((A, B), [])
        rev = dirs.get((B, A), []) if (B, A) != (A, B) else []
        rs = fwd + rev
        out_rows.append(dict(
            side_a=A, side_b=B, n=len(rs), n_fwd=len(fwd), n_rev=len(rev),
            subs=sorted({r["subsystem"] for r in rs}),
            ecs=sorted({e for r in rs for e in r["ec_numbers"]}), reactions=rs))
    out_rows.sort(key=lambda x: (-x["n"], fmt(x["side_a"]).lower(), fmt(x["side_b"]).lower()))
    return out_rows


def strip_wp(side):
    return tuple(x for x in side if x not in WATER_PROTON)


n_directed = len(combos)
rows = merge_reverse(combos)
HEADER = ["side_a", "side_b", "pair", "n_reactions", "n_a_to_b", "n_b_to_a", "n_subsystems", "ec_numbers", "subsystems"]


def write_rows(path, rs, extra=None):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER + (extra[0] if extra else []))
        for x in rs:
            w.writerow([" + ".join(x["side_a"]), " + ".join(x["side_b"]),
                        f"({fmt(x['side_a'])}, {fmt(x['side_b'])})", x["n"], x["n_fwd"], x["n_rev"],
                        len(x["subs"]), "; ".join(x["ecs"]), "; ".join(x["subs"])] + (extra[1](x) if extra else []))


write_rows(out / "cofactor_combinations.csv", rows)

# collapsed view without H2O / H+
core_directed = collections.defaultdict(list)
for (C, P), rs in combos.items():
    core_directed[(strip_wp(C), strip_wp(P))].extend(rs)
core_rows = merge_reverse(core_directed)
for x in core_rows:
    x["variants"] = len({(side(r["cofactors"]["consumed"]), side(r["cofactors"]["produced"])) for r in x["reactions"]})
write_rows(out / "cofactor_combinations_without_H2O_H+.csv", core_rows,
           (["n_variants_with_H2O_H+"], lambda x: [x["variants"]]))

# ---------- cross-subsystem metabolites ----------
nodes = {m["node_id"]: m for m in d["metabolite_nodes"]}
links = list(d["metabolite_links"]) + [dict(l, inter_pathway_class=False) for l in d["unassigned_links"]
                                       if l.get("start_node_id") in nodes and l.get("end_node_id") in nodes]
link_rows = []
for l in links:
    a, b = nodes[l["start_node_id"]], nodes[l["end_node_id"]]
    link_rows.append(dict(
        name_a=a["display_name"], name_b=b["display_name"],
        lid_a=a.get("brenda_ligand_id"), lid_b=b.get("brenda_ligand_id"),
        brenda_name=a["brenda_name"].strip(),
        sub_a=a["subsystem"], sub_b=b["subsystem"], cls_a=a["node_class"], cls_b=b["node_class"],
        node_a=a["node_id"], node_b=b["node_id"], inter=l.get("inter_pathway_class", True)))
with open(out / "cross_subsystem_links.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metabolite", "brenda_ligand_id", "name_end", "brenda_ligand_id_end", "subsystem_start",
                "subsystem_end", "node_start", "node_end", "same_ligand_id", "inter_pw_link_class"])
    for x in sorted(link_rows, key=lambda x: (x["name_a"].lower(), x["sub_a"], x["sub_b"])):
        w.writerow([x["name_a"], x["lid_a"], x["name_b"], x["lid_b"], x["sub_a"], x["sub_b"],
                    x["node_a"], x["node_b"], x["lid_a"] == x["lid_b"], x["inter"]])

# group per compound (ligand ID of start node; mismatches flagged)
per_met = collections.defaultdict(lambda: dict(names=set(), subs=set(), n=0, edges=[]))
for x in link_rows:
    g = per_met[x["lid_a"]]
    g["names"] |= {x["name_a"], x["name_b"]}
    g["brenda"] = x["brenda_name"]
    g["subs"] |= {x["sub_a"], x["sub_b"]}
    g["n"] += 1
    g["edges"].append(x)
met_rows = sorted(per_met.items(), key=lambda kv: (-len(kv[1]["subs"]), sorted(kv[1]["names"])[0].lower()))

# primary metabolites shared by several subsystems (with or without a drawn link)
prim = collections.defaultdict(set)
prim_name = collections.defaultdict(collections.Counter)
for m in d["metabolite_nodes"]:
    if m["node_class"] in ("nodeHM", "nodeSM") and m.get("brenda_ligand_id") is not None:
        prim[m["brenda_ligand_id"]].add(m["subsystem"])
        prim_name[m["brenda_ligand_id"]][m["display_name"]] += 1
linked_lids = {x[k] for x in link_rows for k in ("lid_a", "lid_b")}
shared = sorted(((lid, s) for lid, s in prim.items() if len(s) > 1), key=lambda t: (-len(t[1]), prim_name[t[0]].most_common(1)[0][0].lower()))
with open(out / "primary_metabolites_in_multiple_subsystems.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metabolite", "brenda_ligand_id", "n_subsystems", "has_cross_subsystem_link", "subsystems"])
    for lid, s in shared:
        w.writerow([prim_name[lid].most_common(1)[0][0], lid, len(s), lid in linked_lids, "; ".join(sorted(s))])

# ---------- markdown report ----------
L = []
L.append("# BRENDA map: cofactor combinations and cross-subsystem metabolites\n")
L.append(f"Generated by `scripts/brenda_map_cofactors.py` from `{src.name}`.\n")
L.append("Primary metabolites are the main nodes (`nodeHM`/`nodeSM`). Secondary metabolites are the cofactor nodes "
         "(`nodeCM`), which the map splits into consumed and produced by arrow direction. Compounds are matched by BRENDA "
         "ligand ID and shown under their most frequent label. H₂O and H⁺ are kept. A combination and its "
         "reverse, e.g. (NAD+, NADH) and (NADH, NAD+), are one entry, shown as `(A, B)` in the more frequent "
         "direction; `A→B` and `B→A` count the reactions drawn each way. ∅ means that side is empty.\n")
L.append(f"- Excluded from the combinations: {', '.join(sorted(EXCLUDED))}.")
L.append(f"- Reactions: {len(d['reactions'])}. {len(d['reactions']) - no_cof - only_excluded} have at least one remaining cofactor, "
         f"{only_excluded} have only excluded cofactors and {no_cof} have none.")
L.append(f"- Distinct combinations: **{len(rows)}** after merging reverse pairs ({n_directed} counting each direction separately), "
         f"or {len(core_rows)} once H₂O and H⁺ are also removed "
         "(see `cofactor_combinations_without_H2O_H+.csv`).")
L.append(f"- Distinct cofactor compounds: {len(canon)}.\n")
if aliases:
    L.append("Labels merged because they share a ligand ID: " + "; ".join(
        f"{canon[l]} ← {', '.join(a for a in v if a != canon[l])}" for l, v in sorted(aliases.items(), key=lambda t: canon[t[0]].lower())) + "\n")
    L.append("Labels kept separate because their ligand ID belongs to another compound: "
             + ", ".join(f"{n} (ID {i})" for i, n in sorted(MISLINKED)) + ".\n")
L.append("## 1. Cofactor combinations (all, including H₂O and H⁺)\n")
L.append("| # | Combination (A, B) | Reactions | A→B | B→A | Subsystems |")
L.append("|---|---|---|---|---|---|")
for i, x in enumerate(rows, 1):
    L.append(f"| {i} | ({fmt(x['side_a'])}, {fmt(x['side_b'])}) | {x['n']} | {x['n_fwd']} | {x['n_rev']} | {len(x['subs'])} |")
L.append("\n## 2. Metabolites connected across subsystems\n")
L.append(f"{len(link_rows)} metabolite–metabolite links in total: {sum(x['inter'] for x in link_rows)} with the `inter-pw-link` class, "
         f"plus {sum(not x['inter'] for x in link_rows)} without it (node1999 → node2003). Together they involve "
         f"{len(per_met)} compounds. The table is sorted by the number of subsystems connected.\n")
L.append("| Metabolite | Ligand ID | Links | Subsystems connected |")
L.append("|---|---|---|---|")
for lid, g in met_rows:
    names = sorted(g["names"], key=str.lower)
    L.append(f"| {' / '.join(names)} | {lid} | {g['n']} | {'; '.join(sub_label.get(s, s) for s in sorted(g['subs']))} |")
mism = [x for x in link_rows if x["lid_a"] != x["lid_b"]]
L.append(f"\n**Links joining different ligand IDs ({len(mism)}):**\n")
for x in mism:
    L.append(f"- {x['name_a']} (ID {x['lid_a']}, {sub_label.get(x['sub_a'], x['sub_a'])}) ↔ {x['name_b']} (ID {x['lid_b']}, {sub_label.get(x['sub_b'], x['sub_b'])})")
L.append(f"\nFor comparison, {len(shared)} primary compounds (by ligand ID) appear in two or more subsystems, but only "
         f"{sum(lid in linked_lids for lid, _ in shared)} of them have a drawn cross-subsystem link. See "
         "`primary_metabolites_in_multiple_subsystems.csv`.\n")
(out / "cofactor_and_crosslink_report.md").write_text("\n".join(L) + "\n")
print(len(rows), len(core_rows), len(link_rows), len(per_met), len(shared), len(mism))
