"""Verify the Roche <-> BRENDA compound matching with database identifiers and add
KEGG, ChEBI and BiGG IDs to roche_brenda_merged.csv.

Each side is resolved to KEGG Compound IDs independently, by name:
  * Roche: the Roche compound name.
  * BRENDA: all labels and BRENDA names of the ligand ID on the map.
Names are matched exactly (case, spaces, hyphens and commas ignored) against all KEGG
compound names and synonyms, in tiers: full name; name without / inside a trailing
parenthesis; abbreviations expanded (-P -> phosphate, -PP -> diphosphate,
P2 -> bisphosphate, THF -> tetrahydrofolate, Glc/Fru/Rib/... -> sugar names); and
last, without an L-/D- prefix. Names that KEGG does not know are given a KEGG name by
hand in MANUAL_NAMES (the KEGG ID is still looked up from that name), and where a name
hits several KEGG entries the intended one is chosen in PREFERRED.

ChEBI IDs are the KEGG -> ChEBI cross-references from KEGG (EBI was not reachable).
BiGG IDs are the universal BiGG metabolites whose database links name the KEGG ID; if
there are none, those linking one of the ChEBI IDs (bigg_via = "ChEBI").

A row is "confirmed" when the Roche and BRENDA KEGG IDs overlap, "different KEGG ID"
when both resolve but do not overlap, and "unresolved" when a side has no KEGG ID.

Usage: python scripts/verify_roche_brenda_ids.py [cache_dir] [--linked-only]
Downloads (once, into cache_dir): rest.kegg.jp/list/compound, rest.kegg.jp/conv/chebi/compound,
bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt. These are not committed.
"""
import collections
import csv
import json
import re
import sys
import urllib.request
from pathlib import Path

BASE = Path("BRENDAmap/cofactor_analysis")
args = [a for a in sys.argv[1:] if not a.startswith("--")]
CACHE = Path(args[0] if args else Path.home() / ".cache/ecmapper/compound_ids")
SUFFIX = "_linked" if "--linked-only" in sys.argv else ""
URLS = {
    "kegg_compound.txt": "https://rest.kegg.jp/list/compound",
    "kegg_chebi.txt": "https://rest.kegg.jp/conv/chebi/compound",
    "bigg_mets.txt": "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
}
CACHE.mkdir(parents=True, exist_ok=True)
for fn, url in URLS.items():
    if not (CACHE / fn).exists():
        print("downloading", url)
        urllib.request.urlretrieve(url, CACHE / fn)

# Names KEGG does not list under the Roche/BRENDA wording -> the KEGG name meant.
MANUAL_NAMES = {
    # Roche
    "PROTEINS": "Protein",
    "2'-DEOXYRIBOSE 1-P": "2-Deoxy-D-ribose 1-phosphate",
    "(3-PHOSPHATIDYL)-1D-myo-INOSITOL": "1-Phosphatidyl-1D-myo-inositol",
    "5'P-RIBOSYL-N-FORMYLGLYCINAMIDE (FGAR)": "5'-Phosphoribosyl-N-formylglycinamide",
    "1,3-P2-D-GLYCERATE": "3-Phospho-D-glyceroyl phosphate",
    "3-P-D-GLYCERATE": "3-Phospho-D-glycerate",
    "ADP-D-GLUCOSE": "ADP-glucose",
    "2-(alpha-HYDROXYETHYL)-ThPP": "2-(alpha-Hydroxyethyl)thiamine diphosphate",
    "ISOCAPROIC ALDEHYDE": "4-Methylpentanal",
    "L-4-HYDROXYPROLINE": "trans-4-Hydroxy-L-proline",
    "GLYCYL-tRNA": "Glycyl-tRNA(Gly)",
    "L-GLUTAMYL-tRNA": "L-Glutamyl-tRNA(Glu)",
    "L-HISTIDINYL-tRNA": "L-Histidyl-tRNA(His)",
    "L-ISOLEUCYL-tRNA": "L-Isoleucyl-tRNA(Ile)",
    "L-PHENYLALANYL-tRNA": "L-Phenylalanyl-tRNA(Phe)",
    "L-TRYPTOPHANYL-tRNA": "L-Tryptophanyl-tRNA(Trp)",
    "L-TYROSYL-tRNA": "L-Tyrosyl-tRNA(Tyr)",
    "L-VALYL-tRNA": "L-Valyl-tRNA(Val)",
    "N-FORMYL-L-METHIONYL-tRNA": "N-Formylmethionyl-tRNA",
    # BRENDA ligand IDs
    17973: "UDP-N-acetyl-alpha-D-glucosamine",
    6523: "D-Fructose 1,6-bisphosphate",
    43607: "6-Phospho-D-gluconate",
    203: "L-Phenylalanine",
    191248: "trans-Dec-2-enoyl-[acp]",
    86950: "2-Hydroxy-2,4-pentadienoate",
    5068: "2-Hydroxymuconate",
    5073: "4-Hydroxy-2-oxopentanoate",
    2125: "3-Phospho-D-glyceroyl phosphate",
    81591: "2,5-Diamino-6-(5-phospho-D-ribosylamino)pyrimidin-4(3H)-one",
    33002: "2-Deoxy-D-ribose 1-phosphate",
    41650: "beta-Phenylalanoyl-CoA",
    5824: "(1,4-alpha-D-Glucosyl)n",
    10035: "(1,4-alpha-D-Glucosyl)n+1",
    5970: "D-Glucose 1-phosphate",
    914: "tRNA(Cys)",
}
# Names that hit several KEGG entries -> the entry meant (Roche name or BRENDA ligand ID).
PREFERRED = {
    "D-RIBOSE": "C00121", "D-FRUCTOSE": "C00095", 117: "C00095",
    "PHENYLACETATE": "C07086", 776: "C07086",          # not phenyl acetate (ester), C00548
    194: "C00197", 1708: "C00642", 156: "C00090", 1614: "C03576", 621: "C00110",
    506: "C00182", 9: "C00007",
    328: "C00479",   # ligand 328 is propanal; one map label says "propanol"
}


def norm(s):
    s = s.strip().lower().replace("α", "alpha").replace("β", "beta").replace("–", "-")
    charge = s[-1] if len(s) > 1 and s[-1] in "+-" and s[-2] not in " -" else ""
    s = s[:-1] if charge else s
    return re.sub(r"[\s\-_,]", "", s) + ({"-": "⁻", "+": "⁺"}.get(charge, ""))


kegg = collections.defaultdict(set)
kname = {}
for line in open(CACHE / "kegg_compound.txt"):
    c, n = line.rstrip("\n").split("\t")
    names = [x.strip() for x in n.split(";")]
    kname[c] = names[0]
    for x in names:
        kegg[norm(x)].add(c)
chebi = collections.defaultdict(set)
for line in open(CACHE / "kegg_chebi.txt"):
    c, h = line.split()
    chebi[c.split(":")[1]].add("CHEBI:" + h.split(":")[1])
bigg = collections.defaultdict(set)
bigg_chebi = collections.defaultdict(set)
for r in csv.DictReader(open(CACHE / "bigg_mets.txt"), delimiter="\t"):
    for k in re.findall(r"kegg\.compound/(C\d{5})", r["database_links"] or ""):
        bigg[k].add(r["universal_bigg_id"])
    for h in re.findall(r"chebi/(CHEBI:\d+)", r["database_links"] or ""):
        bigg_chebi[h].add(r["universal_bigg_id"])

ABBR = [(r"-PP\b", " diphosphate"), (r" PP\b", " diphosphate"), (r"\bPPi\b", "diphosphate"),
        (r"-P2\b", " bisphosphate"), (r" P2\b", " bisphosphate"), (r"-P\b", " phosphate"),
        (r" P\b", " phosphate"), (r"\bTHF\b", "tetrahydrofolate"), (r"\bGlc\b", "glucose"),
        (r"\bFru\b", "fructose"), (r"\bRib\b", "ribose"), (r"\bRul\b", "ribulose"),
        (r"\bXul\b", "xylulose"), (r"\bGal\b", "galactose")]


def tiers(n):
    t1 = [n]
    m = re.match(r"^(.*?)\s*\(([^()]*)\)\s*$", n)
    t2 = [m.group(1), m.group(2)] if m else []
    t3 = []
    for v in t1 + t2:
        w = v
        for a, b in ABBR:
            w = re.sub(a, b, w)
        if w not in t1 + t2:
            t3.append(w)
    t4 = [x[2:] for x in t1 + t2 + t3 if x[:2] in ("L-", "D-")]
    return [t1, t2, t3, t4]


def resolve(key, names):
    """-> (set of KEGG IDs, method)"""
    if key in MANUAL_NAMES:
        hits = kegg.get(norm(MANUAL_NAMES[key]), set())
        assert hits, f"manual KEGG name not found: {MANUAL_NAMES[key]}"
        hits, method = hits, "manual name"
    else:
        hits, method = set(), "unresolved"
        for k in range(4):
            for n in names:
                for v in tiers(n)[k]:
                    hits |= kegg.get(norm(v), set())
            if hits:
                method = ["exact name", "name part", "abbreviation expanded", "without L-/D-"][k]
                break
    if key in PREFERRED:
        assert PREFERRED[key] in hits, (key, hits)
        return {PREFERRED[key]}, method + ", chosen"
    if len(hits) > 1:
        method += ", ambiguous"
    return hits, method


def ids_info(ks):
    ks = sorted(ks)
    ch = sorted({h for k in ks for h in chebi[k]}, key=lambda x: int(x[6:]))
    bg = sorted({b for k in ks for b in bigg[k]})
    via = ""
    if not bg:  # fall back to BiGG entries that link one of the ChEBI IDs
        bg = sorted({b for h in ch for b in bigg_chebi[h]})
        via = "ChEBI" if bg else ""
    return dict(kegg="; ".join(ks), kegg_name="; ".join(kname[k] for k in ks),
                chebi="; ".join(ch), bigg="; ".join(bg), bigg_via=via or ("KEGG" if bg else ""))


# BRENDA names per ligand ID from the map
d = json.loads(Path("BRENDAmap/BRENDA_map_reactions.json").read_text())
bnames = collections.defaultdict(set)
for m in d["metabolite_nodes"]:
    if m["node_class"] != "nodeCM":
        bnames[m.get("brenda_ligand_id")] |= {m["display_name"], (m["brenda_name"] or "").strip()} - {""}

rows = list(csv.DictReader(open(BASE / f"roche_brenda_merged{SUFFIX}.csv")))
brenda_kegg = {}
out = []
for r in rows:
    rk, rmeth = (resolve(r["compound"], [r["compound"]]) if r["source"] != "BRENDA only" else (set(), ""))
    bk_all, bmeth = set(), []
    for lid in filter(None, r["brenda_ligand_ids"].split("; ")):
        k, m = resolve(int(lid), sorted(bnames[int(lid)]))
        brenda_kegg[int(lid)] = k
        bk_all |= k
        bmeth.append(f"{lid}: {', '.join(sorted(k)) or '-'} ({m})")
    if r["source"] == "both":
        status = ("unresolved" if not rk or not bk_all else
                  "confirmed" if rk & bk_all and all(brenda_kegg[int(l)] & rk for l in r["brenda_ligand_ids"].split("; ")) else
                  "partly confirmed" if rk & bk_all else "different KEGG ID")
    else:
        status = ""
    ri, bi = ids_info(rk), ids_info(bk_all)
    best = rk or bk_all
    r.update(dict(
        match_check=status,
        kegg=ids_info(best)["kegg"], kegg_name=ids_info(best)["kegg_name"],
        chebi=ids_info(best)["chebi"], bigg=ids_info(best)["bigg"], bigg_via=ids_info(best)["bigg_via"],
        roche_kegg=ri["kegg"], roche_kegg_method=rmeth,
        brenda_kegg=bi["kegg"], brenda_kegg_name=bi["kegg_name"], brenda_kegg_method="; ".join(bmeth)))
    out.append(r)

# Missed matches: Roche-only and BRENDA-only rows sharing a KEGG ID
by_kegg = collections.defaultdict(list)
for r in out:
    for k in filter(None, r["kegg"].split("; ")):
        by_kegg[k].append(r)
missed = []
for k, rs in by_kegg.items():
    ro = [x for x in rs if x["source"] == "Roche only"]
    bo = [x for x in rs if x["source"] == "BRENDA only"]
    if ro and bo:
        missed.append((k, ro, bo))
for k, ro, bo in missed:
    for x in ro + bo:
        x["match_check"] = "possible missed match: " + "; ".join(
            y["compound"] for y in (bo if x in ro else ro)) + f" (KEGG {k})"

cols = ["compound", "source", "n_boxes_roche", "n_subsystems_brenda", "n_total", "match_check",
        "kegg", "kegg_name", "chebi", "bigg", "bigg_via", "roche_kegg", "roche_kegg_method", "brenda_ligand_ids",
        "brenda_labels", "brenda_names", "brenda_kegg", "brenda_kegg_name", "brenda_kegg_method",
        "has_cross_subsystem_link", "grid_squares", "as_structure", "as_name_only", "brenda_subsystems"]
with open(BASE / f"roche_brenda_merged{SUFFIX}_ids.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows({c: r.get(c, "") for c in cols} for r in out)

c = collections.Counter(r["match_check"].split(":")[0] for r in out)
print(dict(c))
print("with KEGG:", sum(bool(r["kegg"]) for r in out), "with ChEBI:", sum(bool(r["chebi"]) for r in out),
      "with BiGG:", sum(bool(r["bigg"]) for r in out), "(via ChEBI:", sum(r["bigg_via"] == "ChEBI" for r in out), ") of", len(out))
for r in out:
    if r["bigg_via"] == "ChEBI" or len(r["bigg"].split("; ")) > 1:
        print("  BiGG", r["bigg_via"], "|", r["compound"], "|", r["kegg_name"], "->", r["bigg"])
for r in out:
    if r["match_check"] not in ("", "confirmed"):
        print(f"[{r['match_check']}] {r['compound']} | Roche {r['roche_kegg']} ({r['roche_kegg_method']}) | BRENDA {r['brenda_kegg_method']}")
amb = [r for r in out if "ambiguous" in r["roche_kegg_method"] + r["brenda_kegg_method"]]
print("ambiguous:")
for r in amb:
    print(" ", r["compound"], "|", r["roche_kegg"], r["roche_kegg_method"], "|", r["brenda_kegg_method"])
