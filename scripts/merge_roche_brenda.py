"""Merge the Roche-chart boxed compounds with the BRENDA-map metabolites that occur
in several subsystems (primary_metabolites_in_multiple_subsystems.csv).

Matching is by a hand-curated map from Roche name to BRENDA ligand IDs (names differ
too much for automatic matching). Several BRENDA IDs can map to one Roche compound:
synonyms with separate IDs (L-malate / (S)-malate, both DHAP IDs) and generic forms
mapped to the only form Roche shows (D-Glc 6-P -> alpha-D-GLUCOSE 6-P). The BRENDA
count is then the number of distinct subsystems over all mapped IDs.
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

base = Path(sys.argv[1] if len(sys.argv) > 1 else "BRENDAmap/cofactor_analysis")
roche = list(csv.DictReader(open(base / "roche_boxed_compounds.csv")))
roche = [r for r in roche if r["compound"].strip()]
brenda = {int(r["brenda_ligand_id"]): r for r in csv.DictReader(open(base / "primary_metabolites_in_multiple_subsystems.csv"))}

# BRENDA reference names per ligand ID (from the map)
d = json.loads((base.parent / "BRENDA_map_reactions.json").read_text())
bnames = {}
for m in d["metabolite_nodes"]:
    if m.get("brenda_ligand_id") in brenda and m.get("brenda_name"):
        bnames.setdefault(m["brenda_ligand_id"], Counter())[m["brenda_name"].strip()] += 1

ROCHE_TO_BRENDA = {
    "L-GLUTAMATE": [41], "2-OXOGLUTARATE": [34], "PYRUVATE": [31], "ACETYL-CoA": [29],
    "ACETATE": [47], "FORMATE": [120], "PRPP": [158], "GLYCINE": [72], "L-GLUTAMINE": [102],
    "SUCCINATE": [58], "TETRAHYDROFOLATE (THF)": [207], "L-ALANINE": [103], "FUMARATE": [170],
    "SUCCINYL-CoA": [224], "D-GLUCOSE": [35, 285], "GLUTATHIONE, REDUCED (GSH)": [44],
    "GLYOXYLATE": [101], "L-ASPARTATE": [97], "OXALOACETATE": [57], "PHOSPHOENOLPYRUVATE": [51],
    "10-FORMYL-THF": [536], "PROPIONYL-CoA": [782], "UREA": [116], "2'-DEOXYRIBOSE 1-P": [33002],
    "5,10-METHYLENE-THF": [401], "5-FORMIMINO-THF": [12390], "ACETALDEHYDE": [90],
    "D-GLYCERALDEHYDE 3-P": [147], "GLUTATHIONE, OXIDIZED (GSSG)": [973], "GLYCEROL": [135],
    "L-CYSTEINE": [74], "L-HOMOCYSTEINE": [305], "L-SERINE": [95], "PHYTYL-PP": [3092],
    "(R)-PANTOTHENATE": [2630], "(S)-METHYLMALONYL-CoA": [2455], "2-OXOBUTYRATE": [390],
    "5-METHYL-THF": [1373], "ACETOACETATE": [419], "ACETYL-P": [358],
    "alpha-D-GLUCOSE 1-P": [107, 5970], "alpha-D-GLUCOSE 6-P": [2944, 105], "BETAINE": [1974],
    "CROTONYL-CoA": [406, 3806], "D-FRUCTOSE": [117], "D-FRUCTOSE 6-P": [87],
    "D-GLUCURONATE": [339], "D-LACTATE": [600], "D-RIBOSE 5-P": [183], "DIMETHYLGLYCINE": [1240],
    "FORMALDEHYDE (HCHO)": [130], "GLYCOLATE": [479], "L-CITRULLINE": [938], "L-THREONINE": [250],
    "MALONATE SEMIALDEHYDE": [2969], "MALONYL-CoA": [76], "PALMITOYL-CoA": [128],
    "PHENYLACETATE": [776], "SARCOSINE": [594], "TARTRONATE SEMIALDEHYDE": [7099],
    "TAURINE": [595], "UBIQUINONE": [410], "UDP-D-GLUCOSE": [64], "(R)-METHYLMALONYL-CoA": [21229],
    "(S)-MALATE": [200, 247], "1,3-P2-D-GLYCERATE": [2125], "3-P-D-GLYCERATE": [201, 194],
    "5,10-METHENYL-THF": [1085], "ADENOSINE 5'-DIPHOSPHATE (ADP)": [13], "ADP-D-GLUCOSE": [8036],
    "beta-D-FRUCTOSE 1,6-P2": [6523], "beta-D-GLUCOSE": [96], "CDP-1,2-DIACYLGLYCEROL": [1248],
    "CHOLINE": [149], "COENZYME A (CoA-SH)": [18], "D-ALANINE": [164], "D-ERYTHROSE 4-P": [500],
    "D-RIBULOSE 5-P": [443], "D-SEDOHEPTULOSE 1,7-P2": [20229], "D-XYLULOSE 5-P": [644],
    "DIACYLGLYCEROL": [2063], "GLYCERONE-P (DIHYDROXYACETONE-P)": [278, 482],
    "GTP (GUANOSINE 5'-TRIPHOSPHATE)": [37], "GUANOSINE 5'-DIPHOSPHATE (GDP)": [53],
    "L-VALINE": [246], "NAD": [7], "NADP": [10], "NICOTINAMIDE": [267], "NICOTINATE": [1029],
    "PALMITATE": [517],
}
names_roche = {r["compound"] for r in roche}
assert set(ROCHE_TO_BRENDA) <= names_roche, set(ROCHE_TO_BRENDA) - names_roche
used = [i for ids in ROCHE_TO_BRENDA.values() for i in ids]
assert len(used) == len(set(used)) and set(used) <= set(brenda), (set(used) - set(brenda))


def brenda_part(ids):
    subs = sorted({s for i in ids for s in brenda[i]["subsystems"].split("; ")})
    return dict(
        brenda_ligand_ids="; ".join(map(str, ids)),
        brenda_labels="; ".join(brenda[i]["metabolite"] for i in ids),
        brenda_names="; ".join(bnames[i].most_common(1)[0][0] for i in ids),
        n_subsystems_brenda=len(subs),
        has_cross_subsystem_link=any(brenda[i]["has_cross_subsystem_link"] == "True" for i in ids),
        brenda_subsystems="; ".join(subs))


EMPTY_B = dict(brenda_ligand_ids="", brenda_labels="", brenda_names="", n_subsystems_brenda=0,
               has_cross_subsystem_link="", brenda_subsystems="")
rows = []
for r in roche:
    ids = ROCHE_TO_BRENDA.get(r["compound"], [])
    b = brenda_part(ids) if ids else EMPTY_B
    rows.append(dict(compound=r["compound"], source="both" if ids else "Roche only",
                     n_boxes_roche=int(r["n_boxes"]), **b, grid_squares=r["grid_squares"],
                     as_structure=r["as_structure"], as_name_only=r["as_name_only"]))
for i in brenda:
    if i in used:
        continue
    b = brenda_part([i])
    rows.append(dict(compound=b["brenda_labels"], source="BRENDA only", n_boxes_roche=0, **b,
                     grid_squares="", as_structure="", as_name_only=""))
for x in rows:
    x["n_total"] = x["n_boxes_roche"] + x["n_subsystems_brenda"]
rows.sort(key=lambda x: (-x["n_total"], -x["n_boxes_roche"], x["compound"].lower()))

cols = ["compound", "source", "n_boxes_roche", "n_subsystems_brenda", "n_total", "brenda_ligand_ids",
        "brenda_labels", "brenda_names", "has_cross_subsystem_link", "grid_squares", "as_structure",
        "as_name_only", "brenda_subsystems"]
with open(base / "roche_brenda_merged.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows({c: x[c] for c in cols} for x in rows)
c = Counter(x["source"] for x in rows)
print(len(rows), dict(c))
