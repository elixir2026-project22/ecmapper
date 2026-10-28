# BRENDA map – analysis summary

Source: `BRENDA_map.svg` (25,378 × 16,085 px). Companion files: `BRENDA_map_reactions.json` (all reactions with layout) and `BRENDA_map_EC_duplicates.xlsx` (detailed tables).

## 1. How the map is encoded

| Element | SVG representation | Count |
|---|---|---|
| Subsystems | `polygon#pw_<name>` plus a `nodePathwayTitle` label node | 193 |
| Reactions (enzymes) | `g.nodeEN`, EC number(s) as text, classes `pw_<subsystem>`, `ec…`, `rID…` | 2,848 |
| Main metabolites | `g.nodeHM` (main) and `g.nodeSM` (small/side) | 2,623 + 195 |
| Cofactors | `g.nodeCM` | 6,257 |
| Links | `path.link` (straight lines, a few 3-point paths) | 12,617 |

Every node belongs to exactly one subsystem through its `pw_` class. Metabolites are drawn once per reaction neighbourhood rather than once per map, so the same compound appears on many nodes.

Links carry no references to the nodes they connect, only coordinates. I assigned each endpoint to the node rectangle it touches; every endpoint fell within 20 map units of a node, and every enzyme node received at least one link.

Arrow semantics:

- **Enzyme ↔ main metabolite** links are double-headed (5,952 of 5,966), so the SVG does not say which main metabolite is the substrate and which the product.
- **Enzyme ↔ cofactor** links have a single arrow, so cofactors can be split into consumed and produced.
- 235 metabolite–metabolite links connect the same compound across subsystems (`inter-pw-link`).

## 2. Reaction export (`BRENDA_map_reactions.json`)

The top-level keys are:

- `notes`: method and coordinate conventions.
- `counts`: totals for each section.
- `subsystems`: polygon, bounding box, title node and member reactions.
- `reactions`: one entry per enzyme node, listed below.
- `metabolite_nodes`: every metabolite node with layout and the reactions that use it.
- `metabolite_links`: the metabolite–metabolite connections.
- `unassigned_links`: one link not attached to an enzyme.

Each reaction entry contains:

- **Identity:** `node_id`, `subsystem`, `ec_label`, `ec_numbers`, `enzyme_name`, `brenda_reaction_ids`.
- **Layout:** `layout` gives the centre, width and height.
- **Main metabolites:** `main_metabolites` lists each metabolite's display name, BRENDA name, ligand ID, node class, link geometry and arrow direction.
- **Cofactors:** `cofactors.consumed`, `cofactors.produced` and `cofactors.undirected`. The undirected list is empty for this map.
- **Duplicate check:** `duplicate_check` gives, for repeated ECs, the relation to the other nodes with the same EC.

All coordinates are in the map frame, inside the group with `translate(30,20) scale(0.08)`.

## 3. EC numbers

The map has 2,077 distinct EC labels. Of these, 2,040 are specific: full ECs or BRENDA provisional ECs of the form `x.x.x.Bn`. The rest are incomplete classes (`x.x.x.-`) or "spontaneous". A further 52 nodes are labelled `?.?.?.?`. Nodes with combined labels (e.g. `3.1.7.11/3.1.7.3`) were split so that each EC is counted individually.

### Repeats within a subsystem

There are 307 cases where an EC appears on more than one node in the same subsystem, spread across 97 subsystems. Of these, 283 involve specific ECs. Comparing the metabolites around each node gives this breakdown for the 283:

| Verdict (most similar pair of nodes) | Cases |
|---|---|
| Genuine duplicate – identical main metabolites and cofactors | 8 |
| Same reaction, cofactor direction reversed | 3 |
| Same conversion, different cofactor (NAD⁺/NADP⁺, FAD/quinone, Fd/acceptor …) | 20 |
| Different reactions, sharing a metabolite (consecutive or branching steps) | 124 |
| Different reactions, no shared metabolite | 128 |

**Genuine duplicates:**

| Subsystem | EC | Reaction |
|---|---|---|
| C4 and CAM carbon fixation | 1.1.1.37 | L-malate ⇄ OAA |
| C4 and CAM carbon fixation | 1.1.1.40 | L-malate ⇄ pyruvate |
| aflatoxin biosynthesis | 1.1.1.352 | 5'-hydroxyaverantin ⇄ 5'-oxoaverantin |
| aflatoxin biosynthesis | 1.14.14.116 | averantin ⇄ 5'-hydroxyaverantin |
| alanine metabolism | 1.2.1.3 | 3-APA ⇄ β-Ala |
| cysteine metabolism | 2.6.1.1 | 3-sulfinoalanine ⇄ 3-sulfinylpyruvate |
| methanogenesis from CO₂ | 1.5.98.1 | methylene-H₄MPT ⇄ methenyl-H₄MPT |
| sulfate reduction | 2.7.7.4 | sulfate ⇄ APS |

The two C4/CAM cases are probably intentional, since the same step is drawn in two cell types.

**Reversed cofactor direction.** 2.6.1.1 and 2.6.1.2 in C4/CAM, and 2.6.1.5 in tyrosine metabolism, are each drawn twice: once with 2-oxoglutarate → glutamate and once with glutamate → 2-oxoglutarate.

**Annotation issues worth checking:**

- **1.13.11.79:** the two nodes differ only in one product label, "unknown substance" versus "dialuric acid". This is effectively a duplicate.
- **5.1.3.14:** one node has H₂O → UDP and the other no cofactor. This looks like the hydrolysing and non-hydrolysing epimerases mixed under one EC.
- **1.1.1.363:** the NAD and NADP variants sit on nodes with different combined labels, 1.1.1.363/1.1.1.388 and 1.1.1.49/1.1.1.363.

**Shared reaction IDs.** Shared BRENDA reaction IDs do not imply duplicates. Of the 17 repeats that share an ID, 8 are genuine duplicates. The others (e.g. 1.1.1.100 in palmitate biosynthesis, 4.4.1.13) use one generic reaction ID for different substrates.

### ECs in several subsystems

229 specific ECs occur in two or more subsystems; 161 of them occur in exactly two. The most widespread are:

- 1.1.1.1 in 8 subsystems.
- 1.2.1.3 and 2.6.1.1 in 6 each.
- 1.1.1.35, 2.6.1.42, 4.2.1.17 and 6.2.1.12 in 5 each.

These are mostly promiscuous dehydrogenases, transaminases, β-oxidation enzymes and CoA ligases. Among non-specific entries, "spontaneous" occurs in 28 subsystems and `1.1.1.-` in 10.

## 4. Metabolite identifiers

**The displayed names are not identifiers.** The usable ID is the BRENDA ligand ID behind each node: the `ligand.php?brenda_ligand_id=N` link, which is identical to the `lidN` class. All metabolite nodes except one carry it, giving 2,344 distinct IDs. BKMS-react offers a "BRENDA Ligand ID" search field, so query by the numeric ID rather than by name. I have not verified how many of these IDs BKMS contains.

**Naming is inconsistent:**

- For 693 IDs the displayed label differs from the BRENDA name shown in the tooltip, often because it is an abbreviation (Pi, Glc, PEP, Suc, GAP …).
- 131 IDs appear under more than one label.
- 40 labels map to more than one ID.

**Likely wrong links**, where the label and the linked ID refer to different compounds:

| Label | Linked ID | BRENDA name of that ID |
|---|---|---|
| NADH (2 nodes) | 12 | H⁺ |
| dimethylallyl-PP | 17 | diphosphate |
| N-methylputrescine | 1931 | S-adenosyl-L-methionine |
| kaempferol | 137 | quercetin |
| tropine | 3848 | tropinone |
| atropine | 4789 | atropinone |
| malonyl-CoA | 811 | malonyl-ACP |
| propanol | 328 | propanal |
| [enzyme]-cysteine | 1410 | ethylene glycol |
| adenosylcobinamide phosphate | 22961 | adenosylcobalamin 5'-phosphate |

This list is heuristic and not exhaustive. The full table is in the spreadsheet's *Likely wrong IDs* sheet.

## 5. Other data-quality notes

- One metabolite node, (2E,8Z,11Z,14Z)-icosatetraenoyl-CoA, has no ligand ID.
- Seven cofactor nodes are not connected to any reaction: four NADH, one NADP⁺, one Val and one oxidized NADPH-hemoprotein reductase.
- One link connects two metabolite nodes without the inter-pathway class (node1999 → node2003).
- 231 enzyme nodes carry no BRENDA reaction ID.
- The curcuminoid biosynthesis title node is tagged `pw_curcumin-biosynthesis`, which doesn't match its polygon `pw_curcuminoid-biosynthesis`.

## 6. Implications for using the map as a template

- **Duplicate check.** Use node ID plus metabolite neighbourhood, not EC number, as the reaction identity. The same EC legitimately represents different steps in most of the repeated cases.
- **Reaction direction.** Main-metabolite direction has to come from an external source, such as BKMS or BRENDA reactions, because the SVG does not encode it.
- **Before matching to BKMS.** Normalise metabolites to their BRENDA ligand IDs and fix the mismatched links listed above.
