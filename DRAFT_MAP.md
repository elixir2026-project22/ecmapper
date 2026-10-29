# Draft metabolic maps from the BRENDA map

`ecmapper.draft_map` turns the matches of a model to the BRENDA map (see [README.md](README.md#finding-model-reactions-on-a-brenda-map))
into a draft map of the model: the subsystems of the BRENDA map that contain model reactions, reduced to the enzyme
nodes with model reactions, completed with model reactions that connect drawn metabolites, condensed and arranged with
few crossings of the inter-pathway links. The subsystem boxes and names are kept. Outputs: an SVG in the style of the
BRENDA map, an Escher map for editing, and box positions for CNApy.

The results for iML1515 are in [BRENDAmap/draft_map_iML1515](BRENDAmap/draft_map_iML1515).

![Draft map of iML1515](BRENDAmap/draft_map_iML1515/figures/draft_map.png)

## Usage

```bash
python -m pip install -e "./ecmapper[draft-map]"   # numpy, scipy, numba
cd ecmapper
# 1. match the model to the map (BiGG's SBML files have no subsystems: take them from the JSON file)
python -m ecmapper.draft_map.prepare iML1515.xml BRENDAmap/BRENDA_map_reactions.json iML1515_prep.pkl \
    --subsystems iML1515.json
# 2. build the draft map
python -m ecmapper.draft_map.run iML1515_prep.pkl BRENDAmap/BRENDA_map_reactions.json BRENDAmap/BRENDA_map.svg OUT \
    --hub-subsystems 6 --iterations 200000 --seed 3
```

`iML1515.xml` and `iML1515.json` are from [BiGG](http://bigg.ucsd.edu/models/iML1515). The first step takes about
20 s with the resource files cached, the second two to three minutes. Options of `run`:

| Option | Default | |
|---|---|---|
| `--min-anchored` | 2 | drawn main metabolites a model reaction needs to be added (step 3) |
| `--no-completion` | | skip step 3 |
| `--compactness` | 0.5 | weight of compactness against keeping the drawing when condensing (step 4) |
| `--hub-subsystems N` | 0 (off) | remove the inter-pathway links of metabolites drawn in at least N subsystems (step 5) |
| `--iterations`, `--seed`, `--seeds` | 20000, 1, 1 | annealing of the arrangement (step 6); with `--seeds k` the best of k runs is kept |
| `--gap` | 80 | space between subsystem boxes |
| `--flips {none,x,xy}` | none | allow mirroring subsystems (left–right, also upside down) |
| `--escher-scale` | 2.5 | scale of the Escher map relative to the BRENDA map |
| `--escher-no-labels` | | no subsystem names in the Escher map |
| `--compare-flips`, `--extra graphviz elk user-order` | | further arrangements for comparison (`map_*.svg`, `stats.json`); Graphviz needs pygraphviz, ELK needs `node` with `elkjs` on `$NODE_PATH` |

Outputs in `OUT`:

| File | |
|---|---|
| `draft_map.svg` | the draft map in the style of the BRENDA map; reaction IDs below the enzyme nodes, model reactions in the tooltips |
| `draft_map.escher.json` | Escher map (see [Escher export](#5-escher-export)) |
| `draft_map.maps` | CNApy box positions `{reaction ID: [x, y]}` for `draft_map.svg` as background |
| `reaction_positions.json` | positions of the enzyme nodes of each drawn model reaction |
| `added_reactions.json` | the model reactions added in step 3, with subsystem, position and drawn metabolites |
| `step1_subsystems.svg`, `step2_pruned.svg`, `step3_completed.svg` | the map after steps 1–3 (original arrangement) |
| `map_*.svg` | the arrangements compared |
| `stats.json` | counts after each step, hub metabolites, metrics of the arrangements, coverage of the model |
| `draft_map.pkl` | the final `DraftMap` |

Enzyme nodes: green: model reaction on the map; cyan: found in the vicinity search; violet, dashed: added from the model
(step 3). Metabolite nodes added in step 3 are dashed, hub metabolites (step 5) have a red border. Inter-pathway links
are dashed grey lines (purple: added in step 3).

`scripts/svg_to_png.js` renders an SVG (or a part of it) as PNG, `scripts/render_escher.js` loads an Escher map with a
model in Escher, saves a screenshot and reports reactions and metabolites missing in the model (both need playwright).

## 1. Method

| Step | What | How |
|---|---|---|
| 0 | Match the model to the map | `prepare`: `MapMatcher` (vicinity search `propagate`); stores the matches, the keys of model and map metabolites, the reactions with stoichiometry, direction, genes and subsystem |
| 1 | Drop subsystems without model reactions | |
| 2 | Remove enzyme nodes without model reactions | also the metabolite nodes that only they use; inter-pathway links through removed metabolite nodes are bridged |
| 3 | Complete the network with model reactions that are not on the map | see 1.1 |
| 4 | Condense each subsystem | order-preserving LP compaction, see 1.2; box and title fitted around the drawing |
| 5 | Disconnect hub metabolites | the inter-pathway links of metabolites drawn in ≥ N subsystems are removed; their nodes stay and are marked |
| 6 | Arrange the subsystems with few crossings of inter-pathway links | simulated annealing over sequence pairs, see 1.3 |
| 7 | Output | SVG, Escher map, CNApy box positions, reaction positions |

The subsystems are arranged last, although the original idea was to arrange them right after dropping the empty
ones: the arrangement needs the final box sizes and the positions of the linked metabolites in the boxes, and both
change with pruning, completion and condensation. Arranging first and condensing afterwards (`--extra user-order`)
ends with 1.6× the crossings (276 vs. 173 for iML1515) and a 23% larger area, because the boxes shrink after the
arrangement was optimised.

Inter-pathway links considered: the BRENDA `inter-pw-link`s (same compound in two subsystems, hidden in the BRENDA
SVG) whose both ends are main metabolites of enzyme nodes with model reactions, plus the links added in step 3.

### 1.1 Completion with model reactions (step 3)

Candidates: model reactions not on the map, without boundary reactions and transporters (as ecmapper defines them),
whose main metabolites (no currency metabolites, carriers such as ACP, thioredoxin, quinones, protons) are in one
compartment. A drawn metabolite node "is" a model metabolite if it is a main metabolite of an enzyme node whose model
reaction contains the metabolite with a common key (KEGG compound, normalised name), as in ecmapper's matching.

- A candidate is added if ≥ 2 (`--min-anchored`) of its main metabolites are drawn, at least one of them by a node
  of the map (not only by a node added for an earlier candidate). Repeated until nothing changes. So reactions close
  gaps between drawn metabolites or extend the map by one step; chains far away from the map are not built (without
  this condition, 205 reactions, mostly lipid chain variants, piled up in purine metabolism for iML1515).
- Subsystem: the one with most of the reaction's drawn main metabolites; ties: the subsystem of the reaction's best
  `related` node (same EC number or reaction ID, other metabolites), the BRENDA subsystem in which most model
  reactions of the same model subsystem are drawn, the nearest subsystem.
- Position: the removed `related` node of the map if there is one in that subsystem, otherwise between the drawn
  metabolites (free place by spiral search). Main metabolites drawn only in other subsystems get a new node, linked
  to the nearest drawn node (an inter-pathway link); metabolites not on the map get a new node; cofactors get
  cofactor nodes. Arrows follow the model's reversibility.
- Variants (same drawn nodes, reaction or new metabolite IDs equal apart from digits, e.g. chain lengths
  `PLIPA2E160pp`, `PLIPA2E161pp`) share one node labelled `PLIPA2E160pp (+1)`.

### 1.2 Condensing a subsystem (step 4)

Alternating x/y 1D compaction solved as an LP (scipy, HiGHS): nodes whose extents overlap in the other axis keep their
order with a minimum gap; nodes closer than 150 units in the other axis keep their order; objective
Σ|pᵢ − s·pᵢ⁰| + 0.5·n·extent, i.e. a uniform shrink by a factor s plus individual moves that close larger holes.
This is the idea of VPSC/PRISM overlap removal (keep the relative layout, remove space). For iML1515, free
compaction (no order band) gave −52% area but +20% crossings of main-metabolite links inside the subsystems; full order
preservation −29% area, +3% crossings; the band of 150 −44% area, +8% crossings.

### 1.3 Arrangement of the subsystems (step 6)

The subsystems are fixed boxes with "ports" (the linked metabolite nodes inside them). An arrangement is encoded as a
sequence pair (Murata et al. 1995, VLSI floorplanning), which always gives an overlap-free, packed placement.
Simulated annealing (numba kernels, 1.5–2.5 minutes for 200,000 moves with 115 subsystems) minimises

`crossings + 0.5·link–box crossings + 0.02·link length/mean box size + 100·area/Σ box area + 200·log(aspect/1.6)²`

with moves: swap two boxes in one sequence, exchange the places of two boxes, optionally mirror a box. It starts from
the original BRENDA arrangement.

## 2. Layout tools

| Tool | Fit for this task | Tried |
|---|---|---|
| **Sequence-pair annealing** (here, Python and numba) | Optimises exactly the wanted objective (crossings of links between fixed drawings with fixed ports, link–box crossings, area); keeps the drawings of the subsystems | yes – used |
| **scipy / HiGHS LP** | Order-preserving compaction (VPSC-like) of nodes and boxes | yes – used |
| Graphviz `neato`/`fdp`/`sfdp` (+ VPSC overlap removal) | Places boxes by stress/forces on the box graph; ignores ports, no crossing objective. `prism` overlap removal needs GTS (not in Ubuntu's build) | yes – 2.3× (neato) to 5.5× (fdp) more crossings |
| ELK (elkjs) `layered` with `FIXED_POS` ports | Sugiyama with real crossing minimisation and fixed port positions, but forces a layered flow; large area. Its orthogonal/polyline edge routing (and libavoid in ELK) is useful for routing inter-pathway links | yes – 4.4–6× more crossings |
| ELK `stress`, `rectpacking` | stress ignores sizes and ports; rectpacking packs boxes without edges | stress – poor |
| OGDF (C++, `ogdf-python`) | Strongest library for crossing minimisation (planarisation, cluster planarity, orthogonal layouts); does not keep given drawings, could arrange the box graph or draw small subsystems from scratch. Needs an OGDF build | no |
| Tulip (used by CBMViz in this repository) | Bundles OGDF layouts and graph hierarchies, scriptable in Python | no |
| adaptagrams (libcola, libvpsc, libavoid) / WebCola | Constraint-based layout with alignment/separation constraints, VPSC overlap removal and libavoid object-avoiding connector routing: best candidate for routing inter-pathway links around boxes and for alignment constraints inside subsystems | no |
| Cytoscape.js `fcose` / Cytoscape (yFiles layouts) | Compound nodes (subsystems) with constraints; yFiles (commercial, free in yEd) has the best grouped hierarchical/orthogonal layouts | no |
| igraph / networkx | Graph metrics, no boxes or ports | no |
| Metabopolis (Wu et al. 2019) | Published method closest to this scheme: pathways as blocks, inter-pathway links minimised, orthogonal routing in channels ("urban map") | no implementation available |
| Escher | Editing target (no automatic layout) | export |

## 3. Results on iML1515

Input: iML1515 from BiGG, matched with vicinity search `propagate`: 714 of 2,375 non-boundary reactions are on 815 of
2,848 enzyme nodes in 115 of 193 subsystems.

### 3.1 Map after each step

| | subsystems | enzyme nodes | main metabolite nodes | cofactor nodes | model reactions | inter-pathway links |
|---|---|---|---|---|---|---|
| BRENDA map | 193 | 2,848 | 2,818 | 6,257 | – | 234 |
| step 1 | 115 | 1,943 | 1,876 | 4,138 | 714 | 161 |
| step 2 | 115 | 815 | 968 | 1,629 | 714 | 117 |
| step 3 | 115 | 1,087 | 1,221 | 1,995 | 1,010 | 212 |
| step 5 (hubs ≥ 6) | 115 | 1,087 | 1,221 | 1,995 | 1,010 | 153 |

Step 3 added 296 of 790 candidate reactions on 272 enzyme nodes (140 at the place of a `related` node, 156 between
their metabolites), 158 nodes of metabolites not on the map and 95 linked nodes of metabolites drawn in other
subsystems. 216 of them connect 2 metabolites that were already drawn, 76 connect 3–5.

Step 4 (sum over the subsystems):

| | before | after |
|---|---|---|
| drawing area | 116.0 M units² | 65.0 M (−44%) |
| overlapping nodes | 7 (from added nodes) | 0 |
| crossings of main-metabolite links inside subsystems | 986 | 1,064 (+8%) |
| crossings of all links (incl. cofactors) | 1,426 | 1,657 (+16%) |

Coverage: 2,712 reactions = 337 boundary + 768 transporters + 1,607 others; the draft map shows 1,010 of the 1,607
(63%; 714 from the map, 296 added). By model subsystem e.g. citric acid cycle 14/14, pentose phosphate pathway 12/12,
glycolysis/gluconeogenesis 20/23, histidine 10/10, cysteine 12/14; glycerophospholipid metabolism 74/246 (chain-length
variants), murein recycling 8/33, murein biosynthesis 0/15, "Transport, Inner Membrane" 1/23.

### 3.2 Hub metabolites

Metabolites drawn in ≥ 6 of the remaining subsystems (subsystems / inter-pathway links): pyruvate 26/8, acetyl-CoA
20/9, GAP 12/4, F6P 8/6, L-serine 8/6, oxaloacetate 8/3, DHAP 8/2, chorismate 7/5, L-glutamate 7/2, acetaldehyde 7/2,
GTP 6/4, fumarate 6/3, succinyl-CoA 6/3, acetate 6/2 — 59 of 212 links.

### 3.3 Arrangements

Inter-pathway link crossings / link–box crossings / area ÷ Σ box area; seed 1, 200,000 moves
([arrangement_experiments.json](BRENDAmap/draft_map_iML1515/arrangement_experiments.json); its `fin_hub*` entries are
the runs with 4 seeds, `x`: with left–right mirroring):

| hubs disconnected | links | original order, compacted | **annealed** | annealed, mirrored (x, y) | arrange, then condense | Graphviz neato | Graphviz fdp | ELK layered | ELK stress |
|---|---|---|---|---|---|---|---|---|---|
| none | 212 | 1321 / 313 / 2.77 | **449 / 184 / 2.20** | 360 / 166 / 2.10 | 706 / 257 / 3.13 | 1863 / 595 / 3.89 | 1625 / 1018 / 2.69 | 2084 / 586 / 3.06 | 1827 / 790 / 3.90 |
| ≥ 8 subsystems (7) | 174 | 984 / 264 / 2.77 | **223 / 117 / 1.89** | 148 / 90 / 1.95 | 360 / 173 / 2.57 | 908 / 530 / 3.24 | 1114 / 779 / 2.84 | 983 / 282 / 3.85 | 1190 / 673 / 4.35 |
| ≥ 6 subsystems (14) | 153 | 880 / 236 / 2.77 | **173 / 82 / 2.07** | 126 / 91 / 2.16 | 276 / 145 / 2.54 | 813 / 422 / 3.10 | 889 / 634 / 2.76 | 1040 / 270 / 2.56 | 806 / 620 / 3.48 |
| ≥ 5 subsystems (21) | 139 | 752 / 216 / 2.77 | **102 / 69 / 1.96** | 90 / 68 / 1.88 | 223 / 114 / 2.35 | 496 / 320 / 3.35 | 780 / 706 / 2.52 | 529 / 193 / 3.37 | 536 / 322 / 4.20 |

The Graphviz columns use Graphviz's default edge lengths; `--extra graphviz` sets them from the box sizes, which gives 406 / 145 / 2.32 (neato) and 950 / 701 / 2.65 (fdp) for hubs ≥ 6. Over 4 seeds (hubs ≥ 6): 148–173
crossings without mirroring, 121–138 with left–right mirroring.

Final map: hubs ≥ 6, no mirroring, seed 3 (lowest cost of 4 seeds): **149 crossings (−83% vs. 880), 87 link–box crossings (−63%),
area 1.89 × Σ box areas (original order: 2.77), 15,900 × 9,350 units**.

| Original arrangement, compacted | Draft map, central metabolism |
|---|---|
| ![original order](BRENDAmap/draft_map_iML1515/figures/map_original_order_compacted.png) | ![central metabolism](BRENDAmap/draft_map_iML1515/figures/draft_map_detail_central.png) |

## 4. Evaluation

What works
- The arrangement is compact and keeps linked subsystems together (pentose phosphate pathway, glycolysis,
  Entner–Doudoroff, gluconeogenesis and TCA cycle are next to each other). Annealing over sequence pairs beats the
  general tools by far because they optimise other things (stress, layers) and cannot use the fixed port positions.
- Disconnecting hubs is the most effective single measure (crossings 449 → 173 with the same annealing) and loses
  little: their nodes stay in every subsystem, marked.
- The subsystems remain recognisable BRENDA drawings, without the unused 58% of enzyme nodes; the arrows show the
  directions of the model.
- Completion adds 296 reactions, mostly where they belong (phospholipid/LPS reactions in lipid A, CDP-diacylglycerol
  and fatty acid subsystems; cofactor reactions; nucleotide salvage in purine/pyrimidine metabolism; e.g. `PYK2 (+2)`,
  pyruvate kinase with other nucleotides, in glycolysis).

Problems
- Added reactions are placed by drawn metabolites, not by biology: phospholipid turnover lands in "lipid A
  biosynthesis" because palmitate/palmitoleate are drawn there; 17 reactions of "Inorganic Ion Transport and
  Metabolism" (Fe³⁺-siderophore reductases, FMN-dependent monooxygenases) in vitamin B12 metabolism; the Fe³⁺-dicitrate
  reactions `FE3DCITexs` (extracellular) and `FE3DCITR5` in the citric acid cycle (citrate is drawn there).
- 58 periplasmic reactions are drawn like cytosolic ones (the BRENDA map has no compartments).
- Added reactions can create long edges inside a subsystem when their drawn metabolites are far apart; condensation
  adds 8% (main) / 16% (all) crossings inside subsystems.
- Generic map metabolites match many model variants (via generic KEGG compounds); the grouping of variants hides
  that, the label shows `(+k)`.
- 87 inter-pathway links still pass through other boxes (straight lines).
- Reactions on several nodes are drawn several times (BRENDA draws some steps in several subsystems); CNApy shows one
  box per reaction (the first node is used).
- Not drawn: 597 non-transport reactions (e.g. many glycerophospholipid chain-length variants, murein biosynthesis,
  where no metabolite of the map matches) and all 768 transporters.

## 5. Escher export

`draft_map.escher.json` (Escher schema 1-0-0, coordinates × 2.5) leaves out the inter-pathway links and the subsystem
boxes; the subsystem names are text labels (`--escher-no-labels` drops them). Each model reaction of an enzyme node is
an Escher reaction with its midmarker at the node (further reactions of the node below it) and multimarkers towards
its substrates and products. Drawn main metabolites are primary nodes, matched cofactors secondary nodes; metabolites
that are not drawn get secondary nodes next to the reaction; protons are not drawn (they stay in the stoichiometry).
There is one Escher node per drawn node and compartment-specific model metabolite.

iML1515: 1,183 Escher reactions (1,010 model reactions; reactions on several nodes are drawn several times),
6,904 nodes, 115 text labels. `escher.validate.validate_map` (escher 1.7.3; schema and consistency of segments,
stoichiometry, genes) passes; Escher 1.7.3 loads it with `iML1515.json` and finds all reactions and metabolites in
the model.

![Escher, central metabolism](BRENDAmap/draft_map_iML1515/figures/escher_detail_central.png)

## 6. Next steps

1. Route inter-pathway links around boxes (libavoid or ELK orthogonal routing) and draw them in the Escher map as
   needed.
2. Completion: compartment handling (restrict to the cytosol or add a "periplasm" subsystem), a target-subsystem score
   that also uses the model subsystem, optional new subsystems for groups of added reactions without a good BRENDA
   subsystem (e.g. murein, membrane lipids).
3. Local refinement inside subsystems after completion (stress with constraints, e.g. WebCola) to shorten the edges of
   added reactions.
4. Weights of the annealing as options; keeping the colour groups (category regions) of the BRENDA map together.
5. A command of the `ecmapper` CLI; runs on the models of [BRENDA_MAP_EVALUATION.md](BRENDA_MAP_EVALUATION.md).
