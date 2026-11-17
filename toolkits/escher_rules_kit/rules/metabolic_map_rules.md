# Metabolic map layout: typed rule set

Recast 2026-11-16 from `metabolic_map_design_principles.md`. Reference numbers [1]-[22] are unchanged; [23]-[32] were added in review.

## How to read a rule

Each rule is one trigger and one action: **When** the trigger holds, **Do** the action.

| Field | Values |
|---|---|
| Type | `motif`, `arrangement`, `encoding` (ID prefix M, A, E) |
| Strength | `must` (breaking it makes the map wrong or unreadable), `should` (default convention), `may` (option), `open` (undecided) |
| Basis | `stated` (a source says so), `inferred` (your reading of canonical maps), `review` (added in this review, not in the original list) |
| Chk | `ok` confirmed against the source text; `fix` confirmed or biology checked, wording corrected; `-` source not opened in this review |

Sources checked in this review (full text or the relevant passage): [7], [8], [9], [10], [12], [13], [15], [17], [18]; [19] for identity only; [4] and [5] only as described in [7]. Everything citing only [1]-[3], [6], [11], [14], [16], [20]-[22] is carried over unchecked.

## Motif rules

How one reaction or one unit is drawn internally.

### Reaction: what one reaction looks like

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| M01 | Every reaction | Draw exactly one substrate-product pair as the primary connection; attach all other participants as side participants. | must | stated | [17] | ok |
| M02 | Choosing the primary connection | Prefer the pair that conserves the carbon / heavy-atom skeleton. Exclude pairs that run through a carrier (CoA, NAD(P)H, ATP ...) whenever a carrier-free pair exists. PDH: pyruvate to acetyl-CoA, not CoA to acetyl-CoA or NAD+ to NADH. | should | stated | [17] | ok |
| M03 | Deciding primary vs. side | Assign the role per reaction in its pathway, never per metabolite. ATP is a side participant in a kinase and primary in purine synthesis; glutamate is primary in its own synthesis and side in a transamination. | must | stated | [7] | ok |
| M04 | No role annotation available | Use a degree threshold only as a fallback and keep the resulting list editable. Degree alone misclassifies central hubs (pyruvate, glutamate, acetyl-CoA). | should | review | [13], [17], [24] | ok |
| M05 | Side participants | Attach as short curved branches at the reaction; keep the backbone straight and visually continuous. | should | stated | [17] | ok |
| M06 | Paired cofactors (ATP/ADP, NAD+/NADH) | Place both members of the pair on the same side of the backbone. | should | stated | [17] | ok |
| M07 | Two antiparallel reactions between the same metabolites (PFK / FBPase) | Draw both along one axis, not as a two-node cycle. | should | stated | [17] | ok |

- **M02** The original list stated the outcome (PDH) but not the selection rule.
- **M03** BioCyc primary/secondary annotation is per pathway. Removing a compound globally makes pathways that end or start at it meaningless.
- **M04** MetDraw, LMME and CBMViz (SMALL_THRSH) are degree-based. MetaCarto reports that connectivity cloning recovers curated pairs only by keeping 1.75-2.86 candidates per reaction (38-61% precision).
- **M06** Source says 'where possible'.

### Shape: what one unit looks like

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| M10 | Unit topology is known | Linear: straight chain. Cyclic: ring. Branched: tree. Compose hierarchically: contract cycles to super-nodes, lay out the rest in layers, then expand the cycles. | should | stated | [4], [5], [6], [17] | ok |
| M11 | Detecting cycles | Search inside a unit, never network-wide. Network-wide detection returns cycles that span several pathways. | must | review | [7] | ok |
| M12 | A unit contains a cycle, even a small one | Draw it as a ring. | should | stated | [4] | - |
| M13 | Placing a ring | Rotate it so the principal input faces the pathway that feeds it. | should | stated | [17] | ok |
| M14 | Ring perimeter | Entries and exits sit on the perimeter; CO2 and redox cofactors point outward. | should | inferred |  | - |
| M15 | Two cycles share a metabolite or reaction (urea cycle / TCA; CBB / photorespiration) | Draw two adjacent rings joined at the shared element instead of one merged cycle. | should | stated | [4] | ok |
| M16 | Iterative module (beta-oxidation, fatty-acid synthesis / elongation) | Spiral, or stacked repeats of one turn. | may | inferred |  | - |
| M17 | Common trunk ending at a branch-point metabolite | Straight trunk, then a fan. Shikimate to chorismate, then Phe / Tyr / Trp. Asp to Lys / Thr / Met, with Ile continuing from Thr. | should | inferred | [11] | fix |
| M18 | Homologous sequences (shared enzymes, or one pathway present in two compartments) | Parallel limbs. Val and Ile run in parallel; Leu branches off 2-oxoisovalerate. Peroxisomal and cytosolic terpenoid branches run in parallel. | should | inferred | [17] | fix |
| M19 | Bypass of trunk steps (ED vs. EMP, glyoxylate shunt, gluconeogenic steps) | Draw beside the trunk. Gluconeogenic steps are counter-direction arrows next to the glycolytic step they bypass. | should | inferred |  | - |

- **M10** Confirmed for [5] and [17]; [4] and [6] not opened. This is where [5], uncited in the original, belongs.
- **M11** Bourqui et al. show this failure on the mouse network (their Fig. 1).
- **M15** Confirmed via the description of [4] in [7]. Examples are mine.
- **M16** In most GEMs every chain length is its own reaction set, so there is no graph cycle to find. Recognise the module by its repeating reaction classes and stack the turns; a true spiral needs pooled metabolites.
- **M17** Branch points made explicit (chorismate; Ile via Thr). Escher's builder works the same way: chain plus perpendicular branches.
- **M18** Original listed Val/Ile/Leu as three siblings; Leu is a branch off the Val limb. Terpenoid example confirmed in [17].

### Orientation and geometry

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| M20 | Orienting a unit | Follow the net / biosynthetic direction (sign of the pFBA flux), not the stored reaction direction and not the reversibility flags. | should | stated | [7], [17] | ok |
| M21 | Flux is zero or unavailable | Fall back to the template direction, then to the surrounding topology. | should | review | [17] | ok |
| M22 | Flux-guided orientation is used | Fix medium and objective, store them with the map, and freeze the result. Orientation is condition-dependent and otherwise violates A30. | must | review |  | - |
| M23 | Routing inside a unit | Long straight axis-aligned runs; minimise crossings and bends; short edges; wide angles. | should | stated | [7], [17], [19], [32] | ok |

- **M20** [7]: valine biosynthesis input could not be identified from reversibility alone.
- **M21** Topology fallback is MetaCarto's; template-first is my addition.
- **M23** Confirmed in [7] and [17]. Force-directed layouts reach 0-2% axis-aligned edges against 62% for MetaCarto, so they cannot satisfy this rule.

## Arrangement rules

What counts as a unit and where units go relative to each other, to compartments and to earlier versions of the map.

### Decomposition: what counts as a unit

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| A01 | Genome-scale model | Decompose into units before layout, in this order: template match, model subsystems, KEGG pathway annotation, modularity clustering on the currency-free graph. | must | stated | [17] | ok |
| A02 | Unit size | 6-60 reactions per unit; merge units of the same functional class up to 120. Hard ceiling of about 200 reactions per Escher map. | should | stated | [12], [17] | ok |
| A03 | Units share nodes and nodes are not duplicated | Choose a priority set that is drawn intact (default: glycolysis plus TCA; the user may name focus pathways). Other units give up the shared nodes. | should | stated | [7] | ok |
| A04 | A metabolite is shared between units | Either duplicate it at the unit boundary (BioCyc) or keep one instance (MetaViz). Decide once, globally. See Conflicts, C1. | open | stated | [1], [7] | ok |
| A05 | Exchange, sink, demand and biomass reactions | Keep them out of pathway units. Exchanges go to the cell boundary. Biomass is one sink glyph whose precursors are cloned; its edges are never routed across the map. | must | review |  | - |

- **A01** Steps 2-4 are MetaCarto's. Template match first (EC number plus main substrates) is the project's own step. Unit boundaries differ between models.
- **A02** The 200 limit is a viewer-performance limit, not a readability finding.
- **A03** E. coli: 658 of 1,140 nodes (compounds plus reactions) lie in several pathways.
- **A05** Missing from the original. CBMViz already aliases the objective's reagents (split_biomass) and can place transport radially.

### Global template: where units go

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| A10 | Arranging units | Specify adjacency explicitly (template or adjacency list). Local motif rules do not add up to a Roche-like global picture. | must | stated | [7] | ok |
| A11 | Function-class template | Biosynthesis on the left, energy / central metabolism in the middle, degradation on the right. Reactions without a pathway go in a grid at the far right. | should | stated | [1], [2], [3] | - |
| A12 | Default flow | Pathways flow top to bottom. | should | stated | [2] | - |
| A13 | Sibling pathways (same ontology class) | Keep them adjacent inside one shaded block. | should | stated | [1], [2] | - |
| A14 | Bow-tie | Inputs on one side, compact core in the middle, biosynthetic fans on the other side. | should | inferred | [21] | - |
| A15 | Canonical neighbours | Glycolysis above TCA. PPP at G6P. Amino-acid families at their precursor: 2-OG (Glu, Gln, Pro, Arg); OAA (Asp, Asn, Lys, Thr, Met, Ile); 3PG (Ser, Gly, Cys); pyruvate (Ala, Val, Leu); PEP + E4P (Phe, Tyr, Trp); R5P / PRPP (His). Nucleotides at PRPP. Cofactors and secondary metabolism at the periphery. | should | inferred |  | fix |
| A16 | Edges between units | Treat units as rectangular blocks and route connecting edges along the gaps. | may | stated | [16] | - |
| A17 | Choosing which units sit next to each other | Lay out the overview graph (unit = node, interface metabolite = edge). Units linked by many interface metabolites are neighbours. | should | review | [24] | ok |
| A18 | Secondary metabolism | Give it a separate overview map instead of crowding the global map. | may | stated | [8], [9] | ok |

- **A12** Overridden by M20 where net flux runs the other way.
- **A14** Compatible with A11 read right to left. KEGG and Roche do not assign sides by function, so the side must be fixed in the template.
- **A15** His family (PRPP) was missing. Fungal Lys comes from 2-OG (aminoadipate route), so the Lys anchor is organism-dependent.
- **A17** Missing from the original; gives A10 something computable for units that the template does not cover.

### Compartments

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| A20 | Cell regions | Regions follow cell architecture: extracellular space, plasma membrane, cytosol; add outer membrane and periplasm for Gram-negatives. | should | stated | [2] | - |
| A21 | Organelles | Draw as nested enclosures inside the cytosol region. | should | stated | [20] | - |
| A22 | Transport reactions | Place subsystems inside their compartment; transport reactions straddle the boundary they cross. Transporters sit in the membrane. | should | stated | [2], [13] | ok |
| A23 | Uptake and secretion | Place on the outer edge, next to the pathway they feed or drain. | should | inferred |  | - |
| A24 | Respiratory chain and ATP synthase | Draw as a band inside the energy-transducing membrane: plasma membrane (bacteria), inner mitochondrial membrane, thylakoid membrane. | should | inferred |  | fix |
| A25 | Shuttles across a membrane (malate-aspartate, carnitine, photorespiration) | Draw as one unit that straddles the boundary, not as separate transport reactions. | should | review |  | - |

- **A22** Confirmed for [13].
- **A24** Original said 'the membrane'; organelle membranes added.
- **A25** Missing from the original. Photorespiration spans chloroplast, peroxisome and mitochondrion.

### Stability, hierarchy, expand / collapse

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| A30 | Any re-layout (new model version, new condition, expand / collapse) | Do not move objects that already have a position. | must | stated | [14], [16] | - |
| A31 | Organism- or model-specific map | Derive it as a subset of the reference layout. Absent reactions are hidden or greyed; positions do not change. | should | stated | [3], [9], [10], [14], [18] | ok |
| A32 | Model too large for one map | Many small maps (per compartment, subsystem or custom) linked by placemarks, not one giant map. | should | stated | [8], [12], [15] | ok |
| A33 | Collapsing a unit | Replace it by one node that keeps its interface metabolites as ports at fixed positions, so that expanding does not move neighbours. | must | review | [24] | - |

- **A31** Confirmed for [10] (933 genomes) and [18].
- **A32** [15] dates from 2021: about 170 hand-arranged maps; the Human-GEM cytosol map needed 5 parts.
- **A33** Missing from the original and central to the planned expand / collapse feature. E04 is the visual counterpart at global scale.

### Organism and condition

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| A40 | Autotroph | The fixation cycle becomes the centre; hexose-phosphate pathways are drawn in the biosynthetic direction. | should | inferred |  | fix |
| A41 | Organism-specific core variant | Give each its own motif; it may replace EMP as the trunk. ED, phosphoketolase, glyoxylate cycle, rTCA, Wood-Ljungdahl, methanogenesis; also 3-hydroxypropionate bicycle, 3HP/4HB and DC/4HB cycles, RuMP and serine cycles (methylotrophs). | should | inferred |  | fix |
| A42 | Lineage-specific periphery | Attach as peripheral units: mycolic acid / arabinogalactan envelope (Mycobacteriales, i.e. Mycobacterium, Corynebacterium, Rhodococcus, Nocardia); TAG storage (Rhodococcus, Streptomyces, Mycobacterium); fermentation branches (yeast ethanol, E. coli mixed-acid); urea cycle, bile acids, cholesterol (mammals); photorespiration (plants). | may | inferred |  | fix |

- **A40** Original said CBB; true only for CBB autotrophs (plants, cyanobacteria, Knallgas bacteria). See A41 for the others.
- **A41** Last five added in review.
- **A42** Original attributed mycolic acids to Actinobacteria as a whole; Streptomyces and Bifidobacterium have none.

## Encoding rules

How things are shown: glyphs, clones, colour, data, legibility.

### Glyphs and zoom

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| E01 | Drawing a reaction | Hyperedge: reactants converge on the reaction marker, then branch to the products. | should | stated | [13], [17] | ok |
| E02 | Zoom tiers | Global: line and filled circle. Overview: arrow and open circle. Detail: box on arrow and open circle. | should | stated | [9] | ok |
| E03 | Labels while zooming | Reveal names progressively (semantic zoom). | should | stated | [3] | - |
| E04 | Global scale | Contract each unbranched segment to a single line. | should | stated | [8] | fix |
| E05 | Reversibility | Arrowheads at both ends for reversible reactions, at the product end only for irreversible ones. | should | review |  | - |
| E06 | Compartment membership | Show it by enclosure (A20-A22); the ID suffix (_c, _m, _e) alone is a fallback. | should | stated | [11], [12] | - |

- **E01** Confirmed for [17] (Escher: metabolite, multimarker, midmarker, multimarker, metabolite).
- **E04** KEGG: 1,534 net-elements and 934 compounds. One drawn line can stand for several net-elements, so 1,534 is not a line count.
- **E05** Missing from the original. Needed because M20 orients reversible reactions by flux.
- **E06** Escher had drawn membranes as future work in 2015.

### Duplication

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| E10 | Currency metabolite | Redraw it at every reaction; never draw it as a shared hub. | must | stated | [13] | ok |
| E11 | Non-currency metabolite | Draw once per compartment. | should | stated | [13] | ok |
| E12 | A non-currency hub is cloned to untangle a unit | Mark every clone and link clones with a dashed line. | should | stated | [13] | ok |
| E13 | Side participants | Draw smaller than backbone metabolites and make them hideable. | may | stated | [11] | - |
| E14 | Any duplication | Each clone hides true connectivity, so duplication must always be visible to the reader. | must | stated | [7] | ok |

- **E10** MetDraw calls these 'minor' metabolites; see Vocabulary.

### Colour and data

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| E20 | Structural colour | Colour encodes functional category. KEGG uses eleven; five of them: blue carbohydrate, purple energy, teal lipid, red nucleotide, orange amino acid. | should | stated | [9] | ok |
| E21 | Data overlay (flux, expression, abundance) | Encode by colour and width / size on top of the fixed layout. Structural colour gives way: switch it off when data is shown. | should | stated | [9], [11], [22] | ok |
| E22 | Gene-reaction rules | Map AND / OR rules onto the reaction to show complexes and isozymes. | may | stated | [11] | - |

- **E20** Eleven categories confirmed; the individual colour values not re-checked.
- **E21** KEGG's 'no colour' mode confirmed in [9].

### Hygiene

| ID | When | Do | Strength | Basis | Src | Chk |
|---|---|---|---|---|---|---|
| E30 | Node overlap | Nodes overlap only for containment. | must | stated | [19] | - |
| E31 | Edge crosses a node | Avoid; where unavoidable the edge is drawn on top. | should | stated | [19] | - |
| E32 | Node labels | At least part of the label lies inside its node. | may | stated | [19] | - |
| E33 | Label size | Check legibility at the target output size before accepting a map. | should | review | [17] | ok |

- **E30** SBGN PD rule; source identity corrected, content not re-checked.
- **E32** Holds for SBGN glyphs. Escher-style maps place labels beside the node, so do not adopt this unchanged.
- **E33** MetaCarto: whole-model maps scaled to a page give a median label of 1.9 pt.

## Premises (not rules)

Findings that justify rules but prescribe nothing themselves.

| ID | Premise | Src | Chk |
|---|---|---|---|
| P1 | Geometric metrics (crossings, bends, spacing) cannot tell whether a reaction sits on the right backbone. Evaluate backbone choice and geometry separately. *bioRxiv preprint, single author, not peer reviewed.* | [17] | ok |
| P2 | The reference maps (KEGG, Roche, WikiPathways, BRENDA) are drawn by hand, so most rules are tacit. *Confirmed for Roche and KEGG in [7]; the WikiPathways and BRENDA part rests on [16], not opened. Michal's own account of the Roche map [23] is the place to replace 'inferred' with 'stated'.* | [7], [16] | ok |
| P3 | Without node duplication only some pathways can be drawn compactly (E. coli: 658 of 1,140 nodes are shared). | [7] | ok |
| P4 | Genome-scale maps have no shared reference coordinate system; layouts vary by drawer. | [14] | - |
| P5 | Metabolism has a bow-tie structure: giant strong component plus IN and OUT sets once currency metabolites are removed; nested inside modules. *'and compartments' in the original could not be confirmed; Zhao et al. describe modules.* | [21] | - |
| P6 | Generic node-edge tools handle metabolic hypergraphs poorly. | [14] | - |
| P7 | Unit boundaries depend on the model's subsystem annotation. | [17] | ok |

## Conflicts between sources

| ID | Conflict | Rules | Suggested resolution |
|---|---|---|---|
| C1 | Duplicate shared metabolites (BioCyc [1]) or not (MetaViz [7]) | A04, E10-E14 | Duplicate currency always; duplicate other metabolites only at unit boundaries, marked as clones. |
| C2 | Fixed downward flow (A12), flux-guided orientation (M20), stable positions (A30) | A12, M20, A30 | Template direction first; flux only for units the template does not cover; then freeze (M22). |
| C3 | Currency by degree (MetDraw, LMME, CBMViz) or by role (BioCyc, MetaCarto) | M03, M04 | Role where annotated or computable; degree as editable fallback. |
| C4 | Sides assigned by function (BioCyc) or not (KEGG, Roche) | A11, A14 | Yours to fix in the template; the sources do not agree. |
| C5 | Unit size 60 / 120 (MetaCarto) or 200 (Escher) | A02 | Different quantities: readability target vs. performance ceiling. Keep both. |

## Vocabulary

One term per concept; the sources use the others.

| Term used here | Meaning | Called elsewhere |
|---|---|---|
| primary connection | the one substrate-product pair that carries the backbone | BioCyc 'primary compounds'; Escher 'primary nodes'; MetaCarto 'primary connection'; KGML substrate / product |
| side participant | any other reactant or product of that reaction | BioCyc 'secondary'; Escher 'secondary metabolite'; MetaCarto 'side participant' |
| currency metabolite | a metabolite that is a side participant in most of its reactions | 'ubiquitous', 'cofactor', 'hub'; MetaCarto 'carrier'; CBMViz 'small species'; MetDraw 'minor metabolite' |
| clone | a repeated drawing of one metabolite | CBMViz / CellDesigner 'alias'; 'duplicate'; SBGN 'clone marker' |
| unit | a set of reactions laid out together by one motif | 'pathway'; model 'subsystem'; KEGG 'module'; MetaCarto 'cluster'; MetaViz / Tulip 'metanode'; Metabopolis 'block' |
| interface metabolite (port) | a metabolite through which a unit connects to another unit | LMME 'interface metabolite'; CBMViz 'boundary species' |
| template | a hand-made reference layout onto which a model is projected | KEGG 'reference map'; ReconMap 'generic map' |

MetDraw's 'major / minor' is a trap: 'minor' means highly connected (currency), and 'major' means drawn once. Neither maps onto primary / side.

## Not yet covered

- Pooled and generic species (lipid classes with R-groups, polymers, acyl-ACP pools): how they are matched to templates and drawn.
- Template recognition itself: minimum fraction of a template that must match, what to do with partial matches and gap reactions, how extra model reactions are attached to a matched template.
- Label placement outside nodes (side, collision handling) for Escher-style maps.
- Stoichiometry, enzyme / EC annotation on the reaction glyph, and off-map connectors between units.
- Exchange formats that constrain what can be encoded: Escher JSON, SBML Layout / Render, SBGN-ML, KGML.
- Requirements work the list does not draw on: Saraiya et al. [30], Siebenhaller et al. [26], and constraint-based layout from an existing KEGG layout [25, 27].
- Automatic or template-based layouts not yet compared: Lambert et al. [28] (successor to MetaViz), MOST-visualization [29], BiGGr template maps [31].

## Changes against the original list

- [19] is Schreiber & Czauderna (2022), not a 'Freivalds/Schreiber-type review'.
- [5] Becker & Rojas was in the reference list but cited nowhere; now cited at M10.
- KEGG global map: 1,534 net-elements, not 1,534 lines (E04).
- Escher's 200-reaction figure is a performance limit (A02).
- 'Force-directed layouts ...; components packed largest-first' merged two unrelated MetaCarto statements. The first is now a note on M23; the second is MetaCarto's own packing step and was dropped as a rule.
- Val / Ile / Leu, Asp family, shikimate trunk: branch points corrected (M17, M18).
- His family added to the neighbour list (A15).
- Mycolic acids narrowed from Actinobacteria to Mycobacteriales (A42).
- 'CBB ring becomes central in autotrophs' narrowed to CBB autotrophs; five further fixation / C1 routes added (A40, A41).
- Bow-tie 'inside modules and compartments' reduced to modules pending a source (P5).
- [21] URL was en.wikipedia.com; primary references given instead.
- [12] cites an RCSB PDB news page for BiGG compartment suffixes; BiGG's own documentation is the appropriate source (not replaced here, not checked).
- MetaCarto details added that the original left out: carrier constraint (M02), topology fallback (M21), merge-to-120 applies within one KEGG BRITE class (A02).
- Several original bullets are premises or trade-offs, not rules; they are kept under Premises and Conflicts.

## References

1. Paley SM, Karp PD (2006). The Pathway Tools cellular overview diagram and Omics Viewer. Nucleic Acids Res 34:3771-3778. https://pmc.ncbi.nlm.nih.gov/articles/PMC1557788
2. Paley S, Billington R, Herson J, Krummenacker M, Karp PD (2021). Pathway Tools visualization of organism-scale metabolic networks. Metabolites 11(2):64. https://pmc.ncbi.nlm.nih.gov/articles/PMC7911265
3. Latendresse M, Karp PD (2011). Web-based metabolic network visualization with a zooming user interface. BMC Bioinformatics 12:176. https://link.springer.com/doi/10.1186/1471-2105-12-176
4. Wegner K, Kummer U (2005). A new dynamical layout algorithm for complex biochemical reaction networks. BMC Bioinformatics 6:212. https://link.springer.com/doi/10.1186/1471-2105-6-212
5. Becker MY, Rojas I (2001). A graph layout algorithm for drawing metabolic pathways. Bioinformatics 17:461-467. doi:10.1093/bioinformatics/17.5.461
6. Karp PD, Paley SM (1994). Automated drawing of metabolic pathways. 3rd Int. Conf. Bioinformatics and Genome Research. https://sri.com/?p=7318
7. Bourqui R, Cottret L, Lacroix V, Auber D, Mary P, Sagot M-F, Jourdan F (2007). Metabolic network visualization eliminating node redundance and preserving metabolic pathways. BMC Syst Biol 1:29. https://link.springer.com/article/10.1186/1752-0509-1-29
8. Okuda S, et al. (2008). KEGG Atlas mapping for global analysis of metabolic pathways. Nucleic Acids Res 36:W423-W426. https://pmc.ncbi.nlm.nih.gov/articles/PMC2447737
9. KEGG documentation: https://www.genome.jp/kegg/document/help_pathway.html ; https://www.genome.jp/kegg/document/help_global.html ; https://www.genome.jp/kegg/kegg1c.html ; https://www.genome.jp/kegg/pathway.html
10. Yamada T, Letunic I, Okuda S, Kanehisa M, Bork P (2011). iPath2.0: interactive pathway explorer. Nucleic Acids Res 39:W412-W415. https://pmc.ncbi.nlm.nih.gov/articles/PMC3125749
11. King ZA, Draeger A, Ebrahim A, Sonnenschein N, Lewis NE, Palsson BO (2015). Escher. PLoS Comput Biol 11:e1004321. doi:10.1371/journal.pcbi.1004321
12. Escher documentation, tips and tricks. https://escher.readthedocs.io/en/1.7.0/tips-and-tricks.html (compartment suffixes: cite BiGG documentation instead of the RCSB PDB news page)
13. Jensen PA, Papin JA (2014). MetDraw. Bioinformatics 30:1327-1328. https://pmc.ncbi.nlm.nih.gov/articles/PMC3998126
14. Noronha A, et al. (2017). ReconMap. Bioinformatics 33:605-607. https://academic.oup.com/bioinformatics/article/33/4/605/2608649 ; Heirendt L, et al. (2019). COBRA Toolbox v3.0. https://arxiv.org/pdf/1710.04038
15. Metabolic Atlas, 'Network layout engine', discussion #539 (2021). https://github.com/orgs/MetabolicAtlas/discussions/539
16. Wu H-Y, Noellenburg M, Sousa FL, Viola I (2019). Metabopolis. BMC Bioinformatics 20:187. https://pmc.ncbi.nlm.nih.gov/articles/PMC6466808
17. Wu T (2026). MetaCarto: biologically faithful automatic layout for genome-scale metabolic maps. bioRxiv preprint. https://www.biorxiv.org/content/10.64898/2026.09.19.752882.full.pdf
18. BioEMMA: automated generation of model-specific Escher-compatible maps from KEGG pathways (2026). Biology 15(17):1493. doi:10.3390/biology15171493
19. Schreiber F, Czauderna T (2022). Design considerations for representing systems biology information with the Systems Biology Graphical Notation. J Integr Bioinform. doi:10.1515/jib-2022-0024. https://pmc.ncbi.nlm.nih.gov/articles/PMC9377698 [corrected]
20. ChloroKB. https://pmc.ncbi.nlm.nih.gov/articles/PMC5462031
21. Ma H-W, Zeng A-P (2003). The connectivity structure, giant strong component and centrality of metabolic networks. Bioinformatics 19:1423-1430. Csete M, Doyle J (2004). Bow ties, metabolism and disease. Trends Biotechnol 22:446-450. Zhao J, et al. (2006). Hierarchical modularity of nested bow-ties in metabolic networks. https://arxiv.org/pdf/q-bio/0605003v1 [primary references substituted; verify page numbers]
22. SRI, Pathway Tools cellular overview and Omics Viewer. https://sri.com/?p=12582
23. Michal G (1998). On representation of metabolic pathways. BioSystems 47:1-7. doi:10.1016/S0303-2647(98)00019-7 [added; not read in this review]
24. Aichem M, et al. (2021). Visual exploration of large metabolic models (LMME). Bioinformatics 37(23):4460-4468. Docs: https://www.cls.uni-konstanz.de/software/lmme/learn/model-decomposition/ [added]
25. Czauderna T, Wybrow M, Marriott K, et al. (2013). Conversion of KEGG metabolic pathways to SBGN maps including automatic layout. BMC Bioinformatics 14:250. doi:10.1186/1471-2105-14-250 [added]
26. Siebenhaller M, et al. (2020). Human-like layout algorithms for signalling hypergraphs: outlining requirements. Brief Bioinform 21(1):62. https://academic.oup.com/bib/article/21/1/62/5115270 [added; verify title]
27. Schreiber F, et al. (2009). A generic algorithm for layout of biological networks. BMC Bioinformatics 10:375 [added]
28. Lambert A, et al. (2011). Pathway preserving representation of metabolic networks. Comput Graph Forum 30(3):1021-1030 [added]
29. Kelley JJ, et al. (2017). MOST-visualization: automated textbook-style maps of genome-scale metabolic networks. Bioinformatics 33(16):2596-2597 [added]
30. Saraiya P, North C, Duca K (2005). Visualizing biological pathways: requirements analysis, systems evaluation and research agenda. Information Visualization 4:1-15 [added]
31. BiGGr: model-specific Escher maps generated from templates or rule sets. https://biggr.org/about [added; not examined]
32. Purchase H (1997). Which aesthetic has the greatest effect on human understanding? Graph Drawing, LNCS 1353:248-261 [added]
