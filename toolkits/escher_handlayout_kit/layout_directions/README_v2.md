# MLL Global Placement Directions v2.0.0

This is a multi-source catalogue of evidence for where metabolic subsystems sit relative to each other on overview maps. Version 1 contained 229 hand-written rules over 298 free-text names. Version 2 rebuilds the catalogue on a controlled vocabulary and keeps every piece of evidence tagged with its source. It does **not** resolve disagreements between maps; it measures them.

## Files

| File | Content |
|---|---|
| `placement_directions_v2.json` | The catalogue (vocabulary, relations, rules, measured observations, consensus) |
| `placement_directions_v2.schema.json` | JSON Schema 2020-12 (enums for directions, relationships, strengths, entity id pattern) |
| `consensus_pairs.csv` | Flat table of the consensus section, for spreadsheets |
| `build_scripts.zip` | `vocab.py`, `legacy_map.py`, `knowledge.py`, `extract.py`, `derive.py`, `make_schema.py` |

## Coordinate frame

The frame has x pointing east and y pointing north. Bearings are measured clockwise from north, and distances are fractions of the canvas diagonal. In every section, `direction` gives the position of the **target relative to the anchor**.

## Vocabulary (`entities`, 402 entries)

| Type | Count | Notes |
|---|---|---|
| category | 12 | KEGG level-2 classes plus `transport_and_shuttles` |
| region | 1 | `central_carbon_metabolism` |
| pathway | 183 | roughly KEGG reference-map granularity; `xrefs.kegg_map` where one exists |
| module | 142 | finer units (gluconeogenesis, oxidative PPP, serine biosynthesis, carnitine shuttle ...) |
| metabolite_hub | 64 | branch-point metabolites with BiGG and KEGG compound ids, used as attachment points |

- **Hierarchy.** Every entity has `part_of`. Rules stated at different granularities can therefore be compared after lifting them to a common level, and the same hierarchy can drive expand/collapse.
- **Aliases.** `aliases` lists the names each source uses for the entity: KEGG, Human-GEM, iJO1366, the Synechocystis FAME model and MLL v1.
- **Process role.** `process_role` is one of central, biosynthesis, degradation, energy, transport or interconversion. It is assigned heuristically and supports role-level conventions such as MetaCyc's.

## Sources

| id | Kind | What was measured |
|---|---|---|
| `kegg_01100` | map | KEGG global map (KGML ko01100, 2012 build, from Biopython test data): 149 pathway-label positions plus hub compounds that occur once |
| `human_gem_maps` | map | 12 Metabolic Atlas compartment maps (Human-GEM 2.1.1, Escher export): median reaction position per Human-GEM subsystem |
| `escher_bigg_maps` | map | 12 Escher example maps (E. coli core/iJO1366, RECON1, iMM904): subsystems from iJO1366 or Human-GEM BiGG cross-references |
| `cbmviz_synechocystis_map` | map | Synechocystis map in your CBMViz repo, with subsystems from the FAME model |
| `cbmviz_generic_map` | map | "Mothermap" in your CBMViz repo; positions estimated from iJO1366 marker metabolites |
| `metacyc_cellular_overview` | literature | documented convention: biosynthesis west, energy middle, degradation east, flow downward |
| `knowledge_textbook` | knowledge | 48 hub-level wall-chart conventions plus 5 cycle layouts |
| `knowledge_biochemistry` | knowledge | 408 orientation-free relations |
| `mll_v1` | knowledge | your 229 rules mapped onto the vocabulary (228 rules after expansion, 20 policies) |

### Coverage gaps

- **Roche/Michal Biochemical Pathways chart.** It is a natural additional source, but I found no machine-readable grid index of it.
- **BRENDA map and KEGG maps 01110/01120.** I didn't include these. The KEGG maps would need KGML from the KEGG API, which wasn't reachable from this environment.

## Sections

### `relations`

These are 408 orientation-free biochemical links. The relationship types are FEEDS (directed material flow), SHARES_HUB, ALTERNATIVE_ROUTE, REVERSE_OF, SAME_FAMILY, COFACTOR_SUPPLY, ENERGY_COUPLING, TRANSPORT and SAME_COMPARTMENT.

- **`via` field.** It names the connecting hub metabolites. This tells a layout engine where to attach a branch, for example at 2-oxoglutarate on the TCA ring rather than anywhere on it.
- **`scope.taxa` field.** It gives the organisms the relation applies to, for filtering.

### `direction_rules`

These are stated (not measured) compass or semantic directions: the textbook hub conventions (`DIR####`) and the legacy rules (`GPR###`, with a suffix letter where a composite v1 name was expanded).

- **Mapping record.** `legacy_index` shows what became of each v1 rule.
- **Reversed rules.** GPR112 and GPR144 carry `flow: target_to_anchor` because their v1 anchor was the product.
- **Reclassified rules.** The cofactor rules mislabelled ENERGY_COUPLED_TO in v1 are now COFACTOR_SUPPLY.

### `observations`

These are 1,952 directions measured from maps. Pairs come from Delaunay neighbours within one canvas. Pairs where one entity is part of the other are skipped, and so are subsystems whose reactions are too dispersed to have a meaningful position.

- **Confidence.** It compares the separation of the two positions with their dispersion.
- **Category-level pairs.** Observations with "category centroids" in `evidence` compare whole categories (for example lipid vs amino acid) within a canvas.

### `absolute_placements`

These record each entity's position within its canvas as fractions plus a nine-cell region (NW, N, NE, W, C, E, SW, S, SE). Use them for "where on the whole map" questions, such as the fact that KEGG puts lipids west and xenobiotics south.

### `consensus`

Each pair of entities gets one record per level:

- **`exact`.** Entities as observed.
- **`pathway`.** Modules and hubs are lifted to their pathway so that sources of different granularity meet.

Each source contributes one averaged vector, so Human-GEM's 12 canvases count once. The record then reports:

| Field | Meaning |
|---|---|
| `agreement_R` | length of the mean unit vector (1 = all sources point the same way) |
| `status` | `single_source`, `agreement` (R ≥ 0.85), `partial` (0.5–0.85) or `conflict` (R < 0.5) |
| `suggested_strength` | strong if ≥3 sources agree, medium if 2 agree, weak if partial or single source, none if conflicting |
| `x_order` / `y_order` | axis-wise votes, usable as ordinal separation constraints even when compass bearings disagree |
| `per_source` | each source's direction, weight and example evidence |

At pathway level there are 1,176 pairs. Of these, 54 show agreement, 93 are partial, 65 conflict and 964 are single-source.

## Examples of what the data says

- **Glycolysis → pyruvate metabolism is *below*.** All 7 sources agree (R = 0.87).
- **Fatty acid biosynthesis sits west/north-west of glycolysis.** KEGG, Human-GEM and the Escher maps agree (R = 0.98).
- **PPP relative to glycolysis is only *partial* (R = 0.56).** Most sources put the PPP to the north, but east vs west varies, and KEGG and the textbook convention put it due east.
- **Propanoate metabolism vs TCA is a genuine *conflict* (R = 0.07).** Each map attaches it at a different side of the cycle. A rule based on the attachment hub (succinyl-CoA) is a better constraint here than a compass direction.
- **Of the v1 rule pairs, 73 now have map evidence.** Sixteen of them agree with the maps and 24 conflict, so it is worth reviewing those 24 before treating them as strong.

## Suggested use in a layout engine

1. **Filter by the organism.** Use `scope.taxa` and drop entities the model does not contain.
2. **Turn hard and strong evidence into separation constraints.** Use the `x_order` and `y_order` votes from consensus, plus your own chosen direction rules. Turn medium and weak evidence into weighted penalties, and ignore `conflict` pairs or replace them with an attachment constraint taken from `relations.via`.
3. **Apply the cycle layouts.** Use `cycle_layouts` and `relations.via` to attach branches at the right intermediate.
4. **Treat cofactor and energy links as visual overlays only.** That means COFACTOR_SUPPLY, ENERGY_COUPLING and the v1 policies.

## Reproducing or extending

The build steps run in this order: `extract.py`, then `derive.py`, then `make_schema.py`.

- **Inputs.** The scripts expect the source files listed in `sources`, and paths are set at the top of each script.
- **Adding a source map.** Supply an extractor that produces `{source, canvas, entity, x, y, n, spread, kind}` points, together with an alias table that maps its subsystem names onto the vocabulary.
