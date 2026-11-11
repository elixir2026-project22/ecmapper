# GUIDE - engine API, layout recipes and how to read the QA metrics

This guide is the technical companion of `PROMPT.md`. Everything here was derived from a real run on the
iML1515 central-metabolism draft map (`specs/iml1515_central_carbon.layout.json` is the worked example).

## 1. Pipeline

```
python run.py inspect --input MAP.json --out analysis/        # facts about the draft map
(write specs/<map>.layout.json by hand-choosing coordinates)    # the creative step (a person, or an LLM: see AGENT.md)
python run.py check-spec --input MAP.json --spec specs/<map>.layout.json   # fast check, no build
python run.py build   --input MAP.json --spec specs/<map>.layout.json --out out/ --curv 0.6
                       -> out/<name>.json, out/<name>_overview.png, out/qa_report.md, out/qa_report.json
python run.py validate --input MAP.json --output out/<name>.json
```

`build` does, in this order: load map -> duplicate shared cofactor nodes (one per reaction) -> run `spec.build(L)`
-> check that every reaction/primary node is placed -> compute geometry (midmarkers, multimarkers, cofactor
positions, Bezier curves) -> place all labels without overlaps -> pack disconnected pieces into free space ->
export JSON -> validate -> compute QA metrics -> render PNG.

## 2. Conventions

- Node and reaction ids are the **strings** used as keys in the input JSON. Use `analysis/analysis.md` to look them up.
- x grows right, y grows down. Escher map units. A primary node has radius 20, a cofactor dot 10.
- Draft-map topology expected by the engine: `substrate -> multimarker(substrate side) -> midmarker ->
  multimarker(product side) -> product`; cofactors hang on the multimarkers. Reactions with other topologies (e.g.
  exchange reactions drawn directly to the midmarker, typical for hand-drawn maps) are reported as unsupported.
- **Primary** nodes (`node_is_primary`) are the backbone: you give them coordinates. **Secondary** nodes (cofactors)
  are placed automatically next to the reaction. Secondary nodes shared by several reactions are duplicated automatically.
- The reaction axis runs from the mean of its primary substrates to the mean of its primary products.
  If one side has **no** primary node (e.g. `ATP -> ADP + Pi + H2O` where only ATP is primary) give `axis=` and `mid=`.
- The direction of a reaction is the sign of the stoichiometry. It can be opposite to the pathway direction; keep it.

## 3. The layout spec (JSON, no code)

A spec is a JSON document (schema: `schema/layout_spec.schema.json`, single source of truth: `escher_handlayout/declspec.py`).
Sections: `zones` (subsystems: name, colour, title, existing label), `plan` (region and reaction list per zone), `splits`
(second node of a metabolite for one reaction), `merges` (hub nodes), `layouts` (blocks with `origin`, `nodes`, `rings`, `reactions`),
`titles`, `small_titles`, `defaults`. Every key and value is described in
`escher_handlayout/agent/prompts/format_reference.md` (the same text a chat model gets); the options of a reaction entry
(`side`, `dist`, `gap`, `step`, `bend`, `shift`, `axis`, `mid`, `curve`, `bow`, `arc`, `arc_nodes`, `center`, `secdist`) are exactly the keyword
arguments of the engine's `Layout.reaction(...)`. `specs/iml1515_central_carbon.layout.json` is a complete worked example,
`specs/layout_template.json` an empty skeleton. Expressions (`{"node":..}`, `{"mean":[..]}`, `{"lerp":[P1,P2,t]}`,
`{"perp":[a,b],"offset":n}`) replace the small arithmetic a Python spec used to contain; regular polygons come from `rings`.

Legacy: `specs/legacy_python/` holds the old Python-code specs (they execute code; `--unsafe-python-spec` is needed).

## 4. Layout recipes (numbers that worked)

- **Linear chain**: primary nodes in a column, 300 apart; cofactors on the right (`side='R', dist=110`); reverse
  reactions of the same two nodes as left lanes (`shift=(-110,0), side='L'`).
- **Branch (aldolase-like)**: `axis=(0,1)`, `mid` slightly below the substrate, the two products 300 below and
  300 apart; an isomerase between the two products goes horizontally.
- **Hourglass (2 substrates -> 2 products, e.g. transketolase, ACOXT)**: two nodes stacked/aligned on each side
  (600-700 apart), midmarker between them. Chain several hourglasses (TKT1 -> TALA -> TKT2). If a metabolite is
  needed at two places use a `splits` entry instead of a long line.
- **Cycle (TCA)**: nodes on a circle (e.g. 9 nodes every 40 degrees, R about 1050, citrate on top, clockwise), every
  ring reaction with `arc=(center, 0)`. Parallel reactions between the same two ring nodes: concentric arcs with
  `arc=(center, +-150)` (inner lanes `side='in'`, outer lanes `side='out'`). Shunts through the interior as ordinary reactions.
- **Isozyme / cofactor variants of one step** (PYK, PYK2, PYK3, PYK4): put them on ONE lane and stack their cofactor
  pairs along it: `shift=(0, +-75..225)`, `curve=False`, `gap=40`; make the distance between the two nodes at least 650.
- **Parallel duplicates** (two MDH): lanes 110-150 apart; cofactors of inner lanes inward, of outer lanes outward.
- **Rows of repeated reactions** (NADH dehydrogenases): one horizontal row per quinone couple, two lanes per row
  (`shift=(0,-45) side='U'` and `shift=(0,+45) side='D'`), rows about 470 apart.
- **Hub with several reactions between two nodes** (ATP <-> ADP): three lanes at -190/0/+190 with
  `gap=85, secdist=430` for the middle one so its cofactors clear the outer lanes.

## 5. Interconnecting subsystems

Goal: one connected map in a biologically plausible flow, connected through the metabolites the subsystems share,
with few crossings. Rearrange the subsystems (move, mirror, rotate) so that the shared metabolites are neighbours.

1. Look at `merge_candidates` in `analysis.md` (same metabolite drawn as several primary nodes). Merge them by rule M
   of `PROMPT.md` (connectivity plan, planarity check) - at least the default hubs: glucose, G6P, F6P, G3P, PEP, pyruvate, acetyl-CoA, acetaldehyde, OAA, 6-phosphogluconate
   (and analogous hubs in other maps: the 3-5 metabolites that connect the big modules).
2. `merges: [{"keep": id, "dups": [ids]}]`: all segments of the duplicates are re-pointed to `keep`; give `keep` one coordinate (only primary nodes of the SAME metabolite can be merged).
3. A metabolite that rule M leaves as a copy keeps its exit label (`-> ...`), which tells the reader where it
   continues. Duplicate reactions drawn in several subsystems become lanes of the host pathway.
4. Never merge cofactors (ATP/ADP/NAD(H)/NADP(H)/CoA/H2O/Pi/CO2/O2/quinones...).
5. Choose the side on which a module hangs so that hubs are close to each other and connecting lines do not cross a
   main chain. Rearrange first; a long line only if no arrangement avoids it.
6. Default flow for central carbon metabolism: glycolysis = central column; pyruvate -> acetyl-CoA below it; TCA ring
   below that; PPP (and ED branching at 6PG) to one side of upper glycolysis; the PEP/pyruvate/OAA node on the side
   between glycolysis and TCA; fermentations beside pyruvate/acetyl-CoA; respiratory chain separate.

## 6. Curvature (`--curv`, 0..1)

For a segment between a metabolite M and a multimarker K (axis `u`, normal `n`, `d=(M-K).n`, `D=(M-K).u`, `s=sign(D)`):
straight if `|d|<25` or `|D|<20`; otherwise Bezier control points (swap b1/b2 for product segments)

```
gentle (0): c_far = M + 0.33 (K-M);                  c_near = K + u s min(0.30|D|, 0.35|d|, 60)
strong (1): c_far = K + u (D + s 0.2|d|);            c_near = K + u s min(0.45|D|, 0.40|d|)
used      : componentwise (1-t)*gentle + t*strong ; lanes k=(0.32+0.10t)*length ; ring handles (0.33+0.03t)*chord
```

0.6 was chosen by the user for the example; 0.75 already caused overlapping curves at shared nodes. If unsure use 0.6.

## 7. What is automatic (do not hand-tune)

- reaction labels: 8 directions x 5 radii around the midmarker, scored against lines, nodes and other labels;
  metabolite and exit labels likewise; white halo only in the PNG.
- packing: disconnected components are moved into free space (grid 40 units, FFT overlap search) minimising the
  bounding box; each moved piece gets a title from `L.small_titles` or a generated one.
- canvas size, label_x/label_y conversion to Escher's start-anchored text.

## 8. QA metrics and what to change

| metric (in `qa_report.md`) | meaning | typical remedy in the spec |
|---|---|---|
| validation errors | structure/biochemistry broken | must be 0; read the message; usually a missing placement or a wrong merge |
| line_crossings | segments of different reactions crossing | see `worst_crossing_pairs`: move one node, flip a `side`, change lane `shift`, split a node instead of a long line |
| labels_on_lines | labels with >=2 line samples inside | widen `dist`/`step`/lane spacing around the named reaction; reduce `curv` |
| labels_on_nodes | label box touches another node | move the node or increase spacing between neighbours |
| label_label_overlaps | two label boxes overlap | spread the two reactions apart (>= 150 units) or change cofactor `side` |
| primary_nodes_closer_than_70 | nodes nearly on top of each other | fix coordinates |
| bbox_area | size of the final rectangle | compact the main component; the packing already removes the pieces' share |

Reference values of the bundled example (91 reactions): crossings 26, labels_on_lines 5, labels_on_nodes 0,
label overlaps 0, close nodes 0, 7444 x 8284 units.

## 9. Known limitations

- Only the draft-map topology is supported (see section 2).
- Label widths are estimated (0.58-0.64 x font size x characters); the real Escher font may differ, so keep clearance.
- The preview PNG is a matplotlib re-implementation, not Escher; the JSON was never loaded into Escher by the tool.
- A cofactor line may cross a neighbouring lane when 4-5 parallel reactions share two nodes.
