# GUIDE - engine API, layout recipes and how to read the QA metrics

This guide is the technical companion of `PROMPT.md`. Everything here was derived from a real run on the
iML1515 central-metabolism draft map (`specs/iml1515_central_carbon.layout.json` is the worked example).

## 1. Pipeline

```
python run.py inspect --input MAP.json --out analysis/        # facts about the draft map + arrangement.md (where the subsystems go)
python run.py directions --input MAP.json [--check SPEC.json]  # the same arrangement brief alone / the directions a spec violates
python run.py start   --input MAP.json --out specs/<map>.layout.json       # a complete spec: the draft's subsystems, arranged and connected (section 3a)
python run.py patch   --input MAP.json --spec specs/<map>.layout.json --json '{...}'   # change a part of the spec (section 3a)
python run.py tidy    --input MAP.json --spec specs/<map>.layout.json   # rows and columns, loose ends and lanes untangled, cofactors in free space (section 8)
(or write specs/<map>.layout.json by hand-choosing coordinates)  # the creative step (a person, or an LLM: see AGENT.md)
python run.py check-spec --input MAP.json --spec specs/<map>.layout.json   # fast check, no build (--partial: unfinished spec)
python run.py build   --input MAP.json --spec specs/<map>.layout.json --out out/ --curv 0.6
                       -> out/<name>.json, out/<name>_overview.png, out/qa_report.md, out/qa_report.json;
                          prints score, targets and what to do next (section 8)
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

## 3a. Start spec and patches (changing a layout instead of writing one)

Writing a complete spec from an empty page is the hardest way to use the kit, in particular for a language model that
is not among the strongest. Two commands make the same work a series of small changes:

- `start` writes a complete, valid spec of the draft map: one zone and one block per subsystem (the title
  above each box of the draft; a component without title is "Module n"), the primary nodes at their draft positions
  relative to the block's `origin`, every reaction with the `side` and `shift` (or `axis`/`mid`) the draft draws.
  Moving a subsystem is then one change of an `origin`. With `--as-is` that is all: the draft as it is, no merges.
  Without it `start` first redraws two shapes (`escher_handlayout/redraw.py`): in a block with a cycle of at least 6
  primary metabolites the longest cycle becomes a circle (the nodes at equal angles, in their order, 300 to 600 units
  apart; the reactions along it become arcs, parallel ones in lanes 150 units apart; the rest of the block keeps its
  place relative to the cycle), and in the other blocks a chain of at least 4 metabolites that runs down or across
  the block and leaves its line by at most 8 % of its length is put on one straight line. Then it does a first
  arrangement (`escher_handlayout/connect.py`): the first block (glycolysis if
  the map has it, else the block that shares most) stays; the others are attached one after the other, the one that
  shares most metabolites with the placed blocks first. A block is moved as a whole. For every metabolite it shares,
  its node is merged into the node of the placed blocks and the block is put where that costs least: near the kept
  node, not on top of what is placed (raster of 50 units, 250 units margin), where the placement catalogue wants the
  subsystem (500 units per penalty point), with few crossings in a straight-line sketch (400 units each). Further
  copies of a metabolite are merged when the two nodes ended up at most 1500 units apart and the sketch gets no new
  crossing; when a block is attached, further metabolites it shares are merged at once (up to 1500 units) if the
  sketch of the block gets shorter or loses crossings by it. A node of a ring is the one that is kept in a merge.
  The others stay copies and are listed in the answer of `start`. Blocks that share only cofactors are left
  to the packer. The result is a connected map (central metabolism of iML1515: 50 crossings, no violated direction,
  score 182 instead of 550 for the draft as it is; E. coli core: score 157 instead of 542, one violated direction and
  one detached piece); it is a first arrangement: blocks are not mirrored, and apart from the ring and the straight
  chain their drawing is the draft's. `node_names` (per block) and `name` (per reaction) say
  what the ids are; they are optional, and `check-spec` reports an id whose name does not fit the map.
- `patch` merges a spec fragment into the spec file, with the rules of the chat agent (blocks by `block`, nodes by id,
  reactions by `id`; a reaction entry replaces the old one completely). A fragment that violates the schema or makes a
  new problem (unknown id, wrong name, a merge of different metabolites ...) is rejected and the file stays as it was.
  A fragment with `merges` also removes the coordinates of the merged `dups`; a merge is undone with
  `{"keep": "604", "remove": true}` or by sending it with fewer `dups`, and the nodes that are copies again get
  coordinates beside the kept node. A node put on a ring loses its coordinates in `nodes`; the title of a subsystem
  moves with the `origin` of its block. The file is rewritten in a compact form
  (one reaction per line) that is easy to read and to edit by hand.

`patch` and `check-spec` also print warnings for specs that are formally right but break a rule (`escher_handlayout/lint.py`):
a merged cofactor, a hub merge that pulls in a reaction with 6 or more primary metabolites (biomass: it keeps its own
copies), a reaction whose metabolites are more than 3000 units apart, primary nodes closer than 70 units.

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

## 5a. Where the subsystems go: the placement catalogue (`layout_directions/`)

`layout_directions/placement_directions_v2.json` records where metabolic subsystems sit relative to each other on published
overview maps (KEGG 01100, Human-GEM, Escher example maps, CBMViz) and in textbook conventions (details: `layout_directions/README_v2.md`).
`escher_handlayout/directions.py` applies it to a map:

1. The subsystems (title labels / zones) are matched to the catalogue's entities: BRENDA map names through
   `brenda_aliases.json`, otherwise entity names and aliases, otherwise word overlap. The map's primary metabolites are matched
   to the catalogue's hub metabolites by BiGG id.
2. `inspect` writes `arrangement.md` (the chat agent shows the same text in its plan turn): for each pair of subsystems the
   direction with the best evidence - stated rules (`direction_rules`) and the consensus of the measured maps (pairs where the
   maps conflict or only one map was measured give no direction) - plus the cycle layouts (TCA clockwise, oxaloacetate at the
   top ...), the attachments (which metabolite joins two subsystems: merge candidates) and textbook directions between hub
   metabolites (glucose above G6P ...). Directions are screen directions: "X is below Y" means X gets the larger y.
   Strength `hard`/`strong`: follow it; `medium`: follow it unless it costs crossings; `weak`: a hint.
3. `build` measures the finished layout against the directions (zone centres and node positions after packing) and reports the
   violated hard/strong/medium ones in `qa_report.md` with a **placement penalty** (5 per hard/strong, 2 per medium
   direction). The chat agent adds the penalty to its score; `directions --check SPEC` gives the same numbers from a spec
   without a build. `--taxa bacteria|plants|animals|fungi|...` drops rules scoped to other organism groups.

The catalogue is evidence from a few maps, not a standard: when a direction costs crossings or long lines, the crossings
count more. A map whose subsystems match no entity (no title labels, unusual names) simply gets no arrangement section.

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
| placement penalty | subsystems (or hub metabolites) not where the placement catalogue puts them (section 5a) | move the block (`origin`) - only if that does not cost crossings; medium directions may stay violated |
| unconnected subsystems | subsystems without a reaction in the main component that draw a metabolite (no cofactor) another piece draws too | move the block next to its partner and merge the metabolite; the build names the merge |
| detached pieces | pieces beside the main component that draw a primary metabolite of it again (cofactors do not count) | merge the metabolite (`merges`), or keep the copy on purpose |

`python run.py tidy` (`escher_handlayout/tidy.py`) repairs what it can of these numbers without touching the design: it
proposes rows and columns (primary metabolites that one reaction joins and that are less than 60 units off the upright
or the level, along chains of such reactions) and a `side` and `dist` for every reaction (where its cofactors and their
branches cross the fewest lines and have the fewest nodes around them), takes the proposals block by block and keeps a
block's only if crossings, labels on lines, labels on nodes, label overlaps, close nodes and the placement penalty
together get better (with long lines, unconnected subsystems and detached pieces, so that no step trades one for the
other), measured as `build` measures them; then it tries the other side and other distances for the reactions that
the metrics still name. Between the two it untangles locally, with the rest of the map pinned: for the reactions that
the metrics name, a metabolite that only this reaction uses is turned around the rest of the reaction (30, 60, 90 and
180 degrees), and a lane (`shift`) is put on the other side of its line unless another lane is there. With `--copies`
a named reaction whose line to a metabolite of three or more reactions is longer than 500 units gets its own copy of
it (`splits`), 400 units from its other metabolites, if that saves more than 3 points. It writes plain values into the
spec. `--block NAME` (several times) limits it to blocks, `--only align` / `untangle` / `cofactors` to one step,
`--seconds N` (90) the time it tries.

`build` turns these numbers into a decision (`escher_handlayout/score.py`), on the screen and in `qa_report.md`:

- the **score** `1000 x validation errors + crossings + 3 x labels on lines + 6 x labels on nodes + 4 x label overlaps
  + 10 x close nodes + 3 x long lines + 3 x crowded reactions + area excess + 40 x unconnected subsystems + 20 x detached pieces + placement
  penalty`; a long line is a reaction with a primary metabolite more than 800 units from its middle, the area excess
  is 1 point per 0.01 million square units of map area per reaction above 0.40; a crowded reaction (3 points each) has
  4 or more other reaction middles within 250 units of its own: a hub without room;
- the **targets**, each `ok` or `MISS`: validation errors 0, crossings at most 0.1 x reactions, labels on lines at most
  max(6, 7 % of the reactions), labels on nodes 0, label overlaps at most 2, close nodes 0, long lines at most max(4, 10 % of the reactions),
  area excess at most 40 (0.80 million square units per reaction), crowded reactions at most max(2, 3 % of the reactions), unconnected subsystems 0,
  no violated hard or strong direction; detached pieces are shown as `look`, because a copy can be deliberate;
- the history of the job: every build is recorded in `builds.json` / `builds.md` beside the `--out` folder (`--history
  FILE`, `--no-history`), a rebuild into the same folder replaces its entry;
- one last line that says what to do (advice on the drawing: the design of the map decides the arrangement, a target
  may stay missed for a reason noted in `DECISIONS.md`, and a build that meets the targets with a higher score than
  another build names that build too): `NOT FINISHED ... Fix first: <the unconnected subsystems with the merge that
  connects each, else the missed target that costs most, its offenders, the remedy>`, `NO IMPROVEMENT` after three changed versions of the spec that were worse than the best one with no improvement between them (a version that ties with the best, such as the way back to it, is not counted), `STOP` after
  `--max-builds` (16) changed versions, `TARGETS MET`.

Packing and unconnected subsystems: the packer moves the small disconnected pieces into free space, but the largest
piece of a subsystem that is not connected yet (and shares a metabolite other than a cofactor with another piece) is a
block of its own and stays where the spec puts it, with the subsystem's title, unless it overlaps something there. So
the arrangement of a start spec survives the build, and the placement check judges the positions the spec gives.

The chat agent (`agent.py`) uses the same score; its own crossing target is 0.35 x reactions. The bundled example
misses the crossing target of 0.1 (26 crossings for 91 reactions) and 5 strong directions: it predates both.

Reference values of the bundled example (91 reactions): crossings 26, labels_on_lines 5, labels_on_nodes 0,
label overlaps 0, close nodes 0, 7444 x 8284 units.

## 9. Known limitations

- Only the draft-map topology is supported (see section 2).
- Label widths are estimated (0.58-0.64 x font size x characters); the real Escher font may differ, so keep clearance.
- The preview PNG is a matplotlib re-implementation, not Escher; the JSON was never loaded into Escher by the tool.
- A cofactor line may cross a neighbouring lane when 4-5 parallel reactions share two nodes.
