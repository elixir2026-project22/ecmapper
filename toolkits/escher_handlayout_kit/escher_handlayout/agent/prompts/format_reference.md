## Reply format: a JSON "spec fragment"

Every reply is ONE JSON object. Allowed top-level keys (all optional, anything else is rejected):

```
{
 "notes": "<=4000 chars of free text (your reasoning; ignored by the program)",
 "defaults": {"bend": 70},
 "zones":  [{"name": "Glycolysis", "color": "#2b5c9e", "title": "Glycolysis", "label": "glycolysis" | null, "bend": 45}],
 "plan":   {"Glycolysis": {"origin": [x, y], "size": [w, h], "reactions": ["492", "501"], "description": "..."}},
 "splits": [{"alias": "@xu5p_copy", "node": "871", "reaction": "909"}],
 "merges": [{"keep": "604", "dups": ["616", "656"]}],
 "layouts": [ BLOCK, ... ],
 "titles": {"Glycolysis": [x, y]},
 "small_titles": [{"reaction": "ATPM", "zone": "Oxidative phosphorylation", "text": "ATP turnover", "reuse_label": null}]
}
```

BLOCK = one piece of the drawing, normally one per subsystem:

```
{"block": "Glycolysis", "zone": "Glycolysis", "origin": [3000, 500],
 "nodes": {"480": [0, 0], "482": [0, 300]},                      // primary metabolite node id -> position RELATIVE to origin
 "rings": [{"name": "TCA ring", "center": [0, 0], "radius": 1050, "nodes": {"135": -90, "136": -50}}],   // angles in degrees, 0 = east, clockwise (y points down)
 "reactions": [ REACTION, ... ]}
```

REACTION (only `id` is required; defaults in brackets):

| key | meaning |
|---|---|
| `id` | reaction id (string, as in the table you are given) |
| `zone` | subsystem (default: the block's zone) |
| `side` | where cofactor dots hang: `"L" "R" "U" "D"` = left/right/up/down on screen [R]; `"out" "in"` = radially outside/inside of a ring reaction |
| `dist` | distance of cofactor dots from the reaction line, 20-800 [150]; use 110 when space is tight, 190 for crowded lanes |
| `gap` | half distance between the two junction points along the line [35]; use 85 for a middle lane with 4+ cofactors per side |
| `step` | spacing between several cofactors on one side [75] |
| `bend` | shift of cofactor dots along the line, gives the curves room (defaults come from `defaults.bend`/zone `bend`) |
| `shift` | lane offset `[dx, dy]` for parallel reactions between the same two nodes, or `{"perp": ["id1","id2"], "offset": n}` = n units perpendicular to the line id1->id2 |
| `axis` | override of the reaction direction: `[dx, dy]` or `{"from": P, "to": P}` (needed for branches and one-sided reactions) |
| `mid` | override of the reaction midpoint: a point P (relative to the block origin when given as `[x, y]`) |
| `curve` | `true`/`false`: force S-shaped lane bows (default: true only when `shift` is non-zero) |
| `bow` | very slight bow for long straight lines, -0.5..0.5 (e.g. 0.05) |
| `arc` | ring reactions: `{"ring": "TCA ring", "offset": 0}` (offset = radial lane offset; +outer, -inner) or `{"center": [x,y], "offset": 0}`; add `"arc_nodes": ["id_a","id_b"]` when the reaction touches other primary nodes besides the two ring nodes |
| `center` | `[x, y]`, only for `side` `"out"/"in"` without `arc` |
| `remove` | `true` in a correction: delete this reaction entry |

A point expression P is `[x, y]` (relative to the block origin), `{"node": "id"}` (global position of that node), `{"mean": ["id", ...]}` or `{"lerp": [P1, P2, t]}`.
Node ids are the ids in the data tables (for nodes created by `splits`, write the alias, e.g. `"@xu5p_copy"`).
`titles` positions are GLOBAL coordinates (baseline-left of the title text).
Coordinates must stay within +-20000. Never invent ids: only ids that appear in the data you are given are valid.

### Corrections / patches
Replies are *merged* into the current spec: blocks by `block`, reactions by `id`, nodes by id (a node listed again is moved),
zones by `name`, merges by `keep`. So a correction only needs the entries that change.
`"remove_nodes": ["id"]` (inside a block) deletes node coordinates, `{"id": "492", "remove": true}` deletes a reaction entry.
