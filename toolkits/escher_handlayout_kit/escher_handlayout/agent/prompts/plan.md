## TASK: PLAN the global arrangement

Below is the draft map as DATA. Decide the subsystems, where each goes on the canvas, which hub metabolites are merged and
which metabolites need a second node. Return ONE JSON fragment with these keys:

- `zones`: one entry per subsystem (start from the PROPOSAL; keep `label` exactly as proposed, it links to an existing title
  label of the map; use null for a new zone). Zone `name` must match the keys used in `plan`.
- `plan`: for EVERY zone an entry with `origin` [x, y] (top-left of its region), `size` [w, h] (generous: leave >= 300 free
  around it), and `reactions` = the ids of its reactions. EVERY reaction id of the map must appear in exactly one zone.
  List the zones in the order you want to lay them out (hubs get their coordinates in the first zone that uses them, so start
  with the zone that forms the central spine).
- `merges`: hub merges (rules in the recipes): merge the metabolites the subsystems share and place the zones so that
  merge partners are neighbours. Use only ids from the MERGE CANDIDATES table.
- `splits`: only if a metabolite is needed twice inside a zone.
- `notes`: your arrangement in two or three sentences.

Do not output `layouts` or `titles` yet; they come in the next turns.

### DATA: map overview
{{OVERVIEW}}

### DATA: proposal (zones from the title labels / connected components)
```json
{{PROPOSAL}}
```

### DATA: merge candidates (same metabolite drawn as several primary nodes; `hub` = in the default hub list)
{{MERGE_TABLE}}
