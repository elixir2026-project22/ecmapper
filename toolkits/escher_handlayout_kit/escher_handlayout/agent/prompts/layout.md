## TASK: LAYOUT of zone "{{ZONE}}" ({{CHUNK_INFO}})

Place the primary nodes and describe every reaction listed below. Reply with ONE JSON fragment:

```
{"layouts": [{"block": "{{BLOCK}}", "zone": "{{ZONE}}", "origin": {{ORIGIN}}, "nodes": {...}, "rings": [...], "reactions": [...]}]{{TITLE_HINT}}, "notes": "..."}
```

- Node coordinates are RELATIVE to `origin`, which is the origin planned for this zone: keep {{ORIGIN}}. Planned region of this
  zone: origin {{ORIGIN}}, size {{SIZE}} (global coordinates). Stay inside it.
- Give coordinates to the primary nodes under "NODES TO PLACE" and nothing else. Nodes under "ALREADY PLACED" belong to other
  zones (hubs): do not place them again, just use their ids in reactions or expressions; their global positions are given.
  To connect to a placed node, put your nodes at a sensible relative position (global position - origin).
- Describe EVERY reaction under "REACTIONS" exactly once, with the options needed (side, shift, axis/mid, arc, ...).
- Occupied space so far (avoid overlaps): {{OCCUPIED}}

### DATA: NODES TO PLACE
{{TO_PLACE}}

### DATA: ALREADY PLACED (global coordinates)
{{PLACED}}

### DATA: REACTIONS  (`[id*]` = primary node id, `@alias` = node created by a split; reversible reactions marked `<->`)
{{REACTIONS}}
