## TASK: IMPROVE the layout (refinement round {{ROUND}} of {{MAX_ROUND}})

The program built the map from the current spec and measured its quality (lower is better). Fix the NAMED problems with a
SMALL correction: send only the entries that change (blocks by `block`, reactions by `id`, nodes by id). Typical remedies:
move one node, flip a `side`, widen `dist`/`step`/lane `shift` of the named reaction, change a lane order, give a metabolite a
second node with `splits`, or drop a merge that forces a long line. Do not rearrange parts that have no problem.

### QUALITY REPORT (DATA)
{{REPORT}}

### CURRENT ENTRIES OF THE REACTIONS NAMED IN THE REPORT (DATA; global node positions in brackets)
```json
{{ENTRIES}}
```

If nothing sensible can be improved, reply `{"notes": "no change"}`.
