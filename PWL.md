# PWL 2.0: text descriptions of metabolic subsystems

PWL (Pathway Layout Language) describes one metabolic subsystem per file in a
compact, line-oriented text. A file records two things:

- **the drawing**: heading, turns, rings, stars, fans and wrapped chains, so a
  layout can be reproduced;
- **the biochemistry**: which compounds are converted, on which side of each
  reaction they are, in which direction it runs, which cofactors it uses, and
  which motifs (cycles, shunts, ...) the subsystem contains.

Version 1.1 came as a separate download. Version 2.0 is part of ecmapper:

| What | Where |
|---|---|
| Parser, validator, motif detection (no dependencies) | `src/ecmapper/pwl/syntax.py`, `motifs.py` |
| Grammar | [`src/ecmapper/pwl/pwl.ebnf`](src/ecmapper/pwl/pwl.ebnf) |
| Generator for the BRENDA map | `src/ecmapper/pwl/brenda.py`, `compounds.py` (needs `pip install .[pwl]`) |
| Check against the BRENDA export | `src/ecmapper/pwl/check.py` |
| The 193 BRENDA subsystems in PWL 2.0 | [`BRENDAmap/pwl/`](BRENDAmap/pwl) |
| Tests | `tests/test_pwl.py` |

```bash
python -m ecmapper.pwl generate BRENDAmap/BRENDA_map_reactions.json BRENDAmap/pwl
```

```bash
python -m ecmapper.pwl parse BRENDAmap/pwl/subsystems/citric-acid-cycle.pwl
```

```bash
python -m ecmapper.pwl motifs my_pathway.pwl
```

`generate` downloads BKMS-react and the KEGG compound list into the ecmapper
cache, as the other ecmapper commands do (`--no-bkms`, `--no-kegg` to do without).
It writes the files and then runs the check, whose result is in
`validation_report.json`.

```python
from ecmapper.pwl import parse_file, detect_motifs

doc = parse_file("BRENDAmap/pwl/subsystems/urea-cycle.pwl")
for r in doc.reactions():
    print(r.rid, r.enzyme_label, [p.token for p in r.substrates or []], "->", [p.token for p in r.products or []])
```

## Example

An excerpt of `BRENDAmap/pwl/subsystems/urea-cycle.pwl` (declarations shortened):

```
pwl 2.0;
subsystem urea-cycle 'urea cycle' {
  flow down;
  compound ADP brenda:13 kegg:C00008;
  compound Arg brenda:123 kegg:C00062;
  compound L_citrulline 'L-citrulline' brenda:938 kegg:C00327;
  compound NH3 brenda:27;
  ...
  port Orn side L to arginine-metabolism;

  ring ccw {
    Orn
      <-L[node4993: 3.5.1.20 | -H2O, +CO2, +NH3 {br=BR19238, ev=cof+eq}]- L_citrulline  // citrullinase
      -[node5843: 6.3.4.5 | ATP>AMP, -Asp, +PPi {br=BR142, ev=cof+eq, match=kegg}]-> L_argininosuccinate  // Argininosuccinate synthase
      -B[node4981: 4.3.2.1 | +fumarate ^R {br=BR3431, ev=cof+eq, match=kegg}]-> Arg;  // argininosuccinate lyase
    Arg -[node4979: 3.5.3.1 | -H2O {br=BR19260, ev=cof+eq}]-> (urea + Orn);  // arginase
    chord (Orn + carbamoyl_phosphate) -[node5838: 2.1.3.3 | +Pi {br=BR43950, ev=cof+eq}]-> L_citrulline;  // ornithine carbamoyltransferase
  }
  ...
  motif cycle { Arg Orn L_citrulline L_argininosuccinate };
}
```

In version 1.1 arginase was written `urea -[3.5.3.1 | -H2O]- (Orn + Arg)`: the
group put the substrate (Arg) and a product (Orn) on one side, and nothing said
which way the reaction runs.

## Language

The grammar is [`pwl.ebnf`](src/ecmapper/pwl/pwl.ebnf). In short:

**File.** `pwl 2.0;` then `subsystem <id> '<title>' { ... }`. Statements end
with `;`. `//` starts a comment.

**Only the layout is required.** PWL is primarily a layout description, so
everything biochemical is optional and can be added or revised later:

- compound and compartment declarations;
- cross-references;
- reaction IDs and enzymes;
- coefficients;
- directions.

A file can start as `A -[]- B -[]- C;` and gain detail step by step. The
smallest file:

```
pwl 2.0;
subsystem demo 'demo' {
  A -[PGI]- B -[]- (C + D);       // a reaction given only by its own ID, and one without any
  direction PGI ->;               // added later, without touching the layout
}
```

**Declarations.**

| Statement | Meaning |
|---|---|
| `flow down;` | Main heading of the drawing (`down`, `up`, `left`, `right`). |
| `compartment side1 'side 1';` | A compartment. |
| `compound Ac_CoA 'Ac-CoA' brenda:29 kegg:C00024;` | A compound: token, label, cross-references `namespace:id`. Optional: undeclared tokens are allowed, and the namespace can be anything (`local:my_id`). |
| `port pyruvate side R to gluconeogenesis;` | The drawn node is linked to another subsystem in that direction. |

**Reactions.** A chain alternates nodes and links:
`A -[r1: 1.1.1.1 | NAD>NADH]-> B -[r2: ...]- C;`

- A **node** is a compound (`pyruvate`, the second drawn copy `pyruvate#2`, in a
  compartment `H@side2`, with an optional coefficient `2*H2O`), a **group** `(A + B)`, or
  the open end `_`. Members of a group are on the **same side** of the reaction.
  The node before the link and the node after it are opposite sides.
- The **link arrows** give the direction:

  | Link | Direction |
  |---|---|
  | `A -[..]- B` | unknown (the sides are known, the direction is not) |
  | `A -[..]-> B` | A to B (reversibility not stated) |
  | `A <-[..]- B` | B to A |
  | `A <-[..]-> B` | reversible (the written order then means nothing) |
  | `A -[..]=> B` | A to B, irreversible |
  | `A <=[..]- B` | B to A, irreversible |

  Every other combination is a syntax error.
- A **fan** `[r1: ... || r2: ...]` holds reactions between the same two nodes.
  A reaction can carry its own direction, which overrides the link arrows:
  `-[r3: -> 2.7.1.11 | ATP>ADP || r4: <- 3.1.3.11 | -H2O, +Pi]-`.
- **A reaction** is `[ID:] [direction] [~][EC[/EC]] [| cofactors] [^L|^R] [{attributes}]`.
  Every part is optional:

  | Body | Meaning |
  |---|---|
  | `[r1: 2.7.1.1 \| ATP>ADP]` | ID and EC number |
  | `[PGI \| ATP>ADP]` | a reaction given only by its ID |
  | `[2.7.1.1]` | EC number, no ID |
  | `[]` | neither |

  - IDs can be any name; no EC number or other primary data source is needed.
  - An ID, where given, is unique in the file. Motifs, effectors and
    `direction` statements can refer only to reactions with an ID.
  - The generator uses the BRENDA node ID.
  - `~?` is a spontaneous reaction.
  - Cofactors are `A>B` (A consumed, B produced), `-A` (consumed) or `+A`
    (produced), with optional coefficient and compartment: `-2*H@side1`.
  - `^L`/`^R` puts the drawn cofactors to the left or right of the link.
  - Attributes are free annotations `{key=value, ...}`. The generator writes:

    | Attribute | Meaning |
    |---|---|
    | `br` | The BRENDA reaction whose equation gave sides and direction. |
    | `ev` | Evidence for the direction, see [Directions](#directions). |
    | `rev` | Evidence for reversibility: `eq` or `map`. |
    | `match` | How drawn metabolites with other names were found in the equation. |

- `react [r8: 1.1.1.1] {A, B, C};` is a reaction whose main metabolites are
  **not** known to be split into two sides (and so has no direction).

**Drawing.** These only describe the picture, not biochemical cycles or hubs:

- `wrap serpentine period 4` after a chain;
- `ring cw|ccw { chains; chord chains; react ...; }`;
- `star HUB arc 210 order cw { spoke; ... }`;
- turn letters `-L[`, `-R[`, `-B[` (90° left, right, back relative to the heading).

**Directions set later.** `direction r1, r2 ->;` sets or changes the direction
of reactions after the layout statements. The direction is relative to the
order in which the reaction's sides are written, and `?` makes it unknown
again. The statements are applied in file order, after everything else, so a
later statement wins. The layout part of the file stays as it is; a
curator, or another source of directions, appends or edits these lines:

```
direction node8565 <->;           // glucose-6-phosphate isomerase is reversible
direction node8581, node11140 =>; // phosphofructokinases are irreversible
direction node2304 ?;             // direction withdrawn
```

**Motifs and regulation.**

```
motif cycle { oxaloacetate citrate cis_aconitate isocitrate ... };
motif substrate_cycle { FBP Fru_6_P } forward { node8421 } backward { node8581 };
motif shunt { isocitrate succinate } bypasses { isocitrate 2OG succinyl_CoA succinate };
motif repeat period=4 count=4 { node3128 node3164 ... };
motif carrier { NAD NADH } forward { ... } backward { ... };
motif transport { H@side1 H@side2 } via { node10039 };
effector ATP inhibits node8581;
effector AMP activates node8581, node11140;
```

The members and roles of the six kinds above are checked: the motif must be
supported by the reactions of the file (`ecmapper.pwl.motifs.check_motif`).
Other motif kinds are accepted unchecked, for hand curation, e.g.
`motif pathway name='upper glycolysis' { ... };`.

## Changes against PWL 1.1, by review point

The [review of PWL 1.1](#review-of-pwl-11) listed nine problems. Version 2.0
addresses points 2, 3, 5, 8 and 9.

### 2. Direction and reversibility

- **Language.** Links and reactions have a direction (`->`, `<-`, `<->`, `=>`,
  `<=`, or none for unknown). A fan can mix directions.
- **Generator.** For 2,277 of the 2,848 reactions (80%) the direction comes
  from the BRENDA equation and the drawn cofactor arrows, with the evidence in
  `ev` (see [Directions](#directions)).
- **Reversibility.** 29 reactions are marked reversible `<->`: BRENDA lists the
  equation in both directions for the EC number (`rev=eq`, 19), or the map draws
  the same reaction twice in opposite directions (`rev=map`, 10).
- **Irreversibility.** Neither source states it, so the generator never writes
  `=>`. The syntax is there for curated files.
- **Chains** are written in the direction of most of their reactions.
  Glycolysis now reads GAP → 1,3-BPG → 3-phosphoglycerate → … → pyruvate, and
  palmitate biosynthesis reads from acetoacetyl-ACP to palmitate.

### 3. Hyperedges: substrate and product sides

- **Groups** `(A + B)` now mean "same side of the reaction"; in 1.1 they meant
  "drawn on the same side of the enzyme".
- **Sides from equations.** 240 reactions have three or more drawn main
  metabolites. For 198 of them the sides come from the BRENDA equation, and for
  **115 of these the geometric split of 1.1 was wrong**: it put a substrate and
  a product into one group. Examples:
  - arginase: now `Arg -> (urea + Orn)`;
  - citrate synthase: now `(Ac_CoA + oxaloacetate) -> citrate`;
  - aldolase: now `FBP -> (DHAP + GAP)`.
- **Unknown sides.** The remaining 42 reactions are written as `react [...] {A, B, C};`
  and make no claim about sides.
- **Same-side pairs.** A drawn pair whose two metabolites are on the same side
  of the equation is written with `_` for the empty side (1 case). Without such
  evidence, two drawn metabolites are taken to be on opposite sides, as drawn.
- **Coefficients.** `2*Fd_ox>2*Fd_red`, `+2*CO2` and `n*H@side2` come from the
  equation where it was matched (point 4, coefficients only; there is no mass
  balance check).

### 5. Compound identity

- **Fixed ligand IDs.** The 1.1 generator read the key `ligand_id`, but the
  export calls it `brenda_ligand_id`. So `node_index.json` had no ligand ID, the
  `_L<ligand>` rule never applied, and the 1.1 check did not notice.
- **One token per compound, on the whole map** (`compounds.json`, 2,308
  compounds). Main metabolites and cofactors share tokens: acetyl-CoA is
  `Ac_CoA` everywhere (1.1 had `acetyl_CoA` 52 times and `Ac_CoA` 50 times), and
  2-oxoglutarate is `2OG`.
- **How compounds are formed** (`CompoundRegistry`). BRENDA names and ligand IDs
  are joined when one stands for the other in most of its occurrences.
  - O2 (ligand IDs 9 and 1180) and CoA (18 and 534) become one compound each
    (32 compounds have several ligand IDs).
  - A drawn metabolite belongs to the compound of its BRENDA name. So 16 drawn
    metabolites with a wrong ligand ID (e.g. NADH with the ID of H+) do not merge
    two compounds; they are listed under `reassigned`.
  - Some display names stay shared by two compounds (31 names, e.g. "DHAP" for
    `DHAP` and `glycerone_phosphate`) because BRENDA gives them different
    ligand IDs and BRENDA names. The two get different tokens.
- **Declarations** carry namespaced cross-references: `compound Arg brenda:123 kegg:C00062;`.
  - 1,793 compounds have a KEGG compound ID: all their names lead to the same
    single entry of KEGG's compound list.
  - The check verifies that every drawn metabolite's ligand ID is declared.
- **Token spelling.** Tokens have no `+`: a trailing charge is dropped
  (`NAD+` → `NAD`, label `'NAD+'`).

### 8. Motifs

`ecmapper.pwl.motifs` detects motifs on the parsed reactions, not on the
drawing, so it also works on hand-written files. Main metabolites are compared
by compound (drawn copies are the same compound). Cycles and shunts do not pass
through compounds that are cofactors elsewhere in the file (acetyl-CoA, NH3, ATP).

| Kind | Definition | Found in the BRENDA map |
|---|---|---|
| `cycle` | Directed cycle of at least 3 compounds, at most half of the steps only reversible | 29: citric acid, glyoxylate, urea, methionine (activated methyl), γ-glutamyl, vitamin K epoxide, sulfate assimilation, … |
| `substrate_cycle` | Different enzymes convert the same compounds in opposite directions, none known to be reversible | 41: PFK/FBPase, PEP/pyruvate (C4), xanthophyll cycle (violaxanthin/antheraxanthin/zeaxanthin), … |
| `shunt` | Route of ≤ 3 steps from u to v bypassing a route ≥ 2 steps longer, disjoint in between; short routes that are converted back are excluded | 41: glyoxylate shunt (isocitrate → succinate bypassing 2OG and succinyl-CoA), … |
| `repeat` | EC numbers repeat with a period of 2–8 steps at least three times along a chain whose directions agree | 1: fatty acid elongation, period 4 |
| `carrier` | Cofactor couple turned over in both directions within the subsystem | 80: NAD+/NADH, ATP/ADP, … |
| `transport` | Reaction with participants in different compartments | 4: H+ translocation in oxidative phosphorylation and photosynthesis |

**Other additions:**

- Compartments are new: `compartment`, `X@side1`, from "H+[side 2]" in the map
  and "Na+/in" in equations.
- Regulation can be written (`effector X inhibits|activates r1, r2;`). The BRENDA
  map has no regulation, so none is generated.
- The ring of the drawing is no longer presented as a cycle. The 1.1 rings were
  chosen by enclosed area and 53 of 103 were probably not cycles; `motif cycle`
  is the biochemical statement.

### 9. Grammar

- **Version line.** `pwl 2.0;` is required, and the major version is checked.
- **Lexing.** Tokens contain no `+`, so `+NAD`, `NAD>NADH` and `(H + A)` are
  unambiguous. Links are matched as a whole (`<-L[ ... ]->`), and the reaction
  body has its own small lexer.
- **Link arrows.** Only the six combinations listed above are valid; `<-[..]->`
  is defined (reversible).
- **Rings** may contain only chains, chords and `react`; `flow`, `port` and
  nested rings are errors.
- **Labels.** Labels are given in `compound` declarations, so the 1.1 rule
  "label on first mention" is gone.
- **Validation.** `parse` reports the following, each with line and column:
  - keywords used as compounds;
  - reaction IDs used twice;
  - motif, effector or `direction` references to unknown reactions;
  - motif references to unknown compounds;
  - cofactor side `^L` without cofactors;
  - copy numbers `#1`;
  - groups with one member;
  - bad coefficients;
  - a direction for a `react` statement (it has no sides).

  Declarations, reaction IDs, enzymes, coefficients and directions are not
  required (see [Language](#language)).
- **Reaction IDs** (`node8579:`) are part of the syntax, because motifs refer
  to reactions. As a result, the check no longer depends on node IDs in
  comments (part of point 6). Comments now hold only the enzyme names.

Not addressed:

- Point 1: drawing and biochemistry are still mixed in the same statements.
- Point 4: there is no mass balance.
- Point 6: provenance other than the IDs is still in comments and attributes.
- Point 7: there is no canonical order of statements.

## Directions

`ecmapper.pwl.brenda.orient` gives a drawn reaction its sides and direction in
these steps:

1. **Find the equations.** Take the BRENDA reaction IDs of the node. Each ID is
   looked up in BKMS-react; a suffix such as `BR1125-a` is dropped if needed.
   This gives BRENDA equations `A + B = C + D`, written in the direction of the
   natural substrates.
2. **Match the drawn metabolites to the equation terms.** Main metabolites and
   cofactors are matched in this order:
   1. by name with stereo descriptors;
   2. without stereo descriptors, if the match is unique;
   3. by a common KEGG compound in KEGG's synonym list (`match=kegg`);
   4. if exactly one drawn main metabolite and one term remain, they are taken
      to be the same compound (`match=single`, e.g. DHAP = "glycerone phosphate");
   5. if several remain, by clearly most similar name (`match=similar`, 35 reactions).

   Sugar ring forms count as the sugar ("beta-D-fructofuranose 6-phosphate" ~
   "D-fructose 6-phosphate"). The best-matching equation is used.
3. **Compare with the drawn cofactor arrows.** A consumed cofactor on the left
   or a produced one on the right supports the equation direction; the opposite
   supports the reverse.

   | `ev` | Meaning | Reactions |
   |---|---|---|
   | `cof+eq` | The drawn cofactors agree with the equation. | 1,766 |
   | `cof-eq` | The map draws the reaction against the equation; the drawn direction is used. | 265 |
   | `eq` | No cofactor to compare; the equation direction is used. | 246 |
   | `conflict` | The cofactors contradict each other; the sides are written but no direction. | 43 |

**Reactions without a direction (571):**

- 231 have no BRENDA reaction ID;
- 238 have IDs without an equation in BKMS-react (mostly `BS…` IDs);
- the rest have an equation whose terms could not all be matched, or the
  `conflict` case.

**The check** (`python -m ecmapper.pwl check`) verifies these directions
independently. For every reaction with a direction it looks up the equation
named in `br` and tests whether the substrates and products are on the sides
that `ev` implies. Current result: 2,277 checked, 0 contradictions.

## Check of the generated files

This check is stricter than `parse`. It compares the files with the BRENDA
export, so it expects what the generator writes: a node ID for every reaction
and declared compounds with their ligand IDs. Hand-written or edited files need
only pass `parse`.

`BRENDAmap/pwl/validation_report.json`:

- 193 files, 0 parse errors;
- every one of the 2,848 reactions appears once, under its node ID;
- EC label, main-metabolite nodes and consumed/produced cofactors match the
  export for every reaction (0 mismatches);
- every drawn metabolite's ligand ID is declared;
- 2,277 directions match their equations;
- all 196 motifs are supported by their reactions.

`tests/test_pwl.py` corrupts generated files and checks that each corruption is
found:

- a changed EC number;
- a changed cofactor;
- a changed ligand ID;
- moved main metabolites;
- a keyword used as compound;
- a reversed direction;
- a deleted reaction.

## Layouts from PWL

`python -m ecmapper.pwl layout` draws PWL subsystems in the style of the BRENDA
map. It uses networkx and Graphviz layout algorithms and the `draft_map`
renderer. With `--brenda-json` it also compares the drawings with the
originals. The code is in `src/ecmapper/pwl/layout.py` and `compare.py`, and the
`draft-map` extra is required (scipy, pygraphviz).

```bash
python -m ecmapper.pwl layout BRENDAmap/pwl OUT \
    --brenda-json BRENDAmap/BRENDA_map_reactions.json --brenda-svg BRENDAmap/BRENDA_map.svg
```

`OUT` gets:

- `layout/<id>.svg`: the chosen layout of each subsystem;
- `compare/<id>.svg`: the original and the chosen layout side by side;
- `metrics.csv` and `summary.md`.

All 193 subsystems take about 3 minutes. [`BRENDAmap/pwl_layout`](BRENDAmap/pwl_layout)
has the summary, the metrics and six examples:

- glycolysis;
- citric acid cycle;
- urea cycle;
- palmitate biosynthesis;
- methionine metabolism;
- isoprenoid biosynthesis.

### Procedure

1. **Converting PWL to a draft map.** `to_draftmap` creates:
   - one enzyme node per reaction;
   - one metabolite node per drawn compound (`token#n`), with the BRENDA node ID
     from `node_index.json`;
   - one cofactor node per cofactor of a reaction.

   Arrows follow the direction; reactions with an unknown direction get
   double-headed arrows, as on the BRENDA map. Box sizes follow the labels, as on
   the map.
2. **Hint layout** (`hints`). The drawing hints are replayed as a turtle walk:
   - **Steps.** Each chain, ring run and spoke starts with the `flow` heading. A
     link turns the heading by its turn letter and takes one step (330 units, the
     median distance of two metabolites of a reaction on the map). A chain is
     walked forwards and backwards from its first node that is already drawn.
     A statement is drawn only once it touches the drawn part: reactions with
     several metabolites, which the generator writes last, often join two chains.
   - **Groups and fans.** Groups `(A + B)` lie across the step. The reactions of
     a fan lie side by side.
   - **Stars.** Spokes are spread over `arc` degrees in `order`, far enough apart
     for their number.
   - **Rings.** A `ring` is the cycle of the drawing, so its metabolites go on a
     circle in the ring's sense, one step apart, starting in the direction of the
     first step. Other metabolites of the ring that react with two or more circle
     metabolites go inside the circle at their centre if those are far apart (a
     chord, e.g. glyoxylate in the citric acid cycle), and outside next to their
     edge if they are neighbours (a triangle, e.g. (3S)-citryl-CoA).
   - **Secondary metabolites.** A main metabolite that is drawn only as a partner
     in groups `(A + B)` of two or more reactions, never as a node of its own,
     is left out of the layout (e.g. acetyl-CoA in the citric acid cycle). As a
     group member it would tie distant reactions together. Unless the drawing
     falls apart without it, it is placed after the layout, next to its enzymes,
     at the free point whose links cross the fewest links (11 metabolites in 10
     subsystems).
   - **Wrap.** In a chain with `wrap serpentine`, a turn back moves to the next
     row.
   - **Taken places.** If a place is taken, the node goes to the free point of a
     half-step grid around it whose link crosses the fewest links drawn so far.
   - **Enzymes** go between the substrate and product sides.
3. **Refinement with Graphviz and networkx:**

   | Layout | What it does |
   |---|---|
   | `hints+neato` | Graphviz `neato` with every node pinned except the ends of crossing links, so crossings are resolved locally and the rest keeps its shape. A round is kept only if it removes crossings. |
   | `hints+spring` | networkx `spring_layout` started from the hint layout (25 iterations); links within 20° of horizontal or vertical are then straightened. |
   | `hints+spring+neato` | Both. |

4. **Baselines** that use only the graph and `flow`:
   - Graphviz `dot` (ranks along the flow), `neato`, `fdp`, `sfdp`;
   - networkx `kamada_kawai_layout` and `spring_layout`;
   - `circle` (`draft_map.relayout.circle_layout`), for subsystems with a cycle of
     at least 5 metabolites.
5. **Finishing every layout:**
   - it is scaled to the BRENDA spacing;
   - overlaps are pushed apart;
   - cofactors go next to their enzyme: on the side given by `^L`/`^R`, consumed
     ones upstream, produced ones downstream;
   - if cofactors overlap other nodes, the drawing is spread by 10%, 20% or 30%
     to make room (the first factor without overlaps is kept);
   - the same cofactor of nearby reactions (nodes at most 200 units apart, e.g. in
     a fan) becomes one node linked to all of them, if that adds no crossing and
     no overlap. Of 6,415 cofactor nodes, 260 are saved.
6. **Choice.** The hint-based layout with the fewest main crossings (then
   overlaps) is drawn. A baseline replaces it only if it has at most half as
   many main crossings and no overlapping nodes. Edges that still pass through
   nodes are routed around them (`draft_map.routes`).

### Comparison with the original drawings

Nodes are matched by BRENDA node ID; cofactors are not compared. The metrics
allow any translation, scale, rotation and reflection:

| Metric | Meaning |
|---|---|
| Procrustes disparity | 0 = same shape |
| directions | Share of enzyme–metabolite links pointing the same of the four directions, best of the 8 symmetries of the square |
| neighbours | Shared 3 nearest neighbours |
| orthogonal | Share of axis-aligned links |

All 193 subsystems (medians; crossings are totals):

| layout | Procrustes | directions | neighbours | orthogonal | main crossings | all crossings | chosen |
|---|---|---|---|---|---|---|---|
| hints | 0.232 | 0.727 | 0.644 | 0.593 | 460 | 1,563 | 84 |
| hints+neato | 0.239 | 0.727 | 0.636 | 0.585 | 334 | 1,336 | 10 |
| hints+spring | 0.322 | 0.605 | 0.667 | 0.462 | 211 | 903 | 59 |
| hints+spring+neato | 0.322 | 0.609 | 0.667 | 0.452 | 172 | 857 | 18 |
| dot | 0.613 | 0.500 | 0.485 | 0.643 | 518 | 3,555 | 1 |
| neato | 0.518 | 0.478 | 0.667 | 0.169 | 117 | 344 | 14 |
| fdp | 0.595 | 0.414 | 0.559 | 0.200 | 450 | 853 | 2 |
| sfdp | 0.535 | 0.458 | 0.619 | 0.250 | 322 | 984 | 0 |
| kamada_kawai | 0.533 | 0.462 | 0.641 | 0.185 | 313 | 857 | 4 |
| spring | 0.528 | 0.474 | 0.667 | 0.195 | 138 | 823 | 0 |
| circle (54 subsystems) | 0.499 | 0.402 | 0.493 | 0.229 | 390 | 1,428 | 1 |
| **chosen** | **0.299** | **0.647** | **0.659** | **0.500** | **75** | **554** | |
| original | 0 | 1 | 1 | 0.889 | 69 | 141 | |

Findings:

- **Similarity.** The hint layouts are much closer to the originals than any
  general-purpose layout: median Procrustes 0.23–0.32, against 0.50–0.61. The
  link directions agree for 73% of the links, against 40–50%.
- **Crossings.** The hint layout alone has many crossings (460). The pinned
  `neato` untangling removes a quarter of them without changing the shape;
  `spring` removes more but loosens the shape.
- **The chosen layouts** have 75 main crossings, fewer than Graphviz `neato`
  (117, the best general-purpose layout) and close to the originals (69). Their
  shape is far more similar (0.30 against 0.52).
- **Orthogonality.** The originals are orthogonal (89% of the links are
  axis-aligned). The hint layouts reach 59%, the general-purpose layouts 17–25%
  (`dot`, at 64%, has 518 crossings).

Effect of the later steps (chosen layouts, totals over all subsystems):

| step | main crossings | all crossings | overlapping nodes | cofactor nodes |
|---|---|---|---|---|
| first version | 114 | 895 | 87 | 6,415 |
| + secondary metabolites placed afterwards | 113 | 907 | 85 | 6,415 |
| + rings on circles | 88 | 875 | 64 | 6,415 |
| + spreading for cofactors | 75 | 573 | 10 | 6,415 |
| + shared cofactors | 75 | 554 | 12 | 6,155 |

The median Procrustes disparity stays at 0.30. Cofactors cause none of the
main crossings, because every layout places enzymes and main metabolites first.
They cause most of the other crossings. In the citric acid cycle these steps
bring the crossings from 4 main and 34 in all down to 3 and 16 (original: 2 and
2). The Procrustes disparity there falls from 0.11 to 0.07.

### Limitations

- **Lost geometry.** Turn letters know only four directions, and the steps of
  one subsystem all have the same length. Diagonal and long links of the original
  are lost. Rings are drawn as circles, including rings drawn as rectangles or
  diamonds (urea cycle); a branch that was drawn far away may be drawn next to
  its neighbours.
- **Spacing.** One step length is used for all subsystems; compactly drawn
  subsystems (glycolysis) come out larger than on the map.
- **Large subsystems with many hubs** (purine, lysine, methionine metabolism)
  keep crossings. The best layouts there are the refined ones, with a Procrustes
  disparity of 0.4–0.5.
- **Cofactors** are placed by rule and do not look for free directions. Their
  crossings (all crossings: 554 against 141) are mostly cofactor links in
  crowded places. Spreading is global: one crowded spot enlarges the whole
  subsystem.

## Limitations

- **Directions** are only as good as BRENDA's natural-substrate equations and
  the map's cofactor arrows. `cof-eq` (265) marks reactions that the map draws
  against BRENDA's direction. Many are probably reversible, but nothing in the
  data says so.
- **Name matching** by similarity (`match=similar`) is a heuristic. All 35
  reactions are marked, and a few pairs are doubtful (cellobiose ↔ cellotriose).
- **Motifs are candidates.**
  - Detection is limited to 3 cycles and 3 shunts per subsystem.
  - Shunts in dense subsystems can be routes a biochemist would not call a shunt.
  - Substrate cycles reflect the drawn directions: e.g. succinyl-CoA synthetase
    is drawn ATP-consuming in the citric acid cycle.
- **KEGG IDs** come from KEGG's compound list (`rest.kegg.jp`, free for
  academic use only); `--no-kegg` leaves them out.

## Review of PWL 1.1

The points the changes refer to:

1. Layout and biochemistry were mixed; the motifs were taken from the drawing.
2. Main reactions had no direction (all `-[..]-`) and no reversibility.
3. Groups `(A + B)` were split by drawing position, mixing substrates and products.
4. There were no stoichiometric coefficients.
5. Compounds had no identifiers: no ligand IDs because of a generator bug, and
   synonym tokens (`Ac_CoA`/`acetyl_CoA`, `FUM`/`fumarate`).
6. Provenance was in comments that the checker depended on.
7. There was no canonical form.
8. The motif vocabulary was missing repeats, shunts, substrate cycles,
   compartments and regulation.
9. Grammar loose ends: `+` in tokens versus `+`/`-` prefixes, `<-[..]->`
   undefined, statements allowed in rings, the label rule not enforced.
