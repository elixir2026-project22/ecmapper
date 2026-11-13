# Escher maps of the BRENDA subsystems with curved reactions

`ecmapper.brenda_escher` turns the drawing of the BRENDA map (`BRENDAmap/BRENDA_map_reactions.json`, the export of
`BRENDA_map.svg`) into Escher maps: one per subsystem and, with `--whole`, one of the whole map.

```bash
cd ecmapper
python -m ecmapper.brenda_escher BRENDAmap/BRENDA_map_reactions.json OUT --pwl BRENDAmap/pwl --png --whole
```

It writes `OUT/<subsystem>.escher.json`, with `--png` a preview of each map (matplotlib, not Escher itself) and
`OUT/quality.json`. `--subsystem ID` (repeatable) selects subsystems, `--no-labels` skips the label placement,
`--scale` (2) and `--radius` (450) change the size and the rounding,
`--no-circles` leaves the cycles as drawn. The module needs no packages; the label
placement is that of `ecmapper.draft_map.escher`, the previews need matplotlib.

## The problem

The BRENDA map draws a reaction as an enzyme box with straight links to its metabolite boxes. Escher draws it as one
line from the substrates through a midmarker to the products. The enzyme box is often not on the line between its
metabolites: of the 2,848 enzyme nodes, 1,324 sit at a corner or beside that line. A direct translation (midmarker
at the enzyme box, straight segments) gives zigzag lines; curves everywhere give meandering ones, where the BRENDA
map has rows and columns.

## Procedure

Positions are those of the map, times 2 (as in the draft maps: Escher's labels and nodes have fixed sizes). Apart
from steps 2 and 3, the main metabolite nodes stay where they are; midmarkers and cofactors move.

1. **Sides.** The PWL description of a reaction (`--pwl`) gives its substrate and product sides, its direction
   (2,277 reactions), coefficients and compound identifiers. Without PWL, or for a reaction without sides there,
   the main metabolites are cut into two sides at the two widest gaps between their directions seen from the enzyme
   box. A reaction without a direction gets the one in which its consumed cofactors are drawn nearer to its
   substrates, else it runs from left to right.
2. **Rows and columns.** Most links of the BRENDA map are horizontal or vertical, but many only nearly: the centres
   of the boxes differ by a few units. The ends of a link that is less than 16 units off the horizontal get the same
   y, along chains of such links, and likewise for vertical links (1,464 enzyme and metabolite nodes move).
3. **Cycles on circles.** The metabolites of a cycle of the drawing (a PWL `ring`) are put on a circle around their
   centre, with their mean distance as radius. The map around is warped with them: a point inside the polygon of
   the cycle keeps its direction and its share of the way to the polygon; a point outside moves as far as the
   polygon in its direction does, less with its distance, and not at all beyond two radii. Rows and columns that
   the warp has bent are straightened again (step 2 with 60 units, the metabolites of the circle stay).
   This is done for cycles of at least 7 metabolites that run once around their centre, are not drawn as a
   rectangle and have no metabolite more than 35% of the radius off the circle: the citric acid cycle,
   photosynthesis (Calvin cycle) and CO2 fixation in Crenarchaeota. `--no-circles` leaves the cycles as drawn.
   - A reaction between two neighbours on the circle is an arc of the circle.
   - Several reactions between the same two metabolites are arcs through both, 180 units apart at their
     midmarkers, on both sides of the circle, in the order of their enzyme boxes.
4. **Trunk.** The trunk of a reaction runs from the main substrate S to the main product P that are drawn most
   nearly opposite each other, seen from the enzyme box c. Further main metabolites join it at the multimarkers.
5. **Midmarker and axis** of the other reactions. a and b are the directions from c to S and P.
   - **Straight:** if the two links are at least 165° apart, the reaction is a straight line from S to P. The
     midmarker is the point of that line nearest to c (1,654 reactions).
   - **Corner:** otherwise the corner at c is rounded off. The two links are cut to the length R (the shorter one,
     at most `--radius`), and the reaction follows the parabola through the two cut points with c as control point.
     The midmarker is its apex, c + R/4 (a + b), the axis its tangent there, along b − a. A U-shaped reaction
     (substrates and products on the same side of the box) moves up to R/2 towards its metabolites. The midmarker
     moves only half as far, or not at all, if it would come closer than 110 units to another midmarker or a main
     metabolite.

   The multimarkers are on the axis, 30 units from the midmarker.
6. **Tangents at the main metabolites.** A segment of a straight reaction runs along its line, one on a circle
   along its arc. Any other leaves its metabolite towards the enzyme box, as the drawn link does. Two segments at a
   metabolite can get a common tangent, so that the pathway runs through it without a kink (at a metabolite with
   two segments if they are at least 75° apart, at one with more for pairs at least 165° apart). Straight and
   upright lines are not bent for it:
   - a segment of a straight reaction or of a circle keeps its direction; the other one takes it if the two are at
     least 135° apart;
   - of two other segments at least 135° apart, a horizontal or vertical one gives the direction; else both turn
     by the same angle;
   - at a corner (less than 135°) both turn only if neither is horizontal or vertical. Otherwise the corner stays
     at the metabolite, as in the BRENDA map.
7. **Main segments** are cubic Bezier curves from the metabolite (tangent of step 6) to the multimarker (tangent:
   the axis). Each handle is 40% of the distance between the two; the one at the multimarker is at most 150 units
   long, so a long link keeps its course and bends only near the reaction. Segments without a bend have no control
   points.
8. **Cofactors** get a secondary node of their own per reaction, beside the multimarker of their side: 50 units
   (+ 40 for each further place) outwards along the axis and 75 units beside it. Of these places on both sides of
   the axis, the one with the fewest nodes and lines around it and along its branch is taken; the side where BRENDA
   draws the cofactor and the places near the multimarker are preferred. The branch leaves the multimarker along
   the axis and bends towards the cofactor.
9. **Labels** go where they cover the fewest lines, nodes and labels: reaction labels at one of 12 places around
   the midmarker, metabolite labels at one of 8 places around their node.

Identifiers: a reaction has its EC label as BiGG ID and the enzyme name as name; a metabolite has its PWL compound
identifier (`Ac_CoA`; without PWL: made from the displayed name) and the BRENDA name. Inter-pathway links and
subsystem boxes are not drawn; the subsystem names are text labels.

## Result

All 193 subsystems (`quality.json`; the counts for the drawing are those of the straight links):

| | BRENDA drawing | Escher maps |
|---|---|---|
| Reactions | 2,848 | 2,848 |
| Corners at enzyme boxes | 1,324 | 0: 1,654 straight reactions, 48 arcs of circles, the others rounded |
| Main segments | 5,966 | 3,310 straight (3,152 exactly horizontal or vertical), 2,656 curved |
| Corners at metabolites with two links (75–170°) | 387 | 328, between straight or upright lines |
| Crossings of main links | 69 | 70 |
| Crossings of all segments | | 93 |
| Metabolite nodes closer than 35 units | | 0 |

`quality.json` gives the counts of the drawing after steps 2 and 3 (`drawing`), which differ slightly from the
table.

## Maps from the PWL descriptions alone

`python -m ecmapper.pwl escher` makes the same kind of map without the drawing of the BRENDA map: from the PWL
description of a subsystem (see [PWL.md](PWL.md)).

```bash
cd ecmapper
python -m ecmapper.pwl escher BRENDAmap/pwl OUT --png
```

It reads `subsystems/*.pwl` (or `*.pwl`) and, if there is one, `node_index.json`, and writes
`OUT/<subsystem>.escher.json`, with `--png` the previews, and `quality.json` (with the layout chosen for each
subsystem). `--subsystem ID` selects subsystems, `--methods hints,hints+neato` the layouts to try, `--no-circles`
draws the reactions of rings like the others. It needs the `draft-map` extra (networkx, scipy, pygraphviz), as
`python -m ecmapper.pwl layout`. In Python: `from_pwl(doc)` gives the reactions for `build`.

1. **Positions** of the enzyme and main metabolite nodes are those of the layout of PWL.md, "Layouts from PWL":
   the drawing hints (flow, turns, rings, stars, fans) are replayed and refined, and the candidate with the fewest
   crossings is taken. Its cofactor positions only say on which side of a reaction a cofactor goes.
2. **Sides, directions, coefficients and identifiers** are those of the description. A reaction without sides
   (`react`) gets the sides of the layout, one without a direction the one of step 1 above.
3. **Then the procedure above from step 2 on.** Rows and columns use a wider tolerance (80 units), because the
   refined layouts move nodes out of the rows and columns of the hints. The layout draws every `ring` as a circle,
   so rings of any size count as circles and their reactions are arcs.

All 193 subsystems (13 minutes, most of it for the layouts): 2,848 reactions, 97 crossings of main links (the
layouts as SVG: 75; the BRENDA drawing: 69), 256 crossings of all segments, 36 circles in 34 subsystems. Nine pairs
of cofactor nodes are closer than 35 units. Purine metabolism alone has 23 of the main crossings.

The maps look like those from the BRENDA drawing where the hint layout is kept (84 subsystems, e.g. glycolysis,
the citric acid cycle). The subsystems that needed a spring or a general-purpose layout (the large ones with hubs,
e.g. arginine or purine metabolism) are wider and have fewer horizontal and vertical lines, as their layouts do.

## Limitations and next steps

- The maps were not loaded in Escher: no Escher installation was at hand. The previews do not draw arrowheads.
- Label positions are chosen among a few places; where BRENDA draws several enzymes between two metabolites,
  some labels still cover lines.
- A cofactor node shared by two reactions in the drawing (159 nodes) is drawn once per reaction.
- Small cycles and cycles drawn as rectangles are not put on circles. Around a circle, lines that were
  horizontal or vertical can be slightly tilted by the warp.
- The draft maps of a model (`ecmapper.draft_map.escher`) still draw straight main lines. Steps 4 to 7 work on
  any enzyme node with its drawn metabolites, so `curve_reactions` can use them; the lanes of several model
  reactions of one enzyme node then need parallel curves.
