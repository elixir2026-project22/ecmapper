# Learnings: laying out metabolic pathways

These are lessons from generating drawings of the 193 BRENDA subsystems from their
PWL encodings and comparing them with the original BRENDA map.

- **Procedure, metrics and full results:** [PWL.md](PWL.md), section "Layouts from PWL".
- **Code:** `src/ecmapper/pwl/layout.py` and `compare.py`.
- **Output:** `BRENDAmap/pwl_layout/`.

The numbers below are totals or medians over all 193 subsystems.

- **Shape distance** is the Procrustes disparity: 0 means the same shape up to
  scaling, rotation and mirroring.
- **Directions** is the share of reaction links that point the same way as on the
  original map.

## 1. A topology-only description is not enough

General-purpose layouts see only the reaction graph. They come out far from the
hand-drawn maps: Graphviz `dot`, `neato`, `fdp`, `sfdp` and networkx `spring`,
`kamada_kawai` give a median shape distance of 0.50–0.61, and only 40–50% of
directions agree.

A few drawing hints recover most of the hand-drawn shape: 0.23–0.32 and 73% of
directions. The hints are:

- a flow direction;
- four turn letters;
- rings, stars, groups and fans;
- serpentine wrapping.

So a pathway description needs a small layer of layout intent on top of the
biochemistry, kept separate from it. This is why PWL treats directions,
coefficients and database IDs as optional: it is primarily a layout description.

## 2. Shape and few crossings pull against each other, so combine methods

| Approach | Result |
|---|---|
| Hint layout alone | Keeps the shape but has many crossings (460 main crossings). |
| Graphviz `neato` with every node pinned except the ends of crossing links | Removes crossings locally and keeps the shape. |
| networkx `spring` started from the hint layout | Removes more crossings but loosens the shape. |
| **Hints, then local refinement; a general-purpose layout only if it has at most half the crossings** | **75 main crossings (originals: 69), shape distance 0.30.** |

## 3. Lay out the backbone first, then add things in a fixed order

The main metabolites and enzymes are laid out alone. The rest is added in this order:

1. secondary main metabolites (e.g. acetyl-CoA), placed where their links cross least;
2. cofactors;
3. sharing of cofactor nodes between nearby reactions;
4. clean-up of overlaps.

Cofactors in the backbone layout change nothing. In the citric acid cycle test,
main-link crossings were identical with and without cofactors.

## 4. Hubs cause tangles, not cofactors

A metabolite drawn as one node shared by several reactions ties distant parts of
the drawing together. Acetyl-CoA in the citric acid cycle is the example: even the
BRENDA original gets both of its crossings from it. Placing such hubs only after
the cycle is laid out, near their enzymes, fixes this locally. Such hubs are rare
across the map: 11 metabolites in 10 subsystems.

## 5. Cycles should be drawn as circles, with a rule for their extra members

Four turn directions can't express diamond-shaped or rounded cycles. Laying each
`ring` out as a circle gives 25% fewer main crossings.

Other metabolites of the ring are placed by where their partners sit on the circle:

- **inside**, if they link distant members (a chord, e.g. glyoxylate);
- **outside**, if they link neighbouring members (a triangle, e.g.
  (3S)-citryl-CoA).

## 6. Cofactors decide whether a drawing looks tidy

Cofactors cause most of the remaining clutter: 554 crossings counting all links,
against 141 in the originals. Two simple steps help a lot:

- **Enlarging the drawing when cofactors overlap other nodes** (by up to 30%).
  This brought all crossings from 875 to 573 and overlaps from 64 to 10.
- **One shared cofactor node for nearby reactions**, applied only where it adds
  no crossing and no overlap. This saved 260 of 6,415 cofactor nodes.

Fixed placement rules (side `^L`/`^R`, consumed upstream, produced downstream) are
the remaining weakness. The obvious next step is cofactor placement that searches
for free directions around each enzyme. Enlarging the whole drawing also makes it
bigger everywhere because of one crowded spot; local spreading would avoid that.

## 7. Measure similarity in several ways

These measures each catch different failures:

- shape distance;
- link directions;
- shared nearest neighbours;
- the share of horizontal or vertical links;
- crossings and overlaps.

One number would hide trade-offs. For example, `dot` has orthogonal links but 518
crossings, and `spring` has few crossings but loses the shape.

## 8. What the encoding still loses

- **Diagonal steps, real link lengths and compactness** are not recorded, because
  PWL uses one step length per subsystem. Compactly drawn subsystems (glycolysis)
  come out larger.
- **Large subsystems with many hubs** (purine, lysine, methionine metabolism) stay
  hard: shape distance 0.4–0.5.
- **Orthogonality:** 89% of the original links are horizontal or vertical,
  against 50–59% generated. Better orthogonal routing would close part of that
  gap.

## In short

A good automatic pathway layout comes from four things together:

- layout intent recorded in the markup;
- a staged layout: backbone first, then hubs and cofactors;
- local, shape-preserving refinement instead of a global layout;
- explicit rules for cycles and for sharing cofactors.
