# Build report

- spec: `specs/iml1515_central_carbon.layout.json` | curvature: 0.6 | packed: True
- reactions: 91 | nodes: 495 -> 510
- area before/after packing (M units^2): 79.7 -> 61.73 | canvas 7750 x 8587
- **validation errors: 0**  warnings: 0


## Layout-quality metrics (lower is better)

- line crossings between different reactions: **26**
- labels sitting on lines: **5**  | labels touching other nodes: **0**  | overlapping label pairs: **0**
- primary nodes closer than 70 units: 0
- bounding box: 7444 x 8284 (61.67 M units^2)

Worst crossing pairs: ATPS4rppxNTP1 (4); PFLxPOR5 (2); CITLxPPCK (2); PPCKxPPS (2); PPCxPPS (2); CITLxPPC (2); CITLxOAADC (2); L_LACD2xL_LACD3 (2); CITLxCS (1); MALDDHxPPCK (1); ME1xPPCK (1); ME2xPPCK (1)

Labels on lines: reaction PPS, metabolite pi, metabolite f6p, exit → Amino acids, exit → Amino acids

Detached pieces that draw a metabolite of the main map again (merge it, or keep the copy on purpose): PGM, GLYCK2 (shares 2pg_c, 3pg_c); FE3DCITR5 (shares cit_c)

## Score and targets

- score 106 = 26 crossings + 3*5 labels on lines + 20*2 detached pieces + 25 placement penalty  (lower is better)
- ok    validation errors: 0 (target: 0)
- MISS  crossings: 26 (target: at most 9)
- ok    labels on lines: 5 (target: at most 6)
- ok    labels on nodes: 0 (target: 0)
- ok    label overlaps: 0 (target: at most 2)
- ok    close nodes: 0 (target: 0)
- MISS  violated hard/strong directions: 5 (target: 0)
- look  detached pieces: 2 (target: 0)

Fix first: crossings (26): ATPS4rpp x NTP1 (4); PFL x POR5 (2); CITL x PPCK (2); PPCK x PPS (2); PPC x PPS (2); CITL x PPC (2); CITL x OAADC (2); L_LACD2 x L_LACD3 (2). Remedy: for the worst pairs: move one node, flip a `side`, change a lane `shift`, or split a node instead of drawing a long line; if many pairs involve two subsystems, move or mirror one block (`origin`).

## Placement of the subsystems (layout_directions catalogue)

- violated directions (hard/strong/medium): **5** of 20 checked | placement penalty: **25**
  - Lactate fermentation should be below Gluconeogenesis (strong; agreement, 7 sources (R 0.87) for their pathways)
  - Entner-Doudoroff pathway should be right of Glycolysis (strong; rule GPR004: ED is an alternative glucose-catabolic route that enters at hexose-pho)
  - Gluconeogenesis should be right of Glycolysis (strong; rule GPR002: Place the reverse-direction pathway beside glycolysis; shared reaction)
  - Pentose phosphate pathway should be above right of Glycolysis (strong; rule GPR003: Glucose-6-phosphate feeds the oxidative PPP and non-oxidative PPP retu)
  - 6pgc_c should be right of g6p_c (strong; textbook convention)

## Packing

- oxidative phosphorylation: ATP synthase / ATP turnover: moved by (-2077, -60)
- Fe(III)-citrate reduction: moved by (-13, 14)
- NADH dehydrogenase (demethylmenaquinone): moved by (-968, -613)
- NADH dehydrogenase (menaquinone): moved by (-1200, -393)
- NADH dehydrogenase (ubiquinone): moved by (289, -1223)
- ED: glycerate branch: moved by (-13, 444)
- Fe(III)-citrate formation (extracellular): moved by (25, -1186)
