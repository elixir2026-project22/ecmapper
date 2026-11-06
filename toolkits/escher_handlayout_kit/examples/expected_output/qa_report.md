# Build report

- spec: `specs/spec_iml1515_central_carbon.py` | curvature: 0.6 | packed: True
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

## Packing

- oxidative phosphorylation: ATP synthase / ATP turnover: moved by (-2077, -60)
- Fe(III)-citrate reduction: moved by (-13, 14)
- NADH dehydrogenase (demethylmenaquinone): moved by (-968, -613)
- NADH dehydrogenase (menaquinone): moved by (-1200, -393)
- NADH dehydrogenase (ubiquinone): moved by (289, -1223)
- ED: glycerate branch: moved by (-13, 444)
- Fe(III)-citrate formation (extracellular): moved by (25, -1186)
