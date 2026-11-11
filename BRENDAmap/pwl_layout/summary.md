# Layouts from PWL: comparison with the BRENDA map

193 subsystems. Median over subsystems (crossings: total). procrustes: lower is more similar; directions, neighbours, orthogonal: higher is more similar; aspect: |log| ratio of aspect ratios. See `ecmapper/pwl/compare.py`.

| method | subsystems | procrustes | directions | neighbours | orthogonal | aspect | main crossings | all crossings | overlapping nodes | best procrustes | chosen |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hints | 193 | 0.232 | 0.727 | 0.644 | 0.593 | 0.39 | 460 | 1563 | 83 | 78 | 84 |
| hints+neato | 193 | 0.239 | 0.727 | 0.636 | 0.585 | 0.41 | 334 | 1336 | 102 | 12 | 10 |
| hints+spring | 193 | 0.322 | 0.605 | 0.667 | 0.462 | 0.42 | 211 | 903 | 46 | 28 | 59 |
| hints+spring+neato | 193 | 0.322 | 0.609 | 0.667 | 0.452 | 0.42 | 172 | 857 | 36 | 3 | 18 |
| dot | 193 | 0.613 | 0.500 | 0.485 | 0.643 | 0.85 | 518 | 3555 | 418 | 18 | 1 |
| neato | 193 | 0.518 | 0.478 | 0.667 | 0.169 | 0.50 | 117 | 344 | 8 | 9 | 14 |
| fdp | 193 | 0.595 | 0.414 | 0.559 | 0.200 | 0.60 | 450 | 853 | 68 | 10 | 2 |
| sfdp | 193 | 0.535 | 0.458 | 0.619 | 0.250 | 0.58 | 322 | 984 | 27 | 9 | 0 |
| kamada_kawai | 193 | 0.533 | 0.462 | 0.641 | 0.185 | 0.50 | 313 | 857 | 94 | 12 | 4 |
| spring | 193 | 0.528 | 0.474 | 0.667 | 0.195 | 0.52 | 138 | 823 | 52 | 13 | 0 |
| circle | 54 | 0.499 | 0.402 | 0.493 | 0.229 | 0.31 | 390 | 1428 | 184 | 1 | 1 |
| chosen | 193 | 0.299 | 0.647 | 0.659 | 0.500 | 0.41 | 75 | 554 | 12 |  |  |

Original drawings: 69 main crossings, 141 crossings in all; median share of orthogonal main links 0.889.
Mean procrustes of the chosen layouts: 0.326.
