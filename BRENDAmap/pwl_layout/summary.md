# Layouts from PWL: comparison with the BRENDA map

193 subsystems. Median over subsystems (crossings: total). procrustes: lower is more similar; directions, neighbours, orthogonal: higher is more similar; aspect: |log| ratio of aspect ratios. See `ecmapper/pwl/compare.py`.

| method | subsystems | procrustes | directions | neighbours | orthogonal | aspect | main crossings | all crossings | overlapping nodes | best procrustes | chosen |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hints | 193 | 0.236 | 0.741 | 0.633 | 0.619 | 0.40 | 543 | 2386 | 288 | 76 | 73 |
| hints+neato | 193 | 0.229 | 0.739 | 0.634 | 0.609 | 0.41 | 329 | 1831 | 248 | 14 | 20 |
| hints+spring | 193 | 0.336 | 0.611 | 0.667 | 0.471 | 0.45 | 241 | 1240 | 129 | 27 | 70 |
| hints+spring+neato | 193 | 0.336 | 0.611 | 0.667 | 0.447 | 0.42 | 173 | 1076 | 125 | 6 | 11 |
| dot | 193 | 0.613 | 0.500 | 0.485 | 0.643 | 0.83 | 517 | 4767 | 740 | 18 | 0 |
| neato | 193 | 0.518 | 0.474 | 0.667 | 0.170 | 0.49 | 120 | 553 | 65 | 9 | 9 |
| fdp | 193 | 0.595 | 0.415 | 0.559 | 0.200 | 0.60 | 457 | 1211 | 169 | 9 | 1 |
| sfdp | 193 | 0.548 | 0.462 | 0.620 | 0.256 | 0.60 | 320 | 1279 | 110 | 8 | 0 |
| kamada_kawai | 193 | 0.533 | 0.462 | 0.643 | 0.185 | 0.48 | 311 | 1147 | 190 | 10 | 6 |
| spring | 193 | 0.528 | 0.472 | 0.667 | 0.191 | 0.50 | 134 | 1192 | 135 | 13 | 3 |
| circle | 56 | 0.490 | 0.395 | 0.517 | 0.231 | 0.36 | 344 | 1646 | 314 | 3 | 0 |
| chosen | 193 | 0.299 | 0.645 | 0.667 | 0.500 | 0.43 | 105 | 882 | 85 |  |  |

Original drawings: 69 main crossings, 141 crossings in all; median share of orthogonal main links 0.889.
Mean procrustes of the chosen layouts: 0.323.
