## Layout recipes (numbers that worked)

Coordinates: x grows to the right, y downwards; a primary node has radius 20. Keep >= 300 between consecutive nodes of a
chain and >= 600 between subsystems. The reaction line runs from the mean of its primary substrates to the mean of its primary
products; cofactor dots (secondary metabolites) are placed automatically beside it.

- Linear chain (glycolysis, fermentations): nodes in a column 300 apart; forward reactions `side "R", dist 110`; the reverse
  reaction of the same two nodes as a lane on the other side: `shift [-110, 0], side "L"`.
- Branch (aldolase-like, one substrate -> two products): `axis [0,1]`, `mid` slightly below the substrate, products 300
  below and 300 apart; an isomerase between the two products runs horizontally.
- Hourglass (two substrates -> two products, e.g. transketolase): two nodes stacked/aligned on each side, 600-700 apart, the
  reaction midpoint falls between them automatically. Chain several hourglasses. If one metabolite is needed in two places use
  a `splits` entry (a second node) instead of a long line.
- Cycle (TCA): put the nodes on a `rings` entry, e.g. 9 nodes every 40 degrees, radius about 1050, citrate on top (-90),
  clockwise; every ring reaction gets `arc {"ring": name, "offset": 0}` and `side "out"`. Parallel reactions between the same two
  ring nodes: concentric lanes `offset +-150`, inner lanes `side "in"`, outer lanes `side "out"`. Shunts through the interior are
  ordinary reactions. If an extra substrate enters from outside add `arc_nodes`.
- Variants of one step with different cofactor pairs (PYK, PYK2, PYK3, PYK4): ONE lane, stack their cofactor pairs along it:
  `shift [0, +-75..225], curve false, gap 40`, with at least 650 between the two nodes.
- Duplicated reactions (two MDH): lanes 110-150 apart; cofactors of inner lanes inward, of outer lanes outward.
- Rows of repeated reactions (NADH dehydrogenases): one horizontal row per quinone couple, two lanes per row
  (`shift [0,-45] side "U"` and `shift [0,45] side "D"`), rows about 470 apart.
- A hub with several reactions between two nodes (ATP <-> ADP): three lanes at -190 / 0 / +190, the middle one with `gap 85`
  and `secdist 430`.
- A reaction whose substrates (or products) are all secondary metabolites needs both `axis` and `mid`.

## Interconnecting subsystems (global arrangement)

One connected map in a biologically plausible flow, linked through FEW shared hub metabolites, with few crossings.
- Merge only important metabolites (default list: glucose, G6P, F6P, G3P, PEP, pyruvate, acetyl-CoA, acetaldehyde, oxaloacetate,
  6-phosphogluconate; analogous hubs in other maps). Never merge cofactors (ATP, ADP, NAD(H), NADP(H), CoA, H2O, Pi, CO2, O2,
  quinones ...). Only primary nodes of the SAME metabolite can be merged. Use the id of the first node as `keep`.
- If a merge would force a line longer than ~2500 units or several new crossings, do not merge it.
- Central carbon metabolism default: glycolysis as the central vertical column; pyruvate -> acetyl-CoA below it; TCA ring below that
  (so glycolysis is above TCA); pentose phosphate pathway on one side of upper glycolysis, with the Entner-Doudoroff pathway
  branching at 6-phosphogluconate; the PEP/pyruvate/oxaloacetate block on the side between glycolysis and TCA; fermentations beside
  pyruvate/acetyl-CoA/acetaldehyde; respiratory chain / ATP hub elsewhere (disconnected pieces are packed automatically).
- Other maps: the largest linear pathway is the central column, cycles become rings after the pathway that feeds them, branches
  hang to the side of the metabolite they start from.
