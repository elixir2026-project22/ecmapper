# Cofactor combinations and cross-subsystem metabolites

This folder holds tables derived from the BRENDA pathway map
(`../BRENDA_map_reactions.json`, extracted from `../BRENDA_map.svg`; see
`../BRENDA_map_findings.md`) and from two tables taken from the Roche
*Biochemical Pathways* chart. Every table here is produced by a script in
`scripts/`. Run the scripts from the repository root.

## Terms

- **Primary metabolites** are the main nodes of the BRENDA map (`nodeHM` and
  `nodeSM`). The SVG draws their links with double-headed arrows, so it does
  not say which one is the substrate and which the product.
- **Secondary metabolites (cofactors)** are the cofactor nodes (`nodeCM`).
  Their arrows show direction, so each one is either *consumed* or *produced*.
- **Compound identity** on the BRENDA side is the BRENDA ligand ID behind each
  node. Different labels with the same ID (2OG and 2-oxoglutarate, Glu and
  L-glutamate, …) are treated as one compound and shown under the most frequent
  label. Six labels sit on an ID that belongs to a different compound, and these
  are kept apart by their label:
  - NADH (ID 12, which is H⁺)
  - dimethylallyl-PP (ID 17)
  - PPi (ID 157)
  - NAD(P)H (ID 5)
  - reduced acceptor (ID 3252)
  - malonyl-CoA (ID 811, which is malonyl-ACP)

## Overview

| File | Produced by | Content |
|---|---|---|
| `cofactor_and_crosslink_report.md` | `brenda_map_cofactors.py` | Readable report: all cofactor combinations and all cross-subsystem metabolites |
| `cofactor_combinations.csv` | `brenda_map_cofactors.py` | Cofactor combinations per reaction, with H₂O and H⁺ |
| `cofactor_combinations_without_H2O_H+.csv` | `brenda_map_cofactors.py` | The same combinations with H₂O and H⁺ also removed |
| `cross_subsystem_links.csv` | `brenda_map_cofactors.py` | Every metabolite–metabolite link between subsystems |
| `primary_metabolites_in_multiple_subsystems.csv` | `brenda_map_cofactors.py` | Primary metabolites occurring in ≥ 2 subsystems |
| `roche_boxed_compounds.csv` | input (Roche chart) | Boxed compounds on the Roche chart with box counts and grid squares |
| `roche_brenda_merged.csv` | `merge_roche_brenda.py` | Roche boxed compounds merged with the multi-subsystem BRENDA metabolites |
| `roche_brenda_merged_ids.csv` | `verify_roche_brenda_ids.py` | The same table with KEGG, ChEBI and BiGG IDs and a match check |
| `roche_brenda_merged_linked.csv`, `roche_brenda_merged_linked_ids.csv` | the same two scripts with `--linked-only` | The same, keeping only BRENDA compounds that have a cross-subsystem link |
| `roche_cofactor_combinations.csv` | input (Roche chart) | Cofactor combinations on the Roche chart, with example grid squares |
| `cofactor_combinations_without_H2O_H+_with_roche.csv` | `merge_roche_cofactors.py` | BRENDA and Roche cofactor combinations in one table |

To regenerate everything:

```bash
python scripts/brenda_map_cofactors.py
python scripts/merge_roche_brenda.py
python scripts/merge_roche_brenda.py --linked-only
python scripts/verify_roche_brenda_ids.py
python scripts/verify_roche_brenda_ids.py --linked-only
python scripts/merge_roche_cofactors.py
```

`verify_roche_brenda_ids.py` needs network access to rest.kegg.jp and
bigg.ucsd.edu the first time it runs. See its section below.

## 1. Cofactor combinations on the BRENDA map

Script: `scripts/brenda_map_cofactors.py`.

For each reaction (enzyme node) the script builds a pair
`(consumed cofactors, produced cofactors)`. Then:

1. **Excluded compounds.** O2, CO2, NH3, NH4+ (ammonium, under its own ligand
   ID) and H2O2 are removed from both sides. The list is the `EXCLUDED` set in
   the script. H₂O and H⁺ are kept in `cofactor_combinations.csv`.
2. **Reverse pairs merged.** A combination and its reverse, e.g. (NAD+, NADH)
   and (NADH, NAD+), are one entry, written as `(A, B)` in the direction drawn
   more often. `n_a_to_b` and `n_b_to_a` count the reactions drawn each way.
3. **Reactions without cofactors.** Reactions with no cofactor at all, and
   reactions whose only cofactors are excluded compounds, are not listed. The
   report gives their numbers.
4. **H₂O and H⁺ removed.** `cofactor_combinations_without_H2O_H+.csv` also
   drops H₂O and H⁺, then merges reverse pairs again. Combinations that become
   identical are combined; `n_variants_with_H2O_H+` says how many combinations
   from step 2 went into each row. Reactions left with no cofactor at all
   appear as `(∅, ∅)`.

`∅` means that side is empty. Each CSV row also lists the EC numbers and the
subsystems involved.

## 2. Metabolites connected across subsystems

Also produced by `scripts/brenda_map_cofactors.py`.

- **`cross_subsystem_links.csv`** has one row per metabolite–metabolite link.
  These are the 234 links with the SVG class `inter-pw-link`, plus the one
  unclassed link node1999 → node2003. The `same_ligand_id` column flags the
  links that join two different ligand IDs. These are mostly synonyms with
  separate IDs, but a few look like mistakes. They are listed in the report.
- **`primary_metabolites_in_multiple_subsystems.csv`** lists every primary
  compound (by ligand ID) that occurs in two or more subsystems, with or without
  a link. Its columns:
  - `n_subsystems`: the number of subsystems the compound occurs in.
  - `has_cross_subsystem_link`: whether any node of the compound is joined by
    a link.

  Cofactor nodes are not counted, so compounds such as SAM and SAH, which the
  map draws only as cofactors, do not appear here.

## 3. Roche boxed compounds merged with BRENDA metabolites

Script: `scripts/merge_roche_brenda.py`. Inputs are `roche_boxed_compounds.csv`
and `primary_metabolites_in_multiple_subsystems.csv`.

- **Matching by hand.** Roche names are matched to BRENDA ligand IDs with the
  table `ROCHE_TO_BRENDA` in the script. The names differ too much between the
  two sources for automatic matching. The matching was checked afterwards
  (section 4).
- **One Roche compound can match several BRENDA IDs.** This happens in two
  cases:
  - The same compound has two BRENDA IDs: L-malate and (S)-malate, the two DHAP
    IDs, crotonyl-CoA and crotonoyl-CoA, and the two 3-phosphoglycerate IDs.
  - An unspecified BRENDA form is matched to the only form the Roche chart
    shows: D-Glc 6-P → alpha-D-GLUCOSE 6-P, Glc 1-P → alpha-D-GLUCOSE 1-P,
    FBP → beta-D-FRUCTOSE 1,6-P2, and alpha-D-Glc → D-GLUCOSE.

  The BRENDA count is then the number of distinct subsystems across all the
  matched IDs.
- **Counts.**
  - `n_boxes_roche` is the box count on the Roche chart.
  - `n_subsystems_brenda` is the number of BRENDA subsystems.
  - `n_total` is their sum. Rows are sorted by `n_total`.
  - `source` is `both`, `Roche only` or `BRENDA only`.
- **`--linked-only`.** With this option only BRENDA compounds with
  `has_cross_subsystem_link = True` are used. Roche compounds whose BRENDA match
  has no link become `Roche only`. Where a Roche compound matched several IDs,
  only the IDs that have a link are kept.

## 4. Database identifiers and match check

Script: `scripts/verify_roche_brenda_ids.py`.

1. **Downloads.** The script fetches three files into a cache directory
   (default `~/.cache/ecmapper/compound_ids`, or the first argument). The
   downloads are not committed.
   - the KEGG compound list (`rest.kegg.jp/list/compound`)
   - KEGG's ChEBI cross-references (`rest.kegg.jp/conv/chebi/compound`)
   - BiGG's metabolite table (`bigg_models_metabolites.txt`)
2. **Each side is resolved to KEGG IDs separately, by name.**
   - Roche side: the Roche name.
   - BRENDA side: all labels and BRENDA names of the ligand ID on the map.

   Names are matched exactly against all KEGG names and synonyms, ignoring case,
   spaces, hyphens and commas. The script tries, in order:
   - the full name;
   - the name without a trailing bracket, or the abbreviation inside it;
   - expanded abbreviations (-P, -PP, P2, THF, Glc, Fru, Rib, …);
   - the name without an L-/D- prefix.

   Names KEGG does not know are given a KEGG name by hand in `MANUAL_NAMES`. The
   ID is still looked up from that name, so a wrong name fails rather than
   matching silently. Where a name hits several KEGG entries, the intended one is
   set in `PREFERRED`. The column `*_kegg_method` records how each ID was found.
3. **ChEBI IDs** come from KEGG's cross-references. EBI was not reachable from
   the environment where the tables were made, so ChEBI is not an independent
   check.
4. **BiGG IDs** are the universal BiGG metabolites whose database links name
   the KEGG ID. If there are none, the script uses BiGG entries that link one of
   the ChEBI IDs; `bigg_via` records which route was used. Some rows have several
   BiGG IDs, either legacy duplicates or questionable links in BiGG (e.g.
   H2O → `oh1`), so check those before using them.
5. **`match_check`** gives the result for each row:
   - `confirmed`: the Roche and BRENDA KEGG IDs agree.
   - `partly confirmed`: only some of the matched BRENDA IDs agree.
   - `different KEGG ID`: both sides resolved, but to different entries. These
     are the generic-to-specific matches from section 3.
   - `possible missed match`: a `Roche only` row and a `BRENDA only` row share a
     KEGG ID. No such row occurred.

## 5. Roche cofactor combinations merged with BRENDA combinations

Script: `scripts/merge_roche_cofactors.py`. Inputs are
`roche_cofactor_combinations.csv` and
`cofactor_combinations_without_H2O_H+.csv`.

- **Same filtering as the BRENDA table.** H₂O, H⁺, O2, CO2, NH3, NH4+ and H2O2
  are removed from the Roche rows too, and reverse pairs are merged.
- **BRENDA synonyms unified first.** Before matching, cofactor labels that name
  the same compound under different ligand IDs or spellings are mapped to one
  label (`BRENDA_SYNONYMS`). Examples: the acceptor/AH2 family, ferredoxin vs.
  ferredoxin [Fe-S] cluster, the ETF spellings, adrenodoxin vs. adrenal
  ferredoxin, UDP-Gal vs. UDP-alpha-D-galactose, and the DMAPP label that sits on
  the diphosphate ID. BRENDA rows that become identical are merged, which takes
  the 447 rows down to 422. The original rows are listed in `merged_brenda_rows`.
  Specific compounds are not merged into generic classes, so UQ/UQH2 stays
  separate from "a quinone"/"a quinol". This step changes only this merged table:
  `cofactor_combinations.csv` and `cofactor_combinations_without_H2O_H+.csv` keep
  the original labels.
- **Roche count.** `n_grid_squares_roche` is the number of grid squares listed
  in `example_grid_squares`. When several Roche rows fall on the same
  combination, their counts are summed. `n_total` is the number of BRENDA
  reactions plus the Roche grid squares, and rows are sorted by it.
- **Alternatives are split.**
  - "X or Y" is paired by position across the row, so
    "ADP or GDP → ATP or GTP" gives ADP → ATP and GDP → GTP.
  - A list in brackets, such as "NTPs (ATP, CTP, GTP, UTP)", gives one
    combination per entry.
  - Each split combination counts the row's grid squares.
  - Generic classes (UDP-sugar, acyl-CoA, NAD(P)+, Acceptor) match only BRENDA's
    own generic entries.
- **Names.** Roche names are translated to BRENDA cofactor labels with
  `ROCHE_NAMES`. One Roche name may stand for several BRENDA labels, e.g.
  reduced ferredoxin is both `Fd_red` and `reduced ferredoxin [Fe-S] cluster`.
- **Exact matches only are merged.** If synonym labels give several exact
  matches, the Roche count goes to the BRENDA row with the most reactions, and
  the others are listed in `also_exact_matches_brenda`.
- **Near matches.** For `Roche only` rows, `near_matches` lists up to five BRENDA
  combinations that contain the Roche combination or are contained in it, in
  either direction, with their reaction counts. These come mostly from the two
  charts splitting primary and cofactor metabolites differently. Example: Roche
  (CoA-SH, NAD+ → succinyl-CoA, NADH) versus BRENDA ((CoA, NAD+), NADH).

## Caveats

- **Counts measure drawn occurrences, not distinct chemistry.** The BRENDA map
  repeats modules across subsystems. For example, the 4-hydroxy-2-oxovalerate
  aldolase / acetaldehyde dehydrogenase pair is drawn in three aromatic
  degradation pathways. Compounds at the start or end of pathways therefore
  score high.
- **Roche counts come from the grid squares listed as examples**, which may not
  be every occurrence on the chart.
- **Hand-made matches and name mappings can contain errors.** They are all
  listed explicitly in the scripts so they can be reviewed.
