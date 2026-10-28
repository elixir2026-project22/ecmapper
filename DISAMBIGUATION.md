# Disambiguation of EC numbers

For most reactions several EC numbers are found (e_coli_core: 39 of 53
reactions). `ecmapper.disambiguation.Disambiguator` selects the ones to use
(`selected`) and records how (`selection`).

## Why several EC numbers are found

1. **Outdated EC numbers.** Sources still list transferred EC numbers
   (4.2.1.52 -> 4.3.3.7, 3.6.3.19 -> 7.5.2.1, ...). BRENDA uses the current ones.
2. **Loose alignments.** A cross-reference points to a reaction that a resource
   assigns to EC numbers of related or generic reactions (BKMS-react aligns
   KEGG R00754, ethanol + NAD+, with the 24 EC numbers of "a primary alcohol +
   NAD+ = an aldehyde + NADH + H+"; ModelSEED lists incomplete or neighbouring
   EC numbers). The reactions of these EC numbers differ from the model
   reaction, which a comparison of metabolites detects.
3. **One reaction, several enzymes.** The reaction itself is catalysed by
   enzymes with different EC numbers: acetate kinase (2.7.2.1) and propionate
   kinase (2.7.2.15) both catalyse ATP + acetate = ADP + acetyl phosphate,
   E. coli has both phosphoglycerate mutases (5.4.2.11 and 5.4.2.12).
   Metabolites can not separate these, only information on the organism or
   the genes of the model can.

## Selection

0. **Current EC numbers.** Transferred EC numbers are replaced by their
   successors, deleted ones dropped (ExPASy ENZYME, `enzyme.dat`: 1,314
   transferred and 178 deleted entries). Incomplete EC numbers that include a
   found complete one (`2.7.1.-` and `2.7.1.165`) are dropped.
1. `annotation` / `notes`: EC numbers that are already in the model are kept
   (made current; if they are all deleted, as they are).
2. `single`: only one EC number was found.
3. `bkms`: all BKMS-react reactions that match the KEGG or MetaCyc reactions
   of the reaction (from the model or found as cross-reference, e.g. in BiGGr)
   have the same EC number. BKMS-react's KEGG alignments are only used for
   this consensus, not as an EC source, because they are too generic (R00771,
   glucose-6-phosphate isomerase, is also aligned with 5.3.1.8 and 2.6.1.16).
4. `metabolites`: the metabolites of the reaction are compared with the KEGG
   reactions of each candidate EC number:
   - KEGG reactions of an EC number: `link/ec/rn` (KEGG source), including
     those of EC numbers transferred to it; their compounds: `link/cpd/rn`.
   - Compounds of the model reaction: every metabolite is mapped to KEGG
     compound IDs by its `kegg.compound` annotation and, via its
     `bigg.metabolite` annotation or its ID without compartment, by the BiGGr
     metabolites download (legacy BiGG IDs like `asn_L` match `asn__L`).
     Protons and metabolites without KEGG compound are ignored.
   - Score: weighted Jaccard index between the metabolites and the compounds of
     a KEGG reaction; a metabolite matches if one of its KEGG compounds is in
     the reaction; currency metabolites (ATP/ADP/AMP, NAD(P)(H), CoA, H2O,
     phosphate, diphosphate, CO2, NH3/NH4+, O2, ...) count 0.5, all others 1.
     Cofactors are kept because they separate e.g. NAD- and NADP-dependent
     enzymes. The score of an EC number is that of its best KEGG reaction. The
     sides of the reactions are not compared (KEGG REST has no bulk download
     of equations).
   - Decision: the EC number with the best score is selected if the score is
     at least `min_score` (0.6) and all other EC numbers score more than
     `min_margin` (0.1) lower. Otherwise, if the best score is at least
     `min_score`, only the EC numbers within `min_margin` of it are kept.
5. Tie-breakers, opt-in (`tie_breakers`, applied in the given order to the
   remaining EC numbers, each selects if one EC number remains):
   - `sources`: the EC numbers supported by most sources;
   - `organism`: the EC numbers of which ENZYME lists a Swiss-Prot protein of
     `organism` (mnemonic, e.g. `ECOLI`).
6. `ambiguous`: the remaining EC numbers are kept.

The scores are in `ReactionECMapping.scores` and in the output tables.

## Evaluation

`ecmapper.evaluate(model, mapper)` (`python -m ecmapper model.xml --evaluate`)
removes the EC numbers of a model, maps it and compares the selected EC numbers
with the removed ones (made current). `consistent`: all selected EC numbers
are in the reference (or in an incomplete reference EC number), `overlap`:
some are, `disjoint`: none is.

### iML1515 (all cross-references of the model, `ec-code` removed)

1,090 reactions have EC numbers, 264 of them several, taken from BiGG's
database links rather than curated per reaction.

| selection | reactions | consistent | overlap | disjoint |
|---|---|---|---|---|
| single | 470 | 470 | | |
| bkms | 10 | 7 | | 3 |
| metabolites | 253 | 243 | | 10 |
| ambiguous | 254 | 181 | 69 | 4 |
| not found | 103 | | | |

Of the reactions with one selected EC number and one reference EC number,
692 of 699 agree. Without the comparison of metabolites 507 reactions stay
ambiguous instead of 254; without current EC numbers 38 instead of 13 of the
selected EC numbers are `disjoint`. With BiGG IDs only (reactions and
metabolites), 862 reactions get EC numbers, 239 are selected by their
metabolites (229 consistent), 227 stay ambiguous.

Most "disjoint" selections are newer or more specific EC numbers than those of
the model: the bkms cases are GART (6.3.1.21, PurT) and the maltose ABC
transporter (7.5.2.1, formerly 3.6.3.19); by metabolites PDH, AKGDH and GLYCL
get the EC numbers of the whole enzyme systems (1.2.1.104, 1.2.1.105,
1.4.1.27, since 2018) instead of those of the subunits, IG3PS gets
imidazole-glycerol-phosphate synthase (4.3.2.10) and ACACT3r/ACACT5r the
thiolase for longer chains (2.3.1.16 instead of 2.3.1.9).

### iAF1260 (BiGG IDs only, EC numbers curated in the model's notes, 2007)

| selection | reactions | consistent | overlap | disjoint |
|---|---|---|---|---|
| single | 315 | 312 | | 3 |
| bkms | 4 | 1 | | 3 |
| metabolites | 180 | 164 | | 16 |
| ambiguous | 155 | 2 | 143 | 10 |
| not found | 303 | | | |

The not found reactions have IDs that are not in BiGGr. Several of the
disagreements correct the 2007 curation, e.g. GLYCLTDy uses NADPH and is
selected as 1.1.1.79 (curated: NAD-dependent 1.1.1.26), GLYCK2 and DHORD2 by
BKMS-react as 2.7.1.165 and 1.3.5.2.

### Tie-breakers

| model | tie-breakers | selected by them | consistent |
|---|---|---|---|
| iML1515 | organism, sources | 132 + 75 | 125 + 73 (96 %) |
| iAF1260 | organism, sources | 83 + 42 | 66 + 33 (79 %) |

Against the independently curated EC numbers of iAF1260, one in five
tie-breaker selections disagrees, so they are not enabled by default. The
better information for the case "one reaction, several enzymes" are the EC
numbers of the model's genes (GPR -> UniProt), which are not used yet.
Tie-breaking by enzyme specificity (fewest KEGG reactions) or by KEGG's "main"
reactions of an EC number was tested as well and agreed less often (43 % and
59 % on iAF1260).
