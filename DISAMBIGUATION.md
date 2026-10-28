# Disambiguation of EC numbers

Several EC numbers are found for most reactions (e_coli_core: 39 of 53
reactions). This document describes the selection that is implemented and a
proposed scheme to select EC numbers by comparing metabolites.

## Implemented selection

`ecmapper.disambiguation.select` sets `selected` and `selection` of every mapping:

1. `annotation` / `notes`: EC numbers that are already in the model are kept.
2. `single`: only one EC number was found. Incomplete EC numbers that are
   included in a found complete one (`2.7.1.-` and `2.7.1.165`) do not count.
3. `bkms`: all BKMS-react reactions that match the KEGG or MetaCyc reactions of
   the reaction (from the model or found as cross-reference, e.g. in BiGGr)
   have the same EC number. Example: GLYCK2, ModelSEED gives 2.7.1.-,
   2.7.1.165 and 2.7.1.31, BKMS-react aligns R08572/GKI-RXN only with 2.7.1.165.
4. `ambiguous`: all EC numbers are kept.

BKMS-react's KEGG alignments are used only for this consensus, not as an EC
source: BKMS aligns KEGG reactions with generic BRENDA reactions (R00754,
ethanol dehydrogenase, with the 24 EC numbers of "a primary alcohol + NAD+ = an
aldehyde + NADH + H+"; R00771, glucose-6-phosphate isomerase, also with 5.3.1.8
and 2.6.1.16) and would make unambiguous reactions ambiguous.

| model (BiGG IDs only) | single | bkms | ambiguous |
|---|---|---|---|
| e_coli_core | 14 | 2 | 37 |
| iAF1260 | 416 | 51 | 464 |

For iAF1260, 16 of the 19 BKMS selections with a curated EC number in the
model agree with it; in two of the three others (DHORD2, GLYCK2) the curated
EC number is the one that is outdated. In 42 of the 51 cases the other
candidates are transferred EC numbers.

## Why several EC numbers are found

1. **Outdated EC numbers.** Sources still list transferred EC numbers
   (4.2.1.52 -> 4.3.3.7, 4.99.1.1 -> 4.98.1.1, ...). BRENDA uses the current ones.
2. **Loose alignments.** A cross-reference points to a reaction that a resource
   assigns to EC numbers of related or generic reactions (BKMS-react rows of
   generic BRENDA reactions; ModelSEED: incomplete or neighbouring EC numbers).
   Here the reactions of the EC numbers differ from the model reaction, which
   a comparison of metabolites can detect.
3. **One reaction, several enzymes.** The reaction itself is catalysed by
   enzymes with different EC numbers, e.g. acetate kinase (2.7.2.1) and
   propionate kinase (2.7.2.15) both catalyse ATP + acetate = ADP + acetyl
   phosphate; KEGG links R00754 (ethanol + NAD+) to 1.1.1.1 and 1.1.1.71.
   Metabolites can not separate these, only information about the organism or
   the genes of the model can.

## Proposed scheme

**Step 0: current EC numbers.** Replace transferred EC numbers by their
successors and drop deleted ones, using the ExPASy ENZYME file
(`ftp.expasy.org/databases/enzyme/enzyme.dat`, 9.6 MB, 1,314 transferred and
178 deleted entries). This runs before all selection rules, so the BKMS
consensus and the comparison only see current EC numbers.

**Step 1: reference reactions per candidate EC number.** From KEGG: the KEGG
reactions of the EC number (`link/ec/rn`, already downloaded) and the
compounds of every KEGG reaction (`link/cpd/rn`, 52,932 links, one download).
KEGG compound IDs do not depend on the protonation state, which makes them a
good common namespace. Rhea (ChEBI participants) can be added later as a
second reference.

**Step 2: compounds of the model reaction.** Every metabolite (compartment
ignored, H+ dropped) is mapped to KEGG compound IDs, from its
`kegg.compound` annotation or, via `bigg.metabolite` or the metabolite ID,
from the BiGGr metabolites download (`bigg.bio/api/v3/download/metabolites`,
1,827 of 2,441 universal metabolites have KEGG IDs). Legacy BiGG IDs
(`asn_L` for `asn__L`) are matched with collapsed underscores.
`seed.compound`, `chebi` and `metanetx.chemical` annotations can be added the
same way via the ModelSEED and BiGGr cross-references.

**Step 3: score.** For each candidate EC number and each of its reference
reactions a weighted Jaccard index between the model's metabolites and the
reference compounds is calculated: a metabolite matches if one of its KEGG IDs
is a compound of the reaction; currency metabolites (ATP/ADP/AMP, NAD(P)(H),
CoA, H2O, phosphate, diphosphate, CO2, NH3, O2, ...) count 0.5, all others 1.
Cofactors are kept because they separate e.g. NAD- and NADP-dependent
enzymes. The score of an EC number is that of its best reference reaction.
(Substrate/product sides are not compared because the KEGG REST API provides
no bulk download of equations.)

**Step 4: decision.** The EC number with the highest score is selected
(`selection = "metabolites"`) if its score is at least 0.6 and at least 0.1
higher than that of the next one. If several EC numbers share the highest
score, the reaction stays ambiguous, but only these EC numbers are kept
(case 3 above). Below 0.6 (generic reference reactions) nothing is changed.

**Step 5 (optional): ties.** Possible tie-breakers, with the agreement measured on iAF1260:
- number of sources that support the EC number: 59 ties resolved, 71 % agree with the curation;
- Swiss-Prot entries of the model organism with the EC number (ENZYME `DR` lines): 7 resolved, 5 agree;
- better: EC numbers of the model's genes (GPR -> UniProt), not available offline.

Given the moderate agreement, ties should stay ambiguous by default and the
tie-breakers be opt-in.

## Prototype results

A prototype of steps 0 to 5 was evaluated on iAF1260, mapped from BiGG IDs only,
with the curated EC numbers in the model's notes (mapped to current EC
numbers) as reference. Of the 957 reactions with a curated EC number, 649 got
candidates:

| step | reactions | agree with the curated EC |
|---|---|---|
| single candidate | 312 | 312 |
| BKMS consensus | 4 | 1 (2 others: curated EC outdated) |
| metabolite comparison | 178 | 163 (92 %) |
| tie (same best score) | 138 | curated EC among the tied ones: 128 |
| score below 0.6 | 17 | (unchanged) |

Several of the 15 disagreements of the comparison look like corrections of the
2007 curation: GLYCLTDy uses NADPH and is selected as 1.1.1.79 (curated:
NAD-dependent 1.1.1.26), LCADi as lactaldehyde dehydrogenase 1.2.1.22
(curated: 1.2.1.21), ALAALAD as D-Ala-D-Ala dipeptidase 3.4.13.22.
Tie-breaking by enzyme specificity (fewest KEGG reactions) or by KEGG's
"main" reactions of an EC number was tested as well and agreed less often
(43 % and 59 %).
