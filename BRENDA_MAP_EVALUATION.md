# Finding model reactions on the BRENDA map: evaluation on curated models

The matching of model reactions to the BRENDA map (see "Finding model
reactions on a BRENDA map" in the README) was developed with iML1515. This
evaluation applies it to further manually curated models of prokaryotes from
[BiGG](https://bigg.ucsd.edu/models). The script is
`scripts/evaluate_map_models.py`; it downloads the models.

## Selection of the models

The successor of BiGG, [bigg.bio](https://bigg.bio/models/), lists 4,881
models. Most are in multi-strain collections built automatically from a
reference model (LactoPanGEM, Ecoli_panGEMs, iPKB_Clostridiaceae, iECx55,
...). 34 models were built individually: 22 of them are prokaryotes other than
*E. coli* K-12. Neither site states curation status, so it was taken from the
publications. The criteria were:

- built or updated by manual, literature-based curation;
- experimental validation, e.g. growth phenotypes or gene essentiality;
- the most recent curated model of an organism;
- only cytosol, periplasm and extracellular space;
- different taxa.

| Model | Organism | Taxon | Publication | Curation and validation | Compartments |
|---|---|---|---|---|---|
| iYO844 | *Bacillus subtilis* 168 | Firmicutes, Bacilli | Oh et al. 2007, J Biol Chem ([PMID 17573341](https://pubmed.ncbi.nlm.nih.gov/17573341/)) | manual, iterative; knockout and substrate growth data | c, e |
| iYS854 | *Staphylococcus aureus* USA300 | Firmicutes, Bacilli | Seif et al. 2019, PLoS Comput Biol ([PMID 30625152](https://pubmed.ncbi.nlm.nih.gov/30625152/)) | manual, with protein structures; successor of iSB619; gene essentiality (RB-TnSeq) | c, e, cell wall (24 metabolites) |
| iCN900 | *Clostridioides difficile* 630 | Firmicutes, Clostridia | Norsigian et al. 2020, npj Syst Biol Appl ([doi](https://doi.org/10.1038/s41540-020-00151-9)) | manual; extends two earlier models | c, e |
| iJN1463 | *Pseudomonas putida* KT2440 | γ-Proteobacteria | Nogales et al. 2020, Environ Microbiol ([PMID 31657101](https://pubmed.ncbi.nlm.nih.gov/31657101/)) | thorough manual update of iJN746; growth screens, gene essentiality | c, p, e |
| iCN718 | *Acinetobacter baumannii* AYE | γ-Proteobacteria | Norsigian et al. 2018, Front Genet ([PMID 29692801](https://pubmed.ncbi.nlm.nih.gov/29692801/)) | manual update of AbyMBEL891; gene essentiality, phenotypes | c, e (one metabolite each in p and r) |
| STM_v1_0 | *Salmonella* Typhimurium LT2 | γ-Proteobacteria, Enterobacterales | Thiele et al. 2011, BMC Syst Biol ([PMID 21244678](https://pubmed.ncbi.nlm.nih.gov/21244678/)) | community consensus of two manual reconstructions | c, p, e |
| iYL1228 | *Klebsiella pneumoniae* MGH 78578 | γ-Proteobacteria, Enterobacterales | Liao et al. 2011, J Bacteriol ([PMID 21296962](https://pubmed.ncbi.nlm.nih.gov/21296962/)) | manual; Biolog phenotype arrays | c, p, e |
| iAF987 | *Geobacter metallireducens* GS-15 | δ-Proteobacteria | Feist et al. 2014, PLoS Comput Biol ([PMID 24762737](https://pubmed.ncbi.nlm.nih.gov/24762737/)) | reconciled earlier model plus manual additions; growth on electron donors | c, p, e |
| iIT341 | *Helicobacter pylori* 26695 | ε-Proteobacteria | Thiele et al. 2005, J Bacteriol ([PMID 16077130](https://pubmed.ncbi.nlm.nih.gov/16077130/)) | manual, with confidence scores; gene essentiality | c, e |
| iEK1008 | *Mycobacterium tuberculosis* H37Rv | Actinobacteria | Kavvas et al. 2018, BMC Syst Biol ([PMID 29499714](https://pubmed.ncbi.nlm.nih.gov/29499714/)) | standardised update of earlier *M. tuberculosis* models; gene essentiality | c, e |
| iAF692 | *Methanosarcina barkeri* Fusaro | Archaea | Feist et al. 2006, Mol Syst Biol | first manually curated archaeal model; later superseded by iMG746, which is not in BiGG | c, e |

Not selected:

| Models | Reason |
|---|---|
| iJB785, iJN678, iSynCJ816 (cyanobacteria) | thylakoid lumen, carboxysome and other internal compartments |
| iJN746, iNJ661, iSB619 | superseded by iJN1463, iEK1008 and iYS854 |
| iPC815 (*Yersinia pestis*) | largely derived from *E. coli* and *Salmonella* models by homology |
| iNF517 (*Lactococcus lactis*) | least documented curation and validation |
| *E. coli* and *Shigella* strain models (iY75_1357, iAPECO1_1312, iSF_1195, iECx55, ...) | derived automatically from iJO1366 |

bigg.bio also provides MEMOTE reports. They were computed on the converted
files, in which for most models all reactions count as blocked, so they were
not used.

## Model files: BiGG instead of bigg.bio

The model files of BiGG (bigg.ucsd.edu) are used. The files of bigg.bio carry
far fewer annotations. Fractions of the reactions and metabolites without
boundary reactions:

| | BiGG | bigg.bio |
|---|---|---|
| Reactions with EC numbers | 38–73% | 0% (iYO844 59%, from its notes) |
| Reactions with KEGG / MetaCyc / SEED reactions | 28–62% / 32–69% / 46–88% | 0% |
| Reactions with Rhea reactions | 32–66% | 0–43% |
| Metabolites with KEGG compounds | 57–90% | 0% (iYS854 66%, iJN1463 58%) |

bigg.bio also gives reactions and metabolites that it could not reconcile with
its universal namespace model-specific IDs, e.g. `__iYL1228__PGI`. ecmapper
reads the bigg.bio annotations and strips these prefixes, but fewer reactions
are still found than with the BiGG files:

| Model | iYO844 | iYS854 | iCN900 | iJN1463 | iCN718 | STM_v1_0 | iYL1228 | iAF987 | iIT341 | iEK1008 | iAF692 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| On map, bigg.bio | 437 | 470 | 535 | 642 | 446 | 595 | 625 | 411 | 281 | 462 | 301 |
| On map, BiGG | 462 | 480 | 554 | 659 | 468 | 641 | 664 | 449 | 293 | 488 | 321 |

## Coverage

Results with full annotation. *Share* is the fraction of the reactions on the
map, after subtracting transporters, which the map does not show. The last
two columns give the reactions and nodes found only by each vicinity mode.

| Model | Reactions | Transporters¹ | On map | Share | Related | Nodes | Subsystems | `neighbours` | `propagate` |
|---|---|---|---|---|---|---|---|---|---|
| iYO844 | 1,021 | 223 | 462 | 58% | 189 | 563 | 94 | 18 / 19 | 19 / 20 |
| iYS854 | 1,210 | 232 | 480 | 49% | 140 | 570 | 99 | 19 / 22 | 19 / 20 |
| iCN900 | 1,097 | 136 | 554 | 58% | 121 | 647 | 95 | 25 / 28 | 22 / 25 |
| iJN1463 | 2,546 | 774 | 659 | 37% | 369 | 757 | 116 | 43 / 48 | 42 / 46 |
| iCN718 | 887 | 113 | 468 | 60% | 138 | 548 | 95 | 20 / 20 | 19 / 18 |
| STM_v1_0 | 2,200 | 700 | 641 | 43% | 423 | 738 | 109 | 40 / 43 | 40 / 42 |
| iYL1228 | 1,971 | 593 | 664 | 48% | 428 | 766 | 113 | 38 / 43 | 39 / 41 |
| iAF987 | 1,188 | 211 | 449 | 46% | 234 | 553 | 100 | 34 / 41 | 34 / 41 |
| iIT341 | 477 | 74 | 293 | 73% | 60 | 354 | 73 | 10 / 11 | 10 / 10 |
| iEK1008 | 1,108 | 112 | 488 | 49% | 160 | 575 | 106 | 19 / 18 | 19 / 18 |
| iAF692 | 620 | 81 | 321 | 60% | 77 | 390 | 84 | 14 / 17 | 13 / 15 |
| iML1515 | 2,375 | 751 | 673 | 41% | 502 | 773 | 115 | 39 / 43 | 41 / 42 |

¹ among the reactions not on the map

## Hold-out of annotations

For each model, the annotation and notes of 30% of the reactions were
removed, three times with different random samples. Removed reactions can
then only be found by their metabolites. The reference is the result with
full annotation. The table counts the removed reactions that are on the map
in the reference:

- **Lost:** not found again on one of their nodes by the regular matching.
- **Recovered:** found again by the vicinity search.
- **Other node:** vicinity matches on nodes that are not among the reference nodes.

| Model | On map (reference) | Lost | Recovered, `neighbours` | Other node | Recovered, `propagate` | Other node |
|---|---|---|---|---|---|---|
| iYO844 | 419 | 45 | 28 | 1 | 26 | 1 |
| iYS854 | 424 | 35 | 21 | 2 | 20 | 0 |
| iCN900 | 507 | 45 | 22 | 0 | 23 | 0 |
| iJN1463 | 610 | 67 | 42 | 2 | 44 | 2 |
| iCN718 | 403 | 36 | 14 | 0 | 14 | 0 |
| STM_v1_0 | 587 | 72 | 37 | 2 | 40 | 2 |
| iYL1228 | 575 | 60 | 32 | 0 | 33 | 0 |
| iAF987 | 398 | 39 | 27 | 2 | 25 | 1 |
| iIT341 | 281 | 37 | 22 | 0 | 22 | 0 |
| iEK1008 | 442 | 53 | 29 | 0 | 29 | 0 |
| iAF692 | 284 | 36 | 17 | 2 | 17 | 2 |
| **Total** | **4,930** | **525** | **291 (55%)** | **11** | **293 (56%)** | **8** |

The vicinity search recovers more than half of the lost matches. For 3–4% of
its matches it chooses another node than the reference. This test only covers
candidates without a common EC number or reaction ID; that path needs a score
of at least 0.6.

## Manual check

Six vicinity matches of each model (`neighbours`, full annotation) were drawn
at random and checked by their reactions, metabolites and EC numbers:

| Verdict | Matches |
|---|---|
| Correct | 55 of 66 (83%) |
| Doubtful: same conversion, different cofactor, electron acceptor or enzyme variant | 7 (11%) |
| Wrong | 4 (6%) |

The doubtful matches include:

- NAD- and NADP-dependent variants of aldehyde, alcohol and
  methylene-tetrahydrofolate dehydrogenase;
- acetate-CoA ligase (ADP-forming) on a propionate-CoA ligase node;
- L-aspartate oxidase on aspartate dehydrogenase;
- precorrin-4 methyltransferase on cobalt-precorrin-4 methyltransferase.

The wrong matches are:

- **Amylomaltase on a maltotriohydrolase node** (two models). The resources
  assign a glycosidase EC number (3.2.1.20) to amylomaltase as well, so the EC
  check does not reject it.
- **A sulfite reductase reaction without EC number on sulfur
  oxygenase/reductase.** They share sulfite, hydrogen sulfide and water.
- **L-lactate dehydrogenase without EC number on D-lactate dehydrogenase.**
  D- and L-forms match each other by design.

One systematic error was not in the sample. Uroporphyrinogen methyltransferase
lands on the node of the next step, precorrin-2 C20-methyltransferase, in 8 of
the 11 models and in iML1515. S-adenosyl-L-methionine and
S-adenosyl-L-homocysteine count as main metabolites, which raises the score to
0.6, and both EC numbers are in sub-subclass 2.1.1.

## Changes made during this evaluation

- **bigg.bio files** (`biggr` annotations, model-specific IDs) can be used,
  see above.
- **Name normalisation.** `alpha-`/`beta-` were removed from all names, which
  merged beta-alanine with alanine, alpha- with beta-isopropylmalate, and
  beta-lysine with lysine. They are now only ignored as anomeric descriptors
  (`alpha-D-glucose`). Anomers and other stereoisomers still match each other,
  but each name also gets a key that keeps its stereo descriptors. The score
  matches each metabolite with at most one metabolite of the other side.
  Before, an anomerase or a glucose transport looked like a conversion of a
  single metabolite and matched e.g. glucose-6-phosphate isomerase or aldose
  1-epimerase.
- **Effect of the name fix.** On the bigg.bio files of the eleven models and
  iML1515, it removed wrong nodes from 42 reactions and added a correct node
  for one (alanyl-tRNA synthetase). The removed nodes belonged to aspartate
  decarboxylases, isopropylmalate dehydratase, lysine 2,3-aminomutase, glucose
  transports and glucose-6-phosphate anomerases.

## Conclusions

- **Coverage:** across the eleven models, 37–73% of the reactions that are not
  transporters are on the map. The vicinity search adds 10 to 43 reactions per
  model.
- **The two modes give nearly the same result.** `propagate` accepts slightly
  fewer matches on other nodes. It needs EC numbers or reaction IDs, because
  only matches through them count as unambiguous starting points.
- **Remaining sources of error:**
  - S-adenosyl-methionine and S-adenosyl-homocysteine count as main
    metabolites;
  - the resources assign EC numbers of other classes to some reactions;
  - D- and L-forms match each other.
