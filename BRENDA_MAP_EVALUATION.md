# Finding model reactions on the BRENDA map: evaluation on curated models

The matching of model reactions to the BRENDA map (see "Finding model
reactions on a BRENDA map" in the README) was developed with iML1515. This
evaluation applies it to further manually curated models of prokaryotes from
[bigg.bio](https://bigg.bio/models/). The script is
`scripts/evaluate_map_models.py`; it downloads the models.

## Selection of the models

bigg.bio lists 4,881 models. Most are in multi-strain collections built
automatically from a reference model (LactoPanGEM, Ecoli_panGEMs,
iPKB_Clostridiaceae, iECx55, ...). 34 models were built individually: 22 of
them are prokaryotes other than *E. coli* K-12. bigg.bio does not state
curation status, so it was taken from the publications. The criteria were:

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
| iAF692 | *Methanosarcina barkeri* Fusaro | Archaea | Feist et al. 2006, Mol Syst Biol | first manually curated archaeal model; later superseded by iMG746, which is not on bigg.bio | c, e |

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

## Model files from bigg.bio

The SBML files of bigg.bio carry fewer annotations than those of the old BiGG
or the publications:

- **Reactions** only have a `biggr` ID, some also a Rhea ID. Some models keep
  their EC numbers in the notes (`PROTEIN_CLASS`).
- **Model-specific IDs:** reactions and metabolites that bigg.bio could not
  reconcile with its universal namespace get IDs with the model as prefix,
  e.g. `__iYL1228__PGI` and `__iYL1228__g6p_c`. In iYL1228, iJN746, iAF692,
  iIT341 and iNJ661 this affects all metabolites and about 90% of the
  reactions. iEK1008 and iCN718 are affected almost as much.

ecmapper now reads `biggr` annotations and strips these prefixes. Otherwise
these models get no EC numbers or KEGG compounds from BiGGr, and only
metabolite names are compared:

| Model | On map, before | On map, with original IDs |
|---|---|---|
| iCN718 | 74 | 447 |
| iYL1228 | 265 | 628 |
| iEK1008 | 213 | 462 |
| iAF692 | 140 | 301 |
| iIT341 | 127 | 282 |

The fully annotated iML1515 file has 676 reactions on the map. The bigg.bio
version of iML1515 has 645 on the map, and 330 instead of 499 `related`.

## Coverage

Results with full annotation. *Share* is the fraction of the reactions on the
map, after subtracting transporters, which the map does not show. The last
two columns give the reactions and nodes found only by each vicinity mode.

| Model | Reactions | Transporters¹ | On map | Share | Related | Nodes | Subsystems | `neighbours` | `propagate` |
|---|---|---|---|---|---|---|---|---|---|
| iYO844 | 1,021 | 225 | 439 | 55% | 229 | 532 | 93 | 29 / 31 | 28 / 29 |
| iYS854 | 1,210 | 233 | 471 | 48% | 103 | 561 | 99 | 26 / 29 | 22 / 24 |
| iCN900 | 1,098 | 136 | 536 | 56% | 104 | 629 | 95 | 38 / 39 | 32 / 35 |
| iJN1463 | 2,548 | 776 | 643 | 36% | 217 | 740 | 117 | 48 / 55 | 44 / 50 |
| iCN718 | 886 | 112 | 447 | 58% | 106 | 526 | 95 | 26 / 28 | 25 / 26 |
| STM_v1_0 | 2,201 | 706 | 598 | 40% | 246 | 692 | 109 | 39 / 42 | 37 / 41 |
| iYL1228 | 1,971 | 598 | 628 | 46% | 253 | 727 | 114 | 39 / 45 | 36 / 40 |
| iAF987 | 1,188 | 212 | 411 | 42% | 101 | 514 | 100 | 33 / 41 | 28 / 36 |
| iIT341 | 477 | 73 | 282 | 70% | 39 | 344 | 74 | 13 / 13 | 13 / 14 |
| iEK1008 | 1,108 | 113 | 462 | 46% | 110 | 556 | 105 | 26 / 25 | 25 / 24 |
| iAF692 | 620 | 81 | 301 | 56% | 58 | 372 | 83 | 25 / 27 | 23 / 24 |
| iML1515 (fully annotated) | 2,375 | 751 | 676 | 42% | 499 | 776 | 116 | 39 / 43 | 41 / 42 |

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
| iYO844 | 396 | 50 | 24 | 1 | 23 | 1 |
| iYS854 | 422 | 33 | 18 | 3 | 17 | 1 |
| iCN900 | 489 | 42 | 23 | 0 | 22 | 0 |
| iJN1463 | 594 | 54 | 29 | 3 | 29 | 3 |
| iCN718 | 394 | 44 | 16 | 3 | 14 | 3 |
| STM_v1_0 | 536 | 54 | 25 | 2 | 27 | 1 |
| iYL1228 | 565 | 64 | 34 | 2 | 36 | 2 |
| iAF987 | 387 | 33 | 16 | 1 | 16 | 1 |
| iIT341 | 260 | 30 | 18 | 0 | 17 | 0 |
| iEK1008 | 428 | 56 | 25 | 0 | 27 | 0 |
| iAF692 | 281 | 32 | 8 | 2 | 8 | 1 |
| **Total** | **4,752** | **492** | **236 (48%)** | **17** | **236 (48%)** | **13** |

The vicinity search recovers about half of the lost matches. For about 6% of
its matches it chooses another node than the reference. This test only covers
candidates without a common EC number or reaction ID; that path needs a score
of at least 0.6.

## Manual check

Six vicinity matches of each model (`neighbours`, full annotation) were drawn
at random and checked by their reactions, metabolites and EC numbers:

| Verdict | Matches |
|---|---|
| Correct | 50 of 66 (76%) |
| Doubtful: same conversion, different cofactor, electron acceptor or enzyme variant | 10 (15%) |
| Wrong | 6 (9%) |

The doubtful matches include acyl-CoA dehydrogenases on acyl-CoA oxidase
nodes, NAD-dependent reactions on NADP-dependent nodes, and L-aspartate
oxidase on aspartate dehydrogenase.

The wrong matches have three causes:

- **Uroporphyrinogen methyltransferase** (three of the six). It lands on the
  node of the next step, precorrin-2 C20-methyltransferase, in 7 of the 11
  models and in iML1515. S-adenosyl-L-methionine and S-adenosyl-L-homocysteine
  count as main metabolites, which raises the score to 0.6, and both EC
  numbers are in sub-subclass 2.1.1.
- **The name normalisation merges isomers.** Removing `beta-` merges
  beta-alanine with alanine. Removing `alpha-`/`D-` merges glucose-6-phosphate
  anomers, so an anomerase looks like a conversion of glucose 6-phosphate and
  matches inositol-3-phosphate synthase. In the beta-alanine case, the model
  reaction in the bigg.bio file is also faulty: its product is octaprenyl
  diphosphate.
- **A phosphatase on a phosphorylase node:** glucose-1-phosphatase without EC
  number on cellobiose phosphorylase. Glucose 1-phosphate, glucose and
  phosphate match.

## Conclusions

- **Coverage:** across the eleven models, 36–70% of the reactions that are
  not transporters are on the map. The vicinity search adds 13 to 48
  reactions per model.
- **The two modes give nearly the same result.** `propagate` accepts slightly
  fewer matches on other nodes. It needs EC numbers or reaction IDs, because
  only matches through them count as unambiguous starting points. Without the
  original IDs, the bigg.bio files had none, and `propagate` found nothing.
- **Possible improvements:**
  - count S-adenosyl-methionine and S-adenosyl-homocysteine (and possibly
    other co-substrates) as currency metabolites;
  - keep `alpha-`/`beta-` where they are not stereo descriptors
    (beta-alanine);
  - reject unlinked candidates whose metabolites match only through anomers.
