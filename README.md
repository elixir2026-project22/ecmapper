# ecmapper

Map the reactions of a constraint-based model (loaded with
[cobrapy](https://github.com/opencobra/cobrapy)) to EC numbers, as a first step
towards linking them to [BRENDA](https://www.brenda-enzymes.org).

For every reaction:

1. EC numbers in the annotation (`ec-code`, also `ec`, `ec-number`, `brenda`) are used.
2. Otherwise EC numbers in COBRA style notes (`EC Number: ...`, `PROTEIN_CLASS: ...`) are used.
3. Otherwise the cross-references in the annotation are looked up in the resources below.
   EC numbers found in a resource for a cross-reference are collected with their
   provenance. If a resource lists further cross-references of the reaction
   (BiGGr does, e.g. `PGI` -> KEGG `R00771`, ModelSEED `rxn00558`, MetaCyc `PGLUCISOM-RXN`),
   these are looked up as well (one step, marked with `via`).

Boundary reactions (exchange, demand, sink) are skipped.

If several EC numbers are found, the ones to use are selected (`selected`,
`selection`). Transferred EC numbers are replaced by their successors and
deleted ones dropped (ExPASy ENZYME). Then EC numbers from the model are kept
(`annotation`, `notes`); a `single` EC number is used; if all BKMS-react
reactions that match the KEGG or MetaCyc reactions of the reaction have the
same EC number, that one is used (`bkms`); otherwise the metabolites of the
reaction are compared with the KEGG reactions of the EC numbers and the
clearly best matching one is used (`metabolites`). Optional tie-breakers
(`sources`, `organism`) can select among the remaining ones, otherwise they
are kept (`ambiguous`). See [DISAMBIGUATION.md](DISAMBIGUATION.md) for details
and an evaluation with iML1515 and iAF1260.

## Resources

| Name | Annotation keys | Data |
|---|---|---|
| `biggr` | `bigg.reaction` | [BiGGr](https://bigg.bio/data_access) reactions download; EC numbers (from the Rhea reference reaction) and cross-references |
| `kegg` | `kegg.reaction` | [KEGG REST](https://rest.kegg.jp/link/ec/rn) `link/ec/rn` (free for academic use only) |
| `modelseed` | `seed.reaction`, `seed` | [ModelSEEDDatabase](https://github.com/ModelSEED/ModelSEEDDatabase) `Unique_ModelSEED_Reaction_ECs.txt` |
| `rhea` | `rhea` | [Rhea](https://www.rhea-db.org) `rhea2ec.tsv`, directed reactions via `rhea-directions.tsv` |
| `bkms` | `metacyc.reaction`, `biocyc` | [BKMS-react](https://bkms.brenda-enzymes.org/download.php) `Reactions_BKMS.tar.gz` (MetaCyc 20.0) |

Reference data for the selection of EC numbers:

| Name | Data | Used for |
|---|---|---|
| `enzyme` | [ExPASy ENZYME](https://ftp.expasy.org/databases/enzyme/enzyme.dat) `enzyme.dat` | transferred and deleted EC numbers, Swiss-Prot organisms per EC number |
| `kegg-compounds` | [KEGG REST](https://rest.kegg.jp/link/cpd/rn) `link/cpd/rn` | compounds of the KEGG reactions of an EC number |
| `biggr-metabolites` | [BiGGr](https://bigg.bio/data_access) metabolites download | KEGG compounds and names of BiGG metabolites |
| `kegg-compound-names` | [KEGG REST](https://rest.kegg.jp/list/compound) `list/compound` | names of KEGG compounds (BRENDA map) |

Each resource is downloaded on first use, completely, and cached in `~/.cache/ecmapper`
(`$ECMAPPER_CACHE` or `$XDG_CACHE_HOME/ecmapper` if set). A resource that can
not be downloaded is skipped with a warning.

### Updating the resource data

The cached files are not updated automatically. Check whether newer versions are
available and download them with

```bash
python -m ecmapper --check-updates     # report only, the cache is not changed
python -m ecmapper --update            # download the files that changed
python -m ecmapper model.xml --update  # update, then map
```

or in Python with `check_updates()` / `update_sources()` (both take source names
or instances and `cache_dir`) or `ECMapper.check_updates()` / `ECMapper.update()`,
which also reload the updated sources. Each file is reported as `current`,
`outdated`, `missing`, `updated`, `downloaded`, `local` (a file passed with
`local_files`, not managed) or `error`.

For every downloaded file, its ETag, Last-Modified date and SHA-256 checksum are
stored in `<file>.meta.json`. The check sends conditional requests, so
unchanged files are not downloaded again. KEGG sends neither ETag nor
Last-Modified; its file (about 240 kB) is downloaded and compared by checksum.
A cached file that does not match its stored checksum (changed or damaged) is
also compared by content. `--refresh` (`refresh=True`, `update(force=True)`)
downloads all files again.

EC numbers are normalised: prefixes like `EC:` are removed and incomplete numbers
are padded (`2.4.1` -> `2.4.1.-`). Preliminary (`1.1.1.n1`) and BRENDA's
provisional numbers (`2.3.1.B43`) are kept, other strings (transporter
classification numbers like `3.A.1.7.-`, letter suffixes like `7.3.2.a`,
`SPONTANEOUS`) are dropped.

## Installation

```bash
python -m pip install ./ecmapper
```

## Usage

```python
import cobra
from ecmapper import ECMapper, summary, to_dataframe

model = cobra.io.read_sbml_model("models/bigg2-ecoli/e_coli_core.xml")
mapper = ECMapper()                     # or ECMapper(["kegg", "rhea"]), see options below
mappings = mapper.map_model(model)      # {reaction ID: ReactionECMapping}

mappings["PGI"].ec_numbers              # ['5.3.1.9'], all EC numbers found
mappings["PGI"].selected                # ['5.3.1.9'], selection: 'single'
mappings["PFK"].selected                # ['2.7.1.11'], selection: 'metabolites'
mappings["PFK"].scores                  # {'2.7.1.105': 0.5, '2.7.1.11': 1.0, '2.7.1.144': 0.2, '2.7.1.146': 0.71}
mappings["PGI"].evidence                # [ECEvidence(ec='5.3.1.9', source='kegg', identifier='R00771', via='bigg.reaction:PGI'), ...]
summary(mappings)                       # counts per origin, selection and source
to_dataframe(mappings)                  # one row per reaction
to_dataframe(mappings, long=True)       # one row per reaction, EC number and evidence
```

Options of `ECMapper`: `use_notes`, `always_query` (also query reactions that
already have EC numbers), `skip_boundary`, `follow_xrefs`, `bkms_consensus`,
`current_ec`, `compare_metabolites`, `min_score` (0.6), `min_margin` (0.1),
`tie_breakers` (`["sources"]`, `["organism"]` or both, in the order given),
`organism` (Swiss-Prot mnemonic, e.g. `"ECOLI"`), `cache_dir`, `refresh`.
`annotate_model(model, mappings)` writes the selected EC numbers found in the
resources to `reaction.annotation["ec-code"]`.

`evaluate(model, mapper)` removes the EC numbers of a model, maps it and
compares the selected EC numbers with the removed ones (one row per reaction,
`evaluation_summary` counts them per selection); `bigg_only=True` keeps only
the BiGG annotations.

Command line:

```bash
python -m ecmapper models/bigg2-ecoli/e_coli_core.xml -o ec_mapping.tsv
```

`--long` writes the long table, `--write-model out.xml` the annotated model,
`--evaluate` (with `--bigg-only`) evaluates the mapping with the model's own EC
numbers; the selection is set with `--no-current-ec`, `--no-bkms-consensus`,
`--no-metabolites`, `--min-score`, `--min-margin`, `--tie-breakers
organism,sources` and `--organism ECOLI`. See `python -m ecmapper --help` for
all options.

### Output

`origin` is `annotation` or `notes` (EC numbers from the model), `resources`
(found in the resources), `none` or `boundary`; `selection` is `annotation`,
`notes`, `single`, `bkms`, `metabolites`, `sources`, `organism`, `ambiguous`,
`none` or `boundary`. The wide table has the columns `reaction_id, origin,
selection, selected, ec_numbers, n_ec, sources, scores, xrefs`, the long table
`reaction_id, origin, selection, ec, selected, score, complete, source,
identifier, via`. `ec_numbers` are the EC numbers as found, `selected` are
current EC numbers.

## Finding model reactions on a BRENDA map

`ecmapper.brenda_map` finds the enzyme nodes of a BRENDA pathway map that
show the reactions of a model. The map is read from its reactions export
(`BRENDA_map_reactions.json`: enzyme nodes with EC numbers, BRENDA reaction
IDs, main metabolites and cofactors with their names, and layout).

```bash
python -m ecmapper model.xml --brenda-map BRENDA_map_reactions.json -o matches.tsv \
    --node-table nodes.tsv --svg BRENDA_map.svg --svg-out BRENDA_map_model.svg
```

```python
from ecmapper import ECMapper
from ecmapper.brenda_map import BrendaMap, MapMatcher, map_matches_to_dataframe

matcher = MapMatcher(BrendaMap.from_json("BRENDA_map_reactions.json"), ECMapper())
matches = matcher.match_model(model)       # {reaction ID: ReactionMapMatch}
matches["PGI"].status                      # 'on map'
matches["PGI"].nodes[0].node.subsystem     # 'glycolysis'
```

For each reaction, candidate nodes are found

- by reaction ID: the BRENDA reaction IDs of a node are aligned in BKMS-react
  with KEGG and MetaCyc reactions, which are compared with those of the model
  reaction (annotation and cross-references in BiGGr);
- by EC number: all EC numbers found for the reaction (model and resources,
  made current); incomplete EC numbers of nodes (`1.1.1.-`) include complete ones;
- by metabolites: at least two common main metabolites (not for transporters).

Each candidate is scored by comparing the metabolites of the reaction with the
main metabolites and cofactors of the node (weighted Jaccard index, cofactors
count half, protons are ignored). Map metabolites are identified only by
names, so metabolites are compared by KEGG compound IDs and by normalised
names: case, punctuation, stereo descriptors (D-/L-, (S)-, alpha-),
bis-/di-phosphate and membrane sides ("H+[side 1]") are ignored, "NAD(P)H"
stands for NADH and NADPH, generic quinones/quinols match ubiquinones and
menaquinones. Map metabolites get KEGG IDs from the names of KEGG compounds
and BiGG metabolites; if a node's BRENDA reaction is aligned with a KEGG
reaction and exactly one metabolite and one compound remain unmatched, they
are taken to be the same (for all nodes with this BRENDA ligand).

A candidate is accepted if its score is at least 0.6 (`--min-map-score`), or
0.8 if it was found only by metabolites (`--min-map-score-unlinked`), and at
most 0.1 lower than the best node of the reaction (`--max-map-gap`): a reaction
drawn in several subsystems is found on all of them, but not on nodes of
variants with other cofactors. A reaction is `on map` if a node is accepted,
`related` if nodes with a common reaction ID or EC number were found but none
matches (e.g. other chain lengths or cofactors), `not on map` otherwise.

Outputs: the matches table (`reaction_id, reaction_name, status, node_id,
subsystem, ec_label, enzyme_name, score, via, vicinity_of, accepted, x, y`;
`--all-candidates` lists all candidates), the node table (one row per map node
with its model reactions) and a copy of the SVG in which nodes with model
reactions are green (cyan if found only in the vicinity of other nodes, see
below) and the best nodes of related reactions yellow; the model reactions are
added to the tooltips of the nodes.

With iML1515, 676 of 2,375 reactions (without boundary reactions) are on the
map, on 776 of its 2,848 enzyme nodes in 116 subsystems; 499 are related and
1,200 not on the map, 788 of them transporters, which the map does not show.
All enzymatic reactions of e_coli_core except the cytochrome bd oxidase are
found.

### Searching the vicinity of found nodes

Optionally, the nodes next to found nodes are searched with lower thresholds:

```bash
python -m ecmapper model.xml --brenda-map BRENDA_map_reactions.json --map-vicinity propagate -o matches.tsv
```

```python
matcher = MapMatcher(brenda_map, ECMapper(), vicinity="propagate")
matches = matcher.match_model(model)
# or afterwards: matcher.search_vicinity(model.reactions, matches, "neighbours")
```

Two enzyme nodes are neighbours if they share a drawn main metabolite;
metabolite nodes connected by a line (e.g. the same metabolite in two
subsystems) count as one, currency metabolites (ATP, NAD, CO2, ...) do not
connect. For a node next to a node of model reaction R, the candidates are the
model reactions (except transporters) that have the shared metabolite, if R
has it, too: the connection on the map then also exists in the model. A
candidate is accepted if

- its score is at least 0.4 (`--min-map-score-vicinity`) if it has a common
  reaction ID or EC number with the node, otherwise at least 0.6
  (`--min-map-score-vicinity-unlinked`), and then only if the EC numbers do not
  contradict each other: if the reaction has complete EC numbers and the node
  has EC numbers, one pair must agree in the first three levels (this rejects
  e.g. phosphatases on kinase nodes with the same main metabolites);
- its score is at most 0.1 (`--max-map-gap`) lower than that of the best
  candidate of the node and that of the best node of the reaction.

There are two modes:

- `neighbours`: the nodes without model reactions next to all found nodes are
  searched once.
- `propagate`: starting only from unambiguous matches (reactions with a single
  accepted node that has a common reaction ID or EC number), all neighbours
  are searched, and the search continues from the newly accepted nodes and the
  found nodes it reaches, until nothing new is found. It finds chains of nodes
  and does not spread errors from reactions found on several nodes (e.g. on
  nodes of D- and L-enzymes, whose metabolite names are the same). Found nodes
  that are not reached are kept.

Nodes accepted this way have `vicinity` in `via`, and the match next to which
they were found in `vicinity_of` (`reaction@node`).

With iML1515, both modes add 47 node matches, 45 of them the same. `neighbours`
puts 39 more reactions on the map, on 43 more nodes; `propagate` puts 41 more
reactions on the map, on 42 more nodes. It starts from 544 of the 846 regular
matches. 34 or 36 of the new reactions were `related` before. Examples are:

- the molybdenum cofactor synthase, and, only with `propagate`, the
  molybdopterin cytidylyltransferase next to it;
- biotin synthase and carnitine-CoA ligase;
- steps of fatty acid synthesis, whose acyl-ACP names partly differ between
  the map and the model;
- pyruvate oxidase, drawn with cytochrome b1 instead of ubiquinone, and
  riboflavin synthase, drawn with an additional metabolite.

Checked by hand, 6 of the 47 matches of `neighbours` are doubtful or wrong.
Three are acyl-CoA dehydrogenases on acyl-CoA oxidase nodes; the resources
assign both EC numbers to them. The others are:

- the uroporphyrinogen methyltransferase on the node of the next methylation
  step;
- an acetate-CoA ligase on a propionate-CoA ligase node;
- a homocysteine methyltransferase on the betaine-homocysteine
  methyltransferase node.

`propagate` also accepts the D-lactaldehyde dehydrogenase on a methylglyoxal
reductase node that produces L-lactaldehyde. e_coli_core gets no new matches.

## Limitations

- The BiGGr namespace does not contain all legacy BiGG reactions, e.g. many
  fatty acid and lipid reactions of iAF1260 are missing; such reactions need
  other annotations.
- BKMS-react contains MetaCyc 20.0, newer MetaCyc reactions are not found. It
  also aligns many EC numbers to a MetaCyc reaction and is the least specific
  resource.
- EC numbers in CBMPy `KeyValueData` annotations are not read by cobrapy and
  therefore not used.
- The comparison of metabolites uses the compounds of KEGG reactions without
  their sides, and only for EC numbers with KEGG reactions.
- If a reaction is catalysed by enzymes with different EC numbers, it stays
  ambiguous unless a tie-breaker is enabled; the genes of the model are not
  used yet.
- BRENDA map metabolites are matched by name; names that are neither KEGG nor
  BiGG names (and can not be inferred from an aligned KEGG reaction) only
  match identical model metabolite names. Generic metabolites of the map
  ("an aldehyde") match no specific metabolite of the model. Synonyms are not
  recognised ("Crotonyl-ACP" on the map, "But-2-enoyl-[acyl-carrier protein]"
  in iML1515), so such nodes are not found even in the vicinity search, which
  needs the shared metabolite in both the map and the model.

## Tests

```bash
python -m pip install -e "./ecmapper[dev]"
cd ecmapper
pytest               # offline tests with excerpts of the resource files
pytest -m network    # download the real files and map e_coli_core
```
