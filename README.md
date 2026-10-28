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
| `biggr-metabolites` | [BiGGr](https://bigg.bio/data_access) metabolites download | KEGG compounds of BiGG metabolites |

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

## Tests

```bash
python -m pip install -e "./ecmapper[dev]"
cd ecmapper
pytest               # offline tests with excerpts of the resource files
pytest -m network    # download the real files and map e_coli_core
```
