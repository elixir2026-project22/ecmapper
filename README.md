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

Boundary reactions (exchange, demand, sink) are skipped. If several EC numbers are
found for a reaction, all of them are kept; disambiguation is left for later.

## Resources

| Name | Annotation keys | Data |
|---|---|---|
| `biggr` | `bigg.reaction` | [BiGGr](https://bigg.bio/data_access) reactions download; EC numbers (from the Rhea reference reaction) and cross-references |
| `kegg` | `kegg.reaction` | [KEGG REST](https://rest.kegg.jp/link/ec/rn) `link/ec/rn` (free for academic use only) |
| `modelseed` | `seed.reaction`, `seed` | [ModelSEEDDatabase](https://github.com/ModelSEED/ModelSEEDDatabase) `Unique_ModelSEED_Reaction_ECs.txt` |
| `rhea` | `rhea` | [Rhea](https://www.rhea-db.org) `rhea2ec.tsv`, directed reactions via `rhea-directions.tsv` |
| `bkms` | `metacyc.reaction`, `biocyc` | [BKMS-react](https://bkms.brenda-enzymes.org/download.php) `Reactions_BKMS.tar.gz` (MetaCyc 20.0) |

Each resource is downloaded once, completely, and cached in `~/.cache/ecmapper`
(`$ECMAPPER_CACHE` or `$XDG_CACHE_HOME/ecmapper` if set). Use `refresh=True` /
`--refresh` to download again. A resource that can not be downloaded is skipped
with a warning.

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

mappings["PGI"].ec_numbers              # ['5.3.1.9']
mappings["PGI"].evidence                # [ECEvidence(ec='5.3.1.9', source='kegg', identifier='R00771', via='bigg.reaction:PGI'), ...]
summary(mappings)                       # counts per origin and source
to_dataframe(mappings)                  # one row per reaction
to_dataframe(mappings, long=True)       # one row per reaction, EC number and evidence
```

Options of `ECMapper`: `use_notes`, `always_query` (also query reactions that
already have EC numbers), `skip_boundary`, `follow_xrefs`, `cache_dir`, `refresh`.
`annotate_model(model, mappings)` writes the EC numbers found in the resources to
`reaction.annotation["ec-code"]`.

Command line:

```bash
python -m ecmapper models/bigg2-ecoli/e_coli_core.xml -o ec_mapping.tsv
```

`--long` writes the long table, `--write-model out.xml` the annotated model;
see `python -m ecmapper --help` for all options.

### Output

`origin` is `annotation` or `notes` (EC numbers from the model), `resources`
(found in the resources), `none` or `boundary`. The wide table has the columns
`reaction_id, origin, ec_numbers, n_ec, sources, xrefs`, the long table
`reaction_id, origin, ec, complete, source, identifier, via`.

## Limitations

- The BiGGr namespace does not contain all legacy BiGG reactions, e.g. many
  fatty acid and lipid reactions of iAF1260 are missing; such reactions need
  other annotations.
- BKMS-react contains MetaCyc 20.0, newer MetaCyc reactions are not found. It
  also aligns many EC numbers to a MetaCyc reaction and is the least specific
  resource.
- EC numbers in CBMPy `KeyValueData` annotations are not read by cobrapy and
  therefore not used.
- EC numbers are not yet checked against BRENDA (transferred or deleted EC numbers).

## Tests

```bash
python -m pip install -e "./ecmapper[dev]"
cd ecmapper
pytest               # offline tests with excerpts of the resource files
pytest -m network    # download the real files and map e_coli_core
```
