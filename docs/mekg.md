# MEKG construction and queries

MEKG separates gene identity, relations, supporting records and quantitative
observations. Source DOI/accession, material, stage, assay and coordinates remain
attached to each record. The current service reads MKMGI scores from a separate
prediction table.

## Data workspace

Builder code is in `mekg/_builders/`; data live in a separate workspace:

```text
mekg-data/
  base/                            # checked parent graph, including mekg.sqlite and source registry
  raw/maizegdb/                    # v5 annotations and v3/v4-to-v5 ID maps
  inputs/parent_predictions.tsv.gz  # historical graph-construction snapshot, not live website scores
  v3/raw/                          # source tables / public expression files
  v31/raw/                         # peaks, chain, cells, TSS and GWAS inputs
  audit_inputs/                    # historical overlap audit files
    model_manifest.json
    source_manifest.tsv
    labels_v9.tsv.gz
    audit_report.json
    feature_manifest.tsv
    existing_data_inventory.json
```

The archived parent and normalized source bundles are inputs. `assemble_base.py`
retains the parent assembly algorithm; run `--help` for its required tables.
No third-party database or annotation dump is bundled with the software.

`mekg/catalog/` retains the public acquisition specifications. The generic
downloader uses the file names, format checks and size limits of those records:

```bash
export MEKG_WORKSPACE=/data/mekg
python mekg/download.py mekg/catalog/data_sources.json
```

Default download destination is `v3/raw/`. Set `MEKG_DOWNLOAD_STAGE=v31` for a
v3.1 specification. Coordinate chains, final Seurat summaries and some other
v3.1 inputs must also be restored from the source-data bundle. The catalogs do
not substitute for a complete versioned input archive.

Install the graph readers in `requirements/mekg.txt`. If rebuilding the final
cell-group summaries from the author's Seurat object, use
`MEKG_WORKSPACE=/data/mekg Rscript mekg/_builders/v31/export_seurat.R` with the
required R packages. Existing archived summaries can be imported directly.

## Build order

```bash
python mekg/build.py public-v3 --workspace /data/mekg
python mekg/build.py v3-semantics --workspace /data/mekg
python mekg/build.py v3-audit --workspace /data/mekg
python mekg/build.py tf-peaks --workspace /data/mekg
python mekg/build.py myb-peaks --workspace /data/mekg
python mekg/build.py physical-contacts --workspace /data/mekg
python mekg/build.py gwas --workspace /data/mekg
python mekg/build.py cells --workspace /data/mekg
python mekg/build.py mutant-rna --workspace /data/mekg
python mekg/build.py tss --workspace /data/mekg
python mekg/build.py tf-mapping-audit --workspace /data/mekg
python mekg/build.py finalize --workspace /data/mekg
python mekg/build.py audit --workspace /data/mekg
```

Use a new workspace when rebuilding. Stage-specific immutable-release checks
are retained. The final database is `v31/release/mekg.sqlite`.

## Evidence classes

The sources cover genetic/functional reports, mutant molecular responses,
TF binding/peak associations, expression/cell context, transcription initiation,
physical contacts, trait associations and functional/homology context.

Peak-to-gene proximity, a PPI, differential expression and a causal genetic
test support different statements. GWAS/QTL associations retain study, trait,
model/year and original statistics. Correlated results are not combined into
independent evidence. Public status alone does not certify independence from
training/label curation. Historical overlap audits retain their original model
scope and do not certify independence for the new model.

Coordinates use explicit native/projected assembly fields. The chain converter
requires a unique, ungapped, whole-interval reciprocal mapping and rejects
ambiguous/gap-crossing intervals.

## Query

```bash
python mekg/query.py --gene YOUR_B73_V5_ID --db /data/mekg/v31/release/mekg.sqlite
```

This is a database query adapter. The web service separately loads current
MKMGI predictions; historical score annotations inside the graph are not the
authority for current gene scores.
