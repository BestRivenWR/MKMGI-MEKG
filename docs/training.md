# Model training and evaluation

## Input contract

Begin with processed B73 v5 gene-level inputs. The final forward tensor has 191
columns. Its temporal encoder receives five measured modalities at eight DAPs
and their observation masks. The remaining inputs encode spatial summaries,
graph statistics, functional features, sequence properties and protein embeddings.
Preprocessing is fitted on eligible training rows within each split.

The input archive must have this layout:

```text
model-inputs/
  data/
    registry.tsv                    # original gene order, eligibility, family/fold assignments
    features.npz                    # main: N_DEV × 92, before graph/background transforms
    graph.npz                       # DEV graph
    background_raw.npz              # gene_ids, sequence25, esm320, availability
    annotation_tokens.npz           # sparse annotation matrix, full gene order
  source_data/matched_candidate_pools.tsv
  reference/
    family80_bg_concat_OOF.tsv.gz
    family80_ExtraTrees_background_OOF.tsv.gz
  universe/
    graph.npz                       # full graph including gene_ids
    sealed/registry.tsv              # non-DEV row mapping; labels are not read in deployment fitting
    sealed/features.npz
    source_data/universe_partition.tsv
```

This is a packaging contract for the existing artifacts; these files are not
included or reconstructed from FASTQs here. Preserve the archived row order,
family groups, folds, label flags, seed columns and control pools.

## Main experiment

```bash
python mkmgi/prepare.py --inputs /data/model-inputs --workspace /data/new-model-run
export MKMGI_WORKSPACE=/data/new-model-run
python mkmgi/freeze.py
python mkmgi/run.py
```

`prepare.py` checks alignment and family isolation before copying into a new
workspace. `freeze.py` hashes its inputs and source. `run.py` executes five outer
folds for each registered condition, reports paired evaluations and fits the
registered main model on DEV for deployment. It can be a substantial GPU job.

Single task:

```bash
python mkmgi/train.py --grouping 80 --fold 0 --condition time_attention
```

The original four inner folds, three seeds and six learning-rate/dropout
configurations are retained. Inner AP chooses settings/epochs. Outer labels
are used for evaluation. The main loss is class-balanced BCE over P and
eligible U; U are not experimentally verified negatives.

Registered conditions are the full model, mean temporal pooling, removal of
explicit DAP encoding, no-message control and a degree-preserving rewired graph.
`report.py` also reads the two archived OOF references listed above. `finalize.py`
needs the completed report and inner-selection outputs; it does not select a
new main model from outer scores.

## Reference methods

```bash
python -c "import sys;sys.path.insert(0,'evaluation');import benchmark_common as b;b.freeze()"
python evaluation/train_baselines.py --job classical --fold 0
python evaluation/train_baselines.py --job neural
python evaluation/run_ppr_float64.py
```

Run classical tasks for folds 0–4. Methods are logistic regression, RBF-SVM,
histogram gradient boosting, MLP, conventional GCN and personalized PageRank.
The PPR wrapper retains the study's float64 transition normalization fix.
`evaluation/config.json` preserves the original grids. Benchmark outputs default
to `$MKMGI_WORKSPACE/benchmarks`; override with `MKMGI_BENCHMARK_WORKSPACE`.

## Input and component ablations

```bash
python evaluation/run_ablation.py --help
```

This retains the published fixed-setting refits and their feature deletion
definitions. Conditions cover all omics, individual modalities, spatial inputs,
functional inputs, sequence/ESM features, graph information, token interactions
and GCN depth. Settings come from the paired main-model inner selections;
these fits are not a new independent tuning experiment. Output defaults to
`$MKMGI_WORKSPACE/ablations`.

## Aligned OOF reporting and querying

```bash
python evaluation/report.py \
  --table MKMGI=/data/new-model-run/predictions/family80_time_attention_OOF.tsv.gz \
  --table LogisticRegression=/data/new-model-run/benchmarks/predictions/family80_LogisticRegression_OOF.tsv.gz \
  --output /data/new-oof-report
python mkmgi/query.py \
  --table /data/new-model-run/predictions/genomewide_time_gcn_consensus.tsv.gz \
  --gene YOUR_B73_V5_ID
```

The OOF reporter requires complete, identically ordered gene tables. It exports
AUROC, AP and Precision@20/50/100, and uses the original paired within-fold family
bootstrap for AUROC/AP comparisons. A fold SD is a computational spread, not a
confidence interval. Benchmark precision measures known-gene recovery in the
specified candidate pool; it does not estimate novel-candidate experimental success.

New workspace paths and source hashes differ from the historical release.
This edition enables reproduction of the algorithms; exact published results
also require the archived data/configuration, environments and frozen outputs.
