# Research agent and web portal

The active conversation service is retained, including DeepSeek integration,
Europe PMC/Crossref/OpenAlex search, public evidence tools, source checking,
GWAS/QTL display, project/PPI/public networks and the current portal frontend.

## Inputs

`agent/setup.py` uses these read-only artifacts:

| Location | Required files |
|---|---|
| Model root | `predictions/genomewide_time_gcn_consensus.tsv.gz`, `predictions/genomewide_mean_gates.npz` |
| KG | Complete MEKG v3.1 `mekg.sqlite` |
| Graph/curation source root | `data/base_graph_tf_supported.npz`; `source_data/ppi_base_edges.tsv.gz`, `regulatory_tf_supported_edges.tsv.gz`, `positives_merged_B73v5.tsv` |

The model table must use the published 44,303-gene schema and three seed columns.
Its priority score, percentile and SD are checked against the stored consensus.
The NPZ provides eight DAP pooling weights, observation flags and relation gates.

## Set up and run

Use a fresh checkout for a new deployment. Install the recorded agent environment
and NumPy for index construction. Setup creates the two local query indexes,
configuration and an artifact-integrity manifest:

```bash
python agent/setup.py --model-root /data/model \
  --kg-database /data/mekg/v31/release/mekg.sqlite \
  --source-root /data/graph-and-curation --share-model-summary
export DEEPSEEK_API_KEY=YOUR_PROVIDER_KEY
python -m uvicorn app:app --app-dir agent --host 127.0.0.1 --port 8788
```

Alternatively pass `--api-key-file /path/to/keyfile`. The API key is supplied at
runtime. `--share-model-summary` enables the current gene's score/rank, coverage
and time-weight summary in DeepSeek requests; raw matrices are not sent by these
tools. Public evidence and the user's question are used for the agent answer.

`agent/config.json` holds deployed paths and the provider model alias. For an
external configuration, set `MEKG_AGENT_CONFIG`. `uvicorn` controls the bind
address/port. A public installation can use its normal HTTPS/reverse-proxy setup.

## Query flow

1. Resolve a B73 v5 ID or recorded symbol/alias.
2. Read current model scores, variability, coverage and temporal pooling weights.
3. Retrieve local evidence and project/PPI edges, then search relevant literature.
4. Select tools according to the user's question and retain citation provenance.
5. Return a conversational answer, supporting evidence and an interactive network.

The graph distinguishes published experimental records from project-derived
candidate relations and inherited PPI. The browser shows trait associations in
the evidence list when available. Evidence checks constrain what the agent can
claim; searches cannot guarantee that every published result will be recovered.

`agent/static/` contains one active portal implementation. The old explanation
agent, older frontend and historical deployment scripts have been removed from
this repository edition. Runtime state, derived indexes and credentials are
excluded from Git tracking.
