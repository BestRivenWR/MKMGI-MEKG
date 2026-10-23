# Research agent and web portal

The final submission portal includes DeepSeek integration,
Europe PMC/Crossref/OpenAlex search, public evidence tools, source checking,
GWAS/QTL display and project/PPI/public networks. It defaults to English and
provides an **EN / 中文** language switch.

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

## Submission interface

Open `http://127.0.0.1:8788/` after startup. Search accepts a B73 v5 ID, a recorded
gene symbol or an alias; the example button is **O2**.

The gene panel shows priority score, rank, percentile, data coverage, temporal
pooling weights and stability. The workspace offers three views:

- **Research chat:** question-driven literature and MEKG retrieval with cited answers.
- **Evidence network:** published relations, project-derived candidates and PPI,
  with selectable edges and source details.
- **Literature & data:** publications and public experimental records; GWAS/QTL
  filters appear when the queried gene has matching records.

The language switch updates interface labels immediately and is remembered in
the browser. New conversations default to English unless the user has selected
Chinese. New answers follow the selected language; existing answers, original
source text and downloadable source records retain their original language.

The root route serves `agent/static/submission.html`, which loads versioned
`submission.js`, `submission.css` and the shared `network.js` renderer. The root
response requests cache revalidation. The older portal files are no longer used.

## Language-aware API

`POST /api/ask` accepts `language: "en"` or `language: "zh"`; omission defaults
to `"en"`. The language is included in the answer-cache identifier. For example:

```bash
curl http://127.0.0.1:8788/api/ask \
  -H 'Content-Type: application/json' \
  -H 'X-MEKG-Request: 1' \
  -d '{"gene_id":"Zm00001eb301570","question":"What is known about O2?","language":"en"}'
```

The response contains a `job_id`. Poll `GET /api/jobs/{job_id}` for progress and
the completed answer. `GET /api/health` reports the portal release, supported
interface languages and default language. Downloads and gene-query endpoints
remain available.

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

The 22 October submission source was synchronized on 23 October 2026. Public
source paths use the existing repository settings; `DEEPSEEK_API_KEY` remains
available as a runtime-only alternative to the protected key file. No credential
value is shipped. Runtime state, derived indexes, biological inputs and raw data
are excluded from Git tracking. Source provenance and code-only verification
are recorded in `provenance.json` and `submission_verification.json`.
