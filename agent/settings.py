"""Load deployment paths independently of source checkout location."""
from pathlib import Path
import json,os,sys

ROOT=Path(__file__).resolve().parent

def defaults():
    return dict(release_id='MKMGI-MEKG',model='deepseek-flash',model_display='DeepSeek',
        api_key_file=str(ROOT/'secrets/deepseek_api_key'),runtime_python=sys.executable,
        scientific_python=sys.executable,host='127.0.0.1',port=8788,daily_run_limit=None,
        hourly_ip_run_limit=None,max_active_jobs=2,max_rounds=4,max_tool_calls=22,max_api_calls=11,
        max_total_tokens=160000,request_deadline_seconds=480,max_evidence_records=None,
        literature_search_seconds=90,literature_cache_days=3,literature_page_size=100,
        agent_context_records_per_page=8,online_sources=['EuropePMC','Crossref','OpenAlex'],
        cache_days=7,send_prediction_values_to_provider=False,current_gene_summary_consent='',
        private_observations_enabled=False,prompt_version='mekg-agent-v3-question-aware-conversation-audited-relations')

def load_config():
    path=Path(os.environ.get('MEKG_AGENT_CONFIG',ROOT/'config.json'))
    c=defaults()
    if path.is_file():c.update(json.loads(path.read_text(encoding='utf-8')))
    else:
        # Enables import/inspection without starting a deployment or opening data.
        model=ROOT.parent/'data/model'
        c.update(model_root=str(model),prediction_file=str(model/'predictions/genomewide_time_gcn_consensus.tsv.gz'),
            kg_database=str(ROOT.parent/'data/mekg/v31/release/mekg.sqlite'),kg_sha256='',
            graph_source=str(ROOT.parent/'data/sources/data/base_graph_tf_supported.npz'),
            source_root=str(ROOT.parent/'data/sources'),
            public_atlas=str(ROOT.parent/'data/sources/source_data/positives_merged_B73v5.tsv'))
    for k in ['model_root','prediction_file','kg_database','api_key_file','graph_source','source_root','public_atlas']:
        c[k]=str(Path(c[k]).expanduser().resolve())
    return c
