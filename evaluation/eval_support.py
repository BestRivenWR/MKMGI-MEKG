"""Workspace/import adapters; the model and metric implementations are shared."""
from pathlib import Path
import hashlib,importlib,importlib.util,json,os,sys,types
import numpy as np

CODE=Path(__file__).resolve().parent
sys.path.insert(0,str(CODE.parent/'mkmgi'))
source_common=importlib.import_module('common')
source_model=importlib.import_module('model')
source_background=importlib.import_module('background')
graphops=source_common.core_module('graph_ops')
SOURCE=source_common.ROOT

def benchmark_config():return json.loads((CODE/'config.json').read_text())
def ablation_config():return json.loads((CODE/'ablation_config.json').read_text())
def benchmark_root():return Path(os.environ.get('MKMGI_BENCHMARK_WORKSPACE',SOURCE/'benchmarks')).resolve()
def ablation_root():return Path(os.environ.get('MKMGI_ABLATION_WORKSPACE',SOURCE/'ablations')).resolve()

def audit_freeze(root,cfg,code):
    source=source_common.verify_freeze(check_sources=True)
    path=root/'protocol/freeze.json';path.parent.mkdir(parents=True,exist_ok=True)
    inputs=dict(source['release_hashes']);inputs.update(source['source_hashes'])
    for p in (SOURCE/'runs').rglob('selection.json'):inputs[str(p)]=source_common.sha(p)
    for p in (SOURCE/'runs').rglob('model.pt'):inputs[str(p)]=source_common.sha(p)
    codepaths=list(CODE.glob('*.py'))+list(CODE.glob('*.json'))
    record={'source_hashes':inputs,'training_hashes':{str(p):source_common.sha(p) for p in codepaths},
            'configuration':cfg}
    if path.exists():
        if json.loads(path.read_text())!=record:raise RuntimeError('Evaluation specification changed')
    else:path.write_text(json.dumps(record,indent=2)+'\n')
    return record

def audit_verify(root):
    record=json.loads((root/'protocol/freeze.json').read_text())
    for section in ['source_hashes','training_hashes']:
        for p,h in record[section].items():
            if source_common.sha(p)!=h:raise RuntimeError('Frozen evaluation input changed: '+p)
    return {'files':len(record['source_hashes']),'changed':[]}

def load_model():
    m,raw,graph=source_common.load()
    bg=source_background.raw_for_indices(m.gene_index.to_numpy(int))
    if not np.array_equal(bg['ids'].astype(str),m.gene_id.to_numpy(str)):raise ValueError('Background row order')
    return m,raw,graph,bg

def matrices_numpy(raw,graph,training,inference,prep,bg,temporal_masks=True):
    gt=graphops.induced(graph,training);gi=graphops.induced(graph,inference)
    assert all(training[e].all() for k,e in gt.items() if k.endswith('_edge_index'))
    rt=np.column_stack([raw,graphops.summary(gt,len(raw))])
    ri=np.column_stack([raw,graphops.summary(gi,len(raw))])
    bx=source_background.apply_background(bg,prep['background'])
    xt=np.column_stack([source_model.apply_preprocess(rt,prep['core']),bx]).astype('float32')
    xi=np.column_stack([source_model.apply_preprocess(ri,prep['core']),bx]).astype('float32')
    if temporal_masks:xt[:,40:80]=raw[:,40:80];xi[:,40:80]=raw[:,40:80]
    assert xt.shape==xi.shape==(len(raw),191) and np.isfinite(xt).all() and np.isfinite(xi).all()
    return xt,xi,gt,gi

pv=types.SimpleNamespace(graphops=graphops,source_model=source_model,source_background=source_background,
    source_common=source_common,sha=source_common.sha,load=load_model,matrices_numpy=matrices_numpy)
