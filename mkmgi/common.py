"""Portable workspace adapter for the published MKMGI implementation."""
from pathlib import Path
import datetime, hashlib, importlib.util, json, os, sys
import numpy as np

sys.dont_write_bytecode = True
CODE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('MKMGI_WORKSPACE', CODE.parent/'runs/mkmgi')).resolve()

def config():
    path=ROOT/'config.json'
    c=json.loads((path if path.is_file() else CODE/'config.json').read_text(encoding='utf-8'))
    c['previous_release']=str(ROOT/'universe')
    c['graph_source']=str(ROOT/'universe/graph.npz')
    return c

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def target(path):
    p=Path(path).resolve()
    if p!=ROOT and ROOT not in p.parents:raise ValueError('Outside model workspace: '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True)
    return p

def save(path,obj):
    p=target(path);tmp=p.with_name(p.name+'.tmp.'+str(os.getpid()))
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(p)

def boolcol(s):
    b=s.astype(str).str.lower().map({'true':True,'false':False,'1':True,'0':False})
    if not b.notna().all():raise ValueError('Invalid Boolean column')
    return b.to_numpy(bool)

def shared_module(name,source):
    alias='_mkmgi_'+source+'_'+name
    if alias not in sys.modules:
        p=CODE/'_shared'/(name+'.py')
        spec=importlib.util.spec_from_file_location(alias,p)
        module=importlib.util.module_from_spec(spec);sys.modules[alias]=module
        spec.loader.exec_module(module)
    return sys.modules[alias]

def core_module(name):return shared_module(name,'core')
def reference_module(name):return shared_module(name,'reference')

def load():
    import pandas as pd
    m=pd.read_csv(ROOT/'data/registry.tsv',sep='\t')
    for k in ['primary','broad','G1','excluded_U','ppr','tf','spatial_observed']:m[k]=boolcol(m[k])
    with np.load(ROOT/'data/features.npz',allow_pickle=False) as z:raw=z['main'].copy()
    with np.load(ROOT/'data/graph.npz',allow_pickle=False) as z:graph={k:z[k].copy() for k in z.files}
    return m,raw,graph

def verify_freeze(check_sources=False):
    f=json.loads((ROOT/'protocol/freeze.json').read_text(encoding='utf-8'))
    for p,h in f['release_hashes'].items():
        if sha(p)!=h:raise RuntimeError('Frozen workspace changed: '+p)
    if check_sources:
        for p,h in f['source_hashes'].items():
            if sha(p)!=h:raise RuntimeError('Source code changed: '+p)
    return f

def finite_random_ap(n,p):
    if not p:return float('nan')
    if n==1:return 1.
    h=np.sum(1./np.arange(1,n+1,dtype=float))
    return float(h/n+(p-1)*(n-h)/(n*(n-1)))
