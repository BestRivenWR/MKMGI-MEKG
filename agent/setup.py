"""Initialize a fresh portal configuration and indexes from frozen local artifacts."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,re,sqlite3,time
from settings import ROOT,defaults

GENE=re.compile(r'^Zm00001eb\d{6}$')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def build_overlay(c):
    out=ROOT/'derived/evidence_index.sqlite'
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(exist_ok=True)
    db=sqlite3.connect(out)
    db.executescript('''
      CREATE TABLE model_edges(edge_id TEXT PRIMARY KEY,source TEXT,target TEXT,relation TEXT,weight REAL,confidence REAL,origin TEXT,source_ref TEXT);
      CREATE INDEX model_source ON model_edges(source,relation);
      CREATE INDEX model_target ON model_edges(target,relation);
      CREATE TABLE curated_links(gene_id TEXT,doi TEXT,source_url TEXT,aliases_json TEXT,phenotype TEXT,limitations TEXT,identity_basis TEXT,PRIMARY KEY(gene_id,doi,source_url));
      CREATE INDEX curated_gene ON curated_links(gene_id);
      CREATE TABLE search_aliases(gene_id TEXT,alias TEXT,source TEXT,PRIMARY KEY(gene_id,alias));
      CREATE INDEX alias_lookup ON search_aliases(alias COLLATE NOCASE);
      CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
    ''')
    src=Path(c['source_root'])/'source_data'
    for relation,name in [('PPI','ppi_base_edges.tsv.gz'),('TF_to_target','regulatory_tf_supported_edges.tsv.gz')]:
        with gzip.open(src/name,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                a,b=r['source'],r['target']
                if not GENE.fullmatch(a) or not GENE.fullmatch(b):raise ValueError('Invalid edge ID')
                if relation=='PPI':a,b=sorted((a,b))
                db.execute('INSERT OR IGNORE INTO model_edges VALUES (?,?,?,?,?,?,?,?)',
                    ('model:'+relation+':'+a+':'+b,a,b,relation,float(r['weight']),float(r['confidence']),
                     'inherited_ppi' if relation=='PPI' else 'model_regulatory_candidate',name))
    with open(c['public_atlas'],encoding='utf-8-sig') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            g=r.get('gene_id_v5','')
            if not GENE.fullmatch(g):continue
            aliases=set(re.split(r'[;|]',r.get('symbols','')+';'+r.get('old_ids','')+';'+r.get('official_symbol','')))-{'','-'}
            for a in aliases:
                if len(a)<=100:db.execute('INSERT OR IGNORE INTO search_aliases VALUES (?,?,?)',(g,a.strip(),'existing_public_literature_curation'))
            links=[(d.rstrip('.,)'), 'https://doi.org/'+d.rstrip('.,)')) for d in set(re.findall(r'10\.\d{4,9}/[^\s;|]+',r.get('primary_dois','')))]
            links += [('',u) for u in set(re.findall(r'https?://[^\s;|]+',r.get('primary_urls',''))) if 'doi.org/' not in u]
            for doi,url in links:
                db.execute('INSERT OR IGNORE INTO curated_links VALUES (?,?,?,?,?,?,?)',
                    (g,doi,url,json.dumps(sorted(aliases),ensure_ascii=False),r.get('phenotype_summary',''),r.get('limitations',''),r.get('identity_basis','')))
    db.execute('INSERT INTO metadata VALUES (?,?)',('source','Frozen graph edges and public literature curation'))
    db.commit();db.close();return out

def build_profiles(c):
    import numpy as np
    source=Path(c['model_root'])/'predictions/genomewide_mean_gates.npz'
    out=ROOT/'derived/model_profiles.sqlite'
    if out.exists():raise FileExistsError(out)
    with np.load(source,allow_pickle=False) as z:
        ids=z['gene_id'];w=z['temporal_pooling_weights'];g=z['relation_gates'];obs=z['temporal_observed'];dap=z['DAP'].tolist()
        if w.shape!=(len(ids),8) or g.shape!=(len(ids),2,4):raise ValueError('Profile shape mismatch')
        if dap!=[9,12,15,18,21,24,27,30]:raise ValueError('Unexpected DAP order')
        if not np.isfinite(w).all() or not np.isfinite(g).all():raise ValueError('Nonfinite profile')
        if not np.allclose(g.sum(2),1,atol=1e-5):raise ValueError('Relation gates do not sum to one')
        if not np.allclose(w.sum(1)[obs.any(1)],1,atol=1e-5) or not np.allclose(w[~obs],0):
            raise ValueError('Temporal missingness/pooling mismatch')
        db=sqlite3.connect(out);db.execute('CREATE TABLE profiles(gene_id TEXT PRIMARY KEY,profile_json TEXT NOT NULL)')
        for i,gene in enumerate(ids):
            payload=dict(DAP=dap,temporal_pooling_weights=w[i].tolist(),temporal_observed=obs[i].tolist(),
                relation_gates=g[i].tolist(),relation_names=['self','PPI','TF_to_target','target_to_TF'],
                semantics='Three-seed mean pooling/gating weights; not causal attribution or direct regulation.')
            db.execute('INSERT INTO profiles VALUES (?,?)',(str(gene),json.dumps(payload)))
        db.commit();db.close()
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--model-root',required=True,type=Path)
    p.add_argument('--kg-database',required=True,type=Path);p.add_argument('--source-root',required=True,type=Path)
    p.add_argument('--share-model-summary',action='store_true');p.add_argument('--api-key-file',type=Path)
    a=p.parse_args()
    if (ROOT/'config.json').exists() or (ROOT/'data_manifest.json').exists() or (ROOT/'derived').exists():
        raise FileExistsError('Initialize in a fresh checkout; existing deployment files are retained')
    c=defaults();c.update(model_root=str(a.model_root.resolve()),kg_database=str(a.kg_database.resolve()),
        source_root=str(a.source_root.resolve()),graph_source=str(a.source_root.resolve()/'data/base_graph_tf_supported.npz'),
        public_atlas=str(a.source_root.resolve()/'source_data/positives_merged_B73v5.tsv'))
    c['prediction_file']=str(a.model_root.resolve()/'predictions/genomewide_time_gcn_consensus.tsv.gz')
    if a.api_key_file:c['api_key_file']=str(a.api_key_file.resolve())
    if a.share_model_summary:c['current_gene_summary_consent']='score_rank_coverage_time_weights_only'
    files=[Path(c[k]) for k in ['prediction_file','kg_database','graph_source','public_atlas']]
    files += [a.source_root.resolve()/'source_data'/n for n in ['ppi_base_edges.tsv.gz','regulatory_tf_supported_edges.tsv.gz']]
    files += [a.model_root.resolve()/'predictions/genomewide_mean_gates.npz']
    for file in files:
        if not file.is_file():raise FileNotFoundError(file)
    before={str(f):sha(f) for f in files};c['kg_sha256']=before[c['kg_database']]
    (ROOT/'config.json').write_text(json.dumps(c,indent=2)+'\n')
    overlay=build_overlay(c);profiles=build_profiles(c)
    after={str(f):sha(f) for f in files}
    if before!=after:raise RuntimeError('Source artifacts changed during setup')
    manifest=dict(release_id=c['release_id'],created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        prediction_sha256=before[c['prediction_file']],kg_sha256=c['kg_sha256'],
        protected_files=[{'path':path,'sha256':h} for path,h in before.items()],
        derived_files=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in [overlay,profiles]])
    (ROOT/'data_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Ready: model/graph indexes and deployment manifest created.')

if __name__=='__main__':main()
