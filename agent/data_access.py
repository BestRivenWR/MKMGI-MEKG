"""Frozen, read-only data adapters. Model values never enter cloud tool outputs."""
import copy,csv,gzip,hashlib,json,math,re,sqlite3,statistics,time
from functools import lru_cache
from pathlib import Path

ROOT=Path(__file__).resolve().parent
from settings import load_config
CONFIG=load_config()
GENE_RE=re.compile(r'^Zm00001eb\d{6}$')
IDENT_RE=re.compile(r'^[A-Za-z0-9_.-]{2,64}$')
MANIFEST=ROOT/'data_manifest.json'

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def connection():
    c=sqlite3.connect(Path(CONFIG['kg_database']).as_uri()+'?mode=ro&immutable=1',uri=True,timeout=5)
    c.row_factory=sqlite3.Row
    return c

@lru_cache(maxsize=1)
def model_rows():
    p=Path(CONFIG['prediction_file']).resolve()
    assert p.parent.parent==Path(CONFIG['model_root']).resolve(),'Prediction must belong to fixed model root'
    rows={};seed_cols=['score_seed20261001','score_seed20261002','score_seed20261003']
    with gzip.open(p,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required=['gene_id','priority_score','priority_percentile','uncertainty_sd','runs']+seed_cols
        if not set(required)<=set(reader.fieldnames or []):raise RuntimeError('prediction_schema')
        for r in reader:
            gene=r['gene_id']
            if not GENE_RE.fullmatch(gene) or gene in rows:raise RuntimeError('prediction_id_duplicate')
            score=float(r['priority_score']);pct=float(r['priority_percentile']);sd=float(r['uncertainty_sd'])
            seeds=[float(r[k]) for k in seed_cols]
            if not all(math.isfinite(v) for v in [score,pct,sd]+seeds) or not 0<=score<=1 or not 0<=pct<=1 or sd<0:raise RuntimeError('prediction_range')
            if int(r['runs'])!=3 or abs(statistics.mean(seeds)-score)>2e-6 or abs(statistics.stdev(seeds)-sd)>2e-6:raise RuntimeError('prediction_consensus_mismatch')
            rows[gene]={'gene_id':gene,'priority_score':score,'priority_percentile':pct,'uncertainty_sd':sd,'runs':3,
                'seed_scores':seeds,'in_development_training':r.get('in_development_training')=='True',
                'in_development_positive':r.get('in_development_positive')=='True',
                'coverage':{'temporal_modalities_observed':int(r.get('temporal_modalities_observed') or 0),
                   'spatial_observed':r.get('spatial_observed')=='True','sequence_observed':r.get('sequence_observed')=='True',
                   'protein_observed':r.get('protein_observed')=='True'},
                'model_release':Path(CONFIG['model_root']).name,'score_semantics':'candidate priority, not phenotype probability'}
    ordered=sorted(rows.values(),key=lambda r:(-r['priority_score'],r['gene_id']))
    # Competition ranks reflect the stored score precision; percentile is untouched.
    prev=None;rank=0
    for i,r in enumerate(ordered,1):
        if r['priority_score']!=prev:rank=i;prev=r['priority_score']
        r['rank']=rank;r['model_gene_count']=len(rows)
    if len(rows)!=44303:raise RuntimeError('prediction_gene_count')
    return rows

@lru_cache(maxsize=1)
def profiles_connection_path():
    p=ROOT/'derived/model_profiles.sqlite'
    return p if p.is_file() else None

def model_profile(gene):
    p=profiles_connection_path()
    if not p:return None
    c=sqlite3.connect(p.as_uri()+'?mode=ro&immutable=1',uri=True)
    try:
        r=c.execute('SELECT profile_json FROM profiles WHERE gene_id=?',(gene,)).fetchone()
        return json.loads(r[0]) if r else None
    finally:c.close()

def resolve_gene(identifier):
    identifier=identifier.strip()
    if not IDENT_RE.fullmatch(identifier):return {'status':'invalid_identifier','candidate_ids':[]}
    rows=model_rows()
    if identifier in rows:return {'status':'ok','requested_gene_id':identifier,'gene_id':identifier,'id_mapping':'exact_B73_v5'}
    c=connection()
    try:
        found=set()
        for table,node_col in [('aliases','canonical_node_id'),('public_aliases','canonical_node_id')]:
            cols={r['name'] for r in c.execute('PRAGMA table_info('+table+')')}
            if node_col not in cols or 'alias_id' not in cols:continue
            for r in c.execute('SELECT '+node_col+' FROM '+table+' WHERE lower(alias_id)=lower(?) LIMIT 12',(identifier,)):
                val=r[0].removeprefix('gene:')
                if GENE_RE.fullmatch(val):found.add(val)
        for r in c.execute('SELECT gene_id FROM public_aliases WHERE lower(alias)=lower(?) LIMIT 20',(identifier,)):
            if GENE_RE.fullmatch(r[0]):found.add(r[0])
        for r in c.execute("SELECT canonical_id FROM nodes WHERE lower(node_type)='gene' AND lower(label)=lower(?)",(identifier,)):
            if GENE_RE.fullmatch(r[0]):found.add(r[0])
        from public_repository import overlay
        extra=overlay()
        try:
            for r in extra.execute('SELECT gene_id FROM search_aliases WHERE alias=? COLLATE NOCASE',(identifier,)):
                if GENE_RE.fullmatch(r[0]):found.add(r[0])
        finally:extra.close()
        if GENE_RE.fullmatch(identifier) and c.execute('SELECT 1 FROM nodes WHERE node_id=?',('gene:'+identifier,)).fetchone():found.add(identifier)
        if len(found)==1:return {'status':'ok','requested_gene_id':identifier,'gene_id':next(iter(found)),'id_mapping':'recorded_alias'}
        return {'status':'ambiguous' if found else 'not_found','requested_gene_id':identifier,'candidate_ids':sorted(found)}
    finally:c.close()

@lru_cache(maxsize=1024)
def node_info(gene):
    c=connection()
    try:
        r=c.execute('SELECT * FROM nodes WHERE node_id=?',('gene:'+gene,)).fetchone()
        aliases=[dict(x) for x in c.execute("SELECT alias_type,alias,source_node_id,mapping_status FROM public_aliases WHERE gene_id=? AND alias_type IN ('gene_symbol','gene_synonym','ensembl_external_gene_name') ORDER BY CASE alias_type WHEN 'gene_symbol' THEN 0 WHEN 'gene_synonym' THEN 1 ELSE 2 END,alias LIMIT 16",(gene,))]
    finally:c.close()
    if not r:return {'gene_id':gene,'symbol':gene,'coordinate_system':'B73 v5','annotation_available':False}
    x=dict(r);attrs=json.loads(x.get('attributes_json') or '{}')
    symbols=[a['alias'] for a in aliases if a['alias_type']=='gene_symbol']
    label=symbols[0] if symbols else x.get('label') or gene
    return {'gene_id':gene,'symbol':label,'coordinate_system':'B73 v5','chromosome':x.get('chromosome'),
      'start':x.get('start'),'end':x.get('end'),'strand':x.get('strand'),'annotation_available':True,
      'aliases':aliases,
      'description':str(attrs.get('description') or attrs.get('function') or attrs.get('annotation') or next((a['alias'] for a in aliases if a['alias_type']=='gene_synonym'),'') or '')[:600]}

def verify_release():
    m=json.loads(MANIFEST.read_text())
    for field,path in [('prediction_sha256',CONFIG['prediction_file']),('kg_sha256',CONFIG['kg_database'])]:
        if sha256(path)!=m[field]:raise RuntimeError('scientific_input_integrity_failed')
    for record in m['protected_files']:
        if sha256(record['path'])!=record['sha256']:raise RuntimeError('protected_model_input_changed')
    for record in m.get('derived_files',[]):
        if sha256(ROOT/record['path'])!=record['sha256']:raise RuntimeError('derived_profile_integrity_failed')
    model_rows()
    return {'model_gene_count':len(model_rows()),'model_sha256':m['prediction_sha256'],'kg_sha256':m['kg_sha256']}

@lru_cache(maxsize=64)
def public_context(gene):
    from public_repository import all_public_evidence,public_paths_all
    from collections import Counter
    raw=all_public_evidence(gene)
    return {'gene_id':gene,'evidence':raw['evidence'],'paths':public_paths_all(gene),
      'public_summary':{'evidence_record_count':raw['total'],'evidence_type_counts':dict(Counter(f['evidence_type'] for f in raw['evidence']))},
      'retrieval':{'total_facts':raw['total'],'record_cap':None,'invalid_source_records':raw['invalid_source_records'],
         'selection':'All local public records, with paginated browser delivery. Online literature is searched separately.'}}

def gene_response(identifier):
    resolved=resolve_gene(identifier)
    if resolved['status']!='ok':return resolved
    gene=resolved['gene_id'];ctx=public_context(gene)
    return {'status':'ok',**resolved,'gene':node_info(gene),'prediction':copy.deepcopy(model_rows().get(gene)),
      'model_profile':model_profile(gene),'evidence':ctx['evidence'][:20],'evidence_total':len(ctx['evidence']),
      'evidence_next_offset':20 if len(ctx['evidence'])>20 else None,
      'paths':ctx['paths'][:20],'paths_total':len(ctx['paths']),'public_summary':ctx['public_summary'],
      'retrieval':ctx['retrieval'],'kg_release':'MEKG public v3.1',
      'boundaries':['优先级不是表型概率；跨种子SD表示计算稳定性。','公共证据解释候选机制，不改变模型分数。',
                    '无记录不等于无功能；公共来源不自动属于独立验证。']}
