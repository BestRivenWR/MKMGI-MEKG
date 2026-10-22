"""Question-aware public research, auditable literature relations, separate overlay.

Only this release's cache is writable. Source KG, models and project data are
read-only. Public API failures remain visible; no claims of exhaustive search.
"""
import concurrent.futures, hashlib, json, re, sqlite3, time
from pathlib import Path
import httpx
from data_access import ROOT, CONFIG, node_info, resolve_gene, GENE_RE
import literature as lit

STORE = ROOT / 'state/research.sqlite'
TOPICS = {
    'regulation': r'调控|靶|结合|上游|下游|互作|regulat|target|binding|interact',
    'phenotype': r'表型|突变|遗传|粒长|粒重|粒宽|敲除|phenotyp|mutant|genetic',
    'starch': r'淀粉|starch|amylose|amylopectin',
    'protein': r'赖氨酸|储藏蛋白|醇溶|zein|lysine|storage protein',
    'model': r'模型|得分|分数|排名|预测|注意力|model|rank|score|attention',
    'recent': r'最新|近期|最近|进展|recent|latest',
}
ENGLISH = {'regulation':'regulation target promoter binding interaction',
           'phenotype':'mutant phenotype endosperm kernel', 'starch':'starch metabolism',
           'protein':'zein lysine storage protein', 'model':'endosperm development',
           'recent':'endosperm kernel'}

def db():
    c = sqlite3.connect(STORE, timeout=20); c.row_factory = sqlite3.Row
    return c

def init():
    c = db()
    c.execute('PRAGMA journal_mode=WAL')
    c.executescript('''CREATE TABLE IF NOT EXISTS query_runs(
      run_key TEXT PRIMARY KEY,gene_id TEXT,question TEXT,queries_json TEXT,
      results_json TEXT,trace_json TEXT,updated REAL);
      CREATE TABLE IF NOT EXISTS relations(
      relation_id TEXT PRIMARY KEY,source_id TEXT,target_id TEXT,record_json TEXT,
      origin TEXT,job_id TEXT,updated REAL);
      CREATE INDEX IF NOT EXISTS relation_source ON relations(source_id);
      CREATE INDEX IF NOT EXISTS relation_target ON relations(target_id);''')
    c.commit(); c.close()
    path = ROOT / 'curation/classic_relations.json'
    if path.is_file():
        for rel in json.loads(path.read_text(encoding='utf-8')):
            save_relation(rel, 'curated_primary_literature', None)

def intent(question):
    tags = [k for k,p in TOPICS.items() if re.search(p,question,re.I)]
    return tags or ['phenotype','regulation']

def query_score(p, question):
    """Rank by the actual question, not merely age, title, or assay keywords."""
    title = p.get('title',''); abstract = p.get('abstract','')
    text = (title + ' ' + abstract).lower()
    score = float(p.get('priority',0)) / 5
    if p.get('retracted'): return -10000
    if lit.classify_study(p).get('is_secondary'): score -= 30
    if p.get('relevance_scope') in ['cross_species_context','food_or_application_context']: score -= 80
    if p.get('identity_match') == 'indexed_search_only': score -= 35
    if p.get('species_check') == 'maize_in_text': score += 10
    if p.get('abstract'): score += 5
    for tag in intent(question):
        pattern = {'regulation':r'regulat|target|promoter|transactiv|interact|transcriptional|chip|emsa|binding',
         'phenotype':r'mutant|phenotyp|genetic|loss.of.function|complementation|knock.?out|kernel|endosperm',
         'starch':r'starch|amylopectin|amylose|branching',
         'protein':r'zein|lysine|storage.protein',
         'model':r'development|endosperm', 'recent':r'endosperm|kernel'}[tag]
        score += min(len(re.findall(pattern, text)),14)*2
        score += min(len(re.findall(pattern, title.lower())),4)*4
    words = re.findall(r'[A-Za-z][A-Za-z0-9_-]{2,}',question)
    for word in words:
        if word.lower() in text: score += 4
    if 'recent' in intent(question): score += max(0,p.get('year',0)-2015)*2
    # Genomics studies often report experiments in the abstract, not the title.
    if re.search(r'chip.{0,20}(?:seq|sequence)',text) and re.search(r'rna.{0,10}seq|transcriptom',text):score+=12
    return score

def all_papers(gene):
    c=lit.db()
    try:return [lit.classify_study(json.loads(r[0])) for r in c.execute('SELECT record_json FROM papers WHERE gene_id=?',(gene,))]
    finally:c.close()

def ranked_papers(gene,question,offset=0,page_size=12):
    items=sorted(all_papers(gene),key=lambda p:(-query_score(p,question),-p.get('year',0),p['paper_key']))
    return {'papers':[{**p,'question_relevance_score':round(query_score(p,question),2)} for p in items[offset:offset+page_size]],
            'total':len(items),'next_offset':offset+page_size if offset+page_size<len(items) else None}

def default_queries(gene,question):
    terms=lit.search_terms(gene)
    identities=[a for a in terms if re.search(r'[A-Za-z]',a) and not a.startswith(('Zm000','GRMZM','AC0'))]
    # Long descriptive aliases can be poor literal search identities.
    identities=sorted(identities,key=lambda a:(len(a.split())>3,not bool(re.search(r'opaque|amylose|ereb',a,re.I)),len(a)))
    alias=identities[0] if identities else gene
    topics=' '.join(ENGLISH[t] for t in intent(question) if t!='recent')
    return [f'maize "{alias}" {topics}', f'maize "{alias}" endosperm']

def _search_one(source,gene,query):
    until=time.monotonic(); params={}; results=[]
    with httpx.Client(timeout=httpx.Timeout(18,connect=7),headers={'User-Agent':lit.UA},trust_env=False) as client:
        if source=='EuropePMC':
            j=client.get(lit.API[source],params={'query':query,'format':'json','resultType':'core','pageSize':40,'sort':'RELEVANCE'}).raise_for_status().json()
            results=[lit.unpack_epmc(x) for x in j.get('resultList',{}).get('result',[])];hits=j.get('hitCount',0)
        elif source=='Crossref':
            j=client.get(lit.API[source],params={'query.bibliographic':query,'rows':40,'filter':'type:journal-article'}).raise_for_status().json()['message']
            results=[lit.unpack_crossref(x) for x in j.get('items',[])];hits=j.get('total-results',0)
        else:
            j=client.get(lit.API[source],params={'search':query,'per-page':40}).raise_for_status().json()
            results=[lit.unpack_openalex(x) for x in j.get('results',[])];hits=j.get('meta',{}).get('count',0)
    keys=[];terms=lit.search_terms(gene)
    for r in results:
        clean=lit.annotate(r,gene,terms)
        if clean:
            lit.save_paper(clean);keys.append(clean['paper_key'])
    return {'source':source,'query':query,'status':'ok','indexed_hits':hits,'accepted':len(keys),
            'elapsed_seconds':round(time.monotonic()-until,2)},keys

def search(gene,question,queries=None,refresh=False,initial=False):
    queries=queries or default_queries(gene,question)
    queries=list(dict.fromkeys(q.strip() for q in queries if isinstance(q,str) and 3<len(q.strip())<=350))[:4]
    key=hashlib.sha256(json.dumps([gene,question,queries,'question-relevance-v1'],ensure_ascii=False).encode()).hexdigest()
    c=db()
    try:
        old=c.execute('SELECT * FROM query_runs WHERE run_key=?',(key,)).fetchone()
        if initial and not refresh and not old:
            old=c.execute('SELECT * FROM query_runs WHERE gene_id=? AND question=? AND updated>? ORDER BY updated DESC LIMIT 1',
                (gene,question,time.time()-86400)).fetchone()
    finally:c.close()
    if old and not refresh and time.time()-old['updated']<86400:
        return {'queries':json.loads(old['queries_json']),'trace':json.loads(old['trace_json']),'cached':True,**ranked_papers(gene,question)}
    trace=[];keys=[]
    # Question searches are bounded work slices, not an evidence-count cap.
    # Persistent general discovery retains cursors for all remaining results.
    jobs=[(src,q) for q in queries for src in CONFIG['online_sources']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        pending={pool.submit(_search_one,s,gene,q):(s,q) for s,q in jobs}
        for f in concurrent.futures.as_completed(pending):
            source,query=pending[f]
            try:t,k=f.result();trace.append(t);keys+=k
            except Exception as exc:
                code=str(exc.response.status_code) if isinstance(exc,httpx.HTTPStatusError) else type(exc).__name__
                trace.append({'source':source,'query':query,'status':'failed','error':code})
    c=db()
    try:
        c.execute('INSERT OR REPLACE INTO query_runs VALUES (?,?,?,?,?,?,?)',
           (key,gene,question,json.dumps(queries,ensure_ascii=False),json.dumps(sorted(set(keys))),json.dumps(trace,ensure_ascii=False),time.time()));c.commit()
    finally:c.close()
    return {'queries':queries,'trace':trace,'cached':False,**ranked_papers(gene,question)}

def lookup_doi(gene,doi):
    doi=lit.normal_doi(doi)
    if not re.fullmatch(r'10\.\d{4,9}/[^\s<>]{2,200}',doi):raise ValueError('doi_invalid')
    records=[]
    with httpx.Client(timeout=18,trust_env=False,headers={'User-Agent':lit.UA}) as client:
        try:
            j=client.get(lit.API['EuropePMC'],params={'query':'DOI:"'+doi+'"','format':'json','resultType':'core'}).raise_for_status().json()
            records += [lit.unpack_epmc(x) for x in j.get('resultList',{}).get('result',[])]
        except httpx.HTTPError:pass
        if not records:
            from urllib.parse import quote
            j=client.get(lit.API['Crossref']+'/'+quote(doi,safe='')).raise_for_status().json()
            records.append(lit.unpack_crossref(j['message']))
    for r in records:
        clean=lit.annotate(r,gene,lit.search_terms(gene))
        if clean:lit.save_paper(clean)
    return read(gene,'doi:'+doi,'')

def read(gene,key,question):
    p=lit.read_paper(gene,key)
    paragraphs=p.get('passages',[])
    terms=lit.search_terms(gene)
    def value(s):
        return sum(len(re.findall(pattern,s,re.I)) for pattern in [
           r'we (?:show|found|demonstrat|identif)|our results|these results|in conclusion',
           r'chip|emsa|reporter|transactiv|binding|mutant|complement|luciferase',
           r'zein|starch|kernel|endosperm'])+len(lit.matches(s,terms))*2
    # Keep study-result passages before background paragraphs; complete available
    # matching text remains in the local cache, not silently called full audit.
    p['selected_passages']=sorted(paragraphs,key=value,reverse=True)[:6]
    p['available_matching_passages']=len(paragraphs)
    p['question']=question
    return p

def source_text(fact):
    d=fact.get('details',{})
    return '\n'.join([d.get('abstract') or '',*(d.get('selected_passages') or []),*(d.get('passages') or [])])

def clean_quote(s):return re.sub(r'\s+',' ',str(s or '')).strip()

def resolve_entity(symbol):
    if not isinstance(symbol,str) or not 1<=len(symbol)<=80 or re.search(r'[<>\r\n]',symbol):raise ValueError('entity_invalid')
    found=resolve_gene(symbol) if re.fullmatch(r'[A-Za-z0-9_.-]{2,64}',symbol) else {'status':'not_found'}
    if found['status']=='ok':
        return {'id':found['gene_id'],'label':node_info(found['gene_id'])['symbol'],'mapping_status':'recorded_B73_v5_mapping','requested_symbol':symbol}
    # Unmapped target symbols are legitimate paper entities, not invented v5 IDs.
    if GENE_RE.fullmatch(symbol):raise ValueError('unverified_v5_identifier')
    return {'id':'symbol:'+symbol,'label':symbol,'mapping_status':'ambiguous' if found['status']=='ambiguous' else 'v5_mapping_pending',
            'candidate_ids':found.get('candidate_ids',[])}

def relation_facts(gene):
    c=db()
    try:
        rows=[json.loads(r[0]) for r in c.execute('SELECT record_json FROM relations WHERE source_id=? OR target_id=? ORDER BY CASE origin WHEN \'curated_primary_literature\' THEN 0 ELSE 1 END,relation_id',(gene,gene))]
        seen=set();out=[]
        for r in rows:
            key=(r['source'],r['target'],r['relation'],r['doi'])
            if key in seen:continue
            seen.add(key);out.append(r)
        return out
    finally:c.close()

def save_relation(record,origin,job_id):
    rel=dict(record)
    a=resolve_entity(rel['source_symbol']);b=resolve_entity(rel['target_symbol'])
    rel['source']=a['id'];rel['target']=b['id'];rel['source_entity']=a;rel['target_entity']=b
    rel['origin']=origin;rel['job_id']=job_id
    fingerprint=[a['id'],b['id'],rel['relation'],rel['doi'],origin]
    rel['id']='research:'+hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()[:24]
    c=db()
    try:
        c.execute('INSERT OR REPLACE INTO relations VALUES (?,?,?,?,?,?,?)',
           (rel['id'],a['id'],b['id'],json.dumps(rel,ensure_ascii=False),origin,job_id,time.time()));c.commit()
    finally:c.close()
    return rel

def network_relations(gene):
    out=[]
    for r in relation_facts(gene):
        curated=r['origin']=='curated_primary_literature'
        experimental=curated and r['support_level'] in ['binding_and_functional','promoter_functional','binding_only','protein_interaction']
        out.append({'id':r['id'],'source':r['source'],'target':r['target'],'source_entity':r['source_entity'],'target_entity':r['target_entity'],
          'predicate':r['relation'],'directed':r['relation'] not in ['PROTEIN_INTERACTION','COOPERATIVE_REGULATION'],
          'layer':'classic_experimental' if experimental else 'literature_reported','origin':r['origin'],
          'public_experimental':experimental,'used_by_current_model':False,'weight':None,'confidence':None,
          'support':r['finding'],'support_level':r['support_level'],'assays':r.get('assays',[]),
          'source_urls':['https://doi.org/'+r['doi']],'evidence_ids':r.get('citations',[]),'evidence_types':['PRIMARY_LITERATURE_RELATION'],
          'stage':r.get('stage'),'tissue':r.get('tissue'),'genotype':r.get('genotype'),
          'attributes':{'study_title':r.get('title'),'support_quote':r.get('quote'),'audit':r.get('audit'),
           'limitations':r.get('limitations'),'not_independent_validation':True,'mapping':[r['source_entity'],r['target_entity']]}})
    return out

def evidence_fact(rel,gene):
    return {'evidence_id':rel['id'],'subject_node_id':'gene:'+rel['source'] if GENE_RE.fullmatch(rel['source']) else rel['source'],
      'object_node_id':'gene:'+rel['target'] if GENE_RE.fullmatch(rel['target']) else rel['target'],
      'subject_label':rel['source_symbol'],'object_label':rel['target_symbol'],'evidence_type':'PRIMARY_LITERATURE_RELATION',
      'origin':rel['origin'],'predicate':rel['relation'],'source_url':'https://doi.org/'+rel['doi'],
      'focus_role':'subject' if rel['source']==gene else 'object','category':'published_literature',
      'details':{k:rel.get(k) for k in ['finding','support_level','assays','doi','title','quote','stage','genotype','tissue','limitations','audit']},
      'caveat':'原始研究报告的关系；具体实验支持和材料条件见原文，不等同于本项目独立验证。'}
