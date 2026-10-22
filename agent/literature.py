"""Resumable scholarly discovery. Fixed trusted APIs, provenance and no total-hit cap.

Search hits are NOT automatically verified gene function. Metadata/abstract/full
text reading levels, gene identity, retractions and source failures stay explicit.
Only public identifiers/aliases/DOIs are sent to scholarly APIs, never predictions.
"""
import collections,hashlib,html,json,re,sqlite3,threading,time,xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import quote,urlencode,urlsplit,urljoin
import httpx
from functools import lru_cache
from data_access import ROOT,CONFIG,node_info
from public_repository import aliases,curated_links
STORE=ROOT/'state/literature.sqlite';LOCKS=collections.defaultdict(threading.Lock)
UA='MEKG-ResearchEvidence/2.0 (public maize gene literature; scholarly research)'
API={'EuropePMC':'https://www.ebi.ac.uk/europepmc/webservices/rest/search',
     'Crossref':'https://api.crossref.org/works','OpenAlex':'https://api.openalex.org/works'}

def db():
    STORE.parent.mkdir(exist_ok=True,mode=0o700)
    c=sqlite3.connect(STORE,timeout=20);c.row_factory=sqlite3.Row;return c

def init():
    c=db();c.execute('PRAGMA journal_mode=WAL')
    c.executescript('''CREATE TABLE IF NOT EXISTS papers(gene_id TEXT,paper_key TEXT,record_json TEXT,priority REAL,year INTEGER,updated REAL,PRIMARY KEY(gene_id,paper_key));
    CREATE INDEX IF NOT EXISTS paper_sort ON papers(gene_id,priority DESC,year);
    CREATE TABLE IF NOT EXISTS searches(gene_id TEXT PRIMARY KEY,status TEXT,updated REAL,plan_json TEXT,trace_json TEXT,aliases_json TEXT);
    CREATE TABLE IF NOT EXISTS read_cache(paper_key TEXT PRIMARY KEY,result_json TEXT,updated REAL);
    CREATE TABLE IF NOT EXISTS settings(key TEXT PRIMARY KEY,value TEXT);''')
    version=c.execute("SELECT value FROM settings WHERE key='relevance_policy'").fetchone()
    c.execute("UPDATE searches SET status='partial' WHERE status='running'");c.commit();c.close()
    if not version or version[0]!='maize-primary-first-v2':
        rerank_cached();c=db()
        c.execute('INSERT OR REPLACE INTO settings VALUES (?,?)',('relevance_policy','maize-primary-first-v2'));c.commit();c.close()

def text(value):return re.sub(r'\s+',' ',html.unescape(re.sub(r'<[^>]+>',' ',str(value or '')))).strip()
def normal_doi(x):return re.sub(r'^https?://(?:dx\.)?doi\.org/','',str(x or ''),flags=re.I).strip().lower().rstrip('.,;')
def paper_key(r):
    if r.get('doi'):return 'doi:'+normal_doi(r['doi'])
    if r.get('pmid'):return 'pmid:'+str(r['pmid'])
    return 'title:'+hashlib.sha256(re.sub(r'\W','',r.get('title','').lower()).encode()).hexdigest()[:24]

@lru_cache(maxsize=1024)
def search_terms(gene):
    values=[]
    for a in aliases(gene):
        a=re.sub(r'\.\d+$','',a).strip()
        if a in ['-','unknown'] or a.isdigit() or len(a)>100:continue
        # Exclude generic functional descriptions as query identities.
        if len(a.split())>7:continue
        if not re.search(r'\d',a) and re.search(r'\b(protein|domain|family|kinase|hypothetical|uncharacterized|putative|repeat)\b',a,re.I):continue
        values.append(a)
        m=re.fullmatch(r'opaque(?: endosperm)?\s*-?\s*(\d+)',a,re.I)
        if m:values += ['opaque'+m[1],'opaque-'+m[1],'opaque '+m[1]]
        if re.fullmatch(r'[A-Za-z]{2,}\d+',a) and not a.startswith(('Zm','GRM','AC','AF')):values.append('Zm'+a)
    values.append(gene)
    return sorted(set(values),key=lambda s:(s!=gene,s.lower()))

def make_plan(gene,terms):
    safe=[x.replace('"',' ').replace('\\',' ') for x in terms]
    groups=[safe[i:i+18] for i in range(0,len(safe),18)]
    plan=[]
    for row in curated_links(gene):
        if row['doi']:
            plan += [{'source':'Crossref','mode':'doi','query':row['doi'],'cursor':'*','done':False,'seed':row},
                     {'source':'EuropePMC','mode':'doi','query':'DOI:"'+row['doi']+'"','cursor':'*','done':False,'seed':row}]
    for group in groups:
        q='('+' OR '.join('"'+a+'"' for a in group)+') AND (maize OR "Zea mays")'
        plan += [{'source':'EuropePMC','mode':'search','query':q,'cursor':'*','done':False},
                 {'source':'OpenAlex','mode':'search','query':q,'cursor':'*','done':False}]
    # Crossref bibliographic search is broad: exact post-filtering is mandatory.
    # All cursor pages remain resumable; raw database counts are not gene counts.
    symbols=[x for x in safe if len(x)>3 and not re.match(r'^(Zm000|GRMZM|AC\d|LOC\d)',x)]
    if symbols:
        plan.append({'source':'Crossref','mode':'search','query':'maize '+symbols[0],'cursor':'*','done':False})
    return plan

def unpack_epmc(x):
    types=x.get('pubTypeList',{}).get('pubType',[])
    journal=x.get('journalInfo',{}).get('journal',{}).get('title','')
    return {'title':text(x.get('title')),'doi':normal_doi(x.get('doi')),'pmid':x.get('pmid') or (x.get('id') if x.get('source')=='MED' else None),
      'pmcid':x.get('pmcid'),'year':int(x.get('pubYear') or 0),'authors':x.get('authorString',''),'journal':journal,
      'abstract':text(x.get('abstractText')),'publication_types':types,'citation_count':int(x.get('citedByCount') or 0),
      'retracted':x.get('isRetracted')=='Y' or any('retract' in t.lower() for t in types),
      'fulltext_available':bool(x.get('pmcid')),'source_ids':{'EuropePMC':str(x.get('id',''))},'sources':['EuropePMC'],
      'related_primary_dois':list(dict.fromkeys(normal_doi(m) for item in x.get('commentCorrectionList',{}).get('commentCorrection',[])
         if item.get('type')=='Comment on' for m in re.findall(r'10\.\d{4,9}/[^\s<>]+',item.get('reference',''))))}

def classify_study(record):
    """Apply source publication types also to old cached records, without edits."""
    r=dict(record);types=' '.join(r.get('publication_types') or []).lower()
    r['is_review']=bool(r.get('is_review') or 'review' in types)
    r['is_commentary']=bool(re.search(r'\bcomment\b|editorial|news|perspective',types))
    r['is_secondary']=r['is_review'] or r['is_commentary']
    return r

def unpack_crossref(x):
    pub=x.get('published',x.get('issued',{})).get('date-parts',[[0]])[0]
    updates=x.get('update-to',[]);assertions=x.get('assertion',[])
    flags=json.dumps([updates,assertions,x.get('relation',{})],ensure_ascii=False).lower()
    return {'title':text((x.get('title') or [''])[0]),'doi':normal_doi(x.get('DOI')),'year':int(pub[0] or 0),
      'authors':', '.join((a.get('given','')+' '+a.get('family','')).strip() for a in x.get('author',[])),
      'journal':(x.get('container-title') or [''])[0],'abstract':text(x.get('abstract')),
      'publication_types':[x.get('type','journal-article')],'citation_count':int(x.get('is-referenced-by-count') or 0),
      'retracted':False,'retraction_check':'metadata_not_complete','update_metadata':updates,
      'update_notice':bool(updates) or 'retract' in flags,'fulltext_available':False,
      'source_ids':{'Crossref':normal_doi(x.get('DOI'))},'sources':['Crossref']}

def unpack_openalex(x):
    words={p:w for w,positions in (x.get('abstract_inverted_index') or {}).items() for p in positions}
    ids=x.get('ids',{});pmid=str(ids.get('pmid','')).rstrip('/').rsplit('/',1)[-1] or None
    location=x.get('primary_location') or {};journal=(location.get('source') or {}).get('display_name','')
    return {'title':text(x.get('display_name')),'doi':normal_doi(x.get('doi')),'pmid':pmid,
      'year':int(x.get('publication_year') or 0),'authors':', '.join((a.get('author') or {}).get('display_name','') for a in x.get('authorships',[])),
      'journal':journal,'abstract':' '.join(words[p] for p in sorted(words)),
      'publication_types':[x.get('type','article')],'citation_count':int(x.get('cited_by_count') or 0),
      'retracted':bool(x.get('is_retracted')),'fulltext_available':bool(x.get('open_access',{}).get('is_oa')),
      'source_ids':{'OpenAlex':x.get('id')},'sources':['OpenAlex']}

def matches(s,terms):
    return [a for a in terms if re.search(r'(?<![A-Za-z0-9])'+re.escape(a)+r'(?![A-Za-z0-9])',s,re.I)]

def annotate(r,gene,terms,seed=None):
    if not r.get('title'):return None
    r['title']=text(r['title']);r['abstract']=text(r.get('abstract'))
    tm=matches(r['title'],terms);am=matches(r['abstract'],terms)
    r['matched_aliases']=sorted(set(tm+am));r['curated_gene_link']=bool(seed)
    r['identity_match']='curated_publication_link' if seed else 'title_and_abstract' if tm and am else 'title' if tm else 'abstract' if am else 'indexed_search_only'
    blob=r['title']+' '+r['abstract'];species=bool(re.search(r'\bmaize\b|\bZea mays\b|\bcorn\b',blob,re.I))
    r['species_check']='maize_in_text' if species else 'curated_maize_link_needs_original_check' if seed else 'species_needs_original_check'
    other_in_title=bool(re.search(r'\b(sorghum|rice|arabidopsis|wheat|barley|mouse|mice|human|rat|chicken)\b',r['title'],re.I))
    maize_in_title=bool(re.search(r'\b(maize|corn|Zea mays)\b',r['title'],re.I))
    food_application=bool(re.search(r'\b(mice|mouse|rat|diabetes|gut barrier|microbiota|digestibility|glycemic)\b',r['title'],re.I))
    r['relevance_scope']='food_or_application_context' if food_application and not seed else 'cross_species_context' if other_in_title and not maize_in_title and not seed else 'maize_gene_research' if (species and (tm or am)) or seed else 'identity_needs_reading'
    ptypes=' '.join(r.get('publication_types',[])).lower()
    r['is_review']='review' in ptypes or bool(re.search(r'\breview\b',r['title'],re.I))
    r['is_preprint']=bool(re.search(r'preprint|posted-content',ptypes))
    tags=[]
    for name,pattern in [('遗传扰动',r'\bmutant|knock.?out|CRISPR|RNAi|transgenic|complementation|loss.of.function'),
                         ('结合实验',r'ChIP|CUT.?Tag|EMSA|DNA.binding|DAP.seq|yeast.one.hybrid'),
                         ('报告实验',r'luciferase|reporter assay|dual.luciferase'),
                         ('生化实验',r'enzyme activity|enzymatic|biochemical|immunoprecipitation|pull.down'),
                         ('表达/组学',r'RNA.seq|transcriptom|ATAC|spatial|single.cell|single.nucleus')]:
        if re.search(pattern,blob,re.I):tags.append(name)
    r['experiment_clues']=tags
    related=bool(tm or am or seed)
    r['relevance_status']='gene_mentioned' if related else 'indexed_match_needs_reading'
    r['reading_level']='abstract' if r['abstract'] else 'metadata'
    r['gene_id']=gene;r['paper_key']=paper_key(r)
    r['source_url']='https://doi.org/'+r['doi'] if r.get('doi') else 'https://pubmed.ncbi.nlm.nih.gov/'+str(r['pmid'])+'/' if r.get('pmid') else next(iter(r.get('source_ids',{}).values()),'')
    if seed:
        r['curation_note']=text(seed.get('phenotype'))
        r['curation_limitations']=text(seed.get('limitations'))
        r['identity_basis']=text(seed.get('identity_basis'))
        if re.search(r'撤稿|retract',r['curation_limitations']+' '+r['curation_note'],re.I):r['retracted']=True;r['retraction_check']='existing_curation_warning'
    r['priority']=(-1000 if r.get('retracted') else 0)+(35 if seed else 0)+(25 if tm else 10 if am else 0)+(12 if tags and not r['is_review'] else 0)+(5 if species else 0)+(4 if r.get('year',9999)<2010 else 0)+min(r.get('citation_count',0),500)/100
    if r['relevance_scope'] in ['cross_species_context','food_or_application_context']:r['priority']-=55
    if r['identity_match']=='indexed_search_only':r['priority']-=40
    if maize_in_title and (tm or am):r['priority']+=15
    if r['relevance_scope']=='maize_gene_research' and re.search(r'\b(kernel|endosperm|seed|starch|zein|grain)\b',r['title'],re.I):r['priority']+=10
    if r['is_review']:r['priority']-=15
    if r['is_preprint']:r['priority']-=20
    # Crossref's bibliographic pool is broad; never display an unrelated hit as gene evidence.
    if r['sources']==['Crossref'] and not related:return None
    r['interpretation']='检索命中不是遗传或直接调控认证；需核对基因主体、实验条件和原文。'
    return r

def merge(old,new):
    for k in ['doi','pmid','pmcid','title','journal','authors','curation_note','curation_limitations','identity_basis']:
        if not old.get(k) and new.get(k):old[k]=new[k]
    if len(new.get('abstract',''))>len(old.get('abstract','')):old['abstract']=new['abstract'];old['reading_level']=new['reading_level']
    old['sources']=sorted(set(old.get('sources',[])+new.get('sources',[])));old.setdefault('source_ids',{}).update(new.get('source_ids',{}))
    old['matched_aliases']=sorted(set(old.get('matched_aliases',[])+new.get('matched_aliases',[])))
    old['experiment_clues']=sorted(set(old.get('experiment_clues',[])+new.get('experiment_clues',[])))
    old['publication_types']=sorted(set(old.get('publication_types',[])+new.get('publication_types',[])))
    old['related_primary_dois']=sorted(set(old.get('related_primary_dois',[])+new.get('related_primary_dois',[])))
    for k in ['retracted','fulltext_available','curated_gene_link','update_notice']:old[k]=bool(old.get(k) or new.get(k))
    old['citation_count']=max(old.get('citation_count',0),new.get('citation_count',0));old['priority']=max(old.get('priority',0),new.get('priority',0))
    if old['retracted']:old['priority']=-1000
    if new.get('identity_match')!='indexed_search_only':old['identity_match']=new['identity_match'];old['relevance_status']=new['relevance_status']
    return old

def save_paper(r):
    c=db()
    try:
        old=c.execute('SELECT record_json FROM papers WHERE gene_id=? AND paper_key=?',(r['gene_id'],r['paper_key'])).fetchone()
        if old:r=merge(json.loads(old[0]),r)
        else:
            # Fallback duplicate matching when one source omitted DOI.
            nt=re.sub(r'\W','',r['title'].lower())
            for candidate in c.execute('SELECT paper_key,record_json FROM papers WHERE gene_id=?',(r['gene_id'],)):
                x=json.loads(candidate[1]);xt=re.sub(r'\W','',x['title'].lower())
                if nt==xt and nt:
                    r=merge(x,r);r['paper_key']=candidate[0];break
        seed={'phenotype':r.get('curation_note'),'limitations':r.get('curation_limitations'),'identity_basis':r.get('identity_basis')} if r.get('curated_gene_link') else None
        r=annotate(r,r['gene_id'],search_terms(r['gene_id']),seed)
        if r:c.execute('INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?)',(r['gene_id'],r['paper_key'],json.dumps(r,ensure_ascii=False),r['priority'],r.get('year',0),time.time()));c.commit()
    finally:c.close()

def rerank_cached():
    """Update display-only relevance labels in our own new literature cache."""
    c=db()
    try:rows=[json.loads(r[0]) for r in c.execute('SELECT record_json FROM papers')]
    finally:c.close()
    for r in rows:save_paper(r)

def fetch_page(client,p):
    source=p['source'];size=CONFIG['literature_page_size']
    if source=='Crossref' and p['mode']=='doi':
        j=client.get(API[source]+'/'+quote(p['query'],safe='')).raise_for_status().json();return [unpack_crossref(j['message'])],None,1
    if source=='EuropePMC':
        params={'query':p['query'],'format':'json','resultType':'core','pageSize':size,'cursorMark':p['cursor']}
        j=client.get(API[source],params=params).raise_for_status().json();items=j.get('resultList',{}).get('result',[])
        cursor=j.get('nextCursorMark')
        if not items or cursor==p['cursor'] or len(items)<size:cursor=None
        return [unpack_epmc(x) for x in items],cursor,int(j.get('hitCount',0))
    if source=='Crossref':
        j=client.get(API[source],params={'query.bibliographic':p['query'],'rows':size,'cursor':p['cursor'],'filter':'type:journal-article'}).raise_for_status().json()['message']
        items=j.get('items',[]);cursor=j.get('next-cursor') if len(items)==size else None
        if cursor==p['cursor']:cursor=None
        return [unpack_crossref(x) for x in items],cursor,int(j.get('total-results',0))
    j=client.get(API[source],params={'search':p['query'],'per-page':size,'cursor':p['cursor']}).raise_for_status().json()
    items=j.get('results',[]);cursor=j.get('meta',{}).get('next_cursor') if len(items)==size else None
    if cursor==p['cursor']:cursor=None
    return [unpack_openalex(x) for x in items],cursor,int(j.get('meta',{}).get('count',0))

def update_search(gene,status,plan,trace,terms):
    c=db()
    try:c.execute('INSERT OR REPLACE INTO searches VALUES (?,?,?,?,?,?)',(gene,status,time.time(),json.dumps(plan,ensure_ascii=False),json.dumps(trace,ensure_ascii=False),json.dumps(terms,ensure_ascii=False)));c.commit()
    finally:c.close()

def status(gene):
    c=db()
    try:
        r=c.execute('SELECT * FROM searches WHERE gene_id=?',(gene,)).fetchone();n=c.execute('SELECT count(*) FROM papers WHERE gene_id=?',(gene,)).fetchone()[0]
        if not r:return {'gene_id':gene,'status':'not_started','paper_count':n,'queries':[],'trace':[]}
        plan=json.loads(r['plan_json']);trace=json.loads(r['trace_json'])
        return {'gene_id':gene,'status':r['status'],'updated':r['updated'],'paper_count':n,'aliases':json.loads(r['aliases_json']),
           'queries':[{k:v for k,v in p.items() if k not in ['seed','cursor']} for p in plan],'trace':trace,
           'can_continue':any(not p['done'] for p in plan),'sources':CONFIG['online_sources'],
           'coverage_note':'联合数据库与已有策展链接，不保证穷尽全部已发表文献；摘要/全文可读性与索引覆盖不同。'}
    finally:c.close()

def run_search(gene,refresh=False):
    """One bounded work slice. Continue resumes cursors, not an arbitrary hit cap."""
    if not LOCKS[gene].acquire(blocking=False):return status(gene)
    try:
        c=db();old=c.execute('SELECT * FROM searches WHERE gene_id=?',(gene,)).fetchone();c.close()
        terms=search_terms(gene)
        if old and not refresh:
            if old['status']=='complete' and time.time()-old['updated']<CONFIG['literature_cache_days']*86400:return status(gene)
            plan=json.loads(old['plan_json']);trace=json.loads(old['trace_json'])
        else:plan=make_plan(gene,terms);trace=[]
        if old and old['status']=='complete' and not refresh:plan=make_plan(gene,terms)
        update_search(gene,'running',plan,trace,terms)
        until=time.monotonic()+CONFIG['literature_search_seconds'];failed=set()
        with httpx.Client(headers={'User-Agent':UA},timeout=httpx.Timeout(15,connect=8),follow_redirects=False,trust_env=False) as client:
            while time.monotonic()<until:
                pending=[p for i,p in enumerate(plan) if not p['done'] and i not in failed]
                if not pending:break
                # Breadth-first across sources before consuming later cursor pages.
                p=min(pending,key=lambda x:x.get('pages',0));idx=plan.index(p)
                try:
                    records,cursor,hits=fetch_page(client,p);kept=0
                    for r in records:
                        clean=annotate(r,gene,terms,p.get('seed'))
                        if clean:save_paper(clean);kept+=1
                    p['pages']=p.get('pages',0)+1;p['indexed_hits']=hits;p['records_seen']=p.get('records_seen',0)+len(records);p['accepted']=p.get('accepted',0)+kept
                    p['cursor']=cursor;p['done']=cursor is None;p.pop('error',None)
                    trace.append({'source':p['source'],'query':p['query'],'status':'ok','records_seen':len(records),'accepted':kept,'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
                except (httpx.HTTPError,ValueError,KeyError,TypeError,AttributeError,IndexError) as exc:
                    code=exc.response.status_code if isinstance(exc,httpx.HTTPStatusError) else type(exc).__name__
                    p['error']=str(code);failed.add(idx)
                    trace.append({'source':p['source'],'query':p['query'],'status':'failed','error':str(code),'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
                update_search(gene,'running',plan,trace,terms)
                time.sleep(.35)  # Respect public APIs; never a user access quota.
        final='complete' if all(p['done'] for p in plan) else 'partial'
        update_search(gene,final,plan,trace,terms);return status(gene)
    finally:LOCKS[gene].release()

def papers_page(gene,offset=0,page_size=20,filter_by='all'):
    c=db()
    try:items=[json.loads(r[0]) for r in c.execute('SELECT record_json FROM papers WHERE gene_id=? ORDER BY priority DESC,year,paper_key',(gene,))]
    finally:c.close()
    items=[classify_study(x) for x in items]
    if filter_by=='experimental':items=[x for x in items if x.get('experiment_clues') and not x.get('is_secondary') and not x.get('retracted')]
    elif filter_by=='classic':items=[x for x in items if 0<x.get('year',0)<=2010 and not x.get('retracted')]
    elif filter_by=='review':items=[x for x in items if x.get('is_secondary')]
    elif filter_by=='warning':items=[x for x in items if x.get('retracted') or x.get('update_notice')]
    return {'gene_id':gene,'papers':items[offset:offset+page_size],'total':len(items),'offset':offset,
       'next_offset':offset+page_size if offset+page_size<len(items) else None,'search':status(gene)}

def get_paper(gene,key):
    c=db()
    try:r=c.execute('SELECT record_json FROM papers WHERE gene_id=? AND paper_key=?',(gene,key)).fetchone()
    finally:c.close()
    if not r:raise ValueError('paper_not_in_gene_search')
    return classify_study(json.loads(r[0]))

def fingerprint(gene):
    c=db()
    try:rows=[r[0] for r in c.execute('SELECT record_json FROM papers WHERE gene_id=? ORDER BY paper_key',(gene,))]
    finally:c.close()
    return hashlib.sha256('\n'.join(rows).encode()).hexdigest()

PUBLISHER_HOSTS={'doi.org','academic.oup.com','oup.silverchair-cdn.com','pmc.ncbi.nlm.nih.gov',
  'europepmc.org','www.pnas.org','www.nature.com','onlinelibrary.wiley.com','www.cell.com',
  'link.springer.com','www.sciencedirect.com','journals.plos.org','www.frontiersin.org',
  'www.tandfonline.com','academic.oup.com.cn'}

def publisher_passages(record,terms):
    """Ordinary public HTTPS only; no CAPTCHA/paywall/proxy circumvention.

    DOI comes from a verified indexed record, not a model-supplied arbitrary URL.
    Every redirect is allowlisted before fetching; internal URLs are impossible.
    """
    from bs4 import BeautifulSoup
    doi=record.get('doi')
    if not doi:return [],'no_indexed_doi'
    url='https://doi.org/'+quote(doi,safe='/():._-')
    with httpx.Client(timeout=httpx.Timeout(16,connect=6),headers={'User-Agent':UA},trust_env=False,follow_redirects=False) as client:
        for _ in range(6):
            u=urlsplit(url)
            if u.scheme!='https' or u.hostname not in PUBLISHER_HOSTS or u.port not in [None,443] or u.username or u.password:
                return [],'publisher_redirect_not_allowed'
            response=client.get(url)
            if response.status_code in [301,302,303,307,308]:
                location=response.headers.get('location')
                if not location:return [],'publisher_redirect_missing'
                url=urljoin(url,location);continue
            response.raise_for_status()
            if len(response.content)>8_000_000:return [],'publisher_document_too_large'
            if 'html' not in response.headers.get('content-type',''):return [],'publisher_not_html'
            soup=BeautifulSoup(response.content,'html.parser')
            canonical=soup.find('meta',attrs={'name':'citation_doi'})
            html_doi=normal_doi(canonical.get('content')) if canonical else None
            if html_doi and html_doi!=doi:return [],'publisher_doi_mismatch'
            if not html_doi and doi.lower() not in response.text.lower():return [],'publisher_identity_not_confirmed'
            headings=' '.join(h.get_text(' ',strip=True) for h in soup.find_all(['h2','h3']))
            if not re.search(r'\b(results|discussion|materials and methods|experimental procedures)\b',headings,re.I):
                return [],'publisher_abstract_or_access_page_only'
            body=soup.find('article') or soup.find('main') or soup
            paras=[]
            for p in body.find_all('p'):
                if any(re.search(r'referenc|ref-list|citation|supplementary',str(a.get('id',''))+' '+' '.join(a.get('class',[])),re.I) for a in p.parents if getattr(a,'attrs',None)):continue
                s=text(p.get_text(' ',strip=True))
                if len(s)>70 and matches(s,terms):paras.append(s)
            return list(dict.fromkeys(paras)),None
    return [],'publisher_redirect_limit'

def read_paper(gene,key):
    r=get_paper(gene,key);c=db()
    cache_key=gene+'|'+key+'|primary-reader-v3'
    try:cached=c.execute('SELECT result_json FROM read_cache WHERE paper_key=? AND updated>?',(cache_key,time.time()-7*86400)).fetchone()
    finally:c.close()
    if cached:return classify_study({**r,**json.loads(cached[0])})
    terms=search_terms(gene);paras=[];level=r.get('reading_level','metadata');error=None
    if r.get('pmcid') and re.fullmatch(r'PMC\d+',r['pmcid']):
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+r['pmcid']+'/fullTextXML'
        try:
            with httpx.Client(headers={'User-Agent':UA},timeout=20,trust_env=False) as client:
                response=client.get(url).raise_for_status()
            # No arbitrary URL fetching / publisher scraping / access-control bypass.
            if len(response.content)>12_000_000:raise ValueError('document_too_large')
            root=ET.fromstring(response.content)
            body=root.find('body')
            if body is not None:
                for p in body.iter('p'):
                    s=text(' '.join(p.itertext()))
                    if matches(s,terms):paras.append(s)
                level='open_fulltext_gene_passages'
        except (httpx.HTTPError,ET.ParseError,ValueError) as exc:error=type(exc).__name__
    if not paras:
        try:
            paras,publisher_error=publisher_passages(r,terms)
            if paras:level='open_fulltext_gene_passages';error=None
            else:error=publisher_error or error
        except (httpx.HTTPError,ValueError) as exc:
            error='publisher_'+str(exc.response.status_code) if isinstance(exc,httpx.HTTPStatusError) else type(exc).__name__
    result={'reading_level':level,'passages':paras,'fulltext_error':error,
       'fulltext_scope':'Gene-matching body passages only, not complete original-paper or supplementary-data certification.'}
    # Gene-dependent passages must never be cached under a paper-only key.
    c=db()
    try:
        c.execute('INSERT OR REPLACE INTO read_cache VALUES (?,?,?)',(cache_key,json.dumps(result,ensure_ascii=False),time.time()));c.commit()
    finally:c.close()
    return {**r,**result}

def as_fact(r):
    """A publication record; experiment clues are not a fabricated regulation edge."""
    r=classify_study(r)
    return {'evidence_id':'literature:'+r['paper_key'],'evidence_type':'PUBLISHED_RESEARCH' if not r.get('is_secondary') else 'PUBLISHED_REVIEW',
       'subject_node_id':'gene:'+r['gene_id'],'subject_label':node_info(r['gene_id'])['symbol'],
       'predicate':'DISCUSSED_IN_PUBLICATION','object_node_id':'publication:'+r['paper_key'],'object_label':r['title'],
       'source_url':r['source_url'],'category':'published_literature','origin':'online_literature','focus_role':'subject','independence_status':'not_certified',
       'details':{k:r.get(k) for k in ['paper_key','title','doi','pmid','year','journal','authors','abstract','reading_level','identity_match','species_check','experiment_clues','is_review','is_commentary','is_secondary','related_primary_dois','is_preprint','retracted','curation_note','curation_limitations','interpretation','passages']},
       'caveat':'论文命中/摘要读取不等于实验结论认证；同源或同名需核对；撤稿记录不能支持正向结论。'}
