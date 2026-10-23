"""Intranet portal with bounded public-evidence agent and immutable model lookup."""
import asyncio,hashlib,json,os,secrets,sqlite3,time
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal
from urllib.parse import urlsplit
from fastapi import FastAPI,Request,HTTPException
from fastapi.responses import FileResponse,JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel,Field,ConfigDict
from data_access import ROOT,CONFIG,MANIFEST,gene_response,resolve_gene,verify_release,node_info,sha256
from conversation_agent import run_agent

STATE=ROOT/'state/agent.sqlite';TASKS=set();WORKER=asyncio.Semaphore(1);RELEASE_INFO={}
LITERATURE_TASKS={};LITERATURE_WORKER=asyncio.Semaphore(2)

def db():
    c=sqlite3.connect(STATE,timeout=10);c.row_factory=sqlite3.Row
    return c

def init_state():
    STATE.parent.mkdir(exist_ok=True,mode=0o700)
    c=db()
    c.executescript('''CREATE TABLE IF NOT EXISTS jobs(
      job_id TEXT PRIMARY KEY, cache_key TEXT NOT NULL, gene_id TEXT NOT NULL,
      question TEXT NOT NULL, ip_hash TEXT NOT NULL, created REAL NOT NULL, updated REAL NOT NULL,
      status TEXT NOT NULL, events_json TEXT NOT NULL, result_json TEXT, error_code TEXT);
      CREATE INDEX IF NOT EXISTS jobs_cache ON jobs(cache_key,status);
      CREATE INDEX IF NOT EXISTS jobs_time ON jobs(created);
      CREATE TABLE IF NOT EXISTS conversations(conversation_id TEXT PRIMARY KEY,gene_id TEXT,created REAL,updated REAL);
      CREATE TABLE IF NOT EXISTS exchanges(conversation_id TEXT,job_id TEXT UNIQUE,turn_index INTEGER,
          PRIMARY KEY(conversation_id,turn_index));''')
    c.execute("UPDATE jobs SET status='interrupted',error_code='service_restarted',updated=? WHERE status IN ('queued','running')",(time.time(),))
    c.commit();c.close();os.chmod(STATE,0o600)

@asynccontextmanager
async def lifespan(app):
    global RELEASE_INFO
    os.umask(0o077)
    RELEASE_INFO=await asyncio.to_thread(verify_release)
    init_state()
    from literature import init as init_literature
    await asyncio.to_thread(init_literature)
    from research import init as init_research
    await asyncio.to_thread(init_research)
    yield
    for task in list(TASKS):task.cancel()
    for task in list(LITERATURE_TASKS.values()):task.cancel()
    if TASKS:await asyncio.gather(*TASKS,return_exceptions=True)

app=FastAPI(title='MKMGI × MEKG',docs_url=None,redoc_url=None,openapi_url=None,lifespan=lifespan)
app.mount('/static',StaticFiles(directory=ROOT/'static'),name='static')

@app.middleware('http')
async def security_headers(request,call_next):
    length=request.headers.get('content-length','0')
    try:large=int(length)>16000
    except ValueError:large=True
    if large:return JSONResponse({'detail':'request_too_large'},status_code=413)
    response=await call_next(request)
    response.headers['Content-Security-Policy']="default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'"
    response.headers['X-Content-Type-Options']='nosniff'
    response.headers['Referrer-Policy']='no-referrer'
    response.headers['X-Frame-Options']='DENY'
    response.headers['Permissions-Policy']='camera=(), microphone=(), geolocation=()'
    if request.url.path.startswith('/api/'):response.headers['Cache-Control']='no-store'
    return response

@app.get('/')
async def index():return FileResponse(ROOT/'static/submission.html',headers={'Cache-Control':'no-cache'})

@app.get('/api/health')
async def health():
    return {'status':'ready','release_id':CONFIG['release_id'],'model_gene_count':RELEASE_INFO.get('model_gene_count'),
       'model_frozen':True,'kg_read_only':True,'provider':'DeepSeek','agent_enabled':True,'endotwin_enabled':False,
       'online_literature_enabled':True,'network_provenance_enabled':True,'user_access_quotas':False,'evidence_total_cap':None,
       'conversation_enabled':True,'question_aware_retrieval':True,'claim_level_source_audit':True,
       'portal_release':'submission_20261022','interface_languages':['en','zh'],'default_language':'en'}

@app.get('/api/info')
async def info():
    return {'title':'MKMGI × MEKG','coordinate_system':'B73 v5','model_gene_count':44303,'deployment_seeds':3,
      'agent_provider':'DeepSeek','model_release':Path(CONFIG['model_root']).name,'kg_release':'MEKG public v3.1',
      'score_semantics':'候选基因优先级，不是表型概率或新候选的实验成功率。',
      'privacy':'提问、最近对话与选中的公共文献会发送给DeepSeek；经授权，仅附加当前基因分数、排名、覆盖和时间权重摘要。不发送原始组学、其他基因预测或未发表实验数据。',
      'limits':{'daily_requests':CONFIG['daily_run_limit'],'hourly_requests_per_client':CONFIG['hourly_ip_run_limit'],
                'evidence_total_cap':None,'pagination_is_not_total_cap':True,
                'operational_safeguards':'Background timeouts, bounded model context/concurrency and source API rate protections retained.'},
      'literature_sources':CONFIG['online_sources'],'literature_coverage':'Multi-source discovery with recorded aliases; not guaranteed exhaustive.',
      'integrity':RELEASE_INFO,'endotwin_enabled':False}

async def lookup(identifier):
    if len(identifier)>64:raise HTTPException(400,'invalid_identifier')
    result=await asyncio.to_thread(gene_response,identifier)
    if result['status']=='ambiguous':
        result['candidates']=[node_info(g) for g in result['candidate_ids']]
        return JSONResponse(result,status_code=409)
    if result['status']!='ok':raise HTTPException(404 if result['status']=='not_found' else 400,result['status'])
    return result

@app.get('/api/gene/{identifier}')
async def gene(identifier:str):return await lookup(identifier)

class Ask(BaseModel):
    model_config=ConfigDict(extra='forbid')
    gene_id:str=Field(min_length=2,max_length=64)
    question:str=Field(min_length=2,max_length=3000)
    conversation_id:str|None=Field(default=None,min_length=20,max_length=100)
    language:Literal['en','zh']='en'

def authorize_post(request):
    if request.headers.get('x-mekg-request')!='1':raise HTTPException(403,'request_header_required')
    origin=request.headers.get('origin')
    if origin:
        u=urlsplit(origin)
        if u.netloc.lower()!=request.headers.get('host','').lower() or u.scheme!=request.url.scheme:
            raise HTTPException(403,'cross_origin_request_denied')
    # No CORS; require a browser same-origin request or an intentional local API client.

def reserve_job(gene,question,ip,conversation_id=None,language='en'):
    stamp=time.time()
    from literature import fingerprint
    manifest=json.loads(MANIFEST.read_text())
    key=hashlib.sha256(json.dumps([gene,question,manifest['prediction_sha256'],manifest['kg_sha256'],
      CONFIG['prompt_version'],CONFIG['model'],fingerprint(gene),conversation_id,language,
      *[sha256(ROOT/name) for name in ['conversation_agent.py','research_tools.py','research.py','public_retrieval.py','public_repository.py','literature.py']]],ensure_ascii=False).encode()).hexdigest()
    iphash=hashlib.sha256((CONFIG['release_id']+'|'+ip).encode()).hexdigest()
    c=db()
    try:
        c.execute('BEGIN IMMEDIATE')
        if conversation_id:
            row=c.execute('SELECT * FROM conversations WHERE conversation_id=?',(conversation_id,)).fetchone()
            if not row or row['gene_id']!=gene:raise HTTPException(404,'conversation_not_found')
            pending_turn=c.execute("SELECT j.job_id,j.question,j.status FROM exchanges x JOIN jobs j ON j.job_id=x.job_id WHERE x.conversation_id=? AND j.status IN ('queued','running')",(conversation_id,)).fetchone()
            if pending_turn:
                if pending_turn['question']==question:
                    c.rollback();return {'job_id':pending_turn['job_id'],'conversation_id':conversation_id,'status':pending_turn['status'],'cached':False}
                raise HTTPException(409,'conversation_busy')
        else:
            conversation_id=secrets.token_urlsafe(32)
            c.execute('INSERT INTO conversations VALUES (?,?,?,?)',(conversation_id,gene,stamp,stamp))
        turn=c.execute('SELECT count(*) FROM exchanges WHERE conversation_id=?',(conversation_id,)).fetchone()[0]
        key=hashlib.sha256((key+'|'+conversation_id+'|'+str(turn)).encode()).hexdigest()
        cached=c.execute("SELECT job_id FROM jobs WHERE cache_key=? AND status='complete' AND updated>? ORDER BY updated DESC LIMIT 1",(key,stamp-CONFIG['cache_days']*86400)).fetchone()
        if cached:c.rollback();return {'job_id':cached['job_id'],'cached':True,'status':'complete'}
        pending=c.execute("SELECT job_id,status FROM jobs WHERE cache_key=? AND status IN ('queued','running') LIMIT 1",(key,)).fetchone()
        if pending:c.rollback();return {'job_id':pending['job_id'],'cached':False,'status':pending['status']}
        daily=c.execute('SELECT count(*) FROM jobs WHERE created>?',(stamp-86400,)).fetchone()[0]
        hourly=c.execute('SELECT count(*) FROM jobs WHERE ip_hash=? AND created>?',(iphash,stamp-3600)).fetchone()[0]
        active=c.execute("SELECT count(*) FROM jobs WHERE status IN ('queued','running')").fetchone()[0]
        if CONFIG.get('daily_run_limit') is not None and daily>=CONFIG['daily_run_limit']:raise HTTPException(429,'daily_agent_limit')
        if CONFIG.get('hourly_ip_run_limit') is not None and hourly>=CONFIG['hourly_ip_run_limit']:raise HTTPException(429,'hourly_agent_limit')
        if active>=CONFIG['max_active_jobs']:raise HTTPException(429,'agent_busy')
        token=secrets.token_urlsafe(24)
        c.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?,?)',(token,key,gene,question,iphash,stamp,stamp,'queued','[]',None,None))
        c.execute('INSERT INTO exchanges VALUES (?,?,?)',(conversation_id,token,turn))
        c.execute('UPDATE conversations SET updated=? WHERE conversation_id=?',(stamp,conversation_id))
        c.commit();return {'job_id':token,'conversation_id':conversation_id,'cached':False,'status':'queued'}
    finally:c.close()

def change_job(token,**values):
    allowed={'status','events_json','result_json','error_code','updated'}
    assert set(values)<=allowed
    c=db()
    try:
        c.execute('UPDATE jobs SET '+','.join(k+'=?' for k in values)+' WHERE job_id=?',list(values.values())+[token]);c.commit()
    finally:c.close()

def conversation_history(conversation_id):
    c=db()
    try:
        rows=c.execute('''SELECT j.question,j.result_json FROM exchanges x JOIN jobs j ON j.job_id=x.job_id
            WHERE x.conversation_id=? AND j.status='complete' ORDER BY x.turn_index DESC LIMIT 5''',(conversation_id,)).fetchall()
        out=[]
        for r in reversed(rows):
            result=json.loads(r['result_json']);answer=result['answer']
            out.append({'question':r['question'],'answer_summary':answer['summary']['text'],
               'answer_paragraphs':[u['text'] for s in answer['sections'] for u in s['items']][:12],
               'paper_keys':result.get('audit',{}).get('read_paper_keys',[])})
        return out
    finally:c.close()

async def work(token,gene,question,conversation_id,language='en'):
    events=[]
    async def event(item):
        # Display only operation summaries; never provider reasoning or credentials.
        safe={k:item[k] for k in ['stage','message','tool','status'] if k in item}
        events.append(safe)
        await asyncio.to_thread(change_job,token,events_json=json.dumps(events,ensure_ascii=False),updated=time.time())
    try:
        async with WORKER:
            await asyncio.to_thread(change_job,token,status='running',updated=time.time())
            history=await asyncio.to_thread(conversation_history,conversation_id)
            answer=await asyncio.wait_for(run_agent(gene,question,event,history,token,language=language),CONFIG['request_deadline_seconds'])
            answer['conversation_id']=conversation_id
            frozen=await asyncio.to_thread(gene_response,gene)
            answer['gene_result']=frozen
            answer['provenance']={'release_id':CONFIG['release_id'],**RELEASE_INFO,'generated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
               'model_root':CONFIG['model_root'],'prediction_file':CONFIG['prediction_file'],'kg_database':CONFIG['kg_database']}
            await asyncio.to_thread(change_job,token,status='complete',result_json=json.dumps(answer,ensure_ascii=False,allow_nan=False),updated=time.time())
    except asyncio.CancelledError:
        await asyncio.to_thread(change_job,token,status='interrupted',error_code='service_restarted',updated=time.time());raise
    except Exception as exc:
        import traceback
        stack=[{'file':Path(f.filename).name,'line':f.lineno,'function':f.name} for f in traceback.extract_tb(exc.__traceback__)]
        print(json.dumps({'event':'agent_exception','type':type(exc).__name__,'stack':stack}),flush=True)
        code=str(exc) if isinstance(exc,RuntimeError) and str(exc).replace('_','').isalnum() else 'agent_unavailable'
        if isinstance(exc,asyncio.TimeoutError):code='agent_timeout'
        await asyncio.to_thread(change_job,token,status='failed',error_code=code,updated=time.time())

@app.post('/api/ask')
async def ask(body:Ask,request:Request):
    authorize_post(request)
    question=body.question.strip()
    if len(question)<2:raise HTTPException(400,'question_too_short')
    resolved=await asyncio.to_thread(resolve_gene,body.gene_id)
    if resolved['status']!='ok':raise HTTPException(400,resolved['status'])
    job=await asyncio.to_thread(reserve_job,resolved['gene_id'],question,request.client.host if request.client else 'unknown',body.conversation_id,body.language)
    if not job['cached'] and job['status']=='queued':
        # A duplicate queued request may already own a task; single-process server plus DB reservation.
        if not any(getattr(t,'job_id',None)==job['job_id'] for t in TASKS):
            task=asyncio.create_task(work(job['job_id'],resolved['gene_id'],question,job['conversation_id'],body.language));task.job_id=job['job_id']
            TASKS.add(task);task.add_done_callback(TASKS.discard)
    return job

@app.get('/api/jobs/{token}')
async def job_status(token:str):
    if len(token)>100:raise HTTPException(404,'job_not_found')
    c=db()
    try:r=c.execute('SELECT * FROM jobs WHERE job_id=?',(token,)).fetchone()
    finally:c.close()
    if not r:raise HTTPException(404,'job_not_found')
    result={'job_id':token,'status':r['status'],'events':json.loads(r['events_json']),'error_code':r['error_code']}
    if r['status']=='complete':result['result']=json.loads(r['result_json'])
    return result

@app.get('/api/download/gene/{identifier}')
async def download(identifier:str):
    result=await lookup(identifier)
    if isinstance(result,JSONResponse):return result
    result['provenance']={'release_id':CONFIG['release_id'],**RELEASE_INFO,
      'prediction_file':CONFIG['prediction_file'],'kg_database':CONFIG['kg_database']}
    return JSONResponse(result,headers={'Content-Disposition':'attachment; filename="'+result['gene_id']+'_MEKG.json"'})

async def canonical(identifier):
    r=await asyncio.to_thread(resolve_gene,identifier)
    if r['status']!='ok':raise HTTPException(400 if r['status']=='ambiguous' else 404,r['status'])
    return r['gene_id']

def check_page(offset,page_size):
    if offset<0 or not 1<=page_size<=500:raise HTTPException(400,'invalid_pagination')

@app.get('/api/evidence/{identifier}')
async def evidence(identifier:str,offset:int=0,page_size:int=40,focus:str='all'):
    check_page(offset,page_size)
    if focus not in ['all','genetics','molecular','regulation','function','traits']:raise HTTPException(400,'invalid_focus')
    from public_repository import evidence_page
    return await asyncio.to_thread(evidence_page,await canonical(identifier),offset,page_size,focus)

@app.get('/api/network/{identifier}')
async def network(identifier:str,offset:int=0,page_size:int=60,layers:str='all',direction:str='all'):
    check_page(offset,page_size)
    from network import network_page,LEGEND
    if layers!='all' and not set(layers.split(','))<=set(LEGEND):raise HTTPException(400,'invalid_layer')
    if direction not in ['all','outgoing','incoming']:raise HTTPException(400,'invalid_direction')
    return await asyncio.to_thread(network_page,await canonical(identifier),offset,page_size,layers,direction)

@app.get('/api/literature/{identifier}')
async def literature(identifier:str,offset:int=0,page_size:int=20,filter:str='all'):
    check_page(offset,page_size)
    if filter not in ['all','experimental','classic','review','warning']:raise HTTPException(400,'invalid_filter')
    from literature import papers_page
    return await asyncio.to_thread(papers_page,await canonical(identifier),offset,page_size,filter)

@app.get('/api/research/{identifier}')
async def research_overview(identifier:str):
    gene=await canonical(identifier)
    from research import relation_facts,ranked_papers
    rels=await asyncio.to_thread(relation_facts,gene)
    papers=await asyncio.to_thread(ranked_papers,gene,'这个基因有哪些已发表的遗传功能和调控关系？',0,6)
    return {'gene_id':gene,'relations':rels,'papers':papers['papers'],'paper_total':papers['total'],
       'relation_scope':'Curated primary research plus source-checked author reports; distinct from model-input edges.'}

@app.get('/api/conversation/{conversation_id}')
async def conversation(conversation_id:str):
    if not re_token(conversation_id):raise HTTPException(404,'conversation_not_found')
    c=db()
    try:
        row=c.execute('SELECT gene_id FROM conversations WHERE conversation_id=?',(conversation_id,)).fetchone()
        if not row:raise HTTPException(404,'conversation_not_found')
        jobs=c.execute('SELECT j.job_id,j.question,j.status FROM exchanges x JOIN jobs j ON j.job_id=x.job_id WHERE x.conversation_id=? ORDER BY x.turn_index',(conversation_id,)).fetchall()
        return {'conversation_id':conversation_id,'gene_id':row['gene_id'],'turns':[dict(r) for r in jobs]}
    finally:c.close()

def re_token(value):
    import re
    return bool(re.fullmatch(r'[A-Za-z0-9_-]{20,100}',value))

class SearchLiterature(BaseModel):
    model_config=ConfigDict(extra='forbid')
    refresh:bool=False

async def literature_work(gene,refresh):
    from literature import run_search
    async with LITERATURE_WORKER:
        await asyncio.to_thread(run_search,gene,refresh)

@app.post('/api/literature/{identifier}/search')
async def start_literature(identifier:str,body:SearchLiterature,request:Request):
    authorize_post(request);gene=await canonical(identifier)
    from literature import status
    active=LITERATURE_TASKS.get(gene)
    if active and not active.done():return {'gene_id':gene,'status':'running'}
    # Concurrent-resource guard, not per-user/day access quota.
    if sum(not x.done() for x in LITERATURE_TASKS.values())>=CONFIG['max_active_jobs']:
        raise HTTPException(429,'literature_busy')
    task=asyncio.create_task(literature_work(gene,body.refresh));LITERATURE_TASKS[gene]=task
    def finished(t):
        if LITERATURE_TASKS.get(gene) is t:LITERATURE_TASKS.pop(gene,None)
        if not t.cancelled():t.exception()  # Consume exceptions; no source text or secrets in logs.
    task.add_done_callback(finished)
    return {'gene_id':gene,'status':'queued','search':await asyncio.to_thread(status,gene)}

@app.get('/api/download/evidence/{identifier}')
async def download_evidence(identifier:str):
    from public_repository import all_public_evidence
    gene=await canonical(identifier);r=await asyncio.to_thread(all_public_evidence,gene)
    return JSONResponse(r,headers={'Content-Disposition':'attachment; filename="'+gene+'_all_public_evidence.json"'})

@app.get('/api/download/network/{identifier}')
async def download_network(identifier:str):
    from network import full_network
    gene=await canonical(identifier);r=await asyncio.to_thread(full_network,gene)
    return JSONResponse(r,headers={'Content-Disposition':'attachment; filename="'+gene+'_full_network.json"'})

@app.get('/api/download/literature/{identifier}')
async def download_literature(identifier:str):
    from literature import papers_page,status
    gene=await canonical(identifier);s=await asyncio.to_thread(status,gene)
    r=await asyncio.to_thread(papers_page,gene,0,max(s['paper_count'],1))
    return JSONResponse(r,headers={'Content-Disposition':'attachment; filename="'+gene+'_literature.json"'})

@app.get('/api/download/report/{token}')
async def report_download(token:str):
    r=await job_status(token)
    if r['status']!='complete':raise HTTPException(409,'report_not_complete')
    return JSONResponse(r['result'],headers={'Content-Disposition':'attachment; filename="'+r['result']['gene_id']+'_MEKG_agent.json"'})
