'use strict';
const $=id=>document.getElementById(id);
const node=(tag,text,cls)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=String(text);if(cls)n.className=cls;return n;};
const wipe=n=>n.replaceChildren();
const pause=ms=>new Promise(r=>setTimeout(r,ms));
let current=null,serial=0,chatSerial=0,conversationId=null,networkOffset=null,literatureOffset=null,evidenceOffset=null;
let literatureFilter='all',evidenceFilter='all',activeTab='chat',networkEpoch=0,sourceEpoch=0;
let networkLayers=new Set(['classic_experimental','literature_reported','public_experimental','ppi','model_regulatory','project_inferred']);
const errorNames={not_found:'未找到这个ID或已记录别名。请尝试B73 v5编号或旧ID。',invalid_identifier:'请输入合法的基因ID或symbol。',agent_busy:'当前有研究任务正在处理，请稍后再试。',conversation_busy:'本轮尚未结束，请等回答完成后继续追问。',conversation_not_found:'这段对话已失效，请点击“新对话”重新开始。',agent_timeout:'本轮检索或解释超时，尚未生成完整回答。可缩小问题或指定一篇论文后重试。',answer_validation_failed:'本轮回答未通过格式与引用检查，没有发布未核实内容。',agent_unavailable:'本轮研究暂不可用；模型结果、文献和网络仍可查询。',provider_unavailable:'DeepSeek暂时不可用，请稍后重试。',literature_busy:'联网检索任务正在处理中，可先浏览已找到的研究。',service_restarted:'服务重启中断了这轮回答，请重新提交。'};
function errorText(code){return errorNames[code]||(String(code).startsWith('provider_')?'DeepSeek本次未完成回答，请稍后重试。':'请求未完成，请稍后重试。');}
async function request(url,options){const r=await fetch(url,options);let data;try{data=await r.json();}catch{throw Error('agent_unavailable');}if(!r.ok){const e=Error(typeof data.detail==='string'?data.detail:'agent_unavailable');e.data=data;e.status=r.status;throw e;}return data;}
const post=(url,body)=>request(url,{method:'POST',headers:{'Content-Type':'application/json','X-MEKG-Request':'1'},body:JSON.stringify(body)});
function link(url,label,cls){const a=node('a',label,cls);try{const u=new URL(url);if(u.protocol!=='https:'&&!url.startsWith('/api/'))throw Error();a.href=url;a.target='_blank';a.rel='noopener noreferrer';}catch{a.removeAttribute('href');}return a;}
function sourceTitle(f){return f.origin==='current_gene_model_summary'?'当前基因 MKMGI 模型摘要':f.details?.title||f.details?.study_title||f.object_label||f.evidence_type;}
function shortQuote(text){return String(text||'').split(/\s+/).slice(0,24).join(' ');}
function showSource(f){
 const host=$('source-content');wipe(host);const d=f.details||{};
 host.append(node('h4',sourceTitle(f)));
 if(f.origin==='current_gene_model_summary'){
  host.append(node('p','当前基因的冻结模型摘要。它支持分数、排名和数据覆盖的描述，不是遗传或因果证据。'));
  host.append(node('pre',JSON.stringify(d,null,2),'record-json'));
 }else{
  host.append(node('p',[d.year,d.journal,f.genotype,f.development_stage,f.tissue_or_compartment].filter(Boolean).join(' · '),'muted'));
  const reading=d.reading_level==='open_fulltext_gene_passages'?'已读取可用开放原文的基因相关段落':d.reading_level==='abstract'?'已读取摘要；未宣称完成逐图审计':d.reading_level==='metadata'?'仅元数据，不能单独支持实验结论':f.evidence_type==='PRIMARY_LITERATURE_RELATION'?'原始研究关系记录，具体支持层级见下方':'本地公共证据记录';
  host.append(node('p',reading));
  if(d.finding)host.append(node('p',d.finding));
  if(d.assays?.length)host.append(node('p','实验：'+(Array.isArray(d.assays)?d.assays.join('、'):d.assays)));
  if(d.support_level)host.append(node('p','支持层级：'+supportName(d.support_level)));
  if(d.quote)host.append(node('p',shortQuote(d.quote),'source-quote'));
  if(d.retracted)host.append(node('p','撤稿记录，不可作为正向机制支持。','error'));
  if(f.caveat)host.append(node('p',f.caveat,'muted'));
  if(d.limitations)host.append(node('p',d.limitations,'muted'));
  if(f.source_url)host.append(link(f.source_url,'打开原始来源 ↗'));
  const detail=node('details');detail.append(node('summary','证据记录字段'));const compact={evidence_id:f.evidence_id,type:f.evidence_type,predicate:f.predicate,subject:f.subject_node_id,object:f.object_node_id,details:{...d}};
  for(const key of ['abstract','passages','selected_passages'])delete compact.details[key];
  detail.append(node('pre',JSON.stringify(compact,null,2),'record-json'));host.append(detail);
 }
 $('source-dialog').showModal();
}
function supportName(k){return {binding_and_functional:'结合与功能实验支持',promoter_functional:'启动子功能实验，结合证据需另核对',binding_only:'结合证据，不单独证明调控',perturbation_response:'扰动后的响应，非直接调控认证',high_throughput_candidate:'高通量候选靶标',protein_interaction:'蛋白互作，不等于转录调控'}[k]||k;}
function citations(host,refs,map){for(const id of refs||[]){const f=map.get(id);if(!f)continue;const b=node('button',id,'citation');b.type='button';b.title=sourceTitle(f);b.addEventListener('click',()=>showSource(f));host.append(b);}}

async function searchGene(identifier){
 const id=++serial;++chatSerial;current=null;conversationId=null;$('workspace').hidden=true;$('search-button').disabled=true;$('search-status').textContent='正在核对基因身份，读取冻结模型与公共研究…';wipe($('alias-choices'));
 try{
  const data=await request('/api/gene/'+encodeURIComponent(identifier));if(id!==serial)return;
  current=data;renderGene(data);$('workspace').hidden=false;$('search-status').textContent='';resetChat(false);setTab('chat');await restoreConversation(id);
  await Promise.allSettled([loadOverview(id),loadNetwork(true),loadLiterature(true,id),loadEvidence(true,id)]);
  if(id===serial){const u=new URL(location.href);u.searchParams.set('gene',data.gene_id);history.replaceState(null,'',u);}
 }catch(e){if(id!==serial)return;if(e.status===409&&e.data.candidates){$('search-status').textContent='这个别名有多个映射，请核对后选择：';for(const c of e.data.candidates){const b=node('button',c.gene_id+' · '+c.symbol,'subtle-button');b.addEventListener('click',()=>searchGene(c.gene_id));$('alias-choices').append(b);}}else $('search-status').textContent=errorText(e.message);}
 finally{if(id===serial)$('search-button').disabled=false;}
}
function renderGene(r){
 $('gene-symbol').textContent=r.gene.symbol||r.gene_id;$('gene-id').textContent=r.gene_id;$('gene-description').textContent=r.gene.description||'';
 $('gene-location').textContent=[r.gene.chromosome&&r.gene.start!==null?r.gene.chromosome+':'+r.gene.start+'–'+r.gene.end:null,'B73 v5',r.id_mapping==='exact_B73_v5'?'精确ID':'已记录别名映射'].filter(Boolean).join(' · ');
 for(const [id,path] of [['download-gene','gene'],['download-network','network'],['download-literature','literature'],['download-evidence','evidence']])$(id).href='/api/download/'+path+'/'+r.gene_id;
 const p=r.prediction;$('score').textContent=p?p.priority_score.toFixed(4):'未建模';$('percentile').textContent=p?(p.priority_percentile*100).toFixed(1)+' 百分位':'—';$('rank').textContent=p?p.rank.toLocaleString('en-US'):'—';$('rank-total').textContent=p?'/ '+p.model_gene_count.toLocaleString('en-US'):'';$('percentile-fill').style.width=p?p.priority_percentile*100+'%':'0%';
 wipe($('coverage'));if(p)for(const [label,yes] of [['时间组学 '+p.coverage.temporal_modalities_observed+'/5',p.coverage.temporal_modalities_observed>0],['空间',p.coverage.spatial_observed],['序列 / 功能',p.coverage.sequence_observed],['蛋白背景',p.coverage.protein_observed]])$('coverage').append(node('span',(yes?'✓ ':'— ')+label,yes?'':'missing'));
 $('stability').textContent=p?'三种子共识；分数 SD = '+p.uncertainty_sd.toFixed(4)+'。这是计算稳定性，不是置信区间。':'';
 $('training-status').textContent=p?.in_development_positive?'该基因是开发阶段已知正例。其部署分数不能当作新增独立验证。':p?.in_development_training?'该基因进入过开发阶段监督样本，不能用本页声称独立验证。':'全基因组部署输出；新候选的实验命中率尚未知。';
 wipe($('time-bars'));wipe($('relation-gates'));const prof=r.model_profile;
 if(prof){const max=Math.max(...prof.temporal_pooling_weights,.01);prof.DAP.forEach((dap,i)=>{const item=node('div',undefined,'time-item'+(prof.temporal_observed[i]?'':' unavailable'));const bar=node('div',undefined,'bar');bar.style.height=prof.temporal_pooling_weights[i]/max*53+'px';item.append(node('span',prof.temporal_pooling_weights[i].toFixed(2),'weight'),bar,node('span',dap));$('time-bars').append(item);});
  const h=node('div',undefined,'gate-row');for(const text of ['层','自身','PPI','TF→靶','靶→TF'])h.append(node('span',text));$('relation-gates').append(h);prof.relation_gates.forEach((v,i)=>{const row=node('div',undefined,'gate-row');row.append(node('span',i+1));v.forEach(w=>row.append(node('span',w.toFixed(2))));$('relation-gates').append(row);});
 }
 literatureFilter='all';evidenceFilter='all';literatureOffset=null;evidenceOffset=null;networkOffset=null;
 $('evidence-traits-filter')?.remove();
 document.querySelectorAll('[data-filter]').forEach(b=>b.classList.toggle('selected',b.dataset.filter==='all'));document.querySelectorAll('[data-focus]').forEach(b=>b.classList.toggle('selected',b.dataset.focus==='all'));
}
const welcome=$('conversation').firstElementChild.cloneNode(true);
function bindQuestions(host){host.querySelectorAll('[data-question]').forEach(b=>b.addEventListener('click',()=>{$('question').value=b.dataset.question;$('question').focus();}));}
function resetChat(newClick=true){++chatSerial;conversationId=null;$('conversation').replaceChildren(welcome.cloneNode(true));bindQuestions($('conversation'));$('question').value='';$('agent-progress').hidden=true;$('agent-error').hidden=true;$('ask-button').disabled=false;wipe($('tool-trace'));if(newClick){if(current)sessionStorage.removeItem('mekg-conversation:'+current.gene_id);$('question').focus();}}
async function restoreConversation(stamp){
 const gene=current.gene_id,stored=sessionStorage.getItem('mekg-conversation:'+gene),chatStamp=chatSerial;if(!stored)return;
 try{const info=await request('/api/conversation/'+encodeURIComponent(stored));if(stamp!==serial||chatStamp!==chatSerial)return;conversationId=stored;const host=$('conversation');host.querySelector('.welcome')?.remove();
  for(const turn of info.turns){if(stamp!==serial||chatStamp!==chatSerial)return;host.append(node('div',turn.question,'user-message'));if(turn.status==='complete'){const job=await request('/api/jobs/'+turn.job_id);if(stamp!==serial||chatStamp!==chatSerial)return;renderAnswer(job.result,turn.job_id);}else if(['queued','running'].includes(turn.status)){resumePending(turn.job_id,stamp,chatStamp);}}
 }catch{sessionStorage.removeItem('mekg-conversation:'+gene);}
}
async function resumePending(jobId,stamp,chatStamp){
 $('ask-button').disabled=true;$('agent-progress').hidden=false;
 try{for(let i=0;i<360;i++){const state=await request('/api/jobs/'+jobId);if(stamp!==serial||chatStamp!==chatSerial)return;renderEvents(state.events);if(state.status==='complete'){renderAnswer(state.result,jobId);$('agent-progress').hidden=true;return;}if(['failed','interrupted'].includes(state.status))throw Error(state.error_code||'agent_unavailable');await pause(1500);}}catch(e){if(stamp===serial&&chatStamp===chatSerial){$('agent-progress').hidden=true;$('agent-error').textContent=errorText(e.message);$('agent-error').hidden=false;}}finally{if(stamp===serial&&chatStamp===chatSerial)$('ask-button').disabled=false;}
}
function setTab(name){activeTab=name;for(const tab of ['chat','network','sources']){$('panel-'+tab).hidden=tab!==name;$('tab-'+tab).classList.toggle('active',tab===name);$('tab-'+tab).setAttribute('aria-selected',String(tab===name));}if(name==='network')requestAnimationFrame(()=>explorer.draw());}
async function askAgent(){
 if(!current)return;const text=$('question').value.trim();if(text.length<2)return;const stamp=++chatSerial,gene=current.gene_id;
 const chat=$('conversation');chat.querySelector('.welcome')?.remove();chat.append(node('div',text,'user-message'));chat.scrollTop=chat.scrollHeight;$('question').value='';$('ask-button').disabled=true;$('agent-error').hidden=true;$('agent-progress').hidden=false;$('progress-label').textContent='正在理解问题并制定检索方向';wipe($('tool-trace'));
 try{
  const job=await post('/api/ask',{gene_id:gene,question:text,conversation_id:conversationId});if(stamp!==chatSerial)return;conversationId=job.conversation_id;sessionStorage.setItem('mekg-conversation:'+gene,conversationId);
  for(let i=0;i<360;i++){
   const state=await request('/api/jobs/'+job.job_id);if(stamp!==chatSerial)return;renderEvents(state.events);
   if(state.status==='complete'){$('agent-progress').hidden=true;renderAnswer(state.result,job.job_id);await Promise.allSettled([loadOverview(serial),loadNetwork(true)]);return;}
   if(['failed','interrupted'].includes(state.status))throw Error(state.error_code||'agent_unavailable');await pause(1500);
  }
  throw Error('agent_timeout');
 }catch(e){if(stamp===chatSerial){$('agent-error').textContent=errorText(e.message);$('agent-error').hidden=false;$('agent-progress').hidden=true;$('question').value=text;}}
 finally{if(stamp===chatSerial)$('ask-button').disabled=false;}
}
function renderEvents(events){wipe($('tool-trace'));for(const e of events||[])$('tool-trace').append(node('li',e.message));if(events?.length)$('progress-label').textContent=events.at(-1).message;}
function renderAnswer(result,job){
 const a=result.answer,map=new Map(result.evidence.map(f=>[f.citation,f]));const msg=node('article',undefined,'assistant-message');
 const label=node('div',undefined,'assistant-label');label.append(node('span','✦'),node('span','MEKG · 研究回答'));msg.append(label,node('h4',a.title,'answer-title'));
 const lead=node('div',a.summary.text,'answer-summary');citations(lead,a.summary.citations,map);msg.append(lead);
 for(const section of a.sections){const block=node('section',undefined,'answer-section');block.append(node('h5',section.heading));for(const u of section.items){const p=node('p',undefined,'answer-paragraph');const names={hypothesis:'待验证假设',model_context:'模型线索',recommendation:'实验建议'};if(names[u.kind])p.append(node('span',names[u.kind],'paragraph-tag'+(u.kind==='hypothesis'?' hypothesis':'')));p.append(document.createTextNode(u.text));citations(p,u.citations,map);block.append(p);}msg.append(block);}
 const refs=node('details');refs.append(node('summary','本轮引用来源（'+map.size+'项）'));const list=node('div',undefined,'source-list');for(const [id,f] of map){const r=node('div',undefined,'source-mini');const b=node('button',id);b.addEventListener('click',()=>showSource(f));r.append(b);const title=sourceTitle(f);const t=f.source_url?link(f.source_url,title):node('span',title);r.append(t,node('small',f.details?.year||''));list.append(r);}refs.append(list);msg.append(refs);
 const audit=node('details');audit.append(node('summary','查看本轮来源核对与检索记录'));audit.append(node('p','逐段检查引用、主体、实验层级和模型数值。自动核对由大模型完成，不属于独立专家评审。','muted'));
 const sources=node('ul');for(const t of result.search.trace||[])sources.append(node('li',t.source+' · '+(t.status==='ok'?'检索完成':'本轮不可用：'+t.error)+' · '+t.query));audit.append(sources);
 if(result.audit.removed_items.length){audit.append(node('p','以下内容未纳入答案，可查看具体原因：','muted'));const ul=node('ul');for(const x of result.audit.removed_items)ul.append(node('li',(x.item.text||x.item.finding||'关系项')+' — '+x.reason));audit.append(ul);}else audit.append(node('p','本轮展示的段落通过了记录支持性检查。仍建议关键结论回到原始文献核对。','muted'));
 msg.append(audit);const bottom=node('div',undefined,'answer-bottom');bottom.append(node('span','已发表结果、模型线索与假设分开；不把同源来源当独立验证。'));const dl=node('a','回答与证据 JSON ↓','text-link');dl.href='/api/download/report/'+job;bottom.append(dl);msg.append(bottom);
 const follow=node('div',undefined,'followups');for(const q of a.follow_up_questions){const b=node('button',q+' ↗');b.addEventListener('click',()=>{$('question').value=q;$('question').focus();});follow.append(b);}msg.append(follow);$('conversation').append(msg);msg.scrollIntoView({behavior:'smooth',block:'nearest'});
}
async function loadOverview(stamp){
 const data=await request('/api/research/'+current.gene_id);if(stamp!==serial)return;const host=$('classic-findings');wipe(host);
 const curated=data.relations.filter(r=>r.origin==='curated_primary_literature');const rest=data.relations.filter(r=>r.origin!=='curated_primary_literature');const groups=new Map();
 for(const r of [...curated,...rest]){const k=r.doi+'|'+r.finding;if(!groups.has(k))groups.set(k,r);}
 if(groups.size){for(const r of [...groups.values()].slice(0,6)){const row=node('article',undefined,'finding');row.append(node('span',r.origin==='curated_primary_literature'?'原始研究 · 已核对':'论文作者报告 · 自动核对','evidence-badge'+(r.origin==='curated_primary_literature'?'':' auto')));row.append(node('h4',r.source_symbol+(r.relation==='PROTEIN_INTERACTION'?' ↔ ':' → ')+r.target_symbol),node('p',r.finding),node('p',supportName(r.support_level),'muted'),link('https://doi.org/'+r.doi,'原始研究 ↗','text-link'));host.append(row);}}
 else if(data.papers.length){const note=node('p','已有论文检索线索；尚无人工核对的结构化关系。请在对话中指定问题，助手将读取原始研究并核对具体结论。','muted');host.append(note);for(const p of data.papers.slice(0,3)){const row=node('article',undefined,'finding');row.append(node('h4',undefined));row.firstChild.append(link(p.source_url,p.title));row.append(node('p',[p.year,p.journal].filter(Boolean).join(' · '),'muted'));host.append(row);}}
 else host.append(node('p','本地暂无结构化经典关系。可以通过上方提问，按旧ID、symbol或具体论文继续联网追查；无记录不代表没有功能。','muted'));
}

const explorer=new window.MEKGNetwork($('evidence-graph'),{onEdge:showEdge,onNode:showNode});
let latestLegend={};
async function loadNetwork(reset=false){
 if(!current)return;const stamp=serial,epoch=reset?++networkEpoch:networkEpoch;const gene=current.gene_id;if(reset){networkOffset=0;explorer.reset(gene);}
 const layers=[...networkLayers].join(',');if(!layers){$('network-status').textContent='请选择至少一种关系来源。';explorer.host.replaceChildren();return;}
 const direction=$('network-direction').value;const data=await request('/api/network/'+gene+'?offset='+(networkOffset||0)+'&page_size=24&layers='+encodeURIComponent(layers)+'&direction='+direction);
 if(stamp!==serial||epoch!==networkEpoch)return;latestLegend=data.legend;explorer.add(data);for(const layer of Object.keys(data.legend))explorer.filter(layer,networkLayers.has(layer));networkOffset=data.next_offset;
 $('network-status').textContent='当前视图 '+explorer.edges.size+' / '+data.total+' 条关系 · 可继续加载，完整网络可下载。';$('network-more').hidden=networkOffset===null;
 wipe($('network-legend'));for(const [layer,info] of Object.entries(data.legend)){const label=node('label');const input=node('input');input.type='checkbox';input.checked=networkLayers.has(layer);input.dataset.layer=layer;input.addEventListener('change',()=>{if(input.checked)networkLayers.add(layer);else networkLayers.delete(layer);loadNetwork(true).catch(e=>$('network-status').textContent=errorText(e.message));});const line=node('span',undefined,'legend-line'+(info.dash?' dashed':''));line.style.borderColor=info.color;label.append(input,line,document.createTextNode(info.label+' '+(data.layer_totals[layer]||0)));label.title=info.caution;$('network-legend').append(label);}
}
function showEdge(e,nodes){
 const host=$('network-inspector');wipe(host);const a=nodes.get(e.source)?.label||e.source,b=nodes.get(e.target)?.label||e.target;host.append(node('span',latestLegend[e.layer]?.label||e.layer,'evidence-badge'+(['model_regulatory','project_inferred','ppi'].includes(e.layer)?' model':'')),node('h4',a+(e.directed?' → ':' ↔ ')+b));
 host.append(node('p',e.support));if(e.support_level)host.append(node('p',supportName(e.support_level)));if(e.assays?.length)host.append(node('p','实验：'+e.assays.join('、')));
 const dl=node('dl');for(const [key,value] of [['关系类型',e.predicate],['原始研究',e.attributes?.study_title],['材料',e.genotype],['时期',e.stage],['组织',e.tissue],['进入模型',e.used_by_current_model===true?'是：原始输入候选网络，非实验认证':e.used_by_current_model===false?'否：文献解释层':'未逐边审计'],['坐标映射',e.source_entity&&e.target_entity?[e.source_entity.mapping_status,e.target_entity.mapping_status].join(' / '):null]]){if(value!==null&&value!==undefined&&value!==''){dl.append(node('dt',key),node('dd',Array.isArray(value)?value.join(' · '):value));}}host.append(dl);
 if(e.attributes?.limitations)host.append(node('p',e.attributes.limitations,'muted'));if(e.attributes?.support_quote)host.append(node('p',shortQuote(e.attributes.support_quote),'source-quote'));for(const url of e.source_urls||[])host.append(link(url,'原始来源 ↗','text-link'));
 const d=node('details');d.append(node('summary','完整关系字段'),node('pre',JSON.stringify(e,null,2),'record-json'));host.append(d);
}
function showNode(n){
 const host=$('network-inspector');wipe(host);host.append(node('h4',n.label),node('p',n.id,'muted'));if(n.description)host.append(node('p',n.description));
 if(n.id.startsWith('symbol:')){host.append(node('p','原文基因或基因集合，尚未确认唯一v5映射；没有编造ID。','muted'));return;}
 const expand=node('button','展开这个基因的邻居','subtle-button');expand.addEventListener('click',async()=>{expand.disabled=true;try{const data=await request('/api/network/'+n.id+'?page_size=24&layers='+encodeURIComponent([...networkLayers].join(',')));explorer.add(data);$('network-status').textContent='已展开 '+n.label+' 的邻居；本图共加载 '+explorer.edges.size+' 条关系。';}catch{host.append(node('p','该节点邻居暂不可用。','muted'));}finally{expand.disabled=false;}});
 const query=node('button','以此基因开始研究 ↗','subtle-button');query.addEventListener('click',()=>{$('gene-input').value=n.id;searchGene(n.id);});host.append(expand,document.createTextNode(' '),query);
}

function publication(p){const box=node('article',undefined,'publication');const badges=node('div',undefined,'paper-badges');if(p.retracted)badges.append(node('span','撤稿 · 不作正向支持','warning'));if(p.is_review)badges.append(node('span','综述'));if(p.is_preprint)badges.append(node('span','预印本'));if(p.identity_match==='indexed_search_only')badges.append(node('span','身份仍待原文核对'));for(const tag of p.experiment_clues||[])badges.append(node('span',tag));box.append(badges);const h=node('h4');h.append(link(p.source_url,p.title));box.append(h,node('p',[p.year,p.journal].filter(Boolean).join(' · '),'muted'));const authors=String(p.authors||'');box.append(node('p',authors.length>190?authors.slice(0,190)+'…':authors,'muted'));
 const d=node('details');d.append(node('summary','来源、身份与读取范围'));d.append(node('p','命中别名：'+((p.matched_aliases||[]).join('、')||'索引命中，需核对')),node('p','来源：'+(p.sources||[]).join(' / ')),node('p','DOI：'+(p.doi||'未提供')),node('p',p.reading_level==='metadata'?'仅元数据，未取得可读摘要。':'已取得可用摘要，具体实验结论由对话按问题核对。'));box.append(d);return box;}
async function loadLiterature(reset=false,stamp=serial){
 if(!current)return;const epoch=reset?++sourceEpoch:sourceEpoch,gene=current.gene_id;if(reset){literatureOffset=0;wipe($('literature-list'));}
 const data=await request('/api/literature/'+gene+'?offset='+(literatureOffset||0)+'&page_size=8&filter='+literatureFilter);if(stamp!==serial||epoch!==sourceEpoch)return;
 for(const p of data.papers)$('literature-list').append(publication(p));literatureOffset=data.next_offset;$('literature-more').hidden=literatureOffset===null;
 $('literature-status').textContent='本地联合文献索引 '+data.search.paper_count+' 条；问答会按当前问题重新排序、补查和阅读，不局限于此列表。';
 wipe($('literature-audit'));$('literature-audit').append(node('p','别名：'+(data.search.aliases||[]).join('、')));const ul=node('ul');for(const trace of (data.search.trace||[]).slice(-15))ul.append(node('li',trace.source+' · '+(trace.status==='ok'?'可用':'本轮失败：'+trace.error)+' · '+trace.query));$('literature-audit').append(ul);
 if(data.papers.length===0)$('literature-list').append(node('p','该筛选下暂无匹配文献，可继续联网补查或在对话中指定基因旧名和论文。','muted'));
}
async function continueLiterature(){if(!current)return;const stamp=serial,gene=current.gene_id;$('literature-continue').disabled=true;$('literature-status').textContent='正在追查基因别名和经典文献，来源状态会保留…';try{await post('/api/literature/'+gene+'/search',{refresh:false});for(let i=0;i<80;i++){await pause(2000);if(stamp!==serial)return;const data=await request('/api/literature/'+gene+'?page_size=1');if(data.search.status!=='running'){await loadLiterature(true,stamp);return;}}await loadLiterature(true,stamp);}catch(e){if(stamp===serial)$('literature-status').textContent=errorText(e.message);}finally{if(stamp===serial)$('literature-continue').disabled=false;}}
const evidenceTypes={MAIZE_GENETIC_FUNCTION:'遗传功能',CURATED_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE:'经典突变体',BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT:'结合与功能实验',MUTANT_TRANSCRIPTOMIC_RESPONSE:'遗传扰动后的转录响应',MUTANT_CHROMATIN_CHANGE:'遗传扰动后的染色质响应',TF_BINDING_PEAK_GENE_ASSOCIATION:'TF峰–基因关联',COMPUTATIONAL_INTEGRATED_REGULON:'计算整合调控网络',SINGLE_NUCLEUS_EXPRESSION_CONTEXT:'单核表达背景',DATABASE_FUNCTION_DESCRIPTION:'数据库功能描述',PANNZER_COMPUTATIONAL_GO:'计算GO注释'};
Object.assign(evidenceTypes,{GWAS_LOCUS_PROXIMITY:'GWAS · 位点邻近关联',GWAS_PROXIMITY_CONTEXT:'GWAS · 位点邻近关联',GWAS_REPORTED_GENE_TRAIT_ASSOCIATION:'GWAS · 作者报告基因–性状关联',GWAS_MARKER_TRAIT_ASSOCIATION:'GWAS · 标记–性状关联',GRAMENE_CURATED_QTL_TRAIT:'QTL · 区间与性状关联'});
function isTraitEvidence(f){return /GWAS|QTL/i.test(f.evidence_type||'');}
function traitStatistic(value){if(value===null||value===undefined||value==='')return '';const n=Number(value);return Number.isFinite(n)?n.toExponential(3):String(value);}
function appendTraitContext(row,f){
 const d={...(f.details||{}),...(f.trait_context||{})};
 row.append(node('p','关联性状：'+(d.trait||f.object_label||f.object_node_id||'原记录未提供')));
 if(d.title||d.study_title)row.append(node('p',d.title||d.study_title));
 if(d.population)row.append(node('p','群体：'+d.population,'muted'));
 if(d.study_id)row.append(node('p','研究编号：'+d.study_id,'muted'));
 if(d.reported_location)row.append(node('p','原报告位置：'+d.reported_location));
 const p=d.p_value??d.pvalue;if(p!==undefined&&p!==null&&p!=='')row.append(node('p','原报告 P = '+traitStatistic(p)));
 if(d.R2_percent!==undefined&&d.R2_percent!==null&&d.R2_percent!=='')row.append(node('p','原报告解释变异：'+d.R2_percent+'%'));
 const loci=Array.isArray(d.locus_results)?d.locus_results:[];
 if(loci.length){
  const first=loci[0];row.append(node('p','位点：'+(first.chromosome!==undefined?'chr'+first.chromosome+':':'')+(first.position??'位置未提供')+' · 原报告 P = '+traitStatistic(first.pvalue??first.p_value)));
  const records=node('details');records.append(node('summary','全部位点与统计量（'+loci.length+' 条原报告结果）'));
  for(const locus of loci)records.append(node('p',[(locus.chromosome!==undefined?'chr'+locus.chromosome+':':'')+(locus.position??'位置未提供'),'P = '+traitStatistic(locus.pvalue??locus.p_value),locus.model,locus.year].filter(Boolean).join(' · ')));
  row.append(records);
 }else if(d.chromosome!==undefined&&(d.start_position!==undefined||d.start!==undefined)){
  row.append(node('p','原报告区间：chr'+d.chromosome+':'+(d.start_position??d.start)+'–'+(d.stop_position??d.end??'未提供')+(d.coordinate_units?' ('+d.coordinate_units+')':'')));
 }
 if(d.window_bp!==undefined)row.append(node('p','关联规则：基因区间与位点 ±'+(Number(d.window_bp)/1000)+' kb 重叠。','muted'));
 if(f.reference_genome)row.append(node('p','记录参考：'+f.reference_genome,'muted'));
 const caveat=/QTL/i.test(f.evidence_type)?'QTL区间关联不等于因果基因认定；历史遗传图谱坐标不能直接当作B73 v5物理坐标。':f.evidence_type==='GWAS_REPORTED_GENE_TRAIT_ASSOCIATION'?'作者或数据库报告的基因–性状关联，不等于已证明因果；基因ID统一为B73 v5不代表原报告位点已转换。':'位点邻近或标记关联不等于因果基因认定；同研究的模型和年份不算独立研究，P值不合并。';
 row.append(node('p',caveat,'muted'));
 if(/^\d{4,10}$/.test(String(d.PMID||'')))row.append(link('https://pubmed.ncbi.nlm.nih.gov/'+d.PMID+'/','关联研究（PubMed） ↗','text-link'),document.createTextNode(' · '));
}
function evidenceRow(f){const d=f.details||{},row=node('article',undefined,'evidence-row');const trait=isTraitEvidence(f);if(trait)row.dataset.evidenceKind='trait';row.append(node('h4',evidenceTypes[f.evidence_type]||(trait&&/QTL/i.test(f.evidence_type)?'QTL · 区间关联':f.evidence_type)));if(trait)appendTraitContext(row,f);else{const direction=[f.subject_label||f.subject_node_id,f.object_label||f.object_node_id].filter(Boolean).join(' → ');row.append(node('p',direction));const fold=d.log2FC??d.log2FoldChange;if(fold!==undefined)row.append(node('p',(d.histone_mark||d.chromatin_mark||'分子响应')+' · '+(d.comparison||'')+' · log₂FC '+Number(fold).toFixed(2)+(d.FDR!==undefined?' · FDR '+d.FDR:'')));else if(d.phenotype||d.description||d.annotation)row.append(node('p',String(d.phenotype||d.description||d.annotation).slice(0,340)));}row.append(node('p',[f.genotype,f.development_stage,f.tissue_or_compartment].filter(Boolean).join(' · '),'muted'));if(f.category==='molecular_response')row.append(node('p','注意：这条记录的扰动材料不一定是当前基因的突变体；响应不等于直接调控。','muted'));if(f.source_url)row.append(link(f.source_url,trait&&f.source_url.includes('ngdc.cncb.ac.cn')?'GWAS Atlas 来源 ↗':'原始来源 ↗','text-link'));const detail=node('details');detail.append(node('summary','查看记录与支持范围'),node('pre',JSON.stringify(f,null,2),'record-json'));row.append(detail);return row;}
function selectEvidenceFocus(button){evidenceFilter=button.dataset.focus;document.querySelectorAll('[data-focus]').forEach(x=>x.classList.toggle('selected',x===button));loadEvidence(true).catch(e=>$('evidence-status').textContent=errorText(e.message));}
function syncTraitFilter(total){let button=$('evidence-traits-filter');if(!total){button?.remove();return;}if(!button){button=node('button','GWAS / QTL');button.id='evidence-traits-filter';button.type='button';button.dataset.focus='traits';button.addEventListener('click',()=>selectEvidenceFocus(button));$('evidence-filters').append(button);}button.classList.toggle('selected',evidenceFilter==='traits');}
async function loadEvidence(reset=false,stamp=serial){if(!current)return;const gene=current.gene_id,focus=evidenceFilter;if(reset){evidenceOffset=0;wipe($('all-evidence'));}const data=await request('/api/evidence/'+gene+'?offset='+(evidenceOffset||0)+'&page_size=12&focus='+focus);if(stamp!==serial||focus!==evidenceFilter)return;syncTraitFilter(data.trait_total||0);data.evidence.forEach(f=>$('all-evidence').append(evidenceRow(f)));evidenceOffset=data.next_offset;$('evidence-more').hidden=evidenceOffset===null;$('evidence-status').textContent='这个筛选共 '+data.total+' 条记录，可继续加载或完整下载。';if(!data.evidence.length&&focus!=='traits')$('all-evidence').append(node('p','该类别暂未取回记录；不据此推断该基因没有功能。','muted'));}

$('search-form').addEventListener('submit',e=>{e.preventDefault();searchGene($('gene-input').value.trim());});document.querySelectorAll('[data-gene]').forEach(b=>b.addEventListener('click',()=>{$('gene-input').value=b.dataset.gene;searchGene(b.dataset.gene);}));
$('agent-form').addEventListener('submit',e=>{e.preventDefault();askAgent();});$('question').addEventListener('keydown',e=>{if(e.key==='Enter'&&(e.ctrlKey||e.metaKey)&&!$('ask-button').disabled){e.preventDefault();askAgent();}});$('new-chat').addEventListener('click',()=>resetChat());bindQuestions($('conversation'));
document.querySelectorAll('[data-tab]').forEach(b=>b.addEventListener('click',()=>setTab(b.dataset.tab)));document.querySelectorAll('[data-go-tab]').forEach(b=>b.addEventListener('click',()=>setTab(b.dataset.goTab)));
$('close-source').addEventListener('click',()=>$('source-dialog').close());$('source-dialog').addEventListener('click',e=>{if(e.target===$('source-dialog'))$('source-dialog').close();});
$('network-zoom-out').addEventListener('click',()=>explorer.zoom(1/1.2));$('network-zoom-in').addEventListener('click',()=>explorer.zoom(1.2));$('network-reset').addEventListener('click',()=>explorer.fit());$('network-direction').addEventListener('change',()=>loadNetwork(true).catch(e=>$('network-status').textContent=errorText(e.message)));$('network-more').addEventListener('click',()=>loadNetwork(false).catch(e=>$('network-status').textContent=errorText(e.message)));
$('literature-more').addEventListener('click',()=>loadLiterature(false).catch(e=>$('literature-status').textContent=errorText(e.message)));$('literature-continue').addEventListener('click',continueLiterature);document.querySelectorAll('[data-filter]').forEach(b=>b.addEventListener('click',()=>{literatureFilter=b.dataset.filter;document.querySelectorAll('[data-filter]').forEach(x=>x.classList.toggle('selected',x===b));loadLiterature(true).catch(e=>$('literature-status').textContent=errorText(e.message));}));
$('evidence-more').addEventListener('click',()=>loadEvidence(false).catch(e=>$('evidence-status').textContent=errorText(e.message)));document.querySelectorAll('[data-focus]').forEach(b=>b.addEventListener('click',()=>selectEvidenceFocus(b)));
request('/api/info').then(info=>{const n=$('release-detail');n.append(node('p','模型：'+info.model_release+'；图谱：'+info.kg_release+'；坐标：'+info.coordinate_system),node('p',info.privacy),node('p','预测表 SHA256：'+info.integrity.model_sha256),node('p','图谱 SHA256：'+info.integrity.kg_sha256));}).catch(()=>{});
const queryGene=new URLSearchParams(location.search).get('gene');if(queryGene&&/^[A-Za-z0-9_.-]{2,64}$/.test(queryGene)){$('gene-input').value=queryGene;searchGene(queryGene);}
const mobileToggle=node('button','展开数据与时间权重','mobile-model-toggle subtle-button');mobileToggle.type='button';mobileToggle.setAttribute('aria-expanded','false');const modelCard=$('coverage').parentElement;modelCard.insertBefore(mobileToggle,modelCard.querySelector('.sidebar-divider'));mobileToggle.addEventListener('click',()=>{const expanded=modelCard.classList.toggle('expanded');mobileToggle.setAttribute('aria-expanded',String(expanded));mobileToggle.textContent=expanded?'收起数据与时间权重':'展开数据与时间权重';});
