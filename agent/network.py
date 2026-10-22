"""Real model/KG neighbors, provenance-first. Project records never go to LLM."""
import collections,json
from functools import lru_cache
from data_access import connection,node_info
from public_repository import overlay,symbol_map

LEGEND={
 'classic_experimental':{'label':'经典原始研究 · 核对关系','color':'#b58228','dash':False,'caution':'从原始研究核对的具体关系，保留实验、材料和适用范围；不等于本项目独立验证。'},
 'literature_reported':{'label':'原文报告 · 候选 / 响应','color':'#b99564','dash':True,'caution':'含核对摘要的候选/响应关系及agent自动读取的作者报告；保留支持等级，不能全部视为逐靶标功能验证。'},
 'public_experimental':{'label':'公共实验支持','color':'#95752f','dash':False,'caution':'证据类型和实验条件仍需逐条核对；不自动认定当前基因表型。'},
 'public_association':{'label':'公共关联 / 推断','color':'#8775a2','dash':True,'caution':'公共网络、结合峰或物理接触；不能单独认定直接调控。'},
 'ppi':{'label':'PPI 网络','color':'#678797','dash':False,'caution':'继承数据库网络分数；尚未细分experiments-only通道，不能全部视为实测物理互作。'},
 'model_regulatory':{'label':'模型输入 · 候选调控','color':'#519c83','dash':True,'caution':'真实模型使用的输入候选边，不是模型新预测的效应或实验验证；confidence=1是存在标记，不是验证概率。'},
 'project_inferred':{'label':'本项目组学推断','color':'#1d7659','dash':True,'caution':'本项目关系推断 / footprint或峰关联，未标作遗传验证。'},
}
PUBLIC_EXPERIMENTAL={'REGULATES_WITH_EXPERIMENTAL_SUPPORT','PUBLIC_TF_MUTANT_ALTERS_TARGET_EXPRESSION'}

def _gene_node(g):
    if g.startswith('symbol:'):
        return {'id':g,'label':g[7:],'node_type':'publication_entity','description':'原文基因symbol或基因集合；B73 v5映射待确认，不编造基因号。','coordinate_system':None}
    info=node_info(g)
    return {'id':g,'label':info['symbol'],'description':info.get('description',''),'node_type':'gene','coordinate_system':'B73 v5'}

def _gene_nodes(genes):
    """Read network metadata in PK-indexed batches, not N alias-table scans.

    Full gene metadata is still available from /api/gene. This projection only
    supplies graph labels/descriptions; no edge or entity is removed.
    """
    ordered=sorted(genes);native=[g for g in ordered if not g.startswith('symbol:')]
    rows={};symbols=symbol_map();c=connection()
    try:
        for offset in range(0,len(native),400):
            chunk=native[offset:offset+400]
            sql='SELECT node_id,label,attributes_json FROM nodes WHERE node_id IN ('+','.join('?' for _ in chunk)+')'
            for row in c.execute(sql,['gene:'+g for g in chunk]):rows[row['node_id'][5:]]=dict(row)
    finally:c.close()
    result=[]
    for g in ordered:
        if g.startswith('symbol:'):result.append(_gene_node(g));continue
        row=rows.get(g,{})
        attrs=json.loads(row.get('attributes_json') or '{}')
        result.append({'id':g,'label':symbols.get(g) or row.get('label') or g,
          'description':str(attrs.get('description') or attrs.get('function') or attrs.get('annotation') or '')[:600],
          'node_type':'gene','coordinate_system':'B73 v5'})
    return result

@lru_cache(maxsize=64)
def all_network(gene):
    edges=[];c=overlay();kg=connection()
    try:
        for r in c.execute('SELECT * FROM model_edges WHERE source=? OR target=? ORDER BY relation,weight DESC,edge_id',(gene,gene)):
            d=dict(r);layer='ppi' if d['relation']=='PPI' else 'model_regulatory'
            edges.append({'id':d['edge_id'],'source':d['source'],'target':d['target'],'predicate':d['relation'],
              'layer':layer,'directed':layer!='ppi','weight':d['weight'],'confidence':d['confidence'],
              'origin':'frozen_model_graph','used_by_current_model':True,'public_experimental':False,
              'support':LEGEND[layer]['caution'],'source_ref':d['source_ref'],'evidence_ids':[],
              'source_urls':[],'evidence_types':[],'genotype':None,'stage':None,'tissue':None})
        # One joined scan instead of one evidence query per edge (o2 >10k edges).
        joined={}
        for r in kg.execute('''SELECT e.*,v.evidence_id,v.evidence_type,v.source_url,v.genotype,v.development_stage,v.tissue_or_compartment
          FROM edges e LEFT JOIN evidence v ON v.edge_id=e.edge_id
          WHERE (subject_node_id=? OR object_node_id=?)
          AND subject_node_id LIKE 'gene:Zm00001eb%' AND object_node_id LIKE 'gene:Zm00001eb%'
          AND predicate IN ('REGULATES','BINDS_PROMOTER_OF','REGULATES_WITH_EXPERIMENTAL_SUPPORT',
           'PUBLIC_TF_MUTANT_ALTERS_TARGET_EXPRESSION','PUBLIC_INTEGRATED_TF_TARGET','PUBLIC_ENDOSPERM_REGULON_TARGETS',
           'PUBLIC_HICHIP_ASSOCIATED_GENE_PAIR','HAS_BINDING_PEAK_ASSOCIATED_GENE',
           'HAS_REPLICATED_CUTTAG_PEAK_ASSOCIATED_GENE','CANDIDATE_REPRESSIVE_TARGET')
          ORDER BY CASE predicate WHEN 'REGULATES_WITH_EXPERIMENTAL_SUPPORT' THEN 0
             WHEN 'BINDS_PROMOTER_OF' THEN 1 WHEN 'REGULATES' THEN 2 ELSE 3 END,e.edge_id,v.evidence_id''',('gene:'+gene,'gene:'+gene)):
            x=dict(r);eid=x['edge_id']
            if eid not in joined:joined[eid]={'edge':x,'evidence':[]}
            if x.get('evidence_id'):joined[eid]['evidence'].append({k:x.get(k) for k in ['evidence_id','evidence_type','source_url','genotype','development_stage','tissue_or_compartment']})
        for item in joined.values():
            d=item['edge'];ev=item['evidence'];pub=d['independence_class'].startswith('PUBLIC_')
            layer='public_experimental' if pub and d['predicate'] in PUBLIC_EXPERIMENTAL else 'public_association' if pub else 'project_inferred'
            urls=sorted({x['source_url'] for x in ev if (x.get('source_url') or '').startswith('https://')})
            attrs=json.loads(d.get('attributes_json') or '{}')
            edges.append({'id':d['edge_id'],'source':d['subject_node_id'][5:],'target':d['object_node_id'][5:],
              'predicate':d['predicate'],'layer':layer,'directed':bool(d['directed']),'weight':d['confidence_score'],
              'origin':'public_kg' if pub else 'project_kg','used_by_current_model':None,
              'public_experimental':layer=='public_experimental','support':LEGEND[layer]['caution'],
              'evidence_ids':[x['evidence_id'] for x in ev],'source_urls':urls,'evidence_types':sorted({x['evidence_type'] for x in ev}),
              'genotype':sorted({x['genotype'] for x in ev if x.get('genotype')}),
              'stage':sorted({x['development_stage'] for x in ev if x.get('development_stage')}),
              'tissue':sorted({x['tissue_or_compartment'] for x in ev if x.get('tissue_or_compartment')}),
              'attributes':attrs,'evidence_records':[{k:v for k,v in x.items() if k!='attributes_json'} for x in ev]})
    finally:c.close();kg.close()
    priority={'public_experimental':0,'project_inferred':1,'ppi':2,'model_regulatory':3,'public_association':4}
    # Interleave layers so high-degree public regulons do not hide PPI/project data.
    buckets=collections.defaultdict(list)
    for e in edges:buckets[e['layer']].append(e)
    result=[]
    for depth in range(max([len(v) for v in buckets.values()] or [0])):
        for layer in sorted(buckets,key=lambda k:priority[k]):
            if depth<len(buckets[layer]):result.append(buckets[layer][depth])
    return result

def network_page(gene,offset=0,page_size=60,layers='all',direction='all'):
    from research import network_relations
    literature_edges=network_relations(gene)
    items=literature_edges+all_network(gene);totals=dict(collections.Counter(x['layer'] for x in items))
    if layers!='all':
        selected=set(layers.split(','));items=[x for x in items if x['layer'] in selected]
    if direction=='outgoing':items=[x for x in items if not x['directed'] or x['source']==gene]
    elif direction=='incoming':items=[x for x in items if not x['directed'] or x['target']==gene]
    chosen=items[offset:offset+page_size];genes={gene}
    for e in chosen:genes.update([e['source'],e['target']])
    return {'gene_id':gene,'nodes':_gene_nodes(genes),'edges':chosen,'total':len(items),
       'offset':offset,'next_offset':offset+len(chosen) if offset+len(chosen)<len(items) else None,
       'layer_totals':totals,'legend':LEGEND,'graph_kind':'Real incident edges from frozen model graph and MEKG; expandable neighbors, not an invented schematic.',
       'private_project_records_stay_local':True,'record_cap':None,
       'cautions':['计算反向TF边不另画为抑制；图中箭头遵循原始关系方向。','公共实验支持仍需检查原文，不能把上游基因的表型转移给靶基因。',
       '本项目推断与已完成的功能实验不是同一证据层级；本轮未添加未经核实的“本实验验证”关系。']}

def full_network(gene):
    from research import network_relations
    return network_page(gene,0,len(all_network(gene))+len(network_relations(gene)) or 1)
