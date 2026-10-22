"""Read-only agent tools; current-gene cloud projection is deliberately narrow."""
import copy, json, re
from dataclasses import dataclass, field
from data_access import CONFIG, GENE_RE, node_info, resolve_gene, model_rows, model_profile, public_context
from public_repository import all_public_evidence
import literature as lit
import research

@dataclass
class Session:
    gene: str
    question: str
    ledger: dict = field(default_factory=dict)
    keys: dict = field(default_factory=dict)
    executed: dict = field(default_factory=dict)
    outcomes: list = field(default_factory=list)
    tool_calls: int = 0

    def add(self, records):
        result=[]
        for f in records:
            f=copy.deepcopy(f)
            key=(f.get('evidence_id'),f.get('source_url'),f.get('subject_node_id'),f.get('object_node_id'))
            cid=self.keys.get(key)
            if not cid:cid='E'+str(len(self.ledger)+1);self.keys[key]=cid
            f['citation']=cid;self.ledger[cid]=f;result.append(f)
        return result

def current_gene_summary(focal_gene, requested_gene):
    """User consent is ONLY for these fields of the current queried gene.

    Not a general model-export switch. No arrays, seed predictions, relation
    weights, other genes, files, raw observations, or unpublished experiments.
    """
    if requested_gene!=focal_gene:raise ValueError('only_current_gene_summary_authorized')
    if CONFIG.get('current_gene_summary_consent')!='score_rank_coverage_time_weights_only':
        return {'gene_id':focal_gene,'model_summary_shared':False}
    row=model_rows().get(focal_gene)
    if not row:return {'gene_id':focal_gene,'supported_by_model':False}
    profile=model_profile(focal_gene) or {}
    return {'gene_id':focal_gene,
       'priority_score':row['priority_score'],'priority_percentile':row['priority_percentile'],
       'rank':row['rank'],'model_gene_count':row['model_gene_count'],
       'coverage':{k:row['coverage'].get(k) for k in ['temporal_modalities_observed','spatial_observed','sequence_observed','protein_observed']},
       'time_weight_summary':{'DAP':profile.get('DAP',[]),
          'pooling_weights':profile.get('temporal_pooling_weights',[]),
          'observed':profile.get('temporal_observed',[])},
       'semantics':'Candidate priority, NOT phenotype probability; time weights are NOT causal attribution. No new independent validation.'}

def model_fact(gene):
    return {'evidence_id':'frozen-model-summary:'+gene,'evidence_type':'MODEL_CURRENT_GENE_SUMMARY',
       'subject_node_id':'gene:'+gene,'subject_label':node_info(gene)['symbol'],
       'predicate':'HAS_FROZEN_PRIORITY_SUMMARY','origin':'current_gene_model_summary',
       'details':current_gene_summary(gene,gene),'category':'model_context',
       'caveat':'冻结模型输出，不是遗传因果证据；时间注意力只是模型读取权重。'}

def cloud_fact(f, focal_gene):
    out={k:f[k] for k in ['citation','evidence_id','edge_id','evidence_type','subject_node_id','object_node_id','subject_label','object_label',
       'predicate','source_url','genotype','development_stage','tissue_or_compartment','category','focus_role','origin','caveat'] if f.get(k) is not None}
    d=f.get('details',{})
    if f.get('origin')=='current_gene_model_summary':
        # Projecting again at the serialization boundary defeats accidental
        # insertion of arbitrary model data into the ledger.
        gene=f['subject_node_id'].removeprefix('gene:')
        out['details']=current_gene_summary(focal_gene,gene)
    elif f.get('origin')=='online_literature':
        out['details']={k:d.get(k) for k in ['paper_key','title','doi','year','journal','reading_level','identity_match','species_check',
             'is_review','is_commentary','is_secondary','related_primary_dois','is_preprint','retracted','experiment_clues','fulltext_error']}
        out['details']['abstract']=d.get('abstract','')[:4200]
        out['details']['selected_passages']=[p[:1800] for p in (d.get('selected_passages') or [])[:4]]
        out['details']['available_matching_passages']=d.get('available_matching_passages',0)
    elif f.get('evidence_type')=='PRIMARY_LITERATURE_RELATION':
        out['details']=d
    elif f.get('evidence_type')=='MODEL_INPUT_RELATION':
        out['details']={k:d.get(k) for k in ['source','target','relation','used_by_model','support_scope']}
    else:
        keep={'phenotype','description','annotation','directness','assay','method','trait','direction','is_direct',
           'not_direct_regulation','not_phenotype_evidence','target_phenotype_proven','ontology','evidence_code',
           'tissue','stage','species','histone_mark','chromatin_mark','comparison','log2FC','log2FoldChange','FDR',
           'regulator','target','label','mapping_status','tf_gene_id','target_gene_id','gene_role','genotype',
           'experimental_support','native_tf','native_target','functional_reporter','limitations'}
        out['details']={k:v for k,v in d.items() if k in keep and isinstance(v,(str,int,float,bool,list,type(None)))}
    return out

def paper_fact(p):
    f=lit.as_fact(p)
    f['details'].update({k:p.get(k) for k in ['selected_passages','available_matching_passages','fulltext_error']})
    return f

def diverse_knowledge(gene,question,offset=0,size=12,focus='all'):
    facts=all_public_evidence(gene)['evidence']
    if focus=='phenotype':facts=[x for x in facts if x.get('category')=='genetic_or_functional']
    elif focus=='regulation':facts=[x for x in facts if any(t in x.get('evidence_type','') for t in ['BINDING','TF_','CUTTAG'])]
    elif focus=='response':facts=[x for x in facts if x.get('category')=='molecular_response']
    elif focus=='development':facts=[x for x in facts if x.get('development_stage') or x.get('category')=='spatial_context']
    levels={'MAIZE_GENETIC_FUNCTION':0,'CURATED_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE':0,'BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT':1,
       'MUTANT_TRANSCRIPTOMIC_RESPONSE':3,'MUTANT_CHROMATIN_CHANGE':3,'DATABASE_FUNCTION_DESCRIPTION':5}
    ordered=sorted(facts,key=lambda f:(levels.get(f['evidence_type'],6),f.get('evidence_id','')))
    # Round-robin by study+type+focal direction prevents one peak catalogue
    # from occupying all of a small transport page. No facts are discarded.
    groups={}
    for f in ordered:
        key=(f.get('source_url'),f['evidence_type'],f.get('focus_role'))
        groups.setdefault(key,[]).append(f)
    stream=[]
    for i in range(max([len(v) for v in groups.values()] or [0])):
        for group in groups.values():
            if i<len(group):stream.append(group[i])
    return {'records':stream[offset:offset+size],'total':len(stream),
        'next_offset':offset+size if offset+size<len(stream) else None}

def function(name,description,props,required):
    return {'type':'function','function':{'name':name,'description':description,'parameters':{
        'type':'object','properties':props,'required':required,'additionalProperties':False}}}

IDENT={'type':'string','maxLength':64,'description':'Exact B73 v5 ID or recorded maize gene symbol/old ID.'}
TOOLS=[
 function('gene_lookup','Resolve a gene symbol/old ID; ambiguity is not silently guessed.',{'identifier':IDENT},['identifier']),
 function('research_search','Search CURRENT PUBLIC INTERNET literature (EuropePMC/PubMed, Crossref, OpenAlex) by the actual question. Supply focused English queries, not an arbitrary fixed template. All cached papers can be paged. Search hits do not establish claims.',
   {'queries':{'type':'array','maxItems':3,'items':{'type':'string','maxLength':350}},'question':{'type':'string','maxLength':800},
    'offset':{'type':'integer','minimum':0},'refresh':{'type':'boolean'}},[]),
 function('paper_read','Read an exact returned paper key. Get abstract and question-relevant available open full-text result passages; no paywall bypass.',
   {'paper_key':{'type':'string','maxLength':300}},['paper_key']),
 function('source_read','Resolve a specific DOI through public publication APIs and read available text. Useful for user-provided studies or known primary sources. Identity must still match the current gene.',
   {'doi':{'type':'string','maxLength':220}},['doi']),
 function('knowledge_search','Read source-diverse local MEKG records including both upstream and downstream relationships; returns pagination, not a total cap.',
   {'focus':{'type':'string','enum':['all','phenotype','regulation','response','development']},'offset':{'type':'integer','minimum':0}},[]),
 function('model_context','Get ONLY the CURRENT queried gene score, rank, coverage and time-weight summary authorized by the user. Not a phenotype probability.',{},[]),
 function('network_neighbors','Read actual graph neighbors, distinguishing original model-input candidate relationships, PPI and published experimental relations. Model-input edges are NOT learned effect predictions and NOT experimental confirmation.',
   {'layer':{'type':'string','enum':['all','classic_experimental','literature_reported','public_experimental','public_association','ppi','model_regulatory','project_inferred']},
    'offset':{'type':'integer','minimum':0}},[])
]
SCHEMAS={x['function']['name']:x['function']['parameters'] for x in TOOLS}

def execute(name,args,session):
    if name not in SCHEMAS or not isinstance(args,dict):raise ValueError('tool_not_allowed')
    schema=SCHEMAS[name]
    if set(args)-set(schema['properties']) or not set(schema['required'])<=set(args):raise ValueError('tool_arguments_invalid')
    session.tool_calls+=1
    if session.tool_calls>CONFIG['max_tool_calls']:raise ValueError('tool_budget_exhausted')
    sig=json.dumps([name,args],sort_keys=True)
    if sig in session.executed:return {'already_returned':True,'citations':session.executed[sig]}
    if name in ['paper_read','source_read']:
        key=args.get('paper_key') if name=='paper_read' else 'doi:'+lit.normal_doi(args.get('doi'))
        existing=next((f for f in session.ledger.values() if f.get('origin')=='online_literature'
            and f.get('details',{}).get('paper_key')==key and f['details'].get('selected_passages') is not None),None)
        if existing:
            session.executed[sig]=[existing['citation']]
            return {'already_read':True,'evidence':[cloud_fact(existing,session.gene)],
                'instruction':'This source was already read to the available level. If full text is unavailable, report the boundary; do not repeatedly request the same paper.'}
    offset=args.get('offset',0)
    if not isinstance(offset,int) or isinstance(offset,bool) or not 0<=offset<10000000:raise ValueError('offset_invalid')
    records=[];extra={}
    if name=='gene_lookup':return resolve_gene(args['identifier'])
    if name=='model_context':records=[model_fact(session.gene)]
    elif name=='research_search':
        question=args.get('question') or session.question
        if offset:
            data=research.ranked_papers(session.gene,question,offset,12)
        else:data=research.search(session.gene,question,args.get('queries'),args.get('refresh',False))
        records=[paper_fact(p) for p in data['papers']]
        extra={k:v for k,v in data.items() if k!='papers'}
    elif name=='paper_read':records=[paper_fact(research.read(session.gene,args['paper_key'],session.question))]
    elif name=='source_read':records=[paper_fact(research.lookup_doi(session.gene,args['doi']))]
    elif name=='knowledge_search':
        data=diverse_knowledge(session.gene,session.question,offset,12,args.get('focus','all'))
        records=data.pop('records');extra=data
        if offset==0:records=[research.evidence_fact(r,session.gene) for r in research.relation_facts(session.gene)]+records
    elif name=='network_neighbors':
        from network import network_page
        layer=args.get('layer','all')
        data=network_page(session.gene,offset,12,layer)
        # User authorized a model SUMMARY, not arbitrary unpublished project
        # relations. Public input PPI/TF topology is shareable; omit private
        # project_inferred records completely from cloud outputs.
        for e in data['edges']:
            if e['layer'] in ['project_inferred','model_regulatory']:continue
            if e['layer'] in ['classic_experimental','literature_reported']:
                rel=next((r for r in research.relation_facts(session.gene) if r['id']==e['id']),None)
                if rel:records.append(research.evidence_fact(rel,session.gene))
            elif e['origin']=='frozen_model_graph':
                records.append({'evidence_id':e['id'],'evidence_type':'MODEL_INPUT_RELATION','origin':'public_model_input_topology',
                    'subject_node_id':'gene:'+e['source'],'object_node_id':'gene:'+e['target'],
                    'subject_label':node_info(e['source'])['symbol'],'object_label':node_info(e['target'])['symbol'],
                    'predicate':e['predicate'],'category':'model_context',
                    'details':{'source':e['source'],'target':e['target'],'relation':e['predicate'],'used_by_model':True,
                      'support_scope':'input candidate network, not a learned causal/effect prediction; no weights or private observations shared'}})
            else:
                facts=all_public_evidence(session.gene)['evidence']
                records += [f for f in facts if f.get('edge_id')==e['id']][:2]
        extra={k:data[k] for k in ['total','next_offset','layer_totals']}
    selected=session.add(records)
    session.executed[sig]=[x['citation'] for x in selected]
    return {**extra,'evidence':[cloud_fact(f,session.gene) for f in selected]}

def preliminary_check(item,session):
    refs=item.get('citations',[])
    if any(r not in session.ledger for r in refs):return 'unknown_citation'
    kind=item.get('kind')
    if kind not in ['boundary','recommendation','hypothesis','interpretation'] and not refs:return 'factual_claim_needs_citation'
    facts=[session.ledger[r] for r in refs]
    if facts and all(f.get('details',{}).get('retracted') for f in facts):return 'retracted_only'
    if kind=='published_finding':
        readable=any(f['evidence_type']=='PRIMARY_LITERATURE_RELATION' or
          (f.get('origin')=='online_literature' and f['details'].get('reading_level') not in [None,'metadata']
           and f['details'].get('identity_match')!='indexed_search_only' and not f['details'].get('retracted')) for f in facts)
        if not readable:return 'publication_not_read_or_identity_unchecked'
    if kind=='model_context' and not any(f['evidence_type'] in ['MODEL_CURRENT_GENE_SUMMARY','MODEL_INPUT_RELATION'] for f in facts):return 'missing_model_source'
    if re.search(r'保证.*表型|敲除成功率|因果概率为|完全没有功能',item.get('text','')):return 'unsupported_probability_or_absence'
    # Direct regulation in a READ primary paper is NOT rejected just because
    # the old KG has no pre-existing typed edge. Entailment is checked below.
    return None

def validate_relation(proposal,session):
    refs=proposal.get('citations',[])
    if not refs or any(x not in session.ledger for x in refs):return 'relation_missing_citation'
    facts=[session.ledger[x] for x in refs]
    papers=[f for f in facts if f.get('origin')=='online_literature' and not f['details'].get('retracted') and not f['details'].get('is_review') and not f['details'].get('is_secondary')
       and f['details'].get('identity_match')!='indexed_search_only']
    if not papers:return 'relation_requires_primary_reading'
    q=research.clean_quote(proposal.get('quote'))
    if not 15<=len(q)<=500 or not any(q.casefold() in research.clean_quote(research.source_text(f)).casefold() for f in papers):
        return 'supporting_quote_not_found'
    try:a=research.resolve_entity(proposal['source_symbol']);b=research.resolve_entity(proposal['target_symbol'])
    except (ValueError,KeyError):return 'relation_entity_unresolved'
    if session.gene not in [a['id'],b['id']]:return 'relation_does_not_involve_current_gene'
    if proposal.get('support_level') not in ['binding_and_functional','promoter_functional','binding_only','perturbation_response','high_throughput_candidate','protein_interaction']:
        return 'relation_support_level_invalid'
    if proposal.get('relation') not in ['ACTIVATES','REPRESSES','BINDS_PROMOTER','PERTURBATION_RESPONSE','CANDIDATE_DIRECT_TARGET','PROTEIN_INTERACTION','COOPERATIVE_REGULATION']:
        return 'relation_type_invalid'
    return None
