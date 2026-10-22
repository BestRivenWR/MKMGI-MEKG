"""Uncapped local evidence access, paginated transport and immutable sources."""
import json,sqlite3,hashlib
from functools import lru_cache
from pathlib import Path
from data_access import ROOT,CONFIG,connection,node_info
from public_retrieval import public_edge,valid_source,compact_fact
from trait_display import is_trait_evidence,display_context,display_order

ORDER="""CASE v.evidence_type WHEN 'MAIZE_GENETIC_FUNCTION' THEN 0
WHEN 'CURATED_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE' THEN 1
WHEN 'BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT' THEN 2
WHEN 'MUTANT_TRANSCRIPTOMIC_RESPONSE' THEN 3 WHEN 'MUTANT_CHROMATIN_CHANGE' THEN 4
WHEN 'BINDING_CHROMATIN_RNA_INTERSECTION' THEN 5
WHEN 'DATABASE_FUNCTION_DESCRIPTION' THEN 6 ELSE 20 END"""
SQL='''SELECT e.edge_id,e.subject_node_id,e.object_node_id,e.predicate,e.independence_class,
 n.label AS object_label,s.label AS subject_label,v.*
 FROM edges e JOIN evidence v ON v.edge_id=e.edge_id
 JOIN nodes n ON n.node_id=e.object_node_id JOIN nodes s ON s.node_id=e.subject_node_id
 WHERE (e.subject_node_id=? OR e.object_node_id=?) AND e.independence_class LIKE 'PUBLIC_%' '''

def overlay():
    c=sqlite3.connect((ROOT/'derived/evidence_index.sqlite').as_uri()+'?mode=ro&immutable=1',uri=True)
    c.row_factory=sqlite3.Row;return c

@lru_cache(maxsize=2048)
def aliases(gene):
    c=connection();extra=overlay()
    try:
        result={gene}
        result.update(r[0] for r in c.execute('SELECT alias_id FROM aliases WHERE canonical_node_id=?',('gene:'+gene,)))
        result.update(r[0] for r in c.execute("SELECT alias FROM public_aliases WHERE gene_id=? AND alias_type IN ('gene_symbol','gene_synonym','ensembl_external_gene_name')",(gene,)))
        result.update(r[0] for r in extra.execute('SELECT alias FROM search_aliases WHERE gene_id=?',(gene,)))
        return sorted(x for x in result if x and len(x)<=150)
    finally:c.close();extra.close()

def curated_links(gene):
    c=overlay()
    try:return [dict(r) for r in c.execute('SELECT * FROM curated_links WHERE gene_id=? ORDER BY doi,source_url',(gene,))]
    finally:c.close()

@lru_cache(maxsize=1)
def symbol_map():
    c=connection();result={}
    try:
        for r in c.execute("SELECT gene_id,alias FROM public_aliases WHERE alias_type='gene_symbol' ORDER BY gene_id,alias"):
            result.setdefault(r[0],r[1])
        return result
    finally:c.close()

@lru_cache(maxsize=12)
def all_public_evidence(gene):
    """All matching records; no LIMIT, source/type cap, or silent truncation."""
    c=connection();items=[];excluded=0;symbols=symbol_map()
    try:
        for r in c.execute(SQL+' ORDER BY '+ORDER+',v.evidence_id',('gene:'+gene,'gene:'+gene)):
            x=dict(r);x['attributes']=json.loads(x.get('attributes_json') or '{}')
            f=public_edge(x)
            if not f:excluded+=1;continue
            f['edge_id']=x['edge_id'];f['origin']='public_kg';f['focus_role']='subject' if x['subject_node_id']=='gene:'+gene else 'object'
            for k in ['subject','object']:
                n=f.get(k+'_node_id','') or ''
                if n.startswith('gene:'):f[k+'_label']=symbols.get(n[5:],f.get(k+'_label') or n[5:])
            # Preserve typed quantitative/context fields missing from legacy allowlist.
            for k in ['label','detection_fraction','nuclei','detected_nuclei','chromatin_mark','log2FoldChange','FDR','comparison','assay','tissue','stage']:
                val=x['attributes'].get(k)
                if isinstance(val,(str,int,float,bool)) and len(str(val))<800:f['details'][k]=val
            f['citation']='KG:'+f['evidence_id'];items.append(f)
        # Author-final single-nucleus groups: all groups/replicates, not selected examples.
        for r in c.execute('''SELECT e.*,s.group_level,s.replicate,s.label,s.nuclei,s.UMI_total
          FROM v31_cell_expression e JOIN v31_cell_groups s ON e.group_id=s.group_id
          WHERE e.gene_id=? AND s.group_level='cell_type' ORDER BY s.label,s.replicate''',(gene,)):
            d=dict(r);d['detection_fraction']=d['detected_nuclei']/d['nuclei'] if d['nuclei'] else 0
            f=compact_fact('cell:'+gene+':'+d['group_id'],'SINGLE_NUCLEUS_EXPRESSION_CONTEXT',
              'https://doi.org/10.1093/plcell/koaf162',d,gene,genotype='B73',development_stage='19 DAP',tissue_or_compartment=d['label'])
            f.update(origin='public_kg',focus_role='subject',subject_label=node_info(gene)['symbol'],citation='KG:'+f['evidence_id'])
            items.append(f)
    finally:c.close()
    return {'gene_id':gene,'evidence':items,'total':len(items),'invalid_source_records':excluded}

def evidence_page(gene,offset=0,page_size=50,focus='all'):
    src=all_public_evidence(gene);facts=src['evidence']
    trait_total=sum(is_trait_evidence(f) for f in facts)
    if focus=='genetics':facts=[f for f in facts if f['category']=='genetic_or_functional']
    elif focus=='molecular':facts=[f for f in facts if f['category']=='molecular_response']
    elif focus=='regulation':facts=[f for f in facts if any(t in f['evidence_type'] for t in ['BINDING','TF_','CUTTAG'])]
    elif focus=='function':facts=[f for f in facts if any(t in f['evidence_type'] for t in ['ANNOTATION','FUNCTION_DESCRIPTION','ORTHOLOG'])]
    elif focus=='traits':facts=[f for f in facts if is_trait_evidence(f)]
    if focus=='all':facts=display_order(facts)
    page=[dict(f) for f in facts[offset:offset+page_size]]
    # Browser-only fields; leave the cached cloud/agent packet unchanged.
    trait_page=[f for f in page if is_trait_evidence(f)]
    if trait_page:
        c=connection()
        try:
            for f in trait_page:
                row=c.execute('SELECT attributes_json FROM evidence WHERE evidence_id=?',(f['evidence_id'],)).fetchone()
                if row:f['trait_context']=display_context(f,json.loads(row[0] or '{}'))
        finally:c.close()
    return {'gene_id':gene,'evidence':page,'total':len(facts),'offset':offset,
       'trait_total':trait_total,
       'next_offset':offset+len(page) if offset+len(page)<len(facts) else None,'page_size':page_size,
       'retrieval_scope':'All local public evidence accessible by pagination; no total-record cap.',
       'invalid_source_records':src['invalid_source_records']}

def public_paths_all(gene):
    c=connection();paths=[]
    try:
        rows=c.execute('''SELECT r.subject_node_id,r.object_node_id,r.predicate,r.edge_id AS relation_edge_id,
          p.edge_id AS phenotype_edge_id,p.object_node_id AS phenotype_node,n.label AS phenotype_label
          FROM edges r JOIN edges p ON p.subject_node_id=r.subject_node_id JOIN nodes n ON n.node_id=p.object_node_id
          WHERE r.object_node_id=? AND r.independence_class LIKE 'PUBLIC_%' AND p.independence_class LIKE 'PUBLIC_%'
          AND r.predicate IN ('REGULATES_WITH_EXPERIMENTAL_SUPPORT','PUBLIC_INTEGRATED_TF_TARGET')
          AND p.predicate IN ('HAS_GENETIC_PERTURBATION_PHENOTYPE','HAS_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE_EVIDENCE')
          ORDER BY r.edge_id,p.edge_id''',('gene:'+gene,))
        for r in rows:paths.append({**dict(r),'target_phenotype_proven':False,'interpretation':'Regulator phenotype is not target phenotype.'})
        return paths
    finally:c.close()

def exact_edge_facts(gene,edge=None,target=None):
    c=connection();where='';args=['gene:'+gene,'gene:'+gene]
    if edge:where+=' AND e.edge_id=?';args.append(edge)
    if target:where+=' AND (e.subject_node_id=? OR e.object_node_id=?)';args+=['gene:'+target,'gene:'+target]
    try:
        facts=[]
        for r in c.execute(SQL+where+' ORDER BY '+ORDER+',v.evidence_id',args):
            x=dict(r);x['attributes']=json.loads(x.get('attributes_json') or '{}');f=public_edge(x)
            if f:f.update(edge_id=x['edge_id'],focus_role='subject' if x['subject_node_id']=='gene:'+gene else 'object');facts.append(f)
        return facts
    finally:c.close()
