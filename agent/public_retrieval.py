"""Public-only graph retrieval. All scientific inputs are read-only.

No private observations, label tiers, experimental results or model scores are
included in the cloud evidence packet. Every selectable fact retains its source.
"""
import collections,csv,gzip,hashlib,json,re,sqlite3,sys
from pathlib import Path
from functools import lru_cache
from urllib.parse import urlsplit

ROOT=Path(__file__).resolve().parent
from settings import load_config
CONFIG=load_config()
DB=Path(CONFIG['kg_database'])
sys.path.insert(0,str(ROOT.parent/'mekg'))
sys.path.insert(0,str(ROOT.parent/'mekg/_builders/v31'))
sys.dont_write_bytecode=True
from query_v31 import query as legacy_query

PUBLIC_HOSTS=('doi.org','ncbi.nlm.nih.gov','maizegdb.org','gramene.org','ensembl.org',
              'geneontology.org','kegg.jp','arabidopsis.org','uniprot.org','ebi.ac.uk',
              'nature.com','sciencedirect.com','cell.com','springer.com','academic.oup.com',
              'ngdc.cncb.ac.cn')
DETAIL_KEYS={
 'directness','histone_mark','assay','comparison','contrast','log2FC','log2FoldChange',
 'log2fc','FDR','fdr','padj','p_adj','pvalue','p_value','effect','direction','tissue',
 'genotype','stage','development_stage','condition','description','phenotype',
 'evidence_grade','evidence_code','ontology','gene_symbol','regulator','target',
 'tf_gene_id','target_gene_id','ortholog_gene_id','species','method','trait',
 'not_phenotype_evidence','not_direct_regulation','target_phenotype_proven',
 'binding_and_expression_overlap','reporter_support','mutant','mutant_genotype',
 'chromatin_mark','mean_peak_signal','promoter_window','ARGOT_PPV',
 'distance_to_tss','distance_bp','annotation','insertion_types','evidence_layer',
 'regulation_direction','replicate_count','baseMean','log2_fold_change',
 'fdr_adjusted_p','adjusted_p_value','is_direct','treatment','control'
}
STRONG={'MAIZE_GENETIC_FUNCTION','BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT',
        'CURATED_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE'}
MOLECULAR={'BINDING_CHROMATIN_RNA_INTERSECTION','MUTANT_TRANSCRIPTOMIC_RESPONSE',
           'MUTANT_CHROMATIN_CHANGE','RECIPROCAL_CROSS_ALLELE_EXPRESSION'}

def connection():
    c=sqlite3.connect(DB.as_uri()+'?mode=ro&immutable=1',uri=True)
    c.row_factory=sqlite3.Row
    return c

def valid_source(value):
    try:
        u=urlsplit(value or '')
        return u.scheme=='https' and not u.username and any(u.hostname==h or (u.hostname or '').endswith('.'+h) for h in PUBLIC_HOSTS)
    except ValueError:return False

def safe_details(obj):
    """Positive allowlist: arbitrary source attributes never cross the API boundary."""
    if not isinstance(obj,dict):return {}
    result={}
    for k,v in obj.items():
        if k not in DETAIL_KEYS:continue
        if isinstance(v,(int,float,bool)) or v is None:result[k]=v
        elif isinstance(v,str) and len(v)<600 and not re.search(r'(/home/|/tmp/|[A-Z]:\\|sk-[a-zA-Z0-9])',v):result[k]=v
        elif isinstance(v,list) and len(v)<10 and all(isinstance(x,(str,int,float,bool)) for x in v):
            result[k]=[str(x)[:120] for x in v if not str(x).startswith('/')]
    return result

def band(typ):
    if typ in STRONG:return 'genetic_or_functional'
    if typ in MOLECULAR:return 'molecular_response'
    return 'context'

def public_edge(x):
    if not valid_source(x.get('source_url')):return None
    y={k:x.get(k) for k in ['evidence_id','edge_id','subject_node_id','predicate','object_node_id',
       'object_label','subject_label','evidence_type','source_url','development_stage',
       'tissue_or_compartment','genotype','reference_genome','effect_direction','effect_size','experiment_group','gene_role']}
    y['details']=safe_details(x.get('attributes',{}))
    y['category']=band(y['evidence_type'])
    y['independence_status']='not_certified'
    return y

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

@lru_cache(maxsize=1)
def predictions():
    # This compatibility layer handles public evidence only. The new model table
    # is loaded and audited independently in data_access.py; never read KG scores.
    return {}

def diversity(items,per_group=3,maximum=30):
    buckets=collections.OrderedDict()
    seen=set()
    for item in items:
        key=(item['source_url'],item['evidence_type'],item.get('subject_node_id'),item.get('object_node_id'),json.dumps(item.get('details',{}),sort_keys=True))
        if key in seen:continue
        seen.add(key)
        group=(item['source_url'],item['evidence_type'])
        buckets.setdefault(group,[]).append(item)
    out=[]
    for depth in range(per_group):
        for values in buckets.values():
            if depth<len(values):out.append(values[depth])
            if len(out)>=maximum:return out
    return out

def public_paths(gene):
    """Bounded two-hop, with BOTH edges and sources required to be public."""
    c=connection()
    try:
        rows=c.execute('''SELECT r.subject_node_id,r.object_node_id,r.predicate,
          r.edge_id AS relation_edge_id,p.edge_id AS phenotype_edge_id,
          p.object_node_id AS phenotype_node,n.label AS phenotype_label
          FROM edges r JOIN edges p ON p.subject_node_id=r.subject_node_id
          JOIN nodes n ON n.node_id=p.object_node_id
          WHERE r.object_node_id=? AND r.independence_class LIKE 'PUBLIC_%'
          AND p.independence_class LIKE 'PUBLIC_%'
          AND r.predicate IN ('REGULATES_WITH_EXPERIMENTAL_SUPPORT','PUBLIC_INTEGRATED_TF_TARGET')
          AND p.predicate IN ('HAS_GENETIC_PERTURBATION_PHENOTYPE','HAS_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE_EVIDENCE')
          ORDER BY r.edge_id,p.edge_id LIMIT 6''',('gene:'+gene,)).fetchall()
        paths=[];facts=[]
        for row in rows:
            ids=[];local=[]
            for edge in ['relation_edge_id','phenotype_edge_id']:
                evidence=c.execute('''SELECT e.subject_node_id,e.object_node_id,e.predicate,
                  n.label AS object_label,s.label AS subject_label,v.*
                  FROM edges e JOIN evidence v ON v.edge_id=e.edge_id
                  JOIN nodes n ON n.node_id=e.object_node_id JOIN nodes s ON s.node_id=e.subject_node_id
                  WHERE e.edge_id=? ORDER BY v.evidence_id LIMIT 3''',(row[edge],)).fetchall()
                for ev in evidence:
                    x=dict(ev);x['attributes']=json.loads(x.get('attributes_json') or '{}')
                    fact=public_edge(x)
                    if fact:local.append(fact);ids.append(fact['evidence_id']);break
            if len(ids)==2:
                paths.append({**dict(row),'supporting_evidence_ids':ids,'target_phenotype_proven':False,
                              'interpretation':'Regulator phenotype does not establish target phenotype; context-dependent hypothesis.'})
                facts+=local
        return paths,facts
    finally:c.close()

def compact_fact(eid,typ,url,details,subject,**extra):
    return {'evidence_id':eid,'evidence_type':typ,'source_url':url,'details':details,
            'subject_node_id':'gene:'+subject,'category':'context','independence_status':'not_certified',**extra}

@lru_cache(maxsize=64)
def retrieve(gene):
    raw=legacy_query(gene,500,DB)
    if 'gene_id' not in raw:return raw
    gene=raw['gene_id']
    facts=[y for x in raw['public_evidence'] if (y:=public_edge(x))]
    paths,pathfacts=public_paths(gene)
    facts+=pathfacts
    cells=[x for x in raw.get('cell_type_and_cluster_expression',[]) if x['group_level']=='cell_type']
    cellfacts=[]
    if cells:
        # Author-final cell types supersede legacy pre-clustering aggregate summaries.
        facts=[x for x in facts if not (x['evidence_type']=='SINGLE_NUCLEUS_EXPRESSION_CONTEXT'
               and '10.1093/plcell/koaf162' in x['source_url'])]
        groups=collections.defaultdict(list)
        for x in cells:groups[x['label']].append(x)
        for label,rows in sorted(groups.items()):
            total=sum(x['nuclei'] for x in rows);detected=sum(x['detected_nuclei'] for x in rows)
            details={'label':label,'nuclei':total,'detected_nuclei':detected,
              'detection_fraction':detected/total if total else 0,
              'summary_method':'pooled detected nuclei / pooled nuclei across listed biological replicates; not replicate independence proof',
              'replicates':[{k:x.get(k) for k in ['replicate','UMI_sum','detected_nuclei','nuclei','pseudobulk_CPM','detection_fraction']} for x in rows]}
            cellfacts.append(compact_fact('celltype:'+gene+':'+label,'SINGLE_NUCLEUS_EXPRESSION_CONTEXT',rows[0]['source_url'],
              details,gene,genotype='B73',development_stage='19 DAP',
              tissue_or_compartment=label+'; placenta-chalazal is maternal, not endosperm',
              caveat='Author-final cell types; cell-type and cluster summaries are the same cells, not independent evidence.'))
        facts+=cellfacts
    # Source audits supply publication identity, never labels or private training metadata.
    source_urls={a['source_id']:'https://doi.org/'+a['doi'] for a in raw.get('public_source_audit',[]) if a.get('doi')}
    # Fixed identities for catalog additions, validated against release metadata.
    source_urls.update({'tf_peak_catalog_2025':'https://doi.org/10.1016/j.molp.2025.08.009',
                        'myb127_peaks_2026':'https://doi.org/10.1038/s41477-026-02238-3'})
    for x in raw.get('peak_support',{}).get('examples',[]):
        url=source_urls.get(x['source_id'])
        if not valid_source(url):continue
        facts.append(compact_fact('peak:'+x['peak_id']+':'+x['gene_role'],'TF_PEAK_GENE_ASSOCIATION',url,
          {k:x.get(k) for k in ['tf_gene_id','target_gene_id','native_chrom','native_start0','native_end0',
          'v5_chrom','v5_start0','v5_end0','mapping_status','annotation','distance_to_tss','score','gene_role']},gene,
          caveat='Peak association alone is not direct regulation; catalog assay must be checked per experiment. Coordinates 0-based half-open.'))
    tss=raw.get('rampage_TSS_context',{})
    if tss.get('consensus_pair_count'):
        coords=sorted({(x['chromosome'],x['start0'],x['end0']) for x in tss.get('examples',[])})[:5]
        facts.append(compact_fact('rampage:'+gene,'RAMPAGE_TSS_CONTEXT',tss['source_url'],
          {'replicate_supported_pairs':tss['consensus_pair_count'],'unique_example_coordinates':coords,
           'assay':tss['assay'],'canonical_TSS_replaced':False,'strand':'not_available','coordinates':'0-based half-open'},gene))
    for i,x in enumerate(raw.get('kernel_trait_GWAS_context',[])[:3]):
        facts.append(compact_fact('gwas:'+gene+':'+str(i),'GWAS_PROXIMITY_CONTEXT',x['source_url'],safe_details(x['attributes']),gene,
                     object_node_id=x['object_node_id'],caveat='Within 50 kb; not a causal gene assignment; correlated P values not combined.'))
    contacts=raw.get('physical_contact_context',{}).get('author_gene_pairs',[])
    for i,x in enumerate(contacts[:3]):
        facts.append(compact_fact('hichip:'+gene+':'+str(i),'CROSS_TISSUE_HICHIP_CONTEXT','https://doi.org/10.1038/s41467-025-65065-8',
           {'subject':x['subject_node_id'],'object':x['object_node_id'],'tissue':'seedlings, NOT endosperm'},gene,
           caveat='Cross-tissue physical-contact context, not evidence of endosperm regulation.'))
    selected=diversity(facts,3,32)
    # Two-hop paths must reference facts included in the same cloud packet.
    selected_ids={x['evidence_id'] for x in selected}
    for fact in pathfacts+cellfacts:
        if fact['evidence_id'] not in selected_ids and len(selected)<44:
            selected.append(fact);selected_ids.add(fact['evidence_id'])
    paths=[p for p in paths if set(p['supporting_evidence_ids'])<=selected_ids]
    for i,x in enumerate(selected,1):x['citation']='E'+str(i)
    mapping={x['evidence_id']:x['citation'] for x in selected}
    for p in paths:p['citations']=[mapping[e] for e in p['supporting_evidence_ids']]
    return {'gene_id':gene,'requested_gene_id':raw['requested_gene_id'],'coordinate_system':'B73 v5',
      'prediction':None,
      'kg_release':'MEKG public v3.1','kg_sha256':CONFIG['kg_sha256'],'evidence':selected,'paths':paths,
      'public_summary':raw['public_summary'],'retrieval':{'selected_facts':len(selected),
      'legacy_query_record_limit':500,'per_source_type_cap':3,'selection':'typed evidence with source/type diversity; not exhaustive'},
      'boundaries':['优先级是模型排序分数，不是胚乳发育的因果概率。','公共支持不自动等于独立验证；跨材料、时期和组织的证据不能直接合并。',
       '峰邻近、同源注释和物理接触是不同证据层级；没有记录不等于没有功能。','AI 解读只使用下列公共证据，不改变模型分数，也不上传模型结果或未发表实验数据。']}

def cloud_packet(result,intent):
    return {'gene_id':result['gene_id'],'intent':intent,'coordinate_system':'B73 v5',
      'public_evidence':result['evidence'],'public_graph_paths':result['paths'],
      'cautions':result['boundaries'][:3],'retrieval_not_exhaustive':True}
