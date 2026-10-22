"""Read-only query adapter: frozen model lookup and provenance-aware public context."""
from pathlib import Path
import sys,json,sqlite3,argparse
from mekg_paths import workspace
ROOT=workspace()/'v31'
from query_base_v3 import query as parent_query

def query(gene,limit=30,db=None):
    path=Path(db) if db else ROOT/'release/mekg.sqlite'
    d=parent_query(gene,limit,path)
    if 'gene_id' not in d:return d
    c=sqlite3.connect(path.resolve().as_uri()+'?mode=ro&immutable=1',uri=True);c.row_factory=sqlite3.Row
    g=d['gene_id'];node='gene:'+g;lim=max(1,min(limit,500))
    d['kg_release']='MEKG-public-v3.1-20260907'
    d['evidence_database']=str(path)
    d['supported_by_model']=d.get('prediction') is not None
    d['prediction_mode']='frozen_consensus_lookup'
    d['cell_type_and_cluster_expression']=[]
    for r in c.execute('''SELECT e.*,s.group_level,s.replicate,s.label,s.nuclei,s.UMI_total,s.attributes_json
      FROM v31_cell_expression e JOIN v31_cell_groups s ON e.group_id=s.group_id
      WHERE e.gene_id=? ORDER BY s.group_level,s.label,s.replicate''',(g,)):
        x=dict(r);x['context']=json.loads(x.pop('attributes_json'))
        x['pseudobulk_CPM']=x['UMI_sum']/x['UMI_total']*1e6;x['detection_fraction']=x['detected_nuclei']/x['nuclei']
        x['source_url']='https://doi.org/10.1093/plcell/koaf162';d['cell_type_and_cluster_expression'].append(x)
    d['author_cluster_markers']=[]
    for r in c.execute('SELECT * FROM v31_cluster_markers WHERE gene_id=? ORDER BY avg_log2FC DESC',(g,)):
        x=dict(r);x['attributes']=json.loads(x.pop('attributes_json'));d['author_cluster_markers'].append(x)
    d['rampage_TSS_context']={'assay':'B73 15 DAP endosperm, RAMPAGE, 2 replicates',
      'source_url':'https://doi.org/10.1038/s41467-024-55195-w','canonical_TSS_replaced':False,'examples':[]}
    d['rampage_TSS_context']['consensus_pair_count']=c.execute('SELECT COUNT(*) FROM v31_tss_genes WHERE gene_id=?',(g,)).fetchone()[0]
    for r in c.execute('''SELECT a.association,t.*,p.summit0 AS rep1_summit0,q.summit0 AS rep2_summit0
      FROM v31_tss_genes a JOIN v31_tss_consensus t ON a.consensus_id=t.consensus_id
      JOIN v31_tss_peaks p ON p.peak_id=t.peak1_id JOIN v31_tss_peaks q ON q.peak_id=t.peak2_id
      WHERE a.gene_id=? ORDER BY t.chromosome,t.start0 LIMIT ?''',(g,lim)):
        d['rampage_TSS_context']['examples'].append(dict(r))
    d['peak_support']={'counts':[],'examples':[],'coordinate_convention':'0-based half-open; native B73 v4 and optional projected B73 v5 are separate'}
    for role,col in [('regulator','tf_gene_id'),('associated_gene','target_gene_id')]:
        for r in c.execute(f'''SELECT source_id,COUNT(*) AS peak_gene_records,SUM(v5_chrom IS NOT NULL) AS projected_v5_records
          FROM v31_peaks WHERE {col}=? GROUP BY source_id''',(g,)):
            d['peak_support']['counts'].append({'gene_role':role,**dict(r)})
        for r in c.execute(f'SELECT * FROM v31_peaks WHERE {col}=? ORDER BY source_id,peak_id LIMIT ?',(g,lim)):
            x=dict(r);x['gene_role']=role;d['peak_support']['examples'].append(x)
    r=c.execute('SELECT chromosome,start,end FROM nodes WHERE node_id=?',(node,)).fetchone()
    d['physical_contact_context']={'tissue':'seedlings, not endosperm','overlapping_projected_anchors':[]}
    d['physical_contact_context']['author_gene_pairs']=[]
    for row in c.execute('''SELECT subject_node_id,object_node_id,attributes_json FROM edges
      WHERE predicate='PUBLIC_HICHIP_ASSOCIATED_GENE_PAIR' AND (subject_node_id=? OR object_node_id=?)
      ORDER BY edge_id LIMIT ?''',(node,node,lim)):
        x=dict(row);x['attributes']=json.loads(x.pop('attributes_json'));d['physical_contact_context']['author_gene_pairs'].append(x)
    d['kernel_trait_GWAS_context']=[]
    for row in c.execute('''SELECT object_node_id,attributes_json FROM edges
      WHERE subject_node_id=? AND predicate='NEAR_KERNEL_TRAIT_GWAS_LOCUS' ORDER BY object_node_id''',(node,)):
        x=dict(row);x['attributes']=json.loads(x.pop('attributes_json'))
        x['source_url']='https://doi.org/10.1186/s12870-026-08285-4';d['kernel_trait_GWAS_context'].append(x)
    if r and r[0] is not None and r[1] is not None:
        ch=str(r[0]).removeprefix('chr');s=int(r[1])-1;e=int(r[2])
        for row in c.execute('''SELECT a.*,t.histone_mark,t.source_file,t.source_row,t.native_json,t.v5_anchor1_json,t.v5_anchor2_json
          FROM v31_contact_anchors a JOIN v31_contacts t ON a.contact_id=t.contact_id
          WHERE a.chromosome=? AND a.start0<? AND a.end0>? ORDER BY a.start0,a.contact_id LIMIT ?''',(ch,e,s,lim)):
            x=dict(row)
            for key in ['native_json','v5_anchor1_json','v5_anchor2_json']:
                value=x.pop(key);x[key.removesuffix('_json')]=json.loads(value) if value else None
            x['source_url']='https://doi.org/10.1038/s41467-025-65065-8'
            x['gene_anchor_overlap_is_geometric_not_functional']=True
            d['physical_contact_context']['overlapping_projected_anchors'].append(x)
    audit=c.execute('SELECT audit_json FROM v31_gene_training_context WHERE gene_id=?',(g,)).fetchone()
    d['training_overlap_context']=json.loads(audit[0]) if audit else {'modeled_gene':False,'gene_family_status':'not_ascertained'}
    sources=[]
    for row in c.execute('''SELECT DISTINCT a.source_id FROM v3_assertions a JOIN evidence v ON a.evidence_id=v.evidence_id
      JOIN edges e ON e.edge_id=v.edge_id WHERE e.subject_node_id=? OR e.object_node_id=?''',(node,node)):
        sources.append(row[0])
    if d['cell_type_and_cluster_expression']:sources.append('snrna_final_cells_2025')
    sources+= [x['source_id'] for x in d['peak_support']['counts']]
    d['public_source_audit']=[]
    for source in sorted(set(sources)):
        row=c.execute('SELECT audit_json FROM v31_source_audit WHERE source_id=?',(source,)).fetchone()
        if row:d['public_source_audit'].append(json.loads(row[0]))
    d['distinct_primary_publication_count_for_audited_sources']=len({x['doi'] for x in d['public_source_audit'] if x.get('doi')})
    d['boundaries'] += ['Author final single-nucleus cell annotations replace pre-Seurat summaries for cell-level interpretation; they are not spatial coordinates.',
      'Cell-type and cluster summaries use the same nuclei and are not separate experiments. Maternal placenta-chalazal cells are labeled separately.',
      'Peak proximity, cross-tissue HiChIP contact and GWAS proximity are distinct evidence categories; none alone proves causal regulation.',
      'The absence of a direct artifact or DOI match does not certify sample-level or label-level independence.']
    c.close();return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('gene');p.add_argument('--limit',type=int,default=30);p.add_argument('--db');a=p.parse_args()
    print(json.dumps(query(a.gene,a.limit,a.db),ensure_ascii=False,indent=2,allow_nan=False))
