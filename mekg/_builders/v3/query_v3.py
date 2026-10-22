"""Read-only gene query for MEKG v3, with unmodified frozen prediction lookup."""
from __future__ import annotations
import argparse,collections,json,pathlib,sqlite3
from mekg_paths import workspace
ROOT=workspace()/'v3'

def query(gene,limit=30,db=None):
    path=pathlib.Path(db) if db else ROOT/'release/mekg.sqlite'
    c=sqlite3.connect(path.resolve().as_uri()+'?mode=ro&immutable=1',uri=True);c.row_factory=sqlite3.Row
    candidates={r[0] for r in c.execute("SELECT canonical_id FROM nodes WHERE node_type='Gene' AND (canonical_id=? OR lower(label)=lower(?))",(gene,gene))}
    candidates.update(r[0] for r in c.execute('SELECT gene_id FROM public_aliases WHERE lower(alias)=lower(?)',(gene,)))
    candidates.update(r[0].removeprefix('gene:') for r in c.execute('SELECT canonical_node_id FROM aliases WHERE alias_id=?',(gene,)))
    candidates={x for x in candidates if x.startswith('Zm00001eb')}
    if len(candidates)!=1:
        c.close();return {'requested_gene_id':gene,'status':'ambiguous' if candidates else 'not_found','candidate_ids':sorted(candidates)}
    canonical=next(iter(candidates));node='gene:'+canonical
    row=c.execute('SELECT prediction_json FROM v3_model_predictions WHERE gene_id=?',(canonical,)).fetchone()
    pred=json.loads(row[0]) if row else None
    if pred:
        for k in ['priority_score','priority_percentile','uncertainty_sd','draw_sd','seed_sd_mean']:
            if pred.get(k) not in [None,'']:pred[k]=float(pred[k])
        pred['runs']=int(pred['runs']);pred['prediction_semantics']='frozen ranking score, not a causal probability'
    clauses="(e.subject_node_id=? OR e.object_node_id=?) AND e.independence_class LIKE 'PUBLIC_%'"
    counts=c.execute(f'''SELECT v.evidence_type,COUNT(*),COUNT(DISTINCT COALESCE(a.experiment_group,v.source_node_id))
        FROM edges e JOIN evidence v ON v.edge_id=e.edge_id LEFT JOIN v3_assertions a ON a.evidence_id=v.evidence_id
        WHERE {clauses} GROUP BY v.evidence_type''',(node,node)).fetchall()
    sourcecount=c.execute(f'''SELECT COUNT(DISTINCT COALESCE(a.experiment_group,v.source_node_id)) FROM edges e
        JOIN evidence v ON v.edge_id=e.edge_id LEFT JOIN v3_assertions a ON a.evidence_id=v.evidence_id WHERE {clauses}''',(node,node)).fetchone()[0]
    items=c.execute(f'''SELECT e.subject_node_id,e.predicate,e.object_node_id,n.label AS object_label,
      sn.label AS subject_label,v.evidence_id,v.evidence_type,v.source_url,v.development_stage,v.tissue_or_compartment,
      v.genotype,v.reference_genome,v.effect_direction,v.effect_size,v.attributes_json,
      a.experiment_group,a.source_id AS new_source_id
      FROM edges e JOIN evidence v ON v.edge_id=e.edge_id JOIN nodes n ON n.node_id=e.object_node_id
      JOIN nodes sn ON sn.node_id=e.subject_node_id LEFT JOIN v3_assertions a ON a.evidence_id=v.evidence_id
      WHERE {clauses} ORDER BY CASE v.evidence_type
      WHEN 'MAIZE_GENETIC_FUNCTION' THEN 0 WHEN 'BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT' THEN 1
      WHEN 'CURATED_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE' THEN 2 WHEN 'BINDING_CHROMATIN_RNA_INTERSECTION' THEN 3
      WHEN 'MUTANT_TRANSCRIPTOMIC_RESPONSE' THEN 4 WHEN 'MUTANT_CHROMATIN_CHANGE' THEN 5
      WHEN 'RECIPROCAL_CROSS_ALLELE_EXPRESSION' THEN 6 WHEN 'SINGLE_NUCLEUS_EXPRESSION_CONTEXT' THEN 7
      ELSE 20 END,v.evidence_id LIMIT ?''',(node,node,max(1,min(limit,500)))).fetchall()
    evidence=[]
    for row in items:
        d=dict(row);d['attributes']=json.loads(d.pop('attributes_json'));d['gene_role']='subject' if d['subject_node_id']==node else 'object'
        if d['predicate']=='PUBLIC_INTEGRATED_TF_TARGET':
            d['overlaps_legacy_regulon_pair']=bool(c.execute('SELECT 1 FROM v3_legacy_regulon_overlap WHERE subject_node_id=? AND object_node_id=?',(d['subject_node_id'],d['object_node_id'])).fetchone())
        d['independence_status']='not_certified';evidence.append(d)
    obs=[]
    for row in c.execute('''SELECT o.* FROM observations o JOIN v3_added_records a ON a.table_name='observations' AND a.record_id=o.observation_id
      WHERE o.subject_node_id=? ORDER BY o.observation_type,o.observation_id''',(node,)):
        d=dict(row);d['attributes']=json.loads(d.pop('attributes_json'));obs.append(d)
    internal=c.execute("SELECT COUNT(*) FROM edges WHERE (subject_node_id=? OR object_node_id=?) AND independence_class NOT LIKE 'PUBLIC_%'",(node,node)).fetchone()[0]
    paths=[]
    # A regulator's phenotype is contextual support, not evidence that its target has the same phenotype.
    for row in c.execute('''SELECT r.subject_node_id,r.predicate,r.edge_id AS relation_edge_id,
        p.object_node_id AS phenotype_node,p.edge_id AS phenotype_edge_id,n.label AS phenotype_label
        FROM edges r JOIN edges p ON p.subject_node_id=r.subject_node_id JOIN nodes n ON n.node_id=p.object_node_id
        WHERE r.object_node_id=? AND r.predicate IN ('REGULATES_WITH_EXPERIMENTAL_SUPPORT','PUBLIC_INTEGRATED_TF_TARGET')
        AND p.predicate IN ('HAS_GENETIC_PERTURBATION_PHENOTYPE','HAS_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE_EVIDENCE')
        ORDER BY CASE r.predicate WHEN 'REGULATES_WITH_EXPERIMENTAL_SUPPORT' THEN 0 ELSE 1 END,r.edge_id LIMIT 10''',(node,)):
        paths.append({**dict(row),'target_gene_id':canonical,'interpretation':'phenotyped regulator connects to queried gene; candidate mechanism hypothesis only',
          'target_phenotype_proven':False,'cross_context_transfer_requires_validation':True})
    meta=json.loads(c.execute("SELECT value FROM metadata WHERE key='public_v3_release'").fetchone()[0]);c.close()
    return {'requested_gene_id':gene,'gene_id':canonical,'coordinate_system':'B73 v5','kg_release':meta['release_id'],
      'prediction':pred,'prediction_source':meta['model_prediction_source'],'model_scores_modified':False,
      'public_summary':{'evidence_record_count':sum(r[1] for r in counts),'experiment_or_legacy_source_groups':sourcecount,
        'groups_are_not_certified_independent_replicates':True,'evidence_type_counts':{r[0]:r[1] for r in counts}},
      'internal_relation_count':internal,'public_evidence':evidence,'new_quantitative_observations':obs,'candidate_context_paths':paths,
      'boundaries':['Public evidence supplies biological context, not a causal probability or numerical model attribution.',
        'No record was used to change frozen model scores. Absence of evidence is not evidence of no biological effect.',
        'B73 v5 is the gene coordinate system; each record retains experimental genotype and native reference.',
        'A public source is not certified independent solely because it was excluded from this model run.']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('gene');p.add_argument('--limit',type=int,default=30);p.add_argument('--db');a=p.parse_args()
    print(json.dumps(query(a.gene,a.limit,a.db),ensure_ascii=False,indent=2,allow_nan=False))
