"""Release QA and source-overlap audit. Never certifies independence by exclusion alone."""
from __future__ import annotations
import collections,csv,gzip,json,pathlib,random,re,sqlite3
from build_v3 import ROOT,PARENT,PRED,OLDRAW,sha,js
from query_v3 import query

def main():
    release=ROOT/'release';db=release/'mekg.sqlite'
    c=sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1',uri=True)
    checks={}
    checks['sqlite_integrity']=c.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    for col in ['subject_node_id','object_node_id']:
        checks['edge_'+col+'_exists']=c.execute(f'SELECT COUNT(*) FROM edges e LEFT JOIN nodes n ON n.node_id=e.{col} WHERE n.node_id IS NULL').fetchone()[0]==0
    checks['evidence_edge_exists']=c.execute('SELECT COUNT(*) FROM evidence e LEFT JOIN edges r ON r.edge_id=e.edge_id WHERE r.edge_id IS NULL').fetchone()[0]==0
    checks['new_evidence_no_training']=c.execute('SELECT COUNT(*) FROM evidence e JOIN v3_assertions a USING(evidence_id) WHERE e.used_for_training!=0').fetchone()[0]==0
    checks['no_blanket_independent_evidence']=c.execute("SELECT COUNT(*) FROM evidence WHERE independence_class='PUBLIC_INDEPENDENT_SUPPORT'").fetchone()[0]==0
    checks['numeric_observations_finite']=c.execute('SELECT COUNT(*) FROM observations WHERE value IS NOT NULL AND (value>1e308 OR value< -1e308)').fetchone()[0]==0
    checks['new_evidence_has_group']=c.execute("SELECT COUNT(*) FROM v3_assertions WHERE experiment_group IS NULL OR experiment_group=''").fetchone()[0]==0
    checks['new_observations_have_source']=c.execute('''SELECT COUNT(*) FROM observations o JOIN v3_added_records a ON a.table_name='observations' AND a.record_id=o.observation_id
        LEFT JOIN nodes n ON n.node_id=o.source_node_id WHERE n.node_id IS NULL''').fetchone()[0]==0
    checks['snrna_cpm_nonnegative']=c.execute("SELECT COUNT(*) FROM observations WHERE observation_type='SNRNA_PSEUDOBULK_CPM' AND value<0").fetchone()[0]==0
    checks['snrna_fraction_in_0_1']=c.execute("SELECT COUNT(*) FROM observations WHERE observation_type='SNRNA_DETECTION_FRACTION' AND (value<0 OR value>1)").fetchone()[0]==0
    conditions=c.execute("SELECT COUNT(DISTINCT condition_node_id) FROM observations WHERE observation_type='SNRNA_PSEUDOBULK_CPM'").fetchone()[0]
    checks['snrna_three_replicates']=conditions==3
    expected={}
    with gzip.open(PRED,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):expected[r['gene_id']]=r
    actual={r[0]:json.loads(r[1]) for r in c.execute('SELECT * FROM v3_model_predictions')}
    checks['all_model_scores_exactly_preserved']=actual==expected
    before=json.loads((ROOT/'preservation_before.json').read_text())
    checks['old_kg_hash_unchanged']=sha(PARENT/'mekg.sqlite')==before['parent_kg']
    checks['source_prediction_hash_unchanged']=sha(PRED)==before['frozen_model_predictions']
    checks['unknown_gene_not_fabricated']=query('Zm00001eb999999',db=db).get('status')=='not_found'
    # Exact raw source identifiers are used for overlap screening, not matching paper titles.
    old_registry=(PARENT/'public_source_registry.json').read_text()
    old_accessions=set(re.findall(r'\b(?:GSE\d+|GSM\d+|CRA\d+|PRJNA\d+)\b',old_registry))
    overlaps=[]
    for key,src_json in c.execute('SELECT source_id,source_json FROM v3_sources'):
        src=json.loads(src_json);overlap=sorted(set(src['accessions'])&old_accessions)
        overlaps.append({'source_id':key,'known_legacy_public_accession_overlap':overlap,
          'no_overlap_does_not_prove_independence':True,'training_raw_samples_audited':False,
          'label_or_model_selection_overlap_audited':False,'validation_eligibility':'not_certified'})
    counts={t:c.execute('SELECT COUNT(*) FROM '+t).fetchone()[0] for t in ['nodes','edges','evidence','observations']}
    added=dict(c.execute('SELECT table_name,COUNT(*) FROM v3_added_records GROUP BY table_name'))
    source_counts=dict(c.execute('SELECT source_id,COUNT(*) FROM v3_assertions GROUP BY source_id'))
    new_genes={r[0].removeprefix('gene:') for r in c.execute('''SELECT e.subject_node_id FROM edges e JOIN v3_added_records a ON a.table_name='edges' AND a.record_id=e.edge_id WHERE e.subject_node_id LIKE 'gene:Zm00001eb%'
        UNION SELECT e.object_node_id FROM edges e JOIN v3_added_records a ON a.table_name='edges' AND a.record_id=e.edge_id WHERE e.object_node_id LIKE 'gene:Zm00001eb%' ''')}
    covered_model=len(new_genes&expected.keys())
    # Pair overlap between old and new integrative networks: one shared relationship is not a new replication.
    oldpairs=set(c.execute("SELECT subject_node_id,object_node_id FROM edges WHERE predicate='PUBLIC_ENDOSPERM_REGULON_TARGETS'"))
    newpairs=set(c.execute("SELECT subject_node_id,object_node_id FROM edges WHERE predicate='PUBLIC_INTEGRATED_TF_TARGET'"))
    repeat_pairs=len(oldpairs&newpairs)
    source_gene_sets=collections.defaultdict(set)
    for source,sub,obj in c.execute('SELECT a.source_id,e.subject_node_id,e.object_node_id FROM v3_assertions a JOIN evidence v USING(evidence_id) JOIN edges e USING(edge_id)'):
        for node in [sub,obj]:
            if node.startswith('gene:Zm00001eb'):source_gene_sets[source].add(node.removeprefix('gene:'))
    tier_overlap={source:dict(collections.Counter(expected[gene]['tier'] if gene in expected else 'outside_model_universe' for gene in genes))
      for source,genes in source_gene_sets.items()}
    tested=random.Random(20260907).sample(sorted(expected),10)+['Zm00001eb125330']
    samples={g:query(g,limit=12,db=db) for g in tested}
    checks['gene_query_smoke']=all(s.get('gene_id')==g and s.get('prediction') is not None for g,s in samples.items())
    report={'passed':all(checks.values()),'checks':checks,'total_counts':counts,'added_counts':added,'new_source_evidence_counts':source_counts,
      'genes_with_new_edges':len(new_genes),'model_genes_with_new_edges':covered_model,
      'model_gene_denominator':len(expected),'old_new_integrated_regulon_shared_pairs':repeat_pairs,
      'source_overlap_audit':overlaps,'cell_type_annotation_status':'not_loaded; replicate-level filtered single-nucleus pseudobulk only',
      'source_genes_by_existing_model_label_tier':tier_overlap,
      'coverage_limitations':['No Hi-C/MapMan/QTL refresh performed in this release.','No genomewide MYB127 peak-to-gene liftover performed.',
        'No complete raw-sample/label/model-selection independence audit; these are contextual public evidence, not a new independent benchmark.',
        'Imprinting lists include B73 allele groups only; MEG/PEG direction and per-DAP RER were not present in those lists.'],
      'live_website_changed':False,'database_sha256':sha(db)}
    c.close()
    (release/'QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    out=ROOT/'examples';out.mkdir(exist_ok=True)
    for gene,payload in samples.items():(out/(gene+'.json')).write_text(json.dumps(payload,ensure_ascii=False,indent=2))
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
