"""One-time, documented semantic cleanup on this new release only, before QA."""
from __future__ import annotations
import csv,json,pathlib,sqlite3
from datetime import datetime,timezone
from build_v3 import ROOT,PUBLIC,sha,js

def main():
    db=ROOT/'release/mekg.sqlite';c=sqlite3.connect(db)
    c.execute("UPDATE public_aliases SET independence_class=? WHERE independence_class='PUBLIC_INDEPENDENT_SUPPORT'",(PUBLIC,))
    c.execute("UPDATE nodes SET attributes_json=json_set(attributes_json,'$.evidence_layer',?,'$.independence_status','legacy_not_certified') WHERE json_extract(attributes_json,'$.evidence_layer')='PUBLIC_INDEPENDENT_SUPPORT'",(PUBLIC,))
    # Relational pair overlap is explicitly shown; it does not manufacture a replicate.
    c.executescript('''CREATE TEMP TABLE old_regulon_pairs(subject_node_id TEXT,object_node_id TEXT,PRIMARY KEY(subject_node_id,object_node_id));
      INSERT OR IGNORE INTO old_regulon_pairs SELECT subject_node_id,object_node_id FROM edges WHERE predicate='PUBLIC_ENDOSPERM_REGULON_TARGETS';
      CREATE TABLE IF NOT EXISTS v3_legacy_regulon_overlap(subject_node_id TEXT,object_node_id TEXT,PRIMARY KEY(subject_node_id,object_node_id));
      INSERT OR IGNORE INTO v3_legacy_regulon_overlap
      SELECT DISTINCT n.subject_node_id,n.object_node_id FROM edges n JOIN old_regulon_pairs o
      ON n.subject_node_id=o.subject_node_id AND n.object_node_id=o.object_node_id
      WHERE n.predicate='PUBLIC_INTEGRATED_TF_TARGET';''')
    # Annotation is read on demand by the v3 query; no new assertion is counted.
    counts=dict(c.execute('SELECT table_name,COUNT(*) FROM v3_added_records GROUP BY table_name'))
    c.execute("INSERT OR REPLACE INTO metadata VALUES ('v3_added_counts',?)",(js(counts),))
    c.execute("INSERT OR REPLACE INTO metadata VALUES ('v3_public_support_summary_policy',?)",(js('public_gene_support and gene_cards are retained historical v2 snapshots; v3 query computes current counts and reads v3_model_predictions'),))
    parent_release=c.execute("SELECT value FROM metadata WHERE key='release_id'").fetchone()
    c.execute("INSERT OR IGNORE INTO metadata VALUES ('inherited_release_id',?)",(parent_release[0] if parent_release else js(None),))
    c.execute("INSERT OR REPLACE INTO metadata VALUES ('release_id',?)",(js('MEKG-public-v3-20260907'),))
    c.execute("INSERT OR REPLACE INTO metadata VALUES ('created_utc',?)",(js(datetime.now(timezone.utc).isoformat()),))
    row=c.execute("SELECT source_json FROM v3_sources WHERE source_id='mgrn_2025'").fetchone()
    sets=[]
    for name in ['mGRN.csv','mGRN_seed.csv']:
        with open(ROOT/'raw'/name) as f:sets.append({(r['TF'],r['Target']) for r in csv.DictReader(f)})
    src=json.loads(row[0]);src['file_schema_audit']={
      'full_unique_pairs':len(sets[0]),'full_unique_regulator_ids':len({r[0] for r in sets[0]}),'seed_unique_pairs':len(sets[1]),
      'seed_unique_regulator_ids':len({r[0] for r in sets[1]}),'shared_pairs':len(sets[0]&sets[1]),
      'interpretation':'Deposited network regulator count exceeds 513 profiled TFs; assay-level mapping not resolved. These edges are computational context, not direct binding.',
      'full_and_seed_are_not_subset_relations':True}
    c.execute("UPDATE v3_sources SET source_json=? WHERE source_id='mgrn_2025'",(js(src),))
    c.execute("UPDATE nodes SET attributes_json=? WHERE node_id='source:public:v3:mgrn_2025'",(js(src),))
    for key in ['e2f_prc2_2026','emf2b_2026']:
        r=c.execute('SELECT source_json FROM v3_sources WHERE source_id=?',(key,)).fetchone()
        s=json.loads(r[0]);s['derived_edge_calling_rules']={
          'MUTANT_TRANSCRIPTOMIC_RESPONSE':'FDR < 0.05 and absolute log2FC >= 1',
          'MUTANT_CHROMATIN_CHANGE':'FDR < 0.05 only; sign of log2FC; author Category/Group retained separately; not individual causal-target validation'}
        c.execute('UPDATE v3_sources SET source_json=? WHERE source_id=?',(js(s),key))
        c.execute('UPDATE nodes SET attributes_json=? WHERE node_id=?',(js(s),'source:public:v3:'+key))
    c.commit();source_updates={r[0]:json.loads(r[1]) for r in c.execute('SELECT source_id,source_json FROM v3_sources')};c.close()
    mf=ROOT/'release/manifest.json';m=json.loads(mf.read_text());m['database_sha256']=sha(db)
    m['new_sources']=source_updates
    m['retained_legacy_tables']={'gene_cards':'historical v8 snapshot; not current score endpoint','public_gene_support':'historical v2 summary; use query_v3 dynamic evidence counts'}
    mf.write_text(json.dumps(m,ensure_ascii=False,indent=2))
    print('semantic cleanup complete')

if __name__=='__main__':main()
