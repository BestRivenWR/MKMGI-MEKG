"""Finalize query indexes and explicitly bounded training-overlap audit."""
from pathlib import Path
import csv,gzip,json,collections,re,sqlite3
from build_v31 import ROOT,PARENT,old,sha,js,Graph31

def main():
    assert not (ROOT/'release/COMPLETE.json').exists()
    graph=Graph31();c=graph.c
    c.execute('PRAGMA cache_size=-2097152');c.execute('PRAGMA temp_store=MEMORY')
    c.execute('PRAGMA journal_mode=MEMORY')
    done={r[0] for r in c.execute('SELECT module FROM v31_modules')}
    assert done=={'tf_peaks','myb_peaks','hichip','gwas','cells','kil_rna','tss'},done
    # Even a profile with an unresolved TF may have a valid target gene and useful native peaks.
    for gene, in c.execute('''SELECT DISTINCT target_gene_id FROM v31_peaks WHERE target_gene_id IS NOT NULL
      UNION SELECT DISTINCT tf_gene_id FROM v31_peaks WHERE tf_gene_id IS NOT NULL''').fetchall():graph.gene(gene)
    c.executescript('''CREATE TABLE v31_contact_anchors(contact_id TEXT,anchor_number INTEGER,
      chromosome TEXT,start0 INTEGER,end0 INTEGER,PRIMARY KEY(contact_id,anchor_number)) WITHOUT ROWID;
      CREATE INDEX v31_anchor_interval ON v31_contact_anchors(chromosome,start0,end0);
      CREATE TABLE v31_source_audit(source_id TEXT PRIMARY KEY,audit_json TEXT);
      CREATE TABLE v31_gene_training_context(gene_id TEXT PRIMARY KEY,audit_json TEXT);
    ''')
    for cid,a,b in c.execute('SELECT contact_id,v5_anchor1_json,v5_anchor2_json FROM v31_contacts').fetchall():
        for i,val in enumerate([a,b],1):
            if val:c.execute('INSERT INTO v31_contact_anchors VALUES (?,?,?,?,?)',(cid,i,*json.loads(val)))
    # The rebuilt peak modules already carry explicit per-record projection semantics.
    # Do not rewrite millions of otherwise unchanged evidence records on shared storage.
    project=ROOT.parent
    manifest_path=project/'audit_inputs/model_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    source_manifest=project/'audit_inputs/source_manifest.tsv'
    declared=dict(manifest['files'])
    with source_manifest.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):declared[r['path']]=r['sha256']
    label_path=project/'audit_inputs/labels_v9.tsv.gz'
    assert sha(label_path)==next(v for k,v in manifest['files'].items() if k.replace('\\','/').endswith('/labels_v9.tsv.gz')),'Labels changed since frozen manifest'
    with gzip.open(label_path,'rt') as f:labels={r['gene_id']:r for r in csv.DictReader(f,delimiter='\t')}
    p1={g for g,r in labels.items() if r['tier'].startswith('P1')}
    fams={labels[g]['family_id'] for g in p1}
    primary=collections.defaultdict(set)
    for gene,r in labels.items():
        for doi in re.findall(r'10\.\d{4,9}/[^\s;,|]+',r['primary_doi'].lower()):primary[doi].add(gene)
        data={'tier':r['tier'],'family_id':r['family_id'],'is_training_benchmark_P1':gene in p1,
          'shares_family_with_training_benchmark_P1':r['family_id'] in fams,'primary_doi':r['primary_doi'] or None,
          'family_scope':'existing model universe; absent external genes have not been newly assigned eggNOG families',
          'not_an_external_holdout_certificate':True}
        c.execute('INSERT INTO v31_gene_training_context VALUES (?,?)',(gene,js(data)))
    registry=project/'base/public_source_registry.json'
    regtext=registry.read_text();legacy=set(re.findall(r'\b(?:GSE|GSM|CRA|PRJNA)\d+\b',regtext))
    audit_files=[manifest_path,source_manifest,label_path,registry,
      project/'audit_inputs/audit_report.json',
      project/'audit_inputs/feature_manifest.tsv',
      project/'audit_inputs/existing_data_inventory.json']
    training_text='\n'.join(f.read_text() for f in audit_files if f.suffix!='.gz')
    declared_accessions=set(re.findall(r'\b(?:GSE|GSM|CRA|PRJNA)\d+\b',training_text))
    # The legacy registry is not a training manifest: remove it from the accession search above.
    training_text='\n'.join(f.read_text() for f in audit_files if f.suffix!='.gz' and f!=registry)
    declared_accessions=set(re.findall(r'\b(?:GSE|GSM|CRA|PRJNA)\d+\b',training_text))
    sources={r[0]:json.loads(r[1]) for r in c.execute('SELECT source_id,source_json FROM v3_sources')}
    source_genes=collections.defaultdict(set)
    for source,sub,obj in c.execute('''SELECT a.source_id,e.subject_node_id,e.object_node_id FROM v3_assertions a
       JOIN evidence v ON a.evidence_id=v.evidence_id JOIN edges e ON e.edge_id=v.edge_id'''):
        for node in [sub,obj]:
            if node.startswith('gene:Zm00001eb'):source_genes[source].add(node[5:])
    source_genes['snrna_final_cells_2025'].update(r[0] for r in c.execute('SELECT DISTINCT gene_id FROM v31_cell_expression'))
    for r in c.execute('SELECT DISTINCT source_id,target_gene_id FROM v31_peaks WHERE target_gene_id IS NOT NULL'):
        source_genes[r[0]].add(r[1])
    audit=[]
    for source,src in sources.items():
        genes=source_genes[source];doi=src['doi'].lower();matches=[]
        for file,h in src.get('source_file_sha256',{}).items():
            for trainpath,trainhash in declared.items():
                if h==trainhash:matches.append({'public_file':file,'declared_training_file':trainpath,'sha256':h})
        d={'source_id':source,'doi':doi,'accessions':src['accessions'],
          'same_primary_publication_sources':sorted(k for k,v in sources.items() if v['doi'].lower()==doi and k!=source),
          'direct_file_hash_matches_to_declared_training_inputs':matches,
          'accession_matches_in_training_manifests':sorted(set(src['accessions'])&declared_accessions),
          'legacy_public_registry_accession_matches':sorted(set(src['accessions'])&legacy),
          'source_DOI_matches_in_model_labels':sorted(primary.get(doi,set())),
          'label_tier_counts_among_source_genes':dict(collections.Counter(labels[g]['tier'] if g in labels else 'outside_model_universe' for g in genes)),
          'source_genes_sharing_P1_family_in_modeled_universe':sum(g in labels and labels[g]['family_id'] in fams for g in genes),
          'genes_without_available_family_audit':sum(g not in labels for g in genes),
          'direct_manifest_screen_completed':True,'sample_equivalence_and_derived_data_reuse_fully_audited':False,
          'historical_label_literature_fully_audited':False,'historical_model_selection_reuse_fully_audited':False,
          'validation_eligibility':'PUBLIC_SUPPORT_UNAUDITED','certified_independent':False,
          'reason':'Model manifests lack complete raw-sample accessions; only 39/17917 label rows contain primary DOI; no match is not proof of independence.'}
        c.execute('INSERT INTO v31_source_audit VALUES (?,?)',(source,js(d)));audit.append(d)
    payload={'scope':'file-level and available label/family audit; not independent-validation certification',
      'audited_file_sha256':{str(f):sha(f) for f in audit_files},
      'training_artifact_hash_basis':'declared hashes from frozen manifest and original source manifest; no wholesale rehash of all training arrays',
      'label_file_hash_verified_against_frozen_manifest':True,
      'model_label_rows':len(labels),'primary_doi_populated':sum(bool(r['primary_doi']) for r in labels.values()),
      'training_P1_genes':len(p1),'training_P1_families':len(fams),'sources':audit,
      'user_reported_training_spatial_paper':'10.1038/s41467-023-43006-7; accession-level transformation chain remains unverified',
      'independent_support_count':0}
    (ROOT/'release/TRAINING_OVERLAP_AUDIT.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2))
    meta={'release_id':'MEKG-public-v3.1-20260907','parent_database':str(PARENT/'release/mekg.sqlite'),
      'current_cell_context':'snrna_final_cells_2025','historical_pre_Seurat_context':'snrna_19dap_2025',
      'cell_expression_table':'v31_cell_expression','model_prediction_table':'v3_model_predictions',
      'public_source_audit_table':'v31_source_audit','live_website_changed':False}
    c.execute('INSERT INTO metadata VALUES (?,?)',('public_v31_release',js(meta)))
    c.execute("UPDATE metadata SET value=? WHERE key='release_id'",(js(meta['release_id']),))
    c.commit();c.execute('PRAGMA journal_mode=DELETE');c.close();graph.quarantine.close()
    print('FINALIZED',js({'sources':len(sources),'label_DOIs_populated':payload['primary_doi_populated'],'independent_support_count':0}),flush=True)

if __name__=='__main__':main()
