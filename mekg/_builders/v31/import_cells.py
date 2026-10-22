"""Import final-author snRNA context without changing clustering or model outputs."""
from pathlib import Path
import csv,gzip,zipfile,json,collections
from build_v31 import Graph31,ROOT,PARENT,js,num,uid
from table_io import sheets

def cells(g):
    source='snrna_final_cells_2025'
    if source in g.sources:
      # A previous constructor-only failure may leave the source and empty schema committed.
      for table in ['v31_cell_groups','v31_cell_expression','v31_cluster_markers']:
        assert g.c.execute('SELECT COUNT(*) FROM '+table).fetchone()[0]==0,'Refuse to overwrite partial cell data'
      assert g.c.execute('SELECT COUNT(*) FROM v3_assertions WHERE source_id=?',(source,)).fetchone()[0]==0
    else:
      source=g.source('snrna_final_cells_2025','Author final B73 19 DAP single-nucleus cell annotations',
      '10.1093/plcell/koaf162',['GSE270866_merged_seurat.rds.gz','kil_supplement.zip','kil_minimal.txt'],
      ['GSE270866','GSM8353850','GSM8353851','GSM8353852'],genotype='B73',tissue='dissected endosperm with maternal placenta-chalazal contamination',
      stage='19DAP',native_reference='B73 v4',
      normalization='RNA assay raw UMI sums divided by all native RNA features retained in the author RNA assay per group; SCT values not used; exact feature count in export QC',
      limitation='author annotation.predicted labels, not new validation; nuclei not biological replicates; author non-Singlet predictions retained without additional filtering',
      overlap='same experiment as snrna_19dap_2025; final object replaces pre-Seurat context for current cell summaries, not independent corroboration',
      final_nuclei=31506,cluster_count=23,cell_type_annotation_count=5)
    g.c.executescript('''CREATE TABLE IF NOT EXISTS v31_cell_groups(group_id TEXT PRIMARY KEY,group_level TEXT,replicate TEXT,label TEXT,
      nuclei INTEGER,UMI_total INTEGER,condition_node_id TEXT,attributes_json TEXT);
      CREATE TABLE IF NOT EXISTS v31_cell_expression(gene_id TEXT,group_id TEXT,UMI_sum INTEGER,detected_nuclei INTEGER,
      PRIMARY KEY(gene_id,group_id)) WITHOUT ROWID;
      CREATE TABLE IF NOT EXISTS v31_cluster_markers(marker_id TEXT PRIMARY KEY,gene_id TEXT,cluster TEXT,
      avg_log2FC REAL,padj REAL,pct1 REAL,pct2 REAL,source_row INTEGER,attributes_json TEXT);
      CREATE INDEX IF NOT EXISTS v31_markers_gene ON v31_cluster_markers(gene_id);
    ''')
    groups={};composition=collections.defaultdict(collections.Counter)
    with gzip.open(ROOT/'seurat_metadata.tsv.gz','rt') as f:
      for r in csv.DictReader(f,delimiter='\t'):
        composition[r['orig.ident']+'|'+r['seurat_clusters']][r['annotation.predicted']]+=1
    with (ROOT/'derived/seurat_group_stats.tsv').open() as f:
      for r in csv.DictReader(f,delimiter='\t'):
        gid=r['group_id'];n=int(r['nuclei']);total=int(float(r['UMI_total']));groups[gid]=(n,total)
        attrs={**r,'source_id':source,'author_label_field':'annotation.predicted' if r['group_level']=='cell_type' else 'seurat_clusters',
          'annotation_is_author_prediction':True,'not_spatial_localization':True,
          'maternal_tissue':r['label']=='Placenta-chalazal region','same_primary_experiment':'GSE270866_19DAP',
          'group_levels_overlap':'cell_type and cluster are two summaries of the same nuclei; never add them as independent samples'}
        if r['group_level']=='cluster':
          attrs['author_predicted_cell_type_composition']=dict(composition[r['replicate']+'|'+r['label']])
          assert sum(attrs['author_predicted_cell_type_composition'].values())==n
        attrs['group_label']=attrs.pop('label')
        cond=g.condition('snfinal:'+gid,'B73 19 DAP '+r['label']+' '+r['replicate'],**attrs)
        g.c.execute('INSERT INTO v31_cell_groups VALUES (?,?,?,?,?,?,?,?)',(gid,r['group_level'],r['replicate'],r['label'],n,total,cond,js(attrs)))
    assert len(groups)==84, len(groups)
    mapped={};counts=collections.Counter();umi=collections.Counter();batch=[]
    with gzip.open(ROOT/'derived/seurat_group_expression.tsv.gz','rt') as f:
      for r in csv.DictReader(f,delimiter='\t'):
        native=r['native_gene_id'];gid=r['group_id'];u=int(float(r['UMI_sum']));d=int(r['detected_nuclei'])
        n,total=groups[gid];assert 0<=d<=n and u>=d
        umi[gid]+=u;counts[gid]+=1
        if native not in mapped:mapped[native]=g.map(native,source,'RNA_assay:'+native)
        gene=mapped[native]
        if not gene:continue
        g.gene(gene);batch.append((gene,gid,u,d))
        if len(batch)>=25000:
          g.c.executemany('INSERT INTO v31_cell_expression VALUES (?,?,?,?)',batch);batch=[]
    if batch:g.c.executemany('INSERT INTO v31_cell_expression VALUES (?,?,?,?)',batch)
    assert set(counts)==set(groups)
    with (ROOT/'derived/seurat_export_qc.tsv').open() as f:
      qc={r['metric']:float(r['value']) for r in csv.DictReader(f,delimiter='\t')}
    assert all(counts[k]==int(qc['native_RNA_features']) and umi[k]==groups[k][1] for k in groups)
    marker_count=0;enriched=0
    with zipfile.ZipFile(ROOT/'raw/kil_supplement.zip') as z:
      for name,rows in sheets(z.read('Sup Dataset 2.xlsx')):
        cluster=name.split()[-1]
        for line,r in enumerate(rows[1:],2):
          if not r or not r[0]:continue
          gene=g.map(r[0],source,name+':'+str(line))
          if not gene:continue
          fc,p1,p2,padj=map(num,[r[2],r[3],r[4],r[5]])
          assert str(r[6])==cluster
          attrs={'source_file':'Sup Dataset 2.xlsx','sheet':name,'source_row':line,'native_gene_id':r[0],
            'avg_log2FC':fc,'pct1':p1,'pct2':p2,'padj':padj,
            'padj_method':'author Seurat p_val_adj; do not rename BH FDR','reported_zero_pvalue_is_numerical_underflow':padj==0,
            'unit_of_test':'nuclei; not independent biological replicates','author_cluster':'cluster '+cluster}
          mid=uid('marker:v31',source,cluster,gene)
          g.c.execute('INSERT INTO v31_cluster_markers VALUES (?,?,?,?,?,?,?,?,?)',(mid,gene,cluster,fc,padj,p1,p2,line,js(attrs)))
          marker_count+=1
          if fc is not None and fc>0 and padj is not None and padj<0.05:
            node=g.node('cellcluster:v31:GSE270866:'+cluster,'CellCluster','GSE270866 cluster '+cluster,
              {'author_cluster':cluster,'cell_type':None,
               'author_cluster_description':'EAS cluster' if cluster=='19' else None,
               'cluster_and_predicted_cell_type_are_distinct_fields':True,
               'cluster_type_assignment_source':'main text Fig.3D' if cluster=='19' else 'not assigned by this importer',
               'not_spatial_region':True})
            g.assertion(source,g.gene(gene),'AUTHOR_CLUSTER_ENRICHED_GENE',node,cluster+'|'+gene,
              'AUTHOR_SINGLE_NUCLEUS_CLUSTER_MARKER',attrs,stage='19DAP',tissue='dissected endosperm; author cluster '+cluster,
              genotype='B73',reference='B73 v4 gene identity mapped to v5',direction='higher_in_cluster',effect=fc,
              experiment='GSE270866_19DAP');enriched+=1
    g.stats['cell_expression_rows']=g.c.execute('SELECT COUNT(*) FROM v31_cell_expression').fetchone()[0]
    g.finish_module('cells',{'native_features':len(mapped),'uniquely_mapped_features':sum(v is not None for v in mapped.values()),
      'nuclei':31506,'groups':len(groups),'marker_records':marker_count,'positive_enriched_marker_records':enriched,
      'author_prediction_labels_not_new_cell_type_inference':True})

def kil_rna(g):
    source=g.source('kil_doublemutant_rna_2025','KIL1/KIL2 double-mutant EAS transcriptomic response',
      '10.1093/plcell/koaf162',['kil_supplement.zip','kil_minimal.txt'],['GSE270812'],
      genotype='A188 kil1-1 kil2-1 double mutant versus WT',stage='19DAP',tissue='EAS-enriched dissection',native_reference='B73 v4 gene models',
      limitation='author significant-gene list, not complete quantified genome; double mutant cannot attribute effects to either KIL alone; no direct-target claim',
      same_paper_as='snrna_final_cells_2025; different assay and genotype, not a second independent paper')
    cond=g.condition('KIL1_KIL2_RNA_19DAP','A188 kil1-1 kil2-1 versus WT, EAS-enriched, 19 DAP',
      genotype='A188',perturbed_genes_native=['Zm00001d003052','Zm00001d025861'],
      perturbation='double knockout',cannot_attribute_individual_gene_effects=True)
    with zipfile.ZipFile(ROOT/'raw/kil_supplement.zip') as z:
      for name,rows in sheets(z.read('Sup dataset 4.xlsx')):
        expected=-1 if 'downregulated' in name else 1
        for line,r in enumerate(rows[1:],2):
          gene=g.map(r[0],source,name+':'+str(line))
          if not gene:continue
          fc=num(r[2]);padj=num(r[6]);stat=num(r[4]);assert fc*expected>0
          attrs={'source_file':'Sup dataset 4.xlsx','sheet':name,'source_row':line,'native_gene_id':r[0],
            'log2FC':fc,'padj':padj,'baseMean':num(r[1]),'lfcSE':num(r[3]),'pvalue':num(r[5]),'author_stat':stat,
            'direction_basis':'author worksheet title and reported log2FoldChange agree',
            'stat_sign_disagrees_with_log2FC':stat*fc<0,'stat_not_recomputed_or_used_for_direction':True,
            'contrast':'double_mutant/WT as author worksheet labels','direct_target_proven':False}
          attrs['passes_padj_0_05_abs_log2FC_1']=padj is not None and padj<0.05 and abs(fc)>=1
          ev=None
          if attrs['passes_padj_0_05_abs_log2FC_1']:
            ev=g.assertion(source,g.gene(gene),'DIFFERENTIALLY_EXPRESSED_IN',cond,gene,
              'MUTANT_TRANSCRIPTOMIC_RESPONSE',attrs,stage='19DAP',tissue='EAS-enriched dissection',genotype='A188',reference='B73 v4',
              direction='down_in_mutant' if expected<0 else 'up_in_mutant',effect=fc,experiment='GSE270812_KIL_doublemutant_RNA')
          else:g.skipped['author_list_row_not_passing_padj_0_05_abs_log2FC_1']+=1
          g.observation(source,gene,'MUTANT_RNA_LOG2FC',cond,fc,'log2 fold-change',gene,attrs,ev,genotype='A188',reference='B73 v4')
    g.finish_module('kil_rna')

if __name__=='__main__':
    import sys
    module=sys.argv[1];assert module in ['cells','kil_rna']
    g=Graph31();g.c.execute('PRAGMA cache_size=-262144');g.c.execute('PRAGMA temp_store=MEMORY')
    if g.c.execute('SELECT 1 FROM v31_modules WHERE module=?',(module,)).fetchone():
      print('Already completed; no mutation',module);g.c.close();g.quarantine.close()
    else:globals()[module](g)
