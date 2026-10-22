"""Add traceable public evidence to a separate immutable-parent MEKG snapshot.

No label, score, training artifact, legacy KG or web file is modified.
No public source is certified independent without a sample-level overlap audit.
"""
from __future__ import annotations
import argparse, collections, csv, gzip, hashlib, json, math, pathlib, re, shutil, sqlite3, time, zipfile
from table_io import sheets

from mekg_paths import workspace
PROJECT=workspace()
ROOT=PROJECT/'v3'
PARENT=PROJECT/'base'
OLDRAW=PROJECT/'raw/maizegdb'
PRED=PROJECT/'inputs/parent_predictions.tsv.gz'
PUBLIC='PUBLIC_SUPPORT_UNAUDITED'

def js(x):return json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False)
def uid(prefix,*x):return prefix+':'+hashlib.sha256(js(x).encode()).hexdigest()[:24]
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def num(x):
    try:
        y=float(x);return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None
def records(rows,header_index=0):
    header=rows[header_index]
    for n,row in enumerate(rows[header_index+1:],header_index+2):
        yield n,{k:row[i] if i<len(row) else '' for i,k in enumerate(header) if k}

class Graph:
    def __init__(self,db,resume=False):
        self.c=sqlite3.connect(db);self.c.execute('PRAGMA synchronous=NORMAL')
        self.c.execute('PRAGMA journal_mode=MEMORY');self.c.execute('PRAGMA temp_store=MEMORY')
        self.known={r[0] for r in self.c.execute('SELECT node_id FROM nodes')}
        self.model_genes={r[0] for r in self.c.execute('SELECT gene_id FROM gene_cards')}
        self.metadata={};self.mapping={};self.skipped=collections.Counter();self.stats=collections.Counter()
        self.sources={};self.new_edges=set();self.new_evidence=set();self.new_observations=set()
        self.quarantine=gzip.open(ROOT/'mapping_quarantine.jsonl.gz','at' if resume else 'wt',encoding='utf8')
        with gzip.open(OLDRAW/'Zm00001eb.1.fulldata.txt.gz','rt') as f:
            for row in csv.reader(f,delimiter='\t'):
                if len(row)>=14 and re.fullmatch(r'Zm00001eb\d{6}',row[1]):self.metadata[row[1]]=row
        for version in [3,4]:
            raw=collections.defaultdict(set);reverse=collections.defaultdict(set)
            with open(OLDRAW/f'B73v{version}_to_B73v5.tsv') as f:
                for r in csv.reader(f,delimiter='\t'):
                    if len(r)<2:continue
                    for target in r[1].split(','):
                        if target in self.metadata:raw[r[0]].add(target);reverse[target].add(r[0])
            for key,targets in raw.items():
                if len(targets)==1 and len(reverse[next(iter(targets))])==1:self.mapping[key]=next(iter(targets))
        if resume:
            state=json.loads((ROOT/'base_checkpoint.json').read_text())
            self.stats.update(state['stats']);self.skipped.update(state['skipped']);self.sources=state['sources']
            return
        self.c.executescript('''
        CREATE TABLE v3_sources(source_id TEXT PRIMARY KEY, doi TEXT, accession_json TEXT,
          independence_status TEXT, source_json TEXT NOT NULL);
        CREATE TABLE v3_assertions(evidence_id TEXT PRIMARY KEY, experiment_group TEXT NOT NULL,
          evidence_grade TEXT NOT NULL, source_id TEXT NOT NULL, biological_causality_claim INTEGER NOT NULL);
        CREATE TABLE v3_added_records(table_name TEXT,record_id TEXT,PRIMARY KEY(table_name,record_id));
        CREATE TABLE v3_model_predictions(gene_id TEXT PRIMARY KEY, prediction_json TEXT NOT NULL);
        CREATE INDEX idx_v3_experiment ON v3_assertions(experiment_group);
        ''')
    def map(self,raw,source,record):
        native=str(raw).strip();native=re.sub(r'_(?:T|P)\d+$','',native)
        gene=native if native in self.metadata else self.mapping.get(native)
        if not gene:
            self.skipped[source+':unmapped_or_not_bijective']+=1
            self.quarantine.write(js({'source':source,'record':record,'native_id':raw,'reason':'unmapped_or_not_bijective'})+'\n')
        return gene
    def node(self,id,kind,label,attrs=None,**kw):
        if id in self.known:return id
        vals=[id,kind,kw.get('canonical_id',id),label,kw.get('species',''),kw.get('genotype',''),kw.get('reference_genome',''),kw.get('chromosome'),kw.get('start'),kw.get('end'),kw.get('strand'),js(attrs or {})]
        self.c.execute('INSERT INTO nodes VALUES ('+','.join('?'*12)+')',vals);self.known.add(id)
        self.c.execute('INSERT INTO v3_added_records VALUES (?,?)',('nodes',id));self.stats['nodes']+=1
        return id
    def gene(self,g):
        r=self.metadata[g]
        self.node('gene:'+g,'Gene',r[10] if r[10] not in ['','-'] else g,
          {'description':r[12],'coordinate_convention':'1-based closed','assembly_source':'MaizeGDB B73 v5'},
          canonical_id=g,species='Zea mays',genotype='B73',reference_genome='B73 v5',chromosome=r[7],start=int(r[8]),end=int(r[9]))
        return 'gene:'+g
    def source(self,key,title,doi,files,accessions,**context):
        data=dict(source_id=key,title=title,doi=doi,accessions=accessions,files=files,
          source_file_sha256={x:sha(ROOT/'raw'/x) for x in files},
          used_for_training=False,used_for_model_selection=False,used_for_benchmark=False,
          independence_status='not_certified_sample_level_audit_pending',
          evidence_layer=PUBLIC,**context)
        self.sources[key]=data
        self.c.execute('INSERT INTO v3_sources VALUES (?,?,?,?,?)',(key,doi,js(accessions),data['independence_status'],js(data)))
        self.node('source:public:v3:'+key,'PublicDataSource',title,data,canonical_id=key)
        return key
    def assertion(self,source,sub,pred,obj,record,grade,attrs=None,*,stage='',tissue='',genotype='',reference='',direction='not_assessed',effect=None,experiment=None,causal=False,directed=1):
        assert sub in self.known and obj in self.known,(sub,obj)
        source_node='source:public:v3:'+source
        edge=uid('edge:v3',source,sub,pred,obj,record);ev=uid('evidence:v3',edge)
        if edge in self.new_edges:self.stats['duplicate_rows_deduplicated']+=1;return ev
        self.new_edges.add(edge);self.new_evidence.add(ev)
        group=experiment or source
        data={'evidence_layer':PUBLIC,'used_for_training':False,'independence_status':'not_certified_sample_level_audit_pending',
          'experiment_group':group,'evidence_grade':grade,'biological_causality_claim':causal,
          'coordinate_mapping':'gene_identity_only_no_peak_liftover','record':record,**(attrs or {})}
        self.c.execute('INSERT INTO edges VALUES (?,?,?,?,?,?,?,?,?,?,?)',(edge,sub,pred,obj,directed,None,1,0,PUBLIC,source_node,js(data)))
        self.c.execute('INSERT INTO evidence VALUES ('+','.join('?'*18)+')',
          (ev,edge,grade,direction,effect,None,source_node,record,'https://doi.org/'+self.sources[source]['doi'],genotype,reference,stage,tissue,0,PUBLIC,'B73_v5_gene_identity_mapping;see_native_ids',str(ROOT/'raw'),js(data)))
        self.c.execute('INSERT INTO v3_assertions VALUES (?,?,?,?,?)',(ev,group,grade,source,int(causal)))
        self.c.executemany('INSERT INTO v3_added_records VALUES (?,?)',[('edges',edge),('evidence',ev)])
        self.stats['edges']+=1;self.stats['source:'+source]+=1;self.stats['grade:'+grade]+=1
        return ev
    def observation(self,source,gene,kind,condition,value,unit,record,attrs=None,ev=None,genotype='',reference=''):
        if value is None:return
        mod=self.node('modality:v3:'+kind,'OmicsModality',kind)
        oid=uid('obs:v3',source,gene,kind,condition,record)
        if oid in self.new_observations:raise ValueError('duplicate observation '+oid)
        self.new_observations.add(oid)
        data={'evidence_layer':PUBLIC,'source_record_id':record,'independence_status':'not_certified_sample_level_audit_pending',**(attrs or {})}
        self.c.execute('INSERT INTO observations VALUES ('+','.join('?'*15)+')',(oid,self.gene(gene),kind,mod,condition,value,unit,'see attributes_json','source:public:v3:'+source,genotype,reference,0,0,ev,js(data)))
        self.c.execute('INSERT INTO v3_added_records VALUES (?,?)',('observations',oid));self.stats['observations']+=1
    def condition(self,key,label,**attrs):return self.node('condition:v3:'+key,'ExperimentalCondition',label,attrs)
    def checkpoint(self,label):
        self.c.commit();print(label,js(dict(self.stats)),flush=True)

def mgrn(g):
    source=g.source('mgrn_2025','mGRN+ integrated maize TF-target context','10.1016/j.molp.2025.08.009',
      ['mGRN.csv','mGRN_seed.csv','zenodo_mgrn.json'],['CRA022036','GSE201701'],
      native_reference='B73 v4 gene identifiers',assay='integrated TF binding + accessibility + expression',
      genotype='mixed/assay-specific; do not infer from reference',license='CC-BY-4.0 for Zenodo record',
      known_overlap='GSE201701 is present in legacy public registry; tissue and full networks are overlapping derived resources, not independent experiments',
      limitation='not a universally in vivo validated or endosperm-specific regulatory network')
    with open(ROOT/'raw/mGRN_seed.csv') as f:seed={(r['TF'],r['Target']) for r in csv.DictReader(f)}
    full=set()
    with open(ROOT/'raw/mGRN.csv') as f:
        full={(r['TF'],r['Target']) for r in csv.DictReader(f)}
    for n,pair in enumerate(sorted(full|seed),1):
            a=g.map(pair[0],source,n);b=g.map(pair[1],source,n)
            if not a or not b:continue
            g.assertion(source,g.gene(a),'PUBLIC_INTEGRATED_TF_TARGET',g.gene(b),'|'.join(pair),
              'COMPUTATIONAL_INTEGRATED_REGULON',{'native_tf':pair[0],'native_target':pair[1],
              'seed_module':pair in seed,'full_network':pair in full,'source_record_locator':'TF+Target pair in original CSV',
              'direct_in_vivo_binding':False},
              tissue='seed module (not isolated endosperm)' if pair in seed else 'multi-tissue integrated network',
              reference='B73 v4 gene IDs',experiment='mGRN2025_integrated_TFBS_context')
    g.stats['mgrn_full_unique_pairs']=len(full);g.stats['mgrn_seed_unique_pairs']=len(seed)
    g.stats['mgrn_shared_pairs']=len(seed&full);g.stats['mgrn_union_unique_pairs']=len(seed|full)
    g.stats['mgrn_seed_pairs_absent_from_full']=len(seed-full)
    g.checkpoint('mGRN imported')

def myb127(g):
    source=g.source('myb127_2026','MYB127 endosperm genetics and regulation','10.1038/s41477-026-02238-3',
      ['myb127_tables.xlsx','myb127_binding.xlsx','myb127_rna.xlsx','myb127.html'],['CRA028659','CRA030621'],
      native_reference='B73 RefGen v4',genotype='B73 knockout; FLAG overexpression material assay-specific',
      limitation='genomewide native peaks retained as raw only; no unverified v4-to-v5 coordinate conversion')
    myb=g.map('Zm00001d041935',source,'main-text-gene-id');mybn=g.gene(myb)
    ph=g.node('phenotype:v3:myb127_knockout','Phenotype','Reduced kernel weight and nutrient deposition; altered filling-stage aleurone morphology',
      {'perturbation':'two independent CRISPR alleles and allelic test','figure':'Fig.1;Extended Data Fig.1'})
    g.assertion(source,mybn,'HAS_GENETIC_PERTURBATION_PHENOTYPE',ph,'Fig1_EDFig1','MAIZE_GENETIC_FUNCTION',
      {'alleles':['Zmmyb127-1','Zmmyb127-2'],'evidence_location':'Results first subsection and Fig.1; Extended Data Fig.1'},
      genotype='B73',tissue='kernel/endosperm',stage='10/16DAP and mature kernels',experiment='MYB127_CRISPR_genetics',causal=True)
    st=dict(sheets(ROOT/'raw/myb127_tables.xlsx'))
    targets={}
    for n,row in records(st['Supplementary Table 2'],1):
        if not re.fullmatch('Zm00001d\\d{6}',row.get('Gene ID','')):continue
        gene=g.map(row['Gene ID'],source,'Table2:'+str(n))
        if gene:targets[row['Gene name']]=(gene,row['Gene ID'])
    for name,(gene,native) in targets.items():
        positive=name.upper() in ['NKD1','NKD2']
        g.assertion(source,mybn,'REGULATES_WITH_EXPERIMENTAL_SUPPORT',g.gene(gene),'Fig2-6:'+name,
          'BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT',{'native_target':native,'target_symbol':name,
          'evidence_location':'Fig.2 and Supplementary Table2; Fig.3-6',
          'context':'O2 coactivation complex' if positive else 'LUG3/ABI4 repression complex',
          'assays':['CUT&Tag FLAG overexpression','EMSA','reporter','mutant RNA'],
          'effect_scope':'context-specific regulation; not a universal TF sign'},
          stage='16DAP and reporter assay contexts',tissue='endosperm; reporter tissue per paper',
          genotype='B73 knockout; FLAG overexpression and reporter materials per source',reference='B73 v4',
          direction='activation' if positive else 'repression',experiment='MYB127_mechanism_'+name)
    rows=dict(sheets(ROOT/'raw/myb127_rna.xlsx'))['Sheet1']
    condition=g.condition('MYB127_RNA_16DAP','Zmmyb127 versus WT, 16 DAP endosperm',genotype='B73',replicates_per_group=3)
    for n,row in records(rows,1):
        gene=g.map(row['GeneID'],source,'RNA:'+str(n))
        if not gene:continue
        fc=num(row['logFC']);p=num(row['Pvalue'])
        attrs={'native_gene_id':row['GeneID'],'pvalue':p,'pvalue_type':'author column Pvalue; not presented as adjusted FDR',
          'logFC_base':'log2 per paper','contrast':'mutant/WT',
          'replicate_fpkm':{x:num(row[x]) for x in ['WT-1','WT-2','WT-3','Zmmyb127-1','Zmmyb127-2','Zmmyb127-3']}}
        # Keep continuous results without inventing FDR or using unadjusted P for genomewide target claims.
        g.observation(source,gene,'MUTANT_RNA_LOG2FC',condition,fc,'log2 fold-change','RNA:'+str(n),attrs,genotype='B73',reference='B73 v4')
    g.checkpoint('MYB127 imported')

def emf2b(g):
    source=g.source('emf2b_2026','EMF2b-1/PRC2 mutant and binding evidence','10.1186/s13059-026-04049-3',
      ['emf2b_tables.zip','emf2b.html'],['CRA025416','CRA025432','CRA025452','CRA031333','CRA032760'],
      native_reference='B73 v5 gene identifiers',genotype='B73-backcrossed mutant; B104 transgenic overexpression',
      limitation='bulk samples are developing kernels; peak-associated gene is not automatically a functional direct target')
    complexnode=g.node('complex:v3:EMF2b1_PRC2','ProteinComplex','EMF2b-1-containing PRC2',{'not_a_single_gene':True})
    with zipfile.ZipFile(ROOT/'raw/emf2b_tables.zip') as z:
        for name,rows in sheets(z.read('Table S1.xlsx')):
            stage='10DAP' if '10 DAP' in name else '14DAP'
            condition=g.condition('EMF_RNA_'+stage,'emf2b-1–1 versus WT '+stage+' kernels',genotype='B73-backcrossed mutant',tissue='kernels',replicates_per_group=3)
            for n,row in records(rows):
                native=row.get('Gene',row.get('Gene model'));gene=g.map(native,source,name+':'+str(n))
                if not gene:continue
                fc=num(row['logFC']);fdr=num(row['FDR'])
                attrs={'native_gene_id':native,'log2FC':fc,'FDR':fdr,'contrast':'mutant/WT','source_table':'Table S1.xlsx/'+name,'source_row':n}
                ev=None
                if fc is not None and fdr is not None and fdr<0.05 and abs(fc)>=1:
                    ev=g.assertion(source,g.gene(gene),'DIFFERENTIALLY_EXPRESSED_IN',condition,name+':'+str(n),
                      'MUTANT_TRANSCRIPTOMIC_RESPONSE',attrs,stage=stage,tissue='kernels',genotype='B73-backcrossed mutant',reference='B73 v5',
                      direction='up_in_mutant' if fc>0 else 'down_in_mutant',effect=fc,experiment='EMF2b_mutant_RNA_'+stage)
                g.observation(source,gene,'MUTANT_RNA_LOG2FC',condition,fc,'log2 fold-change',name+':'+str(n),attrs,ev,genotype='B73-backcrossed mutant',reference='B73 v5')
        for name,rows in sheets(z.read('Table S6.xlsx')):
            unique=set()
            for n,row in enumerate(rows,1):
                if not row or not re.fullmatch(r'Zm00001eb\d{6}',row[0]):continue
                gene=g.map(row[0],source,'TableS6:'+str(n))
                if not gene or gene in unique:continue
                unique.add(gene)
                g.assertion(source,complexnode,'CANDIDATE_REPRESSIVE_TARGET',g.gene(gene),'TableS6:'+gene,
                  'BINDING_CHROMATIN_RNA_INTERSECTION',{'evidence_location':'Table S6',
                  'assays':['EMF2b-1-HA ChIP','reduced H3K27me3 in mutant','increased RNA in mutant'],
                  'directness':'author candidate intersection; individual genetic causality untested'},
                  stage='14DAP',tissue='developing kernels',genotype='B73-backcrossed mutant/B104-origin transgene',
                  reference='B73 v5',direction='consistent_with_repression',experiment='EMF2b_PRC2_intersection')
            g.stats['emf2b_intersection_unique_genes']=len(unique)
    g.checkpoint('EMF2b imported')

def e2f(g):
    source=g.source('e2f_prc2_2026','FIE2/E2F endosperm chromatin and mutant RNA','10.1038/s41477-026-02386-6',
      ['e2f_tables.xlsx','e2f.html'],[],native_reference='B73 v4 gene identifiers',
      genotype='WT/fie2-c1/e2f12-c1;e2f13-1; accession and exact sample background pending',
      limitation='raw accessions not recovered; all claims anchored to public supplementary tables; no independent-validation certification')
    st=dict(sheets(ROOT/'raw/e2f_tables.xlsx'))
    c=g.condition('fie2_RNA_12DAP','fie2-c1 versus WT, 12 DAP endosperm',stage='12DAP',tissue='endosperm')
    for n,row in records(st['Wang_Supplementary_Table_6'],1):
        native=row.get('Gene ID','');gene=g.map(native,source,'Table6:'+str(n))
        if not gene:continue
        fc=num(row.get('log2FoldChange'));fdr=num(row.get('padj'));attrs={'native_gene_id':native,'log2FC':fc,'FDR':fdr,'contrast':'fie2-c1/WT','source_table':'Table6','source_row':n}
        ev=None
        if fc is not None and fdr is not None and fdr<0.05 and abs(fc)>=1:
            ev=g.assertion(source,g.gene(gene),'DIFFERENTIALLY_EXPRESSED_IN',c,'Table6:'+str(n),
              'MUTANT_TRANSCRIPTOMIC_RESPONSE',attrs,stage='12DAP',tissue='endosperm',reference='B73 v4',
              direction='up_in_mutant' if fc>0 else 'down_in_mutant',effect=fc,experiment='FIE2_RNA_12DAP')
        g.observation(source,gene,'MUTANT_RNA_LOG2FC',c,fc,'log2 fold-change','Table6:'+str(n),attrs,ev,reference='B73 v4')
    # Gene-body chromatin comparison, not a FIE2/E2F binding claim or coordinate liftover.
    for table,mutant in [(5,'fie2-c1'),(10,'e2f12-c1;e2f13-1')]:
        condition=g.condition('E2F_K27_Table'+str(table),mutant+' versus WT H3K27me3 at 12 DAP',stage='12DAP',tissue='endosperm')
        for n,row in records(st['Wang_Supplementary_Table_'+str(table)],1):
            native=row.get('Gene ID',row.get('GeneID',''));gene=g.map(native,source,f'Table{table}:{n}')
            if not gene:continue
            fc=num(row.get('log2FoldChange',row.get('logFC')));fdr=num(row.get('Adj. P. value',row.get('adj.P.Val')))
            attrs={'native_gene_id':native,'log2FC':fc,'FDR':fdr,'category':row.get('Category',row.get('Group')),
              'native_chr':row.get('Chr'),'native_start':row.get('start',row.get('Start')),'native_end':row.get('end',row.get('End')),
              'native_coordinates_display_only':True,'contrast':mutant+'/WT','source_table':table,'source_row':n}
            ev=None
            if fc is not None and fdr is not None and fdr<0.05:
                ev=g.assertion(source,g.gene(gene),'HAS_DIFFERENTIAL_H3K27ME3',condition,f'Table{table}:{n}',
                  'MUTANT_CHROMATIN_CHANGE',attrs,stage='12DAP',tissue='endosperm',reference='B73 v4 native coordinates',
                  direction='gain_in_mutant' if fc>0 else 'loss_in_mutant',effect=fc,experiment=f'E2F_PR C2_K27_Table{table}'.replace(' ',''))
            g.observation(source,gene,'H3K27ME3_LOG2FC',condition,fc,'log2 fold-change',f'Table{table}:{n}',attrs,ev,reference='B73 v4')
    g.checkpoint('E2F/PRC2 imported')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base-only',action='store_true');parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    release=ROOT/'release'
    if args.resume:
        if (release/'manifest.json').exists():raise RuntimeError('Completed release is immutable; do not resume')
        for i in range(3):
            if not (ROOT/'raw'/f'GSM835385{i}_Rep{i+1}_filtered_matrix.h5').exists():raise RuntimeError('snRNA download still incomplete')
        before=json.loads((ROOT/'preservation_before.json').read_text())
        g=Graph(release/'mekg.sqlite',resume=True)
    else:
        if release.exists():raise RuntimeError('release exists: preserve it; select a new release directory before rebuilding')
        release.mkdir()
        before={'parent_kg':sha(PARENT/'mekg.sqlite'),'frozen_model_predictions':sha(PRED)}
        (ROOT/'preservation_before.json').write_text(js(before))
        shutil.copy2(PARENT/'mekg.sqlite',release/'mekg.sqlite')
        g=Graph(release/'mekg.sqlite')
        # Correct the old blanket independence label in the new copy only.
        for table in ['edges','evidence']:
            g.c.execute(f"UPDATE {table} SET independence_class=?,attributes_json=json_set(attributes_json,'$.evidence_layer',?,'$.independence_status','legacy_not_certified') WHERE independence_class='PUBLIC_INDEPENDENT_SUPPORT'",(PUBLIC,PUBLIC))
        g.c.execute("UPDATE observations SET attributes_json=json_set(attributes_json,'$.evidence_layer',?,'$.independence_status','legacy_not_certified') WHERE json_extract(attributes_json,'$.evidence_layer')='PUBLIC_INDEPENDENT_SUPPORT'",(PUBLIC,))
        g.c.execute("UPDATE public_gene_support SET independence_class=?",(PUBLIC,))
        g.c.commit()
        for fn in [mgrn,myb127,emf2b,e2f]:fn(g)
    if args.base_only:
        (ROOT/'base_checkpoint.json').write_text(js({'stats':dict(g.stats),'skipped':dict(g.skipped),'sources':g.sources}))
        g.quarantine.close();g.c.close();print('BASE_COMPLETE_WAITING_FOR_CONTEXT',flush=True);return
    # Further importers are separate functions to keep raw-to-claim handling inspectable.
    from import_context import imprinting,single_nucleus
    imprinting(g);single_nucleus(g)
    with gzip.open(PRED,'rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            g.c.execute('INSERT INTO v3_model_predictions VALUES (?,?)',(row['gene_id'],js(row)))
    g.c.execute("INSERT OR REPLACE INTO metadata VALUES ('public_v3_release',?)",(js({'release_id':'MEKG-public-v3-20260907','model_prediction_source':str(PRED),'legacy_gene_cards':'historical v8 snapshot; use v3_model_predictions for current frozen scores','public_independence_policy':'not certified by used_for_training=0 alone'}),))
    g.quarantine.close();g.checkpoint('All sources imported')
    g.c.execute('PRAGMA journal_mode=DELETE');g.c.close()
    after={'parent_kg':sha(PARENT/'mekg.sqlite'),'frozen_model_predictions':sha(PRED)}
    assert before==after,'Unexpected parent source change'
    manifest={'release_id':'MEKG-public-v3-20260907','parent_release':str(PARENT),'new_release':str(release),
      'parent_and_model_hashes_before':before,'parent_and_model_hashes_after':after,'legacy_files_unchanged':True,
      'added_counts':dict(g.stats),'mapping_exclusions':dict(g.skipped),'new_sources':g.sources,
      'current_model_predictions':str(PRED),'live_website_changed':False,
      'public_evidence_policy':'All additions interpretation-only. Independence certification requires raw-sample/label/selection overlap audit.',
      'database_sha256':sha(release/'mekg.sqlite')}
    (release/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print('FINISHED',js(manifest['added_counts']),flush=True)

if __name__=='__main__':main()
