"""Append audited evidence to a separate v3.1 copy. Model rows stay immutable."""
from pathlib import Path
import sys,json,csv,gzip,sqlite3,shutil,collections,io,zipfile,re,time,hashlib,argparse
from mekg_paths import workspace
ROOT=workspace()/'v31';PARENT=workspace()/'v3'
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'v3'))
import build_v3 as old
from table_io import sheets
from coordinate_tools import StrictChain,IntervalIndex,chrom,self_test
old.ROOT=ROOT
js,uid,sha,num=old.js,old.uid,old.sha,old.num

class Graph31(old.Graph):
    def __init__(self):
        target=ROOT/'release/mekg.sqlite'
        assert not (ROOT/'release/COMPLETE.json').exists(),'Completed release is immutable'
        if not target.exists():
            expected=json.loads((PARENT/'release/QA.json').read_text())['database_sha256']
            assert sha(PARENT/'release/mekg.sqlite')==expected
            shutil.copy2(PARENT/'release/mekg.sqlite',target)
            (ROOT/'parent_preservation.json').write_text(js({'parent':str(PARENT),'sha256':expected,'model_file':str(old.PRED),'model_sha256':sha(old.PRED)}))
        self.c=sqlite3.connect(target);self.c.execute('PRAGMA journal_mode=MEMORY');self.c.execute('PRAGMA synchronous=NORMAL')
        self.c.execute('PRAGMA cache_size=-2097152');self.c.execute('PRAGMA temp_store=MEMORY')
        self.known={r[0] for r in self.c.execute('SELECT node_id FROM nodes')}
        self.model_genes={r[0] for r in self.c.execute('SELECT gene_id FROM v3_model_predictions')}
        self.metadata={};self.mapping={};self.skipped=collections.Counter();self.stats=collections.Counter()
        self.sources={r[0]:json.loads(r[1]) for r in self.c.execute('SELECT source_id,source_json FROM v3_sources')}
        self.new_edges=set();self.new_evidence=set();self.new_observations=set()
        self.quarantine=gzip.open(ROOT/'mapping_quarantine.jsonl.gz','at',encoding='utf8')
        with gzip.open(old.OLDRAW/'Zm00001eb.1.fulldata.txt.gz','rt') as f:
            for row in csv.reader(f,delimiter='\t'):
                if len(row)>=14 and re.fullmatch(r'Zm00001eb\d{6}',row[1]):self.metadata[row[1]]=row
        for version in [3,4]:
            raw=collections.defaultdict(set);rev=collections.defaultdict(set)
            with open(old.OLDRAW/f'B73v{version}_to_B73v5.tsv') as f:
                for r in csv.reader(f,delimiter='\t'):
                    if len(r)<2:continue
                    for dest in r[1].split(','):
                        if dest in self.metadata:raw[r[0]].add(dest);rev[dest].add(r[0])
            for k,v in raw.items():
                if len(v)==1 and len(rev[next(iter(v))])==1:self.mapping[k]=next(iter(v))
        self.alias=collections.defaultdict(set)
        for a,g in self.c.execute('SELECT alias,gene_id FROM public_aliases'):
            if g in self.metadata:self.alias[a.lower()].add(g)
        for g,r in self.metadata.items():
            if r[10] not in ['','-']:self.alias[r[10].lower()].add(g)
        self.c.executescript('''
          CREATE TABLE IF NOT EXISTS v31_modules(module TEXT PRIMARY KEY, report_json TEXT);
          CREATE TABLE IF NOT EXISTS v31_peaks(peak_id TEXT PRIMARY KEY,source_id TEXT,tf_gene_id TEXT,target_gene_id TEXT,
            source_file TEXT,source_row INTEGER,native_chrom TEXT,native_start0 INTEGER,native_end0 INTEGER,
            v5_chrom TEXT,v5_start0 INTEGER,v5_end0 INTEGER,mapping_status TEXT,annotation TEXT,distance_to_tss REAL,score REAL);
          CREATE INDEX IF NOT EXISTS v31_peak_target ON v31_peaks(target_gene_id,source_id);
          CREATE INDEX IF NOT EXISTS v31_peak_tf ON v31_peaks(tf_gene_id,source_id);
          CREATE TABLE IF NOT EXISTS v31_contacts(contact_id TEXT PRIMARY KEY,source_file TEXT,source_row INTEGER,
            histone_mark TEXT,native_json TEXT,v5_anchor1_json TEXT,v5_anchor2_json TEXT,attributes_json TEXT);
          CREATE TABLE IF NOT EXISTS v31_gwas_loci(locus_id TEXT PRIMARY KEY,chromosome TEXT,position INTEGER,
            trait TEXT,year TEXT,model TEXT,pvalue REAL,source_row INTEGER,source_id TEXT);
          CREATE INDEX IF NOT EXISTS v31_gwas_pos ON v31_gwas_loci(chromosome,position);
        ''')
        self.chain=StrictChain(ROOT/'raw/v4_to_v5.chain');assert self_test()
        assert hashlib.md5((ROOT/'raw/v4_to_v5.chain').read_bytes()).hexdigest()=='492eb7f383b66f2cca2c9349be980eaa'
    def resolve_tf(self,name,source):
        native=name.removesuffix('_filter_after.bed.txt')
        if re.fullmatch(r'1d\d{6}',native):native='Zm0000'+native
        if native in self.mapping:return self.mapping[native]
        if native in self.metadata:return native
        options=self.alias.get(native.lower(),set())
        if len(options)==1:return next(iter(options))
        self.skipped['unresolved_tf_profile']+=1
        self.quarantine.write(js({'source':source,'tf_profile':name,'reason':'unknown_or_ambiguous_tf_identity','candidates':sorted(options)})+'\n')
        return None
    def finish_module(self,name,extra=None):
        d={'stats':dict(self.stats),'skipped':dict(self.skipped),**(extra or {})}
        self.c.execute('INSERT INTO v31_modules VALUES (?,?)',(name,js(d)));self.c.commit()
        print('MODULE_DONE',name,js(d),flush=True)
        self.quarantine.close();self.c.execute('PRAGMA journal_mode=DELETE');self.c.close()

def tf_peaks(g):
    source=g.source('tf_peak_catalog_2025','513 TF author peak-to-gene annotations','10.1016/j.molp.2025.08.009',
      ['513TFs_ChIPSeeker_regulatory_filtered.zip','tf_metadata.html','v4_to_v5.chain'],['CRA022036','GSE201701'],
      native_reference='B73 v4',assay='mixed TF binding profiles; per-profile assay/sample assignment not certified',
      limitation='author ChIPseeker gene association is not proven regulation; same primary study as mGRN; do not count as independent corroboration',
      promoter_definition='author +/-5kb TSS; retain author annotation; no expansion from flank_geneIds',
      coordinate_policy='ChIPseeker 1-based closed to BED; unique ungapped whole-interval reciprocal chain mapping')
    count=0;files=[]
    with zipfile.ZipFile(ROOT/'raw/513TFs_ChIPSeeker_regulatory_filtered.zip') as z:
        names=[n for n in z.namelist() if not n.endswith('/')];assert len(names)==513
        for i,name in enumerate(names,1):
            tf=g.resolve_tf(name,source);nrows=0;nmap=0
            with z.open(name) as raw:
                for line,row in enumerate(csv.DictReader(io.TextIOWrapper(raw,encoding='utf8'),delimiter='\t'),2):
                    nrows+=1;count+=1
                    c=chrom(row['seqnames']);s=int(row['start'])-1;e=int(row['end'])
                    assert e-s==int(row['width'])
                    target=g.mapping.get(row['geneId'])
                    if target is None:g.skipped['tf_peak_target_nonbijective']+=1
                    got=g.chain.map(c,s,e);nmap+=int(got is not None)
                    loc=got[:3] if got else (None,None,None)
                    pk=uid('peak:v31',source,name,c,s,e,row['geneId'])
                    g.c.execute('INSERT OR IGNORE INTO v31_peaks VALUES ('+','.join('?'*16)+')',
                      (pk,source,tf,target,name,line,c,s,e,*loc,'unique_ungapped_reciprocal' if got else 'native_only_unmapped_or_gap',row['annotation'],num(row['distanceToTSS']),num(row.get('V5'))))
            files.append({'profile':name,'tf_gene_id':tf,'rows':nrows,'v5_mapped_rows':nmap})
            if i%50==0:print('TF_PROFILES',i,'raw_rows',count,flush=True)
    # Collapsing multiple peaks prevents a peak-rich pair from looking like many independent experiments.
    pairs=g.c.execute('''SELECT tf_gene_id,target_gene_id,COUNT(*),SUM(v5_chrom IS NOT NULL)
       FROM v31_peaks WHERE source_id=? AND tf_gene_id IS NOT NULL AND target_gene_id IS NOT NULL GROUP BY tf_gene_id,target_gene_id''',(source,)).fetchall()
    print('TF_GENE_PAIRS_TO_IMPORT',len(pairs),flush=True)
    for counter,(tf,target,n,nlift) in enumerate(pairs,1):
        g.assertion(source,g.gene(tf),'HAS_BINDING_PEAK_ASSOCIATED_GENE',g.gene(target),tf+'|'+target,
          'TF_BINDING_PEAK_GENE_ASSOCIATION',{'peak_count':n,'v5_mapped_peak_count':nlift,
          'peak_table':'v31_peaks','association_is_geometric':True,'direct_in_vivo_regulation_proven':False,
          'coordinate_mapping':'gene identity plus strict reciprocal ungapped v4-to-v5 peak projection when available',
          'assay_assignment':'per-profile assay/sample assignment pending','same_primary_study_as':'mgrn_2025'},
          reference='B73 v4 native peaks; B73 v5 coordinates only where mapping passed',
          tissue='assay-specific; not assumed endosperm',experiment='mGRN2025_integrated_TFBS_context')
        if counter%100000==0:print('TF_PAIRS_IMPORTED',counter,flush=True)
    (ROOT/'tf_profile_inventory.json').write_text(js(files))
    g.finish_module('tf_peaks',{'profiles':len(files),'resolved_tf_profiles':sum(x['tf_gene_id'] is not None for x in files),'input_rows':count,'mapped_coordinate_rows':sum(x['v5_mapped_rows'] for x in files)})

def gene_indexes(g):
    gtf=g.c.execute("SELECT value FROM metadata WHERE key='b73v5_gtf'").fetchone()[0];path=Path(json.loads(gtf));strands={}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):continue
            x=line.rstrip().split('\t')
            if len(x)!=9:continue
            mt=re.search(r'gene_id "(Zm00001eb\d{6})"',x[8])
            if mt:strands[mt.group(1)]=x[6]
    body=[];prom=[]
    for gene,row in g.metadata.items():
        c=chrom(row[7]);s=int(row[8])-1;e=int(row[9]);body.append((c,s,e,gene))
        st=strands.get(gene)
        if st=='+':prom.append((c,max(0,s-2000),s+1,gene))
        elif st=='-':prom.append((c,e-1,e+2000,gene))
    return IntervalIndex(body),IntervalIndex(prom)

def myb_peaks(g):
    file=PARENT/'raw/myb127_binding.xlsx'
    source=g.source('myb127_peaks_2026','MYB127 replicated CUT&Tag peak context','10.1038/s41477-026-02238-3',
      [str(file),'v4_to_v5.chain'],['CRA028659'],genotype='MYB127-FLAG overexpression background per paper; not endogenous B73 occupancy',
      tissue='endosperm',stage='16DAP',native_reference='B73 v4',
      limitation='replicate-overlap and WT-overlap filters are contextual QC, not differential binding statistics; no input-controlled BAM reanalysis')
    tf=g.map('Zm00001d041935',source,'TF');data={}
    for name,rows in sheets(file):
        rs=[]
        for n,row in enumerate(rows[10:],11):
            if len(row)<9:continue
            c=chrom(row[0]);s=int(row[1])-1;e=int(row[2]);assert e-s==int(row[3])
            rs.append((c,s,e,{'row':n,'name':row[8],'fold_enrichment':num(row[7]),'pileup':num(row[5])}))
        data[name]=rs;g.stats['input:'+name]=len(rs)
    rep1=data['ZmMYB127-Flag_Rep1_peaks'];rep2=IntervalIndex(data['ZmMYB127-Flag_Rep2_peaks'])
    wt=IntervalIndex(data['WT_Rep1_peaks']+data['WT_Rep2_peaks']);body,prom=gene_indexes(g);bygene=collections.defaultdict(list)
    candidates={}
    for c,s,e,rec in rep1:
        for a,b,r2 in rep2.overlap(c,s,e):
            overlap=min(e,b)-max(s,a)
            if overlap/(e-s)<0.5 or overlap/(b-a)<0.5:continue
            g.stats['replicate_overlap_pairs']+=1
            # Require no overlap of either original peak with any WT-called peak.
            if wt.overlap(c,min(s,a),max(e,b)):g.skipped['MYB_WT_peak_overlap']+=1;continue
            lo,hi=max(s,a),min(e,b)
            candidates[(c,lo,hi)]= {'rep1':rec,'rep2':r2,'rule':'reciprocal overlap >=0.5; neither peak overlaps WT peak union'}
    for (c,s,e),attrs in sorted(candidates.items()):
        got=g.chain.map(c,s,e)
        if not got:g.skipped['MYB_consensus_liftover_failed']+=1;continue
        nc,ns,ne=got[:3];pk=uid('peak:v31',source,c,s,e)
        targets=collections.defaultdict(set)
        for a,b,gene in body.overlap(nc,ns,ne):targets[gene].add('gene_body')
        for a,b,gene in prom.overlap(nc,ns,ne):targets[gene].add('promoter_2kb_upstream_TSS')
        for target,kinds in targets.items():
            peakid=uid('peak:v31',pk,target);bygene[target].append(peakid)
            g.c.execute('INSERT INTO v31_peaks VALUES ('+','.join('?'*16)+')',(peakid,source,tf,target,'myb127_binding.xlsx',attrs['rep1']['row'],c,s,e,nc,ns,ne,'unique_ungapped_reciprocal',';'.join(sorted(kinds)),None,None))
        g.stats['MYB_replicated_noWT_mapped_peaks']+=1
    for gene,peaks in bygene.items():
        g.assertion(source,g.gene(tf),'HAS_REPLICATED_CUTTAG_PEAK_ASSOCIATED_GENE',g.gene(gene),gene,
          'REPLICATED_CUTTAG_PEAK_GENE_ASSOCIATION',{'peak_count':len(peaks),'peak_table':'v31_peaks',
          'rule':'reciprocal replicate overlap >=0.5; no WT called-peak overlap; strict v4-v5 mapping; gene body or 2kb upstream TSS',
          'association_is_geometric':True,'functional_regulation_proven':False,'same_primary_study_as':'myb127_2026',
          'coordinate_mapping':'strict reciprocal ungapped v4-to-v5 peak projection'},
          stage='16DAP',tissue='endosperm',genotype='MYB127-FLAG overexpression; see primary study',reference='B73 v5 projected peaks from B73 v4',experiment='MYB127_CUTTag_16DAP')
    g.finish_module('myb_peaks',{'peak_associated_genes':len(bygene)})

def hichip(g):
    source=g.source('hichip_2025','B73 seedling H3K4me3/H3K27ac HiChIP contacts','10.1038/s41467-025-65065-8',
      [f'hic_{n}.xlsx' for n in [4,5,8,9,10,11]]+['hic.html','v4_to_v5.chain'],['PRJNA825068'],
      genotype='B73 (Mo17 comparison retained in original files)',tissue='seedlings',native_reference='B73 AGPv4',
      limitation='cross-tissue physical-contact context, not measured endosperm regulation; gene-pair tables and loop tables are overlapping representations')
    for n,mark in [(4,'H3K4me3'),(5,'H3K27ac'),(8,'H3K4me3'),(9,'H3K27ac')]:
        rows=next(sheets(ROOT/'raw'/f'hic_{n}.xlsx'))[1]
        for line,row in enumerate(rows[2:],3):
            if len(row)<6:continue
            if n in [8,9] and row[10]!='B73-ref':continue
            if n in [4,5] and not (num(row[13]) is not None and num(row[15]) is not None and float(row[13])<0.01 and float(row[15])<0.01):
                g.skipped['hichip_specific_not_significant_both_replicates']+=1;continue
            c1,s1,e1,c2,s2,e2=chrom(row[0]),int(row[1]),int(row[2]),chrom(row[3]),int(row[4]),int(row[5])
            # Uppercase Chr is B73 and lowercase chr is Mo17 in the combined genome tables.
            if n in [4,5] and (not row[0].startswith('Chr') or not row[3].startswith('Chr')):continue
            a=g.chain.map(c1,s1,e1);b=g.chain.map(c2,s2,e2)
            attrs={'native_coordinate_convention':'0-based half-open FitHiChIP bins','cross_tissue':True,
              'author_class':'B73-specific' if n in [4,5] else 'conserved','original_row':row,
              'replicate_filter':'both B73 Q<0.01' if n in [4,5] else 'author conserved call; replicate counts retained'}
            g.c.execute('INSERT INTO v31_contacts VALUES (?,?,?,?,?,?,?,?)',(uid('contact:v31',n,line),f'hic_{n}.xlsx',line,mark,js([c1,s1,e1,c2,s2,e2]),js(a[:3]) if a else None,js(b[:3]) if b else None,js(attrs)))
            g.stats['contact_records']+=1;g.stats['both_contact_anchors_mapped']+=int(a is not None and b is not None)
    for n,mark in [(10,'H3K4me3'),(11,'H3K27ac')]:
        rows=next(sheets(ROOT/'raw'/f'hic_{n}.xlsx'))[1];seen=set()
        for line,row in enumerate(rows[2:],3):
            a=g.map(row[0],source,f'{n}:{line}:a');b=g.map(row[1],source,f'{n}:{line}:b')
            if not a or not b:continue
            pair=tuple(sorted((a,b)))
            if pair in seen:continue
            seen.add(pair)
            g.assertion(source,g.gene(pair[0]),'PUBLIC_HICHIP_ASSOCIATED_GENE_PAIR',g.gene(pair[1]),mark+'|'+'|'.join(pair),
              'HICHIP_PHYSICAL_CONTACT_CONTEXT',{'mark':mark,'native_gene_ids':row[:2],
              'native_B73_gene_distance':num(row[4]),'cross_tissue':True,'not_endosperm_measured':True,
              'functional_regulation_proven':False,'source_file':f'hic_{n}.xlsx','source_row':line,
              'source_scope':'author syn-syn interaction-associated gene pairs; not full interactome'},
              tissue='seedlings',genotype='B73',reference='B73 v4 gene identity mapped to B73 v5',experiment='B73_HiChIP_2025',directed=0)
    g.finish_module('hichip')

def gwas(g):
    source=g.source('kernel_gwas_2026','Kernel-trait GWAS in 187 maize accessions','10.1186/s12870-026-08285-4',
      ['gwas_3.xlsx','gwas.html'],[],genotype='187 diverse maize accessions',native_reference='B73 RefGen v5 per author Methods',
      limitation='association not causation; two models and correlated traits on the same panel are not independent replication; field biological replication limited',
      source_scope='Table S3 all author reported loci; gene proximity within +/-50kb independently indexed')
    traits=['Kernel Length','Kernel Width','Kernel Thickness','100-Kernel Weight','Kernel Volume'];rows=None
    for name,rs in sheets(ROOT/'raw/gwas_3.xlsx'):
        if name=='Table S3':rows=rs;break
    assert rows is not None and rows[4][:3]==['Chromosome','Position','P-value']
    body,_=gene_indexes(g);hits=collections.defaultdict(list);seen=set()
    for line,row in enumerate(rows[5:],6):
        for block in range(20):
            off=block*3
            if off+2>=len(row):continue
            fields=[x.strip() for x in row[off:off+3]]
            missing={'','/','NA','-'}
            if all(x in missing for x in fields):
                g.skipped['GWAS_empty_or_slash_placeholder_blocks']+=1;continue
            assert all(x not in missing for x in fields),('partially missing GWAS record',line,block,fields)
            c=chrom(row[off]);pos=int(float(row[off+1]));p=num(row[off+2])
            assert c in [str(x) for x in range(1,11)] and pos>0 and p is not None and 0<p<=1
            model='MLM_PCA+K' if block<10 else 'MLM_Q+K';trait=traits[(block%10)//2];year='2020' if block%2==0 else '2024'
            key=(c,pos,trait,year,model)
            if key in seen:continue
            seen.add(key);loc=uid('locus:v31',source,*key)
            g.c.execute('INSERT INTO v31_gwas_loci VALUES (?,?,?,?,?,?,?,?,?)',(loc,c,pos,trait,year,model,p,line,source))
            if p>10**(-6.09):g.skipped['GWAS_author_reported_below_nominal_threshold']+=1;continue
            # Physical proximity is a separate inference; retain every eligible gene, not just hand-picked candidates.
            for s,e,gene in body.overlap(c,max(0,pos-1-50000),pos+50000):
                hits[(gene,trait)].append({'locus_id':loc,'chromosome':c,'position':pos,'pvalue':p,'model':model,'year':year})
    for (gene,trait),records in hits.items():
        node=g.node('trait:v31:'+trait.replace(' ','_'),'QuantitativeTrait',trait,{'not_endosperm_specific':True})
        g.assertion(source,g.gene(gene),'NEAR_KERNEL_TRAIT_GWAS_LOCUS',node,gene+'|'+trait,
          'GWAS_LOCUS_PROXIMITY',{'locus_results':records,'window_bp':50000,'window_policy':'gene body overlaps SNP +/-50kb',
          'unique_locus_count':len({(x['chromosome'],x['position']) for x in records}),
          'independent_study_count':1,'no_combined_pvalue':True,'causal_gene_proven':False,
          'method_repeats_are_not_biological_replicates':True},
          tissue='mature kernels',genotype='187 accession association panel',reference='B73 v5 per source',experiment='kernel_GWAS_187_accessions_2026')
    g.finish_module('gwas',{'locus_model_trait_year_records':len(seen),'gene_trait_pairs':len(hits)})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('module',choices=['tf_peaks','myb_peaks','hichip','gwas']);args=ap.parse_args()
    g=Graph31()
    if g.c.execute('SELECT 1 FROM v31_modules WHERE module=?',(args.module,)).fetchone():
        print('Already completed; no mutation',args.module);g.c.close();g.quarantine.close()
    else:globals()[args.module](g)
