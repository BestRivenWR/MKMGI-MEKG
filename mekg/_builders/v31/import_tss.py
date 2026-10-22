"""Experimental RAMPAGE context; not a replacement for canonical TSS annotation."""
import csv,gzip,collections,json
from build_v31 import Graph31,ROOT,js,num,uid,gene_indexes
from coordinate_tools import IntervalIndex,chrom

def main():
    g=Graph31()
    assert not g.c.execute("SELECT 1 FROM v31_modules WHERE module='tss'").fetchone()
    files=[f'GSM804317{5+i}_B73_endosperm_RAMPAGE_Rep{i}_peaks.narrowPeak.gz' for i in [1,2]]
    source=g.source('maizecode_rampage_2024','MaizeCODE B73 endosperm 15 DAP RAMPAGE TSS context',
      '10.1038/s41467-024-55195-w',files+['maizecode_geo.soft.txt','maizecode.html'],
      ['GSE254487','GSM8043176','GSM8043177','PRJNA438108'],genotype='B73 WT',tissue='endosperm',stage='15DAP',
      native_reference='B73 NAM (v5)',
      assay='RAMPAGE capped-transcript peaks called with matched RNA-seq control per GEO',
      limitation='Unstranded narrowPeak files; peak proximity cannot establish transcription direction, causal regulation or uniquely assign a promoter; no canonical gene annotation replaced',
      processing='reciprocal replicate overlap >=0.5 and summit distance <=100 bp; associate consensus intersection with gene body or 2kb upstream annotated TSS')
    g.c.executescript('''CREATE TABLE v31_tss_peaks(peak_id TEXT PRIMARY KEY,replicate INTEGER,chromosome TEXT,start0 INTEGER,end0 INTEGER,
      summit0 INTEGER,signal REAL,minus_log10_p REAL,minus_log10_q REAL,source_row INTEGER,source_file TEXT);
      CREATE TABLE v31_tss_consensus(consensus_id TEXT PRIMARY KEY,chromosome TEXT,start0 INTEGER,end0 INTEGER,peak1_id TEXT,peak2_id TEXT);
      CREATE TABLE v31_tss_genes(gene_id TEXT,consensus_id TEXT,association TEXT,PRIMARY KEY(gene_id,consensus_id)) WITHOUT ROWID;
    ''')
    data={}
    for rep,file in enumerate(files,1):
      items=[]
      with gzip.open(ROOT/'raw'/file,'rt') as f:
        for line,row in enumerate(csv.reader(f,delimiter='\t'),1):
          assert len(row)==10 and row[5]=='.'
          c=chrom(row[0]);s,e,offset=int(row[1]),int(row[2]),int(row[9]);assert 0<=s<e and 0<=offset<e-s
          pk=uid('rampage:v31',file,row[3]);summit=s+offset
          g.c.execute('INSERT INTO v31_tss_peaks VALUES (?,?,?,?,?,?,?,?,?,?,?)',
            (pk,rep,c,s,e,summit,num(row[6]),num(row[7]),num(row[8]),line,file))
          items.append((c,s,e,(pk,summit)))
      data[rep]=items;g.stats['rep'+str(rep)+'_input_peaks']=len(items)
    idx=IntervalIndex(data[2]);body,prom=gene_indexes(g);bygene=collections.defaultdict(list)
    for c,s,e,(p1,sum1) in data[1]:
      for a,b,(p2,sum2) in idx.overlap(c,s,e):
        lo,hi=max(s,a),min(e,b);over=hi-lo
        if over/(e-s)<0.5 or over/(b-a)<0.5 or abs(sum1-sum2)>100:continue
        cid=uid('rampage_consensus:v31',p1,p2)
        g.c.execute('INSERT INTO v31_tss_consensus VALUES (?,?,?,?,?,?)',(cid,c,lo,hi,p1,p2))
        targets=collections.defaultdict(set)
        for a,b,gene in body.overlap(c,lo,hi):targets[gene].add('gene_body')
        for a,b,gene in prom.overlap(c,lo,hi):targets[gene].add('2kb_upstream_annotated_TSS')
        for gene,kinds in targets.items():
          g.c.execute('INSERT INTO v31_tss_genes VALUES (?,?,?)',(gene,cid,';'.join(sorted(kinds))))
          bygene[gene].append(cid)
        g.stats['replicate_consensus_pairs']+=1
    cond=g.condition('MaizeCODE_RAMPAGE_15DAP','B73 endosperm 15 DAP RAMPAGE, 2 replicates',
      genotype='B73',stage='15DAP',assay='RAMPAGE',replicates=2,native_reference='B73 NAM v5')
    for gene,pks in bygene.items():
      g.assertion(source,g.gene(gene),'HAS_REPLICATED_RAMPAGE_PEAK_ASSOCIATION',cond,gene,
        'RAMPAGE_TRANSCRIPTION_INITIATION_CONTEXT',{'consensus_pair_count':len(pks),'peak_table':'v31_tss_consensus',
          'association_table':'v31_tss_genes','coordinate_mapping':'native B73 NAM v5, no liftOver',
          'strand_not_available_in_deposited_peak_file':True,'canonical_TSS_replaced':False,'causal_regulation_proven':False},
          genotype='B73',stage='15DAP',tissue='endosperm',reference='B73 NAM v5',experiment='MaizeCODE_B73_15DAP_RAMPAGE')
    g.finish_module('tss',{'peak_associated_genes':len(bygene)})

if __name__=='__main__':main()
