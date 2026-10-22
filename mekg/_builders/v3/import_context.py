"""Context evidence, kept distinct from genetic causality and spatial coordinates."""
from __future__ import annotations
import csv,gzip,json,pathlib,re
import numpy as np
import h5py
from scipy import sparse

from mekg_paths import workspace
ROOT=workspace()/'v3'

def imprinting(g):
    files=[f'B73_Group{i}.txt' for i in range(1,5)]
    source=g.source('imprinting_2025','B73-allele imprinting patterns in B73/W22 reciprocal endosperm','10.1093/g3journal/jkaf028',
      files+['timepoint_script_forGH.Rmd','imprinting_README.md','imprinting_github_tree.json'],['PRJNA1110585'],
      genotype='B73 x W22 and W22 x B73; B73 allele only',native_reference='B73 v5 plus W22 concatenated reference',
      github_commit='41ecb9e2e277c9bf3a0de10da96eac084bdb6df6',
      limitation='group lists do not encode MEG/PEG direction or per-timepoint RER; these fields remain unavailable; do not transfer to B73 selfed material')
    labels={1:'Imprinted and expressed at all four timepoints',2:'Imprinted whenever sufficiently expressed',
      3:'At least one imprinted timepoint with consistent parental bias',4:'Variable imprinting with higher expression threshold'}
    seen=set()
    for i in range(1,5):
        node=g.node('imprinting:v3:group'+str(i),'ImprintingPattern',labels[i],
          {'author_group':i,'DAP':[11,14,17,21],'parental_direction':None,'evidence_location':'Fig.4 and author Rmd lines 937-982'})
        genes=re.findall(r'Zm00001eb\d{6}',(ROOT/'raw'/f'B73_Group{i}.txt').read_text())
        for raw in genes:
            gene=g.map(raw,source,f'B73_Group{i}')
            if not gene:continue
            if gene in seen:raise ValueError('Imprinting groups overlap: '+gene)
            seen.add(gene)
            g.assertion(source,g.gene(gene),'HAS_ALLELE_IMPRINTING_PATTERN',node,f'B73_Group{i}:{gene}',
              'RECIPROCAL_CROSS_ALLELE_EXPRESSION',{'source_file':f'B73_Group{i}.txt','native_gene_id':raw,'author_group':i,
              'MEG_or_PEG':None,'allele':'B73','cross_directions':['B73 x W22','W22 x B73'],
              'replicates':'3 each except W22 x B73 at 21 DAP has 2','not_universal_selfed_B73_property':True},
              stage='11,14,17,21DAP',tissue='endosperm',genotype='B73/W22 reciprocal crosses',reference='B73 v5 + W22',
              experiment='PRJNA1110585_reciprocal_timecourse')
        g.stats['imprinting_group_'+str(i)]=len(genes)
    g.checkpoint('Imprinting imported')

def single_nucleus(g):
    files=[f'GSM835385{i+0}_Rep{i+1}_filtered_matrix.h5' for i in range(3)]
    existing=[f for f in files if (ROOT/'raw'/f).exists()]
    if len(existing)!=3:
        raise RuntimeError('All three snRNA replicates required before release: '+str(existing))
    source=g.source('snrna_19dap_2025','B73 19 DAP endosperm single-nucleus expression','10.1093/plcell/koaf162',
      files+['GSE270866.soft.gz'],['GSE270866','GSM8353850','GSM8353851','GSM8353852'],
      genotype='B73',native_reference='B73 v4 Cell Ranger reference',
      limitation='Cell Ranger filtered nuclei; author Seurat cell-type labels not included in this release. Nucleus detection is not spatial localization.',
      normalization='sum UMI per replicate / sum UMI over all native features * 1e6; detection=nonzero nuclei / all filtered nuclei',
      biological_replication='3 libraries; each uses pooled endosperms; nuclei are not independent biological replicates')
    for i,filename in enumerate(files,1):
        gsm='GSM835385'+str(i-1)
        with h5py.File(ROOT/'raw'/filename,'r') as f:
            mat=f['matrix'];shape=tuple(int(x) for x in mat['shape'][:]);
            x=sparse.csc_matrix((mat['data'][:],mat['indices'][:],mat['indptr'][:]),shape=shape)
            genes=[v.decode() for v in mat['features/id'][:]]
            total=float(x.sum());sums=np.asarray(x.sum(axis=1)).ravel();det=np.asarray((x>0).sum(axis=1)).ravel()
            assert len(genes)==shape[0] and shape[1]>0 and total>0
            condition=g.condition(gsm,'B73 endosperm 19 DAP, single-nucleus replicate '+str(i),
              accession=gsm,stage='19DAP',tissue='endosperm',nuclei=shape[1],native_features=shape[0],total_UMI=total,
              cell_type=None,spatial_coordinates=False)
            mapped=set()
            for j,native in enumerate(genes):
                gene=g.map(native,source,gsm+':'+native)
                if not gene:continue
                if gene in mapped:raise ValueError('Duplicate mapped snRNA feature: '+gene)
                mapped.add(gene)
                cpm=float(sums[j]/total*1e6);frac=float(det[j]/shape[1])
                attrs={'native_gene_id':native,'UMI_sum':int(sums[j]),'detected_nuclei':int(det[j]),'nuclei_total':shape[1],
                  'replicate':i,'detection_fraction':frac,'normalization_denominator_all_native_UMI':total,'cell_type':None,
                  'not_spatial_localization':True,'filtering':'author Cell Ranger filtered matrix, before Seurat QC'}
                ev=None
                if sums[j]>0:
                    ev=g.assertion(source,g.gene(gene),'DETECTED_IN_SINGLE_NUCLEUS_SAMPLE',condition,gsm+':'+native,
                      'SINGLE_NUCLEUS_EXPRESSION_CONTEXT',attrs,stage='19DAP',tissue='dissected endosperm',
                      genotype='B73',reference='B73 v4',experiment='GSE270866_19DAP')
                g.observation(source,gene,'SNRNA_PSEUDOBULK_CPM',condition,cpm,'CPM',gsm+':'+native,attrs,ev,genotype='B73',reference='B73 v4')
                g.observation(source,gene,'SNRNA_DETECTION_FRACTION',condition,frac,'fraction of filtered nuclei',gsm+':'+native,attrs,ev,genotype='B73',reference='B73 v4')
            g.stats['snrna_rep'+str(i)+'_filtered_nuclei']=shape[1]
            g.checkpoint('snRNA replicate '+str(i)+' imported')
