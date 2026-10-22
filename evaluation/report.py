"""Summarize aligned OOF scores with the original tie-aware metrics/bootstrap."""
from pathlib import Path
import argparse,json,sys
import numpy as np,pandas as pd
from eval_support import source_common,SOURCE

def report(tables,output):
    stats=source_common.core_module('statistics');m,_,_=source_common.load()
    eligible=m.primary.to_numpy(bool)|~m.excluded_U.to_numpy(bool)
    registry=m.loc[eligible].reset_index(drop=True);scores={};rows=[]
    for name,path in tables.items():
        z=pd.read_csv(path,sep='\t').sort_values('row_index').reset_index(drop=True)
        if not np.array_equal(z.gene_id.to_numpy(str),m.gene_id.to_numpy(str)):raise ValueError('OOF alignment: '+name)
        s=z.score.to_numpy(float);scores[name]=s[eligible]
        for fold in range(5):
            take=eligible&m.outer80.eq(fold).to_numpy(bool)
            rows.append({'method':name,'fold':fold,**stats.metrics(m.primary.to_numpy(bool)[take],s[take])})
    output=Path(output)
    if output.exists():raise FileExistsError('Choose a new report directory')
    output.mkdir(parents=True);fd=pd.DataFrame(rows)
    fd.to_csv(output/'per_fold_metrics.tsv',sep='\t',index=False)
    columns=['AP','AUROC','precision_at_20','precision_at_50','precision_at_100']
    summary=fd.groupby('method')[columns].agg(['mean','std'])
    summary.to_csv(output/'summary.tsv',sep='\t')
    boot=stats.bootstrap(registry,scores,80);pairs=[]
    if 'MKMGI' in scores:
        for other in scores:
            if other=='MKMGI':continue
            for i,metric in enumerate(['AP','AUROC']):
                delta=boot['MKMGI']['all'][:,i]-boot[other]['all'][:,i]
                low,high=stats.interval(delta);al,ah=stats.interval(delta,len(scores)-1)
                means=fd.groupby('method')[metric].mean()
                pairs.append({'comparator':other,'metric':metric,'delta':float(means['MKMGI']-means[other]),
                    'conditional95_low':low,'conditional95_high':high,
                    'family_adjusted_low':al,'family_adjusted_high':ah,'family_size':len(scores)-1})
    pd.DataFrame(pairs).to_csv(output/'paired_comparisons.tsv',sep='\t',index=False)
    print(summary.to_string())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table',action='append',required=True,help='Name=/absolute/OOF.tsv.gz')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    tables={}
    for t in a.table:
        name,path=t.split('=',1)
        if name in tables:raise ValueError('Duplicate method name')
        tables[name]=Path(path)
    report(tables,a.output)
