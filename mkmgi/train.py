"""Reuse the audited fitting algorithm, with explicit time masks and truthful artifacts."""
import argparse, json
import numpy as np
import torch
from common import ROOT,config,reference_module,load,target,save,now

_base=reference_module('train')
_original_context=_base.context
_original_fit=_base.neural_fit
directory=_base.directory
report_progress=_base.report_progress

def context(raw,graph,train_mask,infer_mask,loss_mask,condition,seed,key,device,bgraw):
    ctx=_original_context(raw,graph,train_mask,infer_mask,loss_mask,condition,seed,key,device,bgraw)
    # Unlike the flat baseline, token masks must be real observed 0/1, not z-scored flags.
    masks=torch.as_tensor(raw[:,40:80],device=device)
    ctx['xt'][:,40:80]=masks;ctx['xi'][:,40:80]=masks
    return ctx

def neural_fit(*args,**kwargs):
    result=_original_fit(*args,**kwargs)
    if 'background_gates' in result:
        result['temporal_attention']=result['background_gates']
    return result

_base.context=context
_base.neural_fit=neural_fit

def correct_checkpoint(path):
    ck=torch.load(path,map_location='cpu')
    ck.update(architecture='TimeTransformer32_RelGCN64x2',raw_input_dim=191,fused_input_dim=143,
              time_tokens=8,time_width=32,time_heads=4,DAP=config()['DAP'],
              time_pooling='mean' if ck['condition']=='time_mean' else 'learned_attention',
              explicit_DAP_encoding=ck['condition']!='time_no_pos')
    torch.save(ck,target(path))

def train_job(grouping,fold,condition):
    _base.train_job(grouping,fold,condition)
    path=directory(grouping,fold,condition);m,raw,_=load();test=m.outer80.to_numpy(int)==fold
    observed=(raw[test,40:80].reshape(-1,5,8)>.5).any(1)
    for seed in config()['seeds']:
        p=path/('seed'+str(seed));correct_checkpoint(p/'model.pt')
        with np.load(p/'outer_gates.npz') as z:
            ids=z['gene_id'].copy();relation=z['relation_gates'].copy();attention=z['background_gates'].copy() if 'background_gates' in z.files else z['temporal_pooling_weights'].copy()
        assert attention.shape==(int(test.sum()),8) and not attention[~observed].any()
        np.savez_compressed(target(p/'outer_gates.npz'),gene_id=ids,relation_gates=relation,
            temporal_pooling_weights=attention,temporal_observed=observed,DAP=np.array(config()['DAP']))
    rows=json.loads((path/'gate_summary.json').read_text())
    for row in rows:
        if 'mean_outer_background_gates' in row:
            row['mean_outer_temporal_pooling_weights']=row.pop('mean_outer_background_gates')
    save(path/'gate_summary.json',rows)
    done=json.loads((path/'complete.json').read_text())
    done.update(architecture='TimeTransformer32_RelGCN64x2',attention_weights='temporal pooling, not causal or per-omics attribution',artifact_schema_checked=now())
    save(path/'complete.json',done)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--grouping',type=int,choices=[80],required=True)
    p.add_argument('--fold',type=int,choices=list(range(5)),required=True)
    p.add_argument('--condition',choices=config()['conditions80'],required=True)
    a=p.parse_args();train_job(a.grouping,a.fold,a.condition)
