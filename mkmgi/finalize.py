"""Reuse DEV-only refit algorithm; export truthful time-model metadata and weights."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from common import ROOT,config,reference_module,load,target,save,verify_freeze
from train import correct_checkpoint

def main():
    reference_module('finalize').main()
    c=config();path=ROOT/'models/deployment'/c['registered_main']
    for seed in c['seeds']:
        correct_checkpoint(path/('seed'+str(seed))/'model.pt')
    original=ROOT/'predictions/genomewide_background_gcn_consensus.tsv.gz'
    original.replace(target(ROOT/'predictions/genomewide_time_gcn_consensus.tsv.gz'))
    p=ROOT/'predictions/genomewide_mean_gates.npz'
    with np.load(p) as z:
        ids=z['gene_id'].copy();rel=z['relation_gates'].copy();time=z['background_gates'].copy()
    dev,raw,_=load();observed=np.zeros((len(ids),8),bool)
    observed[dev.gene_index.to_numpy(int)]=(raw[:,40:80].reshape(-1,5,8)>.5).any(1)
    old=Path(c['previous_release']);other=pd.read_csv(old/'sealed/registry.tsv',sep='\t',usecols=['gene_id','gene_index'])
    with np.load(old/'sealed/features.npz') as z:
        flags=z['main'][:,40:80]
    assert np.array_equal(ids[other.gene_index.to_numpy(int)],other.gene_id.to_numpy(str))
    observed[other.gene_index.to_numpy(int)]=(flags.reshape(-1,5,8)>.5).any(1)
    assert not time[~observed].any()
    np.savez_compressed(target(p),gene_id=ids,relation_gates=rel,temporal_pooling_weights=time,
                        temporal_observed=observed,DAP=np.array(c['DAP']))
    manifest=json.loads((ROOT/'models/deployment/manifest.json').read_text())
    manifest.update(architecture='TimeTransformer32_RelGCN64x2',time_width=32,time_heads=4,time_layers=1,DAP=c['DAP'],
        temporal_weights='per-gene temporal pooling weights; not causal/per-omics attribution',
        prediction_file='predictions/genomewide_time_gcn_consensus.tsv.gz')
    save(ROOT/'models/deployment/manifest.json',manifest);verify_freeze(check_sources=True)
    print('Time-model checkpoints, genomewide predictions and 8-DAP weights exported; no website changes',flush=True)

if __name__=='__main__':main()
