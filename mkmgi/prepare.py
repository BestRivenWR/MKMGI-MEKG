"""Prepare one NEW model run from an already processed, versioned input bundle."""
from pathlib import Path
import argparse, json, os, shutil

def prepare(inputs,workspace):
    inputs=Path(inputs).resolve();workspace=Path(workspace).resolve()
    if workspace.exists():raise FileExistsError('Choose a new model workspace')
    if inputs==workspace or inputs in workspace.parents:raise ValueError('Separate input and output directories')
    required=['data/registry.tsv','data/features.npz','data/graph.npz','data/background_raw.npz',
        'data/annotation_tokens.npz','source_data/matched_candidate_pools.tsv',
        'reference/family80_bg_concat_OOF.tsv.gz','reference/family80_ExtraTrees_background_OOF.tsv.gz',
        'universe/graph.npz','universe/sealed/registry.tsv','universe/sealed/features.npz',
        'universe/source_data/universe_partition.tsv']
    absent=[n for n in required if not (inputs/n).is_file()]
    if absent:raise FileNotFoundError('Incomplete prepared model bundle: '+', '.join(absent))
    # Validate before creating output; retain the archived family/split definitions.
    import numpy as np,pandas as pd
    m=pd.read_csv(inputs/'data/registry.tsv',sep='\t')
    if not m.gene_id.is_unique or not np.array_equal(m.row_index,np.arange(len(m))):
        raise ValueError('Registry IDs or order are invalid')
    if m.groupby('group80').outer80.nunique().max()!=1:raise ValueError('Family crosses outer folds')
    for fold in range(5):
        p=m.loc[m.outer80.ne(fold)]
        if p.groupby('group80')['inner80_outer'+str(fold)].nunique().max()!=1:
            raise ValueError('Family crosses inner folds')
    with np.load(inputs/'data/features.npz',allow_pickle=False) as z:
        if z['main'].shape!=(len(m),92):raise ValueError('Expected 92 prepared main features')
    workspace.mkdir(parents=True)
    for rel in required:
        dest=workspace/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(inputs/rel,dest)
    shutil.copy2(Path(__file__).with_name('config.json'),workspace/'config.json')
    (workspace/'protocol').mkdir()
    return {'workspace':str(workspace),'genes':len(m),'input_files':len(required)}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',required=True,type=Path)
    p.add_argument('--workspace',required=True,type=Path);a=p.parse_args()
    print(json.dumps(prepare(a.inputs,a.workspace),indent=2))
