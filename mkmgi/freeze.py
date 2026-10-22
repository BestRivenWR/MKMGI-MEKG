"""Freeze the NEW workspace and packaged source code before fitting."""
import json,sys
from common import ROOT,CODE,config,sha,save,now,load

def main():
    if (ROOT/'protocol/freeze.json').exists() or (ROOT/'runs').exists():
        raise RuntimeError('Workspace is already frozen or fitting has started')
    m,raw,g=load()
    if raw.shape!=(len(m),92):raise ValueError('Invalid feature shape')
    paths=[ROOT/'config.json']
    for sub in ['data','reference','source_data','universe']:
        paths += [p for p in (ROOT/sub).rglob('*') if p.is_file()]
    sources=[p for p in CODE.rglob('*.py') if '__pycache__' not in p.parts]+[CODE/'config.json']
    save(ROOT/'protocol/freeze.json',dict(frozen=now(),release=config()['release'],
        release_hashes={str(p):sha(p) for p in paths},source_hashes={str(p):sha(p) for p in sources},
        registered_main=config()['registered_main'],python=sys.version,
        status='reproduction workspace freeze; historical result archive remains separate'))
    print('Frozen model workspace: '+str(ROOT))

if __name__=='__main__':main()
