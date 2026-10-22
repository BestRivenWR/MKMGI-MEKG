"""Small dispatcher for the source-preserving MEKG public-evidence build stages."""
from pathlib import Path
import argparse,os,runpy,sys

HERE=Path(__file__).resolve().parent
STAGES={'public-v3':('v3','build_v3.py',[]),
    'v3-semantics':('v3','finalize_semantics.py',[]),
    'v3-audit':('v3','audit_v3.py',[]),
    'tf-peaks':('v31','build_v31.py',['tf_peaks']),
    'myb-peaks':('v31','build_v31.py',['myb_peaks']),
    'physical-contacts':('v31','build_v31.py',['hichip']),
    'gwas':('v31','build_v31.py',['gwas']),
    'cells':('v31','import_cells.py',['cells']),
    'mutant-rna':('v31','import_cells.py',['kil_rna']),
    'tss':('v31','import_tss.py',[]),
    'tf-mapping-audit':('v31','audit_tf_mapping.py',[]),
    'finalize':('v31','finalize_v31.py',[]),
    'audit':('v31','audit_release.py',[])}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=STAGES)
    p.add_argument('--workspace',required=True,type=Path);a,extra=p.parse_known_args()
    root=a.workspace.resolve();os.environ['MEKG_WORKSPACE']=str(root)
    for sub in ['v3','v31','v3/raw','v31/raw']: (root/sub).mkdir(parents=True,exist_ok=True)
    folder,script,args=STAGES[a.stage]
    code=HERE/'_builders'/folder/script
    sys.path[:0]=[str(HERE),str(code.parent),str(HERE/'_builders/v3')]
    sys.argv=[str(code),*args,*extra]
    runpy.run_path(str(code),run_name='__main__')
