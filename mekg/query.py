"""Query evidence from an existing MEKG v3.1 SQLite snapshot."""
from pathlib import Path
import argparse,json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'_builders/v31'))
from query_v31 import query

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gene',required=True)
    p.add_argument('--db',required=True,type=Path);p.add_argument('--limit',type=int,default=30)
    a=p.parse_args();print(json.dumps(query(a.gene,a.limit,a.db),ensure_ascii=False,indent=2,allow_nan=False))
