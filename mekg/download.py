"""Fetch public MEKG inputs into this new work directory only (stdlib)."""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, pathlib, time, urllib.request

from mekg_paths import workspace
import os
ROOT=workspace()/os.environ.get('MEKG_DOWNLOAD_STAGE','v3')
ROOT.mkdir(parents=True,exist_ok=True)

def fetch(item):
    target = (ROOT / 'raw' / item['filename']).resolve()
    assert ROOT in target.parents
    target.parent.mkdir(parents=True, exist_ok=True)
    record = dict(item)
    try:
        if not target.exists():
            req = urllib.request.Request(item['url'], headers={'User-Agent':'MEKG-public-data-curation/3.0'})
            with urllib.request.urlopen(req, timeout=90) as response, target.with_suffix(target.suffix+'.part').open('wb') as out:
                expected=response.headers.get('Content-Length')
                size = 0
                while True:
                    block = response.read(1024*1024)
                    if not block: break
                    size += len(block)
                    if size > item.get('max_bytes', 300_000_000):
                        raise ValueError('download size exceeds source-specific limit')
                    out.write(block)
                if expected and size!=int(expected):
                    raise ValueError('truncated download: %d of %s bytes'%(size,expected))
            target.with_suffix(target.suffix+'.part').rename(target)
        data = target.read_bytes()
        if item.get('format') in ('xlsx','zip') and not data.startswith(b'PK'):
            raise ValueError('expected ZIP/XLSX container, not downloaded content')
        if item.get('format') == 'json': json.loads(data)
        record.update(status='downloaded', bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), local_path=str(target))
    except Exception as ex:
        record.update(status='failed', error=str(ex))
    print(json.dumps(record, ensure_ascii=False), flush=True)
    return record

def main():
    p=argparse.ArgumentParser();p.add_argument('spec');a=p.parse_args()
    spec=json.loads(pathlib.Path(a.spec).read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(fetch,spec))
    (ROOT / ('acquisition_'+pathlib.Path(a.spec).stem+'.json')).write_text(json.dumps(records,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
