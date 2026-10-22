"""Read-only accounting of peak target ID mapping losses, no database changes."""
from pathlib import Path
import csv,gzip,zipfile,io,collections,json,re
from mekg_paths import workspace
ROOT=workspace()/'v31'
RAW=workspace()/'raw/maizegdb'
known=set()
with gzip.open(RAW/'Zm00001eb.1.fulldata.txt.gz','rt') as f:
    for r in csv.reader(f,delimiter='\t'):
        if len(r)>=14 and re.fullmatch(r'Zm00001eb\d{6}',r[1]):known.add(r[1])
mapping=collections.defaultdict(set);reverse=collections.defaultdict(set)
for version in [3,4]:
    # Match the builder's separate per-version bijection rule.
    raw=collections.defaultdict(set);rev=collections.defaultdict(set)
    with (RAW/f'B73v{version}_to_B73v5.tsv').open() as f:
        for r in csv.reader(f,delimiter='\t'):
            if len(r)<2:continue
            for dest in r[1].split(','):
                if dest in known:raw[r[0]].add(dest);rev[dest].add(r[0])
    for key,values in raw.items():mapping[key].update(values)
    for key,values in rev.items():reverse[(version,key)].update(values)
counts=collections.Counter();unique=collections.defaultdict(set);examples=collections.defaultdict(collections.Counter)
with zipfile.ZipFile(ROOT/'raw/513TFs_ChIPSeeker_regulatory_filtered.zip') as z:
    for name in z.namelist():
        if name.endswith('/'):continue
        with z.open(name) as f:
            for row in csv.DictReader(io.TextIOWrapper(f,encoding='utf8'),delimiter='\t'):
                key=row['geneId'];vals=mapping.get(key,set());v=4 if key.startswith('Zm00001d') else 3
                if key in known:reason='already_v5_not_expected_in_v4_annotation'
                elif not vals:reason='no_curated_target_in_v5_metadata'
                elif len(vals)>1:reason='one_to_many_mapping'
                elif len(reverse[(v,next(iter(vals)))])>1:reason='many_to_one_mapping'
                else:reason='unique_bijective_mapping'
                counts[reason]+=1;unique[reason].add(key)
                if reason!='unique_bijective_mapping':examples[reason][key]+=1
result={'raw_row_counts':dict(counts),'unique_native_ids_by_reason':{k:len(v) for k,v in unique.items()},
        'top_unmapped_native_ids':{k:dict(v.most_common(20)) for k,v in examples.items()},'database_modified':False}
(ROOT/'release/TF_MAPPING_AUDIT.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result),flush=True)
