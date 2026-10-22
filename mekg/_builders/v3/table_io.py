"""Read author supplementary XLSX without modifying files or Excel dependencies."""
from __future__ import annotations
import io, posixpath, zipfile, xml.etree.ElementTree as ET

NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
REL='{http://schemas.openxmlformats.org/package/2006/relationships}'
RID='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'

def sheets(path_or_bytes):
    with zipfile.ZipFile(io.BytesIO(path_or_bytes) if isinstance(path_or_bytes,bytes) else path_or_bytes) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            strings=[''.join(x.itertext()) for x in root.findall(NS+'si')]
        rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        for s in ET.fromstring(z.read('xl/workbook.xml')).findall(NS+'sheets/'+NS+'sheet'):
            target=rels[s.attrib[RID]]
            member=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
            rows=[]
            with z.open(member) as f:
                for _,el in ET.iterparse(f,events=('end',)):
                    if el.tag!=NS+'row': continue
                    cells={}
                    for c in el.findall(NS+'c'):
                        col=0
                        for char in ''.join(x for x in c.attrib.get('r','A1') if x.isalpha()): col=col*26+ord(char.upper())-64
                        value=c.findtext(NS+'v','')
                        if c.attrib.get('t')=='s' and value: value=strings[int(value)]
                        elif c.attrib.get('t')=='inlineStr': value=''.join(c.find(NS+'is').itertext())
                        cells[col]=value
                    rows.append([cells.get(i,'') for i in range(1,max(cells,default=0)+1)])
                    el.clear()
            yield s.attrib['name'],rows

if __name__=='__main__':
    import pathlib,json
    root=pathlib.Path(__file__).resolve().parent/'raw'
    for f in sorted(root.glob('*.xlsx')):
        print('FILE',f.name,flush=True)
        for name,rows in sheets(f):print(json.dumps({'sheet':name,'nrows':len(rows),'top':rows[:4]},ensure_ascii=False)[:4500],flush=True)
    for f in root.glob('*tables.zip'):
        with zipfile.ZipFile(f) as z:
            print('ZIP',f.name,z.namelist(),flush=True)
            for member in z.namelist():
                if member.lower().endswith('.xlsx'):
                    for name,rows in sheets(z.read(member)):print(json.dumps({'file':member,'sheet':name,'nrows':len(rows),'top':rows[:4]},ensure_ascii=False)[:4500],flush=True)
