"""Conservative UCSC chain subset: unique, whole-interval, ungapped mappings.

All intervals here are zero-based, half-open. An interval crossing a chain gap
is deliberately rejected rather than stretching its endpoints across an indel.
"""
from collections import defaultdict
from bisect import bisect_left,bisect_right

def chrom(s):
    s=str(s).strip()
    if s.lower().startswith('chr'):s=s[3:]
    if s.startswith('ctg'):s='B73V4_'+s
    return s

class IntervalIndex:
    def __init__(self,records):
        data=defaultdict(list)
        for c,s,e,item in records:
            assert 0<=s<e
            data[chrom(c)].append((s,e,item))
        self.data={};self.starts={};self.prefix={}
        for c,rs in data.items():
            rs.sort(key=lambda x:(x[0],x[1]));self.data[c]=rs;self.starts[c]=[x[0] for x in rs]
            mx=0;pre=[]
            for s,e,item in rs:mx=max(mx,e);pre.append(mx)
            self.prefix[c]=pre
    def overlap(self,c,s,e):
        c=chrom(c);rs=self.data.get(c,[]);out=[]
        j=bisect_left(self.starts.get(c,[]),e)-1;pre=self.prefix.get(c,[])
        while j>=0 and pre[j]>s:
            a,b,x=rs[j]
            if a<e and b>s:out.append((a,b,x))
            j-=1
        return out

class StrictChain:
    def __init__(self,path):
        blocks=[];reverse=[];self.source_sizes={};self.target_sizes={}
        with open(path) as f:
            for line in f:
                x=line.split()
                if not x:continue
                if x[0]=='chain':
                    _,score,tc,ts,tst,t0,t1,qc,qs,qst,q0,q1,ci=x
                    assert tst=='+'
                    ts=int(ts);qs=int(qs);t=int(t0);q=int(q0)
                    self.source_sizes[tc]=ts;self.target_sizes[qc]=qs
                    continue
                size=int(x[0]);a=q if qst=='+' else qs-q-size;b=a+size
                item=(qc,a,b,qst,ci)
                blocks.append((tc,t,t+size,item));reverse.append((qc,a,b,(tc,t,t+size,qst,ci)))
                if len(x)==3:t+=size+int(x[1]);q+=size+int(x[2])
        self.forward=IntervalIndex(blocks);self.reverse=IntervalIndex(reverse)
    def _map(self,index,c,s,e):
        hits=[]
        for lo,hi,(dest,a,b,orient,ci) in index.overlap(c,s,e):
            if lo<=s and e<=hi:
                if orient=='+':start=a+s-lo;end=a+e-lo
                else:start=b-(e-lo);end=b-(s-lo)
                hits.append((dest,start,end,orient,ci))
        # Reject multiple supporting chain blocks, even if endpoints happen to match.
        return hits[0] if len(hits)==1 else None
    def map(self,c,s,e):
        if not 0<=s<e:return None
        got=self._map(self.forward,c,s,e)
        if got is None:return None
        back=self._map(self.reverse,*got[:3])
        if back is None or back[:3]!=(chrom(c),s,e):return None
        return got

def self_test():
    import tempfile,pathlib
    text='chain 1 1 1000 + 100 120 2 2000 - 300 320 1\n20\n\nchain 1 1 1000 + 200 215 2 2000 + 100 115 2\n5 5 5\n5\n'
    with tempfile.TemporaryDirectory() as td:
        p=pathlib.Path(td)/'test.chain';p.write_text(text);c=StrictChain(p)
        assert c.map('chr1',102,108)==('2',1692,1698,'-','1')
        assert c.map('1',202,204)==('2',102,104,'+','2')
        assert c.map('1',202,213) is None
        assert c.map('1',99,101) is None
        assert c.map('1',100,120)==('2',1680,1700,'-','1')
    return True
