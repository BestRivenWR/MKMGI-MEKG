"""Offline integration checks with synthetic inputs; no GPU, API or research dataset."""
from pathlib import Path
import ast,csv,gzip,importlib.util,json,sys,tempfile,unittest

ROOT=Path(__file__).resolve().parents[1]

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

class OfflineChecks(unittest.TestCase):
    def test_python_and_json(self):
        for p in ROOT.rglob('*.py'):
            if '__pycache__' not in p.parts:ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
        for p in ROOT.rglob('*.json'):
            if 'state' not in p.parts:json.loads(p.read_text(encoding='utf-8'))

    def test_chain_conversion(self):
        c=module('_smoke_coordinates','mekg/_builders/v31/coordinate_tools.py')
        self.assertTrue(c.self_test())

    def test_graph_workspace_paths(self):
        try:import h5py
        except ImportError:self.skipTest('h5py needed for graph context importer')
        import os
        sys.path[:0]=[str(ROOT/'mekg'),str(ROOT/'mekg/_builders/v3')]
        with tempfile.TemporaryDirectory() as td:
            old=os.environ.get('MEKG_WORKSPACE');os.environ['MEKG_WORKSPACE']=td
            try:
                base=module('_smoke_v3_builder','mekg/_builders/v3/build_v3.py')
                context=module('_smoke_v3_context','mekg/_builders/v3/import_context.py')
                self.assertEqual(base.ROOT,Path(td)/'v3');self.assertEqual(context.ROOT,base.ROOT)
                self.assertEqual(base.PARENT,Path(td)/'base')
            finally:
                if old is None:os.environ.pop('MEKG_WORKSPACE',None)
                else:os.environ['MEKG_WORKSPACE']=old

    def test_stored_query_precision(self):
        q=module('_smoke_query','mkmgi/query.py')
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'synthetic.tsv.gz'
            with gzip.open(p,'wt') as f:
                f.write('gene_id\tpriority_score\tpriority_percentile\nZm00001eb000000\t0.123400\t0.5000\n')
            result=q.read_gene(p,'Zm00001eb000000')
            self.assertEqual(result['prediction']['priority_score'],'0.123400')
            self.assertFalse(q.read_gene(p,'Zm00001eb999999')['supported_by_table'])

    def test_trait_provenance(self):
        traits=module('_smoke_traits','agent/trait_display.py')
        loci=[{'position':100,'pvalue':.01,'model':'synthetic_A'},
              {'position':100,'pvalue':.02,'model':'synthetic_B'}]
        got=traits.display_context({'evidence_type':'GWAS_LOCUS_PROXIMITY'},{'locus_results':loci})
        self.assertEqual(got['locus_results'],loci);self.assertNotIn('combined_pvalue',got)

    def test_agent_index_schema(self):
        sys.path.insert(0,str(ROOT/'agent'))
        setup=module('_smoke_setup','agent/setup.py')
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'source_data').mkdir();(root/'deployment').mkdir()
            for name in ['ppi_base_edges.tsv.gz','regulatory_tf_supported_edges.tsv.gz']:
                with gzip.open(root/'source_data'/name,'wt') as f:
                    f.write('source\ttarget\tweight\tconfidence\nZm00001eb000000\tZm00001eb000010\t1.0\t0.9\n')
            atlas=root/'atlas.tsv'
            atlas.write_text('gene_id_v5\tsymbols\tprimary_dois\tphenotype_summary\n'
                'Zm00001eb000000\tSYNTHETIC\t\tsynthetic fixture\n')
            old=setup.ROOT;setup.ROOT=root/'deployment'
            try:path=setup.build_overlay({'source_root':str(root),'public_atlas':str(atlas)})
            finally:setup.ROOT=old
            import sqlite3
            db=sqlite3.connect(path)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM model_edges').fetchone()[0],2)
            self.assertEqual(db.execute('SELECT alias FROM search_aliases').fetchone()[0],'SYNTHETIC');db.close()

    def test_model_forward(self):
        try:import torch
        except ImportError:self.skipTest('PyTorch not installed')
        sys.path.insert(0,str(ROOT/'mkmgi'))
        model=module('model','mkmgi/model.py')
        torch.set_num_threads(1);torch.manual_seed(19)
        net=model.BackgroundGCN('time_attention',0).eval();n=8
        x=torch.randn(n,191);x[:,40:80]=1;x[0,40:80]=0
        ix=torch.arange(n);edges=torch.stack([ix,(ix+1)%n])
        matrix=torch.sparse_coo_tensor(edges,torch.ones(n),(n,n)).coalesce()
        graph=([matrix,matrix,matrix],torch.ones(n,4,dtype=torch.bool))
        score,gates,weights=net(x,graph,return_gates=True)
        self.assertEqual(tuple(score.shape),(n,));self.assertEqual(len(gates),2)
        self.assertTrue(torch.isfinite(score).all());self.assertEqual(float(weights[0].sum()),0)
        score.square().mean().backward()
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in net.parameters() if p.grad is not None))

if __name__=='__main__':unittest.main(verbosity=2)
