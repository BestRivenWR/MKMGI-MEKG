from eval_support import *
import json, os, hashlib, datetime
import numpy as np
CONFIG = benchmark_config()
ROOT = benchmark_root()
PREVIOUS = ROOT / "optional_reference_cache"
def output(path):
    p = Path(path).resolve()
    assert ROOT.resolve() in p.parents, 'Outside new output directory: ' + str(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save_json(path, data):
    p = output(path)
    t = p.with_name(p.name + '.tmp.' + str(os.getpid()))
    t.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    t.replace(p)


def status(worker, stage, **extra):
    record = dict(time=now(), pid=os.getpid(), worker=worker, stage=stage, **extra)
    save_json(ROOT / 'logs' / (worker + '_status.json'), record)
    print(json.dumps(record), flush=True)


def freeze():
    return audit_freeze(ROOT, CONFIG, __file__)


def verify_training():
    return audit_verify(ROOT)


def verify_sources():
    return audit_verify(ROOT)


def load():
    return pv.load()


def context(m, raw, graph, bg, training, inference, key, scaled=True):
    import joblib
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    eligible = m.primary.to_numpy(bool) | ~m.excluded_U.to_numpy(bool)
    loss = training & eligible
    digest = hashlib.sha256(np.packbits(loss).tobytes()).hexdigest()
    cache = ROOT / 'cache' / key
    path = cache / 'preprocessor.joblib'
    if path.exists():
        fitted = joblib.load(path)
        metadata = json.loads((cache / 'metadata.json').read_text())
        assert metadata['training_mask_hash'] == digest
    else:
        old = PREVIOUS / 'cache/baselines' / key / 'preprocessor.joblib'
        if old.exists():
            metadata = json.loads(old.with_name('metadata.json').read_text())
            assert metadata['training_mask_hash'] == digest
            prep = joblib.load(old)
        else:
            gt = graphops.induced(graph, training)
            rt = np.column_stack([raw, graphops.summary(gt, len(raw))])
            prep = dict(core=source_model.fit_preprocess(rt, loss),
                        background=source_background.fit_background(bg, loss))
        xt, _, _, _ = pv.matrices_numpy(raw, graph, training, inference, prep, bg, temporal_masks=False)
        scaler = StandardScaler().fit(xt[loss])
        fitted = dict(prep=prep, scaler=scaler)
        save_json(cache / 'metadata.json', dict(training_mask_hash=digest, training_N=int(loss.sum()),
                  preprocessing='Only eligible training genes; previous fold-fitted base preprocessing reusable'))
        temp = output(path.with_name(path.name + '.tmp.' + str(os.getpid())))
        joblib.dump(fitted, temp)
        temp.replace(output(path))
    xt, xi, gt, gi = pv.matrices_numpy(raw, graph, training, inference, fitted['prep'], bg, temporal_masks=False)
    if scaled:
        xt = fitted['scaler'].transform(xt).astype('float32')
        xi = fitted['scaler'].transform(xi).astype('float32')
    assert xt.shape == xi.shape == (len(m), 191)
    assert np.isfinite(xt).all() and np.isfinite(xi).all()
    return dict(xt=xt, xi=xi, gt=gt, gi=gi, loss=loss, preprocessing=fitted)


def merged_graph(graph, n):
    """Undirected binary union, no self loops; relation labels/weights discarded."""
    import numpy as np
    from scipy import sparse
    e = np.concatenate([graph['ppi_edge_index'], graph['reg_edge_index']], axis=1)
    row = np.r_[e[0], e[1]]
    col = np.r_[e[1], e[0]]
    a = sparse.coo_matrix((np.ones(len(row), dtype=np.float32), (row, col)), shape=(n, n)).tocsr()
    a.sum_duplicates()
    a.data[:] = 1.
    a.setdiag(0.)
    a.eliminate_zeros()
    assert a.diagonal().sum() == 0 and (a != a.T).nnz == 0
    return a


def normalized_tensor(graph, n, device):
    import numpy as np
    import torch
    from scipy import sparse
    a = merged_graph(graph, n) + sparse.eye(n, dtype=np.float32, format='csr')
    d = np.asarray(a.sum(axis=1)).ravel()
    inv = np.power(d, -.5).astype('float32')
    z = (sparse.diags(inv) @ a @ sparse.diags(inv)).tocoo()
    return torch.sparse_coo_tensor(torch.as_tensor(np.stack([z.row, z.col]), dtype=torch.long, device=device),
             torch.as_tensor(z.data, device=device), (n, n)).coalesce()


def write_predictions(m, fold, condition, scores, metadata):
    import numpy as np
    z = m.loc[m.outer80.eq(fold), ['gene_id', 'row_index', 'primary', 'excluded_U', 'ppr', 'tf', 'group80', 'group50']].copy()
    assert len(scores) == len(CONFIG['seeds'])
    assert all(len(s) == len(z) and np.isfinite(s).all() for s in scores)
    for seed, score in zip(CONFIG['seeds'], scores):
        z['score_seed' + str(seed)] = score
    z['score'] = np.mean(scores, axis=0)
    z['uncertainty_sd'] = np.std(scores, axis=0, ddof=1)
    path = ROOT / 'runs/family80' / ('outer' + str(fold)) / condition
    z.to_csv(output(path / 'outer_predictions.tsv.gz'), sep='\t', index=False, compression='gzip')
    save_json(path / 'complete.json', dict(time=now(), fold=fold, condition=condition, **metadata))
