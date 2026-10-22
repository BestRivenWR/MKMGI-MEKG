"""Fold-only background preprocessing, explicit true availability, shared cache."""
import hashlib, os
from pathlib import Path
import numpy as np
import joblib
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from common import ROOT, config, target, core_module

_raw = None

def load_full_raw():
    global _raw
    if _raw is None:
        z = np.load(ROOT / 'data/background_raw.npz')
        _raw = dict(ids=z['gene_ids'].astype(str), seq=z['sequence25'].astype('float32'),
                    protein=z['esm320'].astype('float32'), seq_present=z['sequence_present'].astype(bool),
                    protein_present=z['protein_present'].astype(bool),
                    tokens=sparse.load_npz(ROOT / 'data/annotation_tokens.npz').tocsr())
    return _raw

def raw_for_indices(ix):
    raw = load_full_raw(); ix = np.asarray(ix, int)
    return {k: v[ix] for k, v in raw.items()}

def masked_scaler(x, present, train):
    use = train & present
    if not use.any():
        return dict(mean=np.zeros(x.shape[1], np.float32), scale=np.ones(x.shape[1], np.float32))
    mean = x[use].mean(0, dtype=np.float64)
    scale = x[use].std(0, dtype=np.float64); scale[scale < 1e-6] = 1.
    return dict(mean=mean.astype('float32'), scale=scale.astype('float32'))

def standardize(x, present, scaler):
    z = ((x - scaler['mean']) / scaler['scale']).astype('float32')
    z[~present] = 0.; assert np.isfinite(z).all()
    return z

def annotation_transform(raw, fit):
    n = len(raw['ids']); tok = raw['tokens'][:, fit['keep']]
    present = np.asarray(tok.getnnz(axis=1)).ravel() > 0
    a = np.zeros((n, config()['annotation_components']), np.float32)
    if fit['svd'] is not None:
        d = fit['svd'].transform(tok).astype('float32'); a[:, :d.shape[1]] = d
    return a, present

def fit_background(raw, train):
    c = config(); assert train.dtype == bool and train.any()
    freq = np.asarray((raw['tokens'][train] > 0).sum(axis=0)).ravel()
    keep = np.flatnonzero(freq >= c['annotation_min_df'])
    svd = None
    if len(keep) >= 2:
        nc = min(c['annotation_components'], len(keep) - 1, int(train.sum()) - 1)
        svd = TruncatedSVD(n_components=nc, n_iter=7, random_state=c['preprocessing_seed'])
        svd.fit(raw['tokens'][train][:, keep].astype(np.float64))
    fit = dict(keep=keep, svd=svd, training_mask_hash=hashlib.sha256(np.packbits(train).tobytes()).hexdigest(),
               fitting_rows=int(train.sum()), fitting_ids_hash=hashlib.sha256(raw['ids'][train].tobytes()).hexdigest())
    a, ap = annotation_transform(raw, fit)
    fit['annotation_scaler'] = masked_scaler(a, ap, train)
    fit['sequence_scaler'] = masked_scaler(raw['seq'], raw['seq_present'], train)
    use = train & raw['protein_present']
    assert use.sum() > c['protein_components']
    pca = PCA(n_components=c['protein_components'], svd_solver='full')
    pca.fit(raw['protein'][use].astype(np.float64)); fit['pca'] = pca
    p = pca.transform(raw['protein']).astype('float32')
    fit['protein_scaler'] = masked_scaler(p, raw['protein_present'], train)
    return fit

def apply_background(raw, fit):
    a, ap = annotation_transform(raw, fit)
    a = standardize(a, ap, fit['annotation_scaler'])
    seq = standardize(raw['seq'], raw['seq_present'], fit['sequence_scaler'])
    p = standardize(fit['pca'].transform(raw['protein']).astype('float32'), raw['protein_present'], fit['protein_scaler'])
    masks = np.column_stack([ap, raw['seq_present'], raw['protein_present']]).astype('float32')
    return np.column_stack([a, seq, p, masks]).astype('float32')

def cached_background(raw, train, key):
    digest = hashlib.sha256(np.packbits(train).tobytes() + raw['ids'].tobytes()).hexdigest()
    base = ROOT / 'cache/background' / (key + '_' + digest[:16])
    meta = base.with_suffix('.joblib'); data = base.with_suffix('.npy')
    if meta.exists() and data.exists():
        fit = joblib.load(meta)
        assert fit['training_mask_hash'] == hashlib.sha256(np.packbits(train).tobytes()).hexdigest()
        return fit, np.load(data)
    fit = fit_background(raw, train); x = apply_background(raw, fit)
    # Same deterministic fit can be requested by CPU and GPU lanes; publish both atomically.
    tmpmeta = target(meta).with_name(meta.name + '.tmp.' + str(os.getpid()))
    tmpdata = target(data).with_name(data.name + '.tmp.' + str(os.getpid()))
    joblib.dump(fit, tmpmeta)
    with tmpdata.open('wb') as handle:
        np.save(handle, x)
    tmpdata.replace(data); tmpmeta.replace(meta)
    return fit, x
