"""Strict induced graphs, transparent summaries, role-stratified degree nulls."""
import hashlib, json
import numpy as np
from common import ROOT, config, target, save

def induced(graph, allowed):
    out = {}
    for kind in ['ppi', 'reg']:
        e = graph[kind + '_edge_index']
        ok = allowed[e[0]] & allowed[e[1]]
        out[kind + '_edge_index'] = e[:, ok].copy()
        out[kind + '_base_weight'] = graph[kind + '_base_weight'][ok].copy()
    return out

def summary(graph, n):
    p = graph['ppi_edge_index']; r = graph['reg_edge_index']
    pw = graph['ppi_base_weight']; rw = graph['reg_base_weight']
    pp = np.bincount(p.ravel(), minlength=n)
    ri = np.bincount(r[1], minlength=n); ro = np.bincount(r[0], minlength=n)
    pws = np.bincount(p[0], weights=pw, minlength=n) + np.bincount(p[1], weights=pw, minlength=n)
    ris = np.bincount(r[1], weights=rw, minlength=n); ros = np.bincount(r[0], weights=rw, minlength=n)
    return np.column_stack([np.log1p(pp), np.log1p(ri), np.log1p(ro), np.log1p(pws),
                            np.log1p(ris), np.log1p(ros), (pp + ri + ro) > 0]).astype('float32')

def degree(graph, n):
    p = graph['ppi_edge_index']; r = graph['reg_edge_index']
    return np.column_stack([np.bincount(p.ravel(), minlength=n), np.bincount(r[1], minlength=n),
                            np.bincount(r[0], minlength=n)])

def _swap(e, strata, directed, seed):
    """Preserve role-specific degrees, no self loops or duplicate edges."""
    c = config(); rng = np.random.default_rng(seed); work = e.T.copy()
    pool = {(int(a), int(b)) if directed else (min(int(a), int(b)), max(int(a), int(b))) for a, b in work}
    assert len(pool) == len(work)
    audits = []
    for stratum in np.unique(strata):
        slots = np.flatnonzero(strata == stratum)
        requested = int(len(slots) * c['rewire_swaps_per_edge']); accepted = attempts = 0
        if len(slots) >= 2:
            while accepted < requested and attempts < requested * c['rewire_attempts_per_requested_swap']:
                attempts += 1
                i, j = rng.choice(slots, size=2, replace=False)
                a, b = map(int, work[i]); u, v = map(int, work[j])
                # Endpoints have already been oriented by role, including undirected cross-role edges.
                if len({a, b, u, v}) < 4:
                    continue
                x, y = (a, v), (u, b)
                keyx = x if directed else (min(x), max(x))
                keyy = y if directed else (min(y), max(y))
                if keyx in pool or keyy in pool or keyx == keyy:
                    continue
                oldi = (a, b) if directed else (min(a, b), max(a, b))
                oldj = (u, v) if directed else (min(u, v), max(u, v))
                pool.remove(oldi); pool.remove(oldj); pool.add(keyx); pool.add(keyy)
                work[i] = x; work[j] = y; accepted += 1
        audits.append(dict(stratum=int(stratum), edges=int(len(slots)), requested=requested,
                           accepted=accepted, attempts=attempts))
    if not directed:
        work.sort(axis=1)
    return work.T, audits

def rewired_pair(graph, train_allowed, inference_allowed, seed, cache_key):
    """Training edges stay identical at inference; rewiring never uses phenotype labels."""
    path = ROOT / 'cache' / (cache_key + '_rewire_' + str(seed) + '.npz')
    sig = hashlib.sha256(np.packbits(train_allowed).tobytes() + np.packbits(inference_allowed).tobytes()).hexdigest()
    auditpath = path.with_suffix('.json')
    if path.exists() and auditpath.exists():
        audit = json.loads(auditpath.read_text()); assert audit['mask_hash'] == sig
        with np.load(path) as z:
            g = {k: z[k] for k in z.files}
        return induced(g, train_allowed), g
    real = induced(graph, inference_allowed); n = len(train_allowed); out = {}; audits = {}
    for k, directed in [('ppi', False), ('reg', True)]:
        e = real[k + '_edge_index'].copy()
        roles = np.where(train_allowed, 0, 1)
        if not directed:
            flip = (roles[e[0]] > roles[e[1]])
            e[:, flip] = e[::-1, flip]
        strata = roles[e[0]] * 2 + roles[e[1]]
        new, aud = _swap(e, strata, directed, seed + (0 if k == 'ppi' else 100000))
        out[k + '_edge_index'] = new
        # Confidence summaries are NOT computed from this null; they stay fixed to the real graph.
        out[k + '_base_weight'] = real[k + '_base_weight'].copy()
        original = {tuple(sorted(x)) if not directed else tuple(x) for x in e.T.tolist()}
        retained = sum((tuple(x) in original) for x in new.T.tolist())
        audits[k] = dict(swaps=aud, retained_original_edges=int(retained), edges=int(e.shape[1]),
                         retained_fraction=float(retained / max(1, e.shape[1])))
    assert np.array_equal(degree(real, n), degree(out, n)), 'Null must preserve inference degrees'
    realtrain = induced(real, train_allowed); nulltrain = induced(out, train_allowed)
    assert np.array_equal(degree(realtrain, n), degree(nulltrain, n)), 'Null must preserve training degrees'
    np.savez_compressed(target(path), **out)
    save(auditpath, dict(mask_hash=sig, seed=int(seed), degree_preservation_exact=True,
                        training_edges_locked_at_inference=True, conditions=audits))
    return nulltrain, out

def tensors(graph, n, device):
    import torch
    p = graph['ppi_edge_index']; r = graph['reg_edge_index']
    channels = [np.concatenate([p, p[::-1]], axis=1), r[::-1], r]
    matrices = []; avail = np.ones((n, 4), dtype=bool)
    for k, e in enumerate(channels):
        counts = np.bincount(e[0], minlength=n)
        vals = (1. / np.maximum(1, counts[e[0]])).astype('float32')
        mat = torch.sparse_coo_tensor(torch.as_tensor(e.astype('int64'), device=device),
                                      torch.as_tensor(vals, device=device), (n, n)).coalesce()
        matrices.append(mat); avail[:, k + 1] = counts > 0
    return matrices, torch.as_tensor(avail, device=device)
