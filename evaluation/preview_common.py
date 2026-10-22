from eval_support import *
import json, os, hashlib, datetime
import numpy as np
CONFIG = ablation_config()
ROOT = ablation_root()
def output(path):
    p = Path(path).resolve()
    if ROOT.resolve() not in p.parents:
        raise RuntimeError('Write outside new experiment directory: ' + str(p))
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save_json(path, data):
    p = output(path)
    tmp = p.with_name(p.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    tmp.replace(p)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda: stream.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def freeze():
    return audit_freeze(ROOT, CONFIG, __file__)


def check_source_hashes():
    return audit_verify(ROOT)


def load():
    import numpy as np
    m, raw, graph = source_common.load()
    assert raw.shape == (35275, 92)
    assert int(m.primary.sum()) == 228
    assert int((m.primary | ~m.excluded_U).sum()) == 35017
    assert np.array_equal(m.row_index.to_numpy(int), np.arange(len(m)))
    bg = source_background.raw_for_indices(m.gene_index.to_numpy(int))
    assert np.array_equal(bg['ids'].astype(str), m.gene_id.to_numpy(str))
    return m, raw, graph, bg


def matrices_numpy(raw, graph, training, inference, prep, bg, temporal_masks=True):
    import numpy as np
    gt = graphops.induced(graph, training)
    gi = graphops.induced(graph, inference)
    assert all(training[e].all() for k, e in gt.items() if k.endswith('_edge_index'))
    rt = np.column_stack([raw, graphops.summary(gt, len(raw))])
    ri = np.column_stack([raw, graphops.summary(gi, len(raw))])
    bx = source_background.apply_background(bg, prep['background'])
    xt = np.column_stack([source_model.apply_preprocess(rt, prep['core']), bx]).astype('float32')
    xi = np.column_stack([source_model.apply_preprocess(ri, prep['core']), bx]).astype('float32')
    if temporal_masks:
        xt[:, 40:80] = raw[:, 40:80]
        xi[:, 40:80] = raw[:, 40:80]
    assert xt.shape == (len(raw), 191) and xi.shape == xt.shape
    assert np.isfinite(xt).all() and np.isfinite(xi).all()
    return xt, xi, gt, gi


FEATURE_COLUMNS = {
    'drop_omics': list(range(80)),
    'drop_spatial': list(range(80, 92)),
    'drop_function': list(range(99, 131)) + [188],
    'drop_protein': list(range(131, 188)) + [189, 190],
    'drop_graph_all': list(range(92, 99)),
    'drop_sequence': list(range(131, 156)) + [189],
    'drop_ESM': list(range(156, 188)) + [190],
}
for j, name in enumerate(['RNA', 'ATAC', 'H3K4me3', 'H3K27me3', 'H3K27ac']):
    FEATURE_COLUMNS['drop_' + name] = list(range(j * 8, j * 8 + 8)) + list(range(40 + j * 8, 48 + j * 8))


def delete_features(xt, xi, condition):
    columns = FEATURE_COLUMNS.get(condition, [])
    if columns:
        xt = xt.clone()
        xi = xi.clone()
        xt[:, columns] = 0.
        xi[:, columns] = 0.
        assert not xt[:, columns].any() and not xi[:, columns].any()
    return xt, xi


def empty_graph(n, device):
    import torch
    z = torch.sparse_coo_tensor(torch.empty((2, 0), dtype=torch.long, device=device),
                               torch.empty(0, device=device), (n, n)).coalesce()
    available = torch.zeros((n, 4), dtype=torch.bool, device=device)
    available[:, 0] = True
    return [z, z, z], available


def reference_selection(fold):
    p = SOURCE / 'runs/family80' / ('outer' + str(fold)) / 'time_attention/selection.json'
    return json.loads(p.read_text())
