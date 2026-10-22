"""Documented numerical-only PageRank repair; initial frozen code is preserved."""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
import train_baselines as baseline
from benchmark_common import ROOT, output, save_json, now, merged_graph, verify_training, pv


def transition_float64(graph, n):
    # Cast before reciprocal calculation: float32 division followed by a float64
    # result array still leaves approximately 6e-8 error in each transition row.
    a = merged_graph(graph, n).astype(np.float64)
    degree = np.asarray(a.sum(axis=1)).ravel()
    inv = np.divide(1., degree, out=np.zeros(n, dtype=np.float64), where=degree > 0)
    t = (sparse.diags(inv) @ a).tocsr()
    row_sum = np.asarray(t.sum(axis=1)).ravel()
    assert np.max(np.abs(row_sum[degree > 0] - 1.)) < 1e-12
    return t


if __name__ == '__main__':
    verify_training()
    path = ROOT / 'protocol/ppr_numeric_amendment.json'
    wrapper = Path(__file__).resolve()
    if path.exists():
        old = json.loads(path.read_text())
        assert old['wrapper_sha256'] == pv.sha(wrapper)
    else:
        save_json(path, dict(time=now(), original_freeze_retained=True,
           change='Use float64 adjacency BEFORE row-sum reciprocal in PageRank normalization',
           reason='Mass conservation assertion: sum 1.0000001508012673, row normalization error 5.9138983488082886e-08',
           input_labels_splits_grid_unchanged=True, outer_metrics_not_consulted=True,
           initial_training_source_unchanged=True, wrapper_sha256=pv.sha(wrapper)))
    baseline.transition = transition_float64
    baseline.run_ppr()
