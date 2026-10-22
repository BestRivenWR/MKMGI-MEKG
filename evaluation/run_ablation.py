"""Bounded fixed-config refits, paired with the unmodified current time-GCN."""
import argparse
import gc
import json
import time
import numpy as np
import pandas as pd
import torch
from torch import nn
import joblib
from preview_common import (ROOT, SOURCE, CONFIG, output, now, save_json, freeze,
    load, matrices_numpy, delete_features, empty_graph, graphops,
    source_model, reference_selection, check_source_hashes, FEATURE_COLUMNS)


class TokenOnlyEncoder(nn.Module):
    """Keep projection, DAP, normalization and pooling; remove token interactions."""
    def __init__(self, reference):
        super().__init__()
        self.token = reference.token
        self.position = reference.position
        self.output_norm = reference.output_norm
        self.pool = reference.pool
        self.register_buffer('dap', reference.dap.detach().clone())

    def forward(self, values, observed):
        present = observed.any(-1)
        valid = present.clone()
        empty = ~valid.any(1)
        valid[empty, 0] = True
        values = values.masked_fill(~observed, 0.)
        token = self.token(torch.cat([values, observed.to(values.dtype)], -1))
        token = self.output_norm(token + self.position(self.dap[:, None])[None, :, :])
        weight = torch.softmax(self.pool(token).squeeze(-1).masked_fill(~valid, -1e9), 1)
        weight = weight.masked_fill(empty[:, None], 0.)
        return (token * weight[:, :, None]).sum(1), weight


class AblationGCN(source_model.BackgroundGCN):
    def __init__(self, condition, dropout):
        super().__init__('time_attention', dropout)
        self.ablation = condition
        if condition == 'no_time_transformer':
            self.time = TokenOnlyEncoder(self.time)
        elif condition == 'gcn_1layer':
            self.blocks = nn.ModuleList([self.blocks[0]])


def fit(xt, xi, gt, gi, mask, labels, condition, selection, seed):
    source_model.seed_all(seed)
    model = AblationGCN(condition, selection['params']['dropout']).to(xt.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=selection['params']['lr'],
                                  weight_decay=source_model.config()['weight_decay'])
    ix = torch.as_tensor(np.flatnonzero(mask), device=xt.device)
    yy = torch.as_tensor(labels.astype('float32'), device=xt.device)
    weight = torch.tensor((len(yy) - yy.sum().item()) / yy.sum().item(), device=xt.device)
    criterion = nn.BCEWithLogitsLoss(pos_weight=weight)
    trace = []
    for epoch in range(1, selection['refit_epochs'] + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        logits = model(xt, gt)
        loss = criterion(logits[ix], yy)
        assert torch.isfinite(loss), 'Non-finite training loss'
        loss.backward()
        grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
        optimizer.step()
        trace.append(dict(epoch=epoch, loss=float(loss.detach().cpu()),
                          gradient_norm=float(grad_norm.detach().cpu())))
    model.eval()
    with torch.no_grad():
        score = model(xi, gi).sigmoid().cpu().numpy()
    assert np.isfinite(score).all()
    state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
    result = dict(score=score, state=state, trace=trace,
                  parameter_count=sum(p.numel() for p in model.parameters()))
    del optimizer, model
    gc.collect()
    torch.cuda.empty_cache()
    return result


def progress(stage, **extra):
    record = dict(time=now(), stage=stage, **extra)
    save_json(ROOT / 'logs/ablation_status.json', record)
    print(json.dumps(record), flush=True)


def make_context(fold, m, raw, graph, bg, device):
    outer = m.outer80.to_numpy(int)
    allowed = outer != fold
    eligible = m.primary.to_numpy(bool) | ~m.excluded_U.to_numpy(bool)
    mask = allowed & eligible
    assert not set(m.loc[allowed, 'group80']) & set(m.loc[~allowed, 'group80'])
    base = SOURCE / 'runs/family80' / ('outer' + str(fold)) / 'time_attention'
    prep = joblib.load(base / 'seed20261001/preprocessor.joblib')
    xt, xi, gt, gi = matrices_numpy(raw, graph, allowed, np.ones(len(m), bool), prep, bg)
    return (torch.as_tensor(xt, device=device), torch.as_tensor(xi, device=device),
            graphops.tensors(gt, len(m), device), graphops.tensors(gi, len(m), device), mask, ~allowed)


def smoke(m, raw, graph, bg, device):
    p = ROOT / 'qa/reference_replay.json'
    if p.exists():
        record = json.loads(p.read_text())
        assert record['pass']
        return
    progress('reference_replay')
    xt, xi, gt, gi, mask, test = make_context(0, m, raw, graph, bg, device)
    base = SOURCE / 'runs/family80/outer0/time_attention'
    original = pd.read_csv(base / 'outer_predictions.tsv.gz', sep='\t')
    assert np.array_equal(original.gene_id.to_numpy(str), m.loc[test, 'gene_id'].to_numpy(str))
    selection = reference_selection(0)
    expected = original.score_seed20261001.to_numpy(float)
    model = AblationGCN('reference_replay', selection['params']['dropout']).to(device)
    ck = torch.load(base / 'seed20261001/model.pt', map_location=device)
    model.load_state_dict(ck['state_dict'])
    model.eval()
    with torch.no_grad():
        exact = model(xi, gi).sigmoid().cpu().numpy()[test]
    checkpoint_error = float(np.max(np.abs(exact - expected)))
    assert checkpoint_error < 2e-6, 'Source checkpoint/input replay failed'
    del model, ck
    start = time.time()
    result = fit(xt, xi, gt, gi, mask, m.primary.to_numpy(bool)[mask],
                 'reference_replay', selection, CONFIG['seeds'][0])
    refit_error = float(np.max(np.abs(result['score'][test] - expected)))
    assert refit_error < 2e-5, 'Source training algorithm replay failed'
    save_json(p, dict(time=now(), checkpoint_max_abs_error=checkpoint_error,
                     refit_max_abs_error=refit_error, runtime_seconds=time.time() - start,
                     epochs=selection['refit_epochs'], pass_=True, **{'pass': True}))
    for condition, columns in FEATURE_COLUMNS.items():
        xx, ii = delete_features(xt, xi, condition)
        assert not xx[:, columns].any() and not ii[:, columns].any()
        if condition == 'drop_omics':
            probe = AblationGCN(condition, .1).to(device).eval()
            with torch.no_grad():
                values = xx[:100, :40].reshape(-1, 5, 8).transpose(1, 2)
                observed = xx[:100, 40:80].reshape(-1, 5, 8).transpose(1, 2) > .5
                z, a = probe.time(values, observed)
                assert not z.any() and not a.any(), 'All-missing time pathway must be zero'
            del probe
        del xx, ii
    eg = empty_graph(len(m), device)
    assert eg[1][:, 0].all() and not eg[1][:, 1:].any()
    assert all(x._nnz() == 0 for x in eg[0])
    save_json(ROOT / 'qa/deletion_tests.json', dict(pass_=True, conditions=list(FEATURE_COLUMNS),
              graph_summary_and_availability_removed=True, temporal_all_missing_output_zero=True))
    del xt, xi, gt, gi, result
    gc.collect()
    torch.cuda.empty_cache()


def main(only=None, smoke_only=False):
    freeze()
    torch.set_num_threads(CONFIG['threads'])
    assert torch.cuda.is_available()
    device = torch.device('cuda:0')
    m, raw, graph, bg = load()
    smoke(m, raw, graph, bg, device)
    if smoke_only:
        print('Reference replay and deletion tests passed', flush=True)
        return
    conditions = CONFIG['ablation_conditions'] if only is None else [only]
    started = time.time()
    for fold in range(5):
        selection = reference_selection(fold)
        ctx = None
        for condition in conditions:
            path = ROOT / 'runs/family80' / ('outer' + str(fold)) / condition
            if (path / 'complete.json').exists():
                continue
            if ctx is None:
                ctx = make_context(fold, m, raw, graph, bg, device)
            xt, xi, gt, gi, mask, test = ctx
            xx, ii = delete_features(xt, xi, condition)
            if condition == 'drop_graph_all':
                gtrain = ginfer = empty_graph(len(m), device)
                assert not xx[:, 92:99].any() and not ii[:, 92:99].any()
            else:
                gtrain, ginfer = gt, gi
            selection_note = {**selection, 'used_from_full_model': True,
                              'retuned_for_ablation': False, 'source_outer_fold': fold}
            save_json(path / 'selection.json', selection_note)
            scores = []
            for seed in CONFIG['seeds']:
                predpath = path / ('seed' + str(seed)) / 'outer_scores.npy'
                if predpath.exists():
                    score = np.load(predpath)
                else:
                    progress('refit', fold=fold, condition=condition, seed=seed,
                             epochs=selection['refit_epochs'], elapsed_seconds=time.time() - started)
                    result = fit(xx, ii, gtrain, ginfer, mask, m.primary.to_numpy(bool)[mask],
                                 condition, selection, seed)
                    score = result['score'][test]
                    np.save(output(predpath), score)
                    torch.save(dict(state_dict=result['state'], condition=condition, seed=seed,
                         params=selection['params'], epochs=selection['refit_epochs'],
                         architecture='CurrentTimeGCN_fixed_config_ablation',
                         parameter_count=result['parameter_count']), output(predpath.parent / 'model.pt'))
                    save_json(predpath.parent / 'training_trace.json', result['trace'])
                    del result
                assert len(score) == int(test.sum()) and np.isfinite(score).all()
                scores.append(score)
            cols = ['gene_id', 'row_index', 'primary', 'excluded_U', 'ppr', 'tf', 'group80', 'group50']
            z = m.loc[test, cols].copy()
            for seed, score in zip(CONFIG['seeds'], scores):
                z['score_seed' + str(seed)] = score
            z['score'] = np.mean(scores, axis=0)
            z['uncertainty_sd'] = np.std(scores, axis=0, ddof=1)
            z.to_csv(output(path / 'outer_predictions.tsv.gz'), sep='\t', index=False, compression='gzip')
            save_json(path / 'complete.json', dict(time=now(), fold=fold, condition=condition,
                      seeds=CONFIG['seeds'], test_nodes=int(test.sum()),
                      fixed_full_model_configuration=True))
            del xx, ii, scores
            gc.collect()
            torch.cuda.empty_cache()
        del ctx
        gc.collect()
        torch.cuda.empty_cache()
    integrity = check_source_hashes()
    save_json(ROOT / 'qa/after_ablation_source_hashes.json', integrity)
    progress('complete', conditions=conditions, fits=len(conditions) * 5 * len(CONFIG['seeds']),
             elapsed_seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--condition', choices=CONFIG['ablation_conditions'])
    parser.add_argument('--smoke-only', action='store_true')
    args = parser.parse_args()
    main(args.condition, args.smoke_only)
