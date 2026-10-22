"""Finite nested family-held-out fitting; outer labels never select a model."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import argparse, gc, json, time
import numpy as np
import pandas as pd
import torch
import joblib
from sklearn.metrics import average_precision_score
from common import ROOT, config, save, target, load, verify_freeze, now, core_module
from background import raw_for_indices, cached_background
from model import BackgroundGCN, seed_all, fit_preprocess, apply_preprocess

graphops = core_module('graph_ops')

def directory(grouping, fold, condition):
    return ROOT / 'runs' / ('family' + str(grouping)) / ('outer' + str(fold)) / condition

def report_progress(path, stage, **kw):
    record = dict(updated=now(), pid=os.getpid(), stage=stage, **kw)
    save(path / 'progress.json', record); print(json.dumps(record), flush=True)

def context(raw, graph, train_mask, infer_mask, loss_mask, condition, seed, key, device, bgraw):
    real_train = graphops.induced(graph, train_mask); real_infer = graphops.induced(graph, infer_mask)
    rt = np.column_stack([raw, graphops.summary(real_train, len(raw))])
    ri = np.column_stack([raw, graphops.summary(real_infer, len(raw))])
    fitted = fit_preprocess(rt, loss_mask)
    bgfit, bx = cached_background(bgraw, loss_mask, key)
    xt = np.column_stack([apply_preprocess(rt, fitted), bx]).astype('float32')
    xi = np.column_stack([apply_preprocess(ri, fitted), bx]).astype('float32')
    assert xt.shape[1] == 191 and np.isfinite(xt).all() and np.isfinite(xi).all()
    if condition == 'gcn_rewired':
        gt, gi = graphops.rewired_pair(graph, train_mask, infer_mask, seed, key)
    else:
        gt, gi = real_train, real_infer
    assert all(train_mask[e].all() for k, e in gt.items() if k.endswith('_edge_index'))
    assert all(infer_mask[e].all() for k, e in gi.items() if k.endswith('_edge_index'))
    payload = dict(preprocessor=dict(core=fitted, background=bgfit), input_dim=191, condition=condition)
    if condition == 'ExtraTrees_background':
        payload.update(xt=xt, xi=xi)
    else:
        payload.update(xt=torch.as_tensor(xt, device=device), xi=torch.as_tensor(xi, device=device),
                       gt=graphops.tensors(gt, len(raw), device), gi=graphops.tensors(gi, len(raw), device))
    return payload

def neural_fit(ctx, train_loss_mask, train_labels, params, seed, validation_mask=None,
               validation_labels=None, epochs=None):
    c = config(); seed_all(seed); device = ctx['xt'].device
    model = BackgroundGCN(ctx['condition'], params['dropout']).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=params['lr'], weight_decay=c['weight_decay'])
    ix = torch.as_tensor(np.flatnonzero(train_loss_mask), device=device)
    yy = torch.as_tensor(np.asarray(train_labels, np.float32), device=device)
    assert len(yy) == len(ix) and 0 < float(yy.sum()) < len(yy)
    pos_weight = torch.tensor((len(yy) - yy.sum().item()) / yy.sum().item(), device=device)
    criterion = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    best_score = -np.inf; best_epoch = 0; bad = 0; history = []
    limit = c['max_epochs'] if epochs is None else int(epochs)
    for epoch in range(1, limit + 1):
        model.train(); optimizer.zero_grad(set_to_none=True)
        logits = model(ctx['xt'], ctx['gt']); loss = criterion(logits[ix], yy)
        assert bool(torch.isfinite(loss)), 'Non-finite loss'
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 5.); optimizer.step()
        if validation_mask is not None and epoch % c['check_every'] == 0:
            model.eval()
            with torch.no_grad():
                score = model(ctx['xi'], ctx['gi']).sigmoid().cpu().numpy()[validation_mask]
            ap = float(average_precision_score(validation_labels, score))
            history.append(dict(epoch=epoch, AP=ap, loss=float(loss.detach().cpu())))
            if ap > best_score + 1e-10:
                best_score = ap; best_epoch = epoch; bad = 0
            else:
                bad += 1
            if epoch >= c['minimum_epochs'] and bad >= c['patience_checks']:
                break
    if validation_mask is not None:
        return dict(AP=best_score, best_epoch=best_epoch, last_epoch=epoch, history=history)
    model.eval()
    with torch.no_grad():
        logits, relation_gates, bg_gate = model(ctx['xi'], ctx['gi'], return_gates=True)
        score = logits.sigmoid().cpu().numpy()
        relation_values = np.stack([g.cpu().numpy() for g in relation_gates], axis=1)
    return dict(score=score, state={k: v.detach().cpu() for k, v in model.state_dict().items()},
                gates=relation_values, background_gates=bg_gate.cpu().numpy(), epochs=limit)

def tree_fit(ctx, mask, labels, params, seed):
    from sklearn.ensemble import ExtraTreesClassifier
    c = config(); model = ExtraTreesClassifier(n_estimators=c['ExtraTrees_estimators'],
        class_weight='balanced', random_state=seed, n_jobs=c['threads'], **params)
    model.fit(ctx['xt'][mask], labels)
    return model

def train_job(grouping, fold, condition):
    verify_freeze(); c = config(); torch.set_num_threads(c['threads'])
    assert condition in c['conditions' + str(grouping)]
    m, raw, graph = load(); bgraw = raw_for_indices(m.gene_index.to_numpy(int))
    path = directory(grouping, fold, condition); target(path / 'progress.json')
    if (path / 'complete.json').exists():
        print('Already complete: ' + str(path), flush=True); return
    outer = m['outer' + str(grouping)].to_numpy(int)
    inner = m['inner%d_outer%d' % (grouping, fold)].to_numpy(int)
    y = m.primary.to_numpy(bool); eligible = y | ~m.excluded_U.to_numpy(bool)
    allowed = outer != fold; test = ~allowed
    is_tree = condition == 'ExtraTrees_background'
    device = torch.device('cpu' if is_tree else 'cuda:0')
    assert is_tree or torch.cuda.is_available()
    grid = c['ExtraTrees_grid'] if is_tree else c['grid']
    total = len(grid) * c['inner_folds'] * len(c['seeds']); rows = []; started = time.time()
    report_progress(path, 'inner_tuning', completed_inner=0, total_inner=total, device=str(device))
    for inner_fold in range(c['inner_folds']):
        train_mask = allowed & (inner != inner_fold); val_mask = allowed & (inner == inner_fold)
        loss_mask = train_mask & eligible; val_eval = val_mask & eligible
        group_col = 'group' + str(grouping)
        assert not set(m.loc[train_mask, group_col]) & set(m.loc[val_mask | test, group_col])
        context_seeds = c['seeds'] if condition == 'gcn_rewired' else [c['seeds'][0]]
        key = 'family%d_outer%d_inner%d' % (grouping, fold, inner_fold)
        contexts = {seed: context(raw, graph, train_mask, allowed, loss_mask, condition, seed, key, device, bgraw)
                    for seed in context_seeds}
        for config_index, params in enumerate(grid):
            for seed in c['seeds']:
                ctx = contexts[seed if condition == 'gcn_rewired' else c['seeds'][0]]
                part = path / 'inner' / ('fold%d_config%d_seed%d.json' % (inner_fold, config_index, seed))
                if part.exists():
                    result = json.loads(part.read_text())
                elif is_tree:
                    fit = tree_fit(ctx, loss_mask, y[loss_mask], params, seed)
                    score = fit.predict_proba(ctx['xi'][val_eval])[:, 1]
                    result = dict(AP=float(average_precision_score(y[val_eval], score))); del fit
                else:
                    result = neural_fit(ctx, loss_mask, y[loss_mask], params, seed, val_eval, y[val_eval])
                record = {**result, 'inner_fold': inner_fold, 'config_index': config_index, 'seed': seed, 'params': params}
                save(part, record); rows.append(record)
                report_progress(path, 'inner_tuning', completed_inner=len(rows), total_inner=total,
                                inner_fold=inner_fold, config_index=config_index, seed=seed,
                                elapsed_seconds=round(time.time() - started, 1))
        del contexts; gc.collect()
        if not is_tree:
            torch.cuda.empty_cache()
    means = [float(np.mean([r['AP'] for r in rows if r['config_index'] == i])) for i in range(len(grid))]
    chosen = int(np.argmax(means)); params = grid[chosen]
    epochs = None if is_tree else max(c['minimum_epochs'], int(round(np.median(
        [r['best_epoch'] for r in rows if r['config_index'] == chosen]))))
    selection = dict(chosen_index=chosen, params=params, inner_AP_means=means, refit_epochs=epochs,
                     inner_evaluations=len(rows), selected_by='inner_AP_only')
    save(path / 'selection.json', selection)
    loss_mask = allowed & eligible; infer_allowed = np.ones(len(m), bool); scores = []; gate_records = []
    for seed in c['seeds']:
        report_progress(path, 'outer_refit', seed=seed, completed_inner=len(rows), total_inner=total)
        ctx = context(raw, graph, allowed, infer_allowed, loss_mask, condition, seed,
                      'family%d_outer%d_refit' % (grouping, fold), device, bgraw)
        seedpath = path / ('seed' + str(seed)); target(seedpath / 'metadata.json')
        joblib.dump(ctx['preprocessor'], target(seedpath / 'preprocessor.joblib'))
        if is_tree:
            fit = tree_fit(ctx, loss_mask, y[loss_mask], params, seed)
            score = fit.predict_proba(ctx['xi'])[:, 1]
            joblib.dump(fit, target(seedpath / 'model.joblib'), compress=3); del fit
        else:
            fit = neural_fit(ctx, loss_mask, y[loss_mask], params, seed, epochs=epochs); score = fit['score']
            torch.save(dict(state_dict=fit['state'], input_dim=191, condition=condition,
                params=params, epochs=epochs, seed=seed, architecture='BackgroundRelGCN64x2'), target(seedpath / 'model.pt'))
            np.savez_compressed(target(seedpath / 'outer_gates.npz'), gene_id=m.gene_id.to_numpy(str)[test],
                                relation_gates=fit['gates'][test], background_gates=fit['background_gates'][test])
            gate_records.append(dict(seed=seed, mean_outer_relation_gates=fit['gates'][test].mean(0).tolist(),
                                     mean_outer_background_gates=fit['background_gates'][test].mean(0).tolist()))
            del fit
        scores.append(score[test]); save(seedpath / 'metadata.json', dict(seed=seed, completed=now()))
        del ctx; gc.collect()
        if not is_tree:
            torch.cuda.empty_cache()
    result = m.loc[test, ['gene_id', 'row_index', 'primary', 'excluded_U', 'ppr', 'tf', 'group80', 'group50']].copy()
    for seed, score in zip(c['seeds'], scores):
        result['score_seed' + str(seed)] = score
    result['score'] = np.mean(scores, axis=0); result['uncertainty_sd'] = np.std(scores, axis=0, ddof=1)
    result.to_csv(target(path / 'outer_predictions.tsv.gz'), sep='\t', index=False, compression='gzip')
    save(path / 'gate_summary.json', gate_records)
    save(path / 'complete.json', dict(completed=now(), grouping=grouping, fold=fold, condition=condition,
        outer_predictions=int(test.sum()), selected=selection, elapsed_seconds=round(time.time() - started, 1)))
    report_progress(path, 'complete', completed_inner=len(rows), total_inner=total)

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--grouping', type=int, required=True, choices=[80, 50])
    p.add_argument('--fold', type=int, required=True, choices=list(range(5)))
    p.add_argument('--condition', required=True, choices=config()['conditions80'])
    a = p.parse_args(); train_job(a.grouping, a.fold, a.condition)
