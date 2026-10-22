"""Exact requested baselines, nested selection, no outer-driven model updates."""
import argparse
import gc
import json
import time
import warnings
import joblib
import numpy as np
import torch
from torch import nn
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score
from threadpoolctl import threadpool_limits
from benchmark_common import (ROOT, CONFIG, output, save_json, status, load, context,
    graphops, merged_graph, normalized_tensor, write_predictions, verify_training, source_model)


def classical_factory(name, params, seed, x):
    if name == 'LogisticRegression':
        from sklearn.linear_model import LogisticRegression
        return LogisticRegression(class_weight='balanced', solver='lbfgs', max_iter=5000,
                                  tol=1e-4, random_state=seed, **params)
    if name == 'RBFSVM':
        from sklearn.svm import SVC
        gamma = params['gamma_factor'] / (x.shape[1] * float(x.var()))
        return SVC(kernel='rbf', C=params['C'], gamma=gamma, class_weight='balanced',
                   cache_size=CONFIG['svm_cache_size_MB'], probability=False, random_state=seed)
    assert name == 'GradientBoosting'
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(class_weight='balanced', learning_rate=.1,
        max_iter=CONFIG['gradient_boosting_max_iter'], min_samples_leaf=20,
        early_stopping=False, random_state=seed, **params)


def classical_score(model, x):
    if hasattr(model, 'decision_function'):
        return model.decision_function(x)
    return model.predict_proba(x)[:, 1]


class NeuralBaseline(nn.Module):
    def __init__(self, condition, dropout):
        super().__init__()
        self.condition = condition
        self.dropout = nn.Dropout(dropout)
        self.fc1 = nn.Linear(191, CONFIG['hidden'])
        if condition == 'MLP':
            self.fc2 = nn.Linear(CONFIG['hidden'], CONFIG['hidden'])
            self.head = nn.Linear(CONFIG['hidden'], 1)
        else:
            assert condition == 'GCN'
            self.fc2 = nn.Linear(CONFIG['hidden'], 1)

    def forward(self, x, graph=None):
        if self.condition == 'MLP':
            h = self.dropout(torch.relu(self.fc1(x)))
            h = self.dropout(torch.relu(self.fc2(h)))
            return self.head(h).squeeze(-1)
        # x is precomputed S@X; biases are added after propagation, as in GCN.
        h = self.dropout(torch.relu(self.fc1(x)))
        return self.fc2(torch.sparse.mm(graph, h)).squeeze(-1)


def neural_fit(ctx, labels, name, params, seed, validation=None, epochs=None):
    source_model.seed_all(seed)
    device = ctx['xt'].device
    model = NeuralBaseline(name, params['dropout']).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=params['lr'], weight_decay=CONFIG['weight_decay'])
    ix = torch.as_tensor(np.flatnonzero(ctx['loss']), device=device)
    yy = torch.as_tensor(labels[ctx['loss']].astype('float32'), device=device)
    pw = torch.tensor((len(yy) - yy.sum().item()) / yy.sum().item(), device=device)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pw)
    history = []
    best_ap, best_epoch, bad = -np.inf, 0, 0
    maximum = CONFIG['max_epochs'] if epochs is None else int(epochs)
    for epoch in range(1, maximum + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        logits = model(ctx['xt'], ctx.get('at'))
        loss = criterion(logits[ix], yy)
        assert bool(torch.isfinite(loss))
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 5.)
        optimizer.step()
        if validation is not None and epoch % CONFIG['check_every'] == 0:
            model.eval()
            with torch.no_grad():
                score = model(ctx['xi'], ctx.get('ai')).sigmoid().cpu().numpy()[validation]
            ap = float(average_precision_score(labels[validation], score))
            history.append(dict(epoch=epoch, AP=ap, loss=float(loss.detach().cpu())))
            if ap > best_ap + 1e-10:
                best_ap, best_epoch, bad = ap, epoch, 0
            else:
                bad += 1
            if epoch >= CONFIG['minimum_epochs'] and bad >= CONFIG['patience_checks']:
                break
    if validation is not None:
        return dict(AP=best_ap, best_epoch=best_epoch, last_epoch=epoch, history=history)
    model.eval()
    with torch.no_grad():
        score = model(ctx['xi'], ctx.get('ai')).sigmoid().cpu().numpy()
    return dict(score=score, state_dict={k: v.detach().cpu() for k, v in model.state_dict().items()},
                epochs=maximum)


def tensor_context(ctx, condition):
    device = torch.device('cuda:0')
    n = len(ctx['xt'])
    ctx['xt'] = torch.as_tensor(ctx['xt'], device=device)
    ctx['xi'] = torch.as_tensor(ctx['xi'], device=device)
    if condition == 'GCN':
        ctx['at'] = normalized_tensor(ctx['gt'], n, device)
        ctx['ai'] = normalized_tensor(ctx['gi'], n, device)
        # Raw input propagation is fixed and can be reused across neural fits.
        ctx['xt'] = torch.sparse.mm(ctx['at'], ctx['xt'])
        ctx['xi'] = torch.sparse.mm(ctx['ai'], ctx['xi'])
    return ctx


def ppr(transition, seeds, restart):
    n = transition.shape[0]
    assert len(seeds) == n and seeds.dtype == bool and seeds.any()
    s = seeds.astype(float) / seeds.sum()
    h = s.copy()
    dangling = np.asarray(transition.sum(axis=1)).ravel() == 0.
    for iteration in range(1, CONFIG['ppr_max_iter'] + 1):
        z = restart * s + (1. - restart) * (transition.T @ h + h[dangling].sum() * s)
        delta = float(np.abs(z - h).sum())
        h = z
        if delta < CONFIG['ppr_tolerance']:
            break
    assert delta < CONFIG['ppr_tolerance'], 'PageRank failed to converge'
    assert np.isfinite(h).all() and h.min() >= 0 and abs(h.sum() - 1.) < 1e-7
    return h, dict(iterations=iteration, L1_error=delta, seed_P=int(seeds.sum()))


def transition(graph, n):
    from scipy import sparse
    a = merged_graph(graph, n)
    degree = np.asarray(a.sum(axis=1)).ravel()
    inv = np.divide(1., degree, out=np.zeros(n), where=degree > 0)
    return (sparse.diags(inv) @ a).tocsr()


def run_classical(fold):
    worker = 'classical_outer' + str(fold)
    m, raw, graph, bg = load()
    y = m.primary.to_numpy(bool)
    eligible = y | ~m.excluded_U.to_numpy(bool)
    allowed = m.outer80.to_numpy(int) != fold
    test = ~allowed
    inner = m['inner80_outer' + str(fold)].to_numpy(int)
    start = time.time()
    with threadpool_limits(limits=CONFIG['threads_per_worker']):
        for name in ['LogisticRegression', 'GradientBoosting', 'RBFSVM']:
            path = ROOT / 'runs/family80' / ('outer' + str(fold)) / name
            if (path / 'complete.json').exists():
                continue
            grid = CONFIG[name + '_grid']
            records = []
            for inner_fold in range(4):
                training = allowed & (inner != inner_fold)
                validation = allowed & (inner == inner_fold) & eligible
                assert not set(m.loc[training, 'group80']) & set(m.loc[validation | test, 'group80'])
                ctx = context(m, raw, graph, bg, training, allowed,
                              'outer%d_inner%d' % (fold, inner_fold), scaled=name != 'GradientBoosting')
                trainx = ctx['xt'][ctx['loss']]
                for index, params in enumerate(grid):
                    part = path / 'inner' / ('fold%d_config%d.json' % (inner_fold, index))
                    if part.exists():
                        record = json.loads(part.read_text())
                    else:
                        status(worker, 'inner_fit', condition=name, fold=fold, inner_fold=inner_fold,
                               config_index=index, elapsed_seconds=time.time() - start)
                        fitstart = time.time()
                        fit = classical_factory(name, params, CONFIG['classical_inner_seed'], trainx)
                        fit.fit(trainx, y[ctx['loss']])
                        score = classical_score(fit, ctx['xi'][validation])
                        record = dict(inner_fold=inner_fold, config_index=index, params=params,
                             AP=float(average_precision_score(y[validation], score)), fit_seconds=time.time() - fitstart)
                        save_json(part, record)
                        del fit
                    records.append(record)
                del ctx, trainx
                gc.collect()
            means = [float(np.mean([r['AP'] for r in records if r['config_index'] == i])) for i in range(len(grid))]
            chosen = int(np.argmax(means))
            selection = dict(params=grid[chosen], chosen_index=chosen, inner_AP_means=means,
                   selected_by='inner_AP_only', inner_evaluations=len(records), deterministic=name != 'GradientBoosting')
            save_json(path / 'selection.json', selection)
            ctx = context(m, raw, graph, bg, allowed, np.ones(len(m), bool), 'outer%d_refit' % fold,
                          scaled=name != 'GradientBoosting')
            trainx = ctx['xt'][ctx['loss']]
            trainy = y[ctx['loss']]
            scores = []
            seeds = CONFIG['seeds'] if name == 'GradientBoosting' else [CONFIG['seeds'][0]]
            for seed in seeds:
                scorepath = path / ('seed' + str(seed)) / 'outer_scores.npy'
                if scorepath.exists():
                    score = np.load(scorepath)
                else:
                    status(worker, 'outer_fit', condition=name, fold=fold, seed=seed,
                           elapsed_seconds=time.time() - start)
                    fit = classical_factory(name, grid[chosen], seed, trainx)
                    fit.fit(trainx, trainy)
                    score = classical_score(fit, ctx['xi'][test])
                    np.save(output(scorepath), score)
                    joblib.dump(fit, output(scorepath.with_name('model.joblib')), compress=3)
                    del fit
                scores.append(score)
            if len(scores) == 1:
                scores = scores * len(CONFIG['seeds'])
            write_predictions(m, fold, name, scores, dict(selection=selection,
                    unique_outer_fits=len(seeds), deterministic_columns_repeated=len(seeds) == 1))
            del ctx, trainx, trainy
            gc.collect()
    status(worker, 'complete', elapsed_seconds=time.time() - start)


def run_neural():
    worker = 'neural'
    torch.set_num_threads(CONFIG['threads_per_worker'])
    assert torch.cuda.is_available()
    m, raw, graph, bg = load()
    y = m.primary.to_numpy(bool)
    eligible = y | ~m.excluded_U.to_numpy(bool)
    start = time.time()
    for fold in range(5):
        allowed = m.outer80.to_numpy(int) != fold
        test = ~allowed
        inner = m['inner80_outer' + str(fold)].to_numpy(int)
        for name in ['MLP', 'GCN']:
            path = ROOT / 'runs/family80' / ('outer' + str(fold)) / name
            if (path / 'complete.json').exists():
                continue
            records = []
            grid = CONFIG['neural_grid']
            for inner_fold in range(4):
                training = allowed & (inner != inner_fold)
                validation = allowed & (inner == inner_fold) & eligible
                assert not set(m.loc[training, 'group80']) & set(m.loc[validation | test, 'group80'])
                ctx = tensor_context(context(m, raw, graph, bg, training, allowed,
                         'outer%d_inner%d' % (fold, inner_fold)), name)
                for index, params in enumerate(grid):
                    for seed in CONFIG['seeds']:
                        part = path / 'inner' / ('fold%d_config%d_seed%d.json' % (inner_fold, index, seed))
                        if part.exists():
                            record = json.loads(part.read_text())
                        else:
                            status(worker, 'inner_fit', condition=name, fold=fold, inner_fold=inner_fold,
                                   config_index=index, seed=seed, elapsed_seconds=time.time() - start)
                            result = neural_fit(ctx, y, name, params, seed, validation)
                            record = dict(**result, inner_fold=inner_fold, config_index=index, seed=seed, params=params)
                            save_json(part, record)
                        records.append(record)
                del ctx
                gc.collect()
                torch.cuda.empty_cache()
            means = [float(np.mean([r['AP'] for r in records if r['config_index'] == i])) for i in range(len(grid))]
            chosen = int(np.argmax(means))
            epochs = max(CONFIG['minimum_epochs'], int(round(np.median(
                     [r['best_epoch'] for r in records if r['config_index'] == chosen]))))
            selection = dict(params=grid[chosen], chosen_index=chosen, inner_AP_means=means,
                  refit_epochs=epochs, selected_by='inner_AP_only', inner_evaluations=len(records))
            save_json(path / 'selection.json', selection)
            ctx = tensor_context(context(m, raw, graph, bg, allowed, np.ones(len(m), bool),
                    'outer%d_refit' % fold), name)
            scores = []
            for seed in CONFIG['seeds']:
                scorepath = path / ('seed' + str(seed)) / 'outer_scores.npy'
                if scorepath.exists():
                    score = np.load(scorepath)
                else:
                    status(worker, 'outer_fit', condition=name, fold=fold, seed=seed,
                           epochs=epochs, elapsed_seconds=time.time() - start)
                    result = neural_fit(ctx, y, name, grid[chosen], seed, epochs=epochs)
                    score = result['score'][test]
                    np.save(output(scorepath), score)
                    torch.save(dict(**{k: v for k, v in result.items() if k != 'score'},
                         condition=name, params=grid[chosen], seed=seed, architecture='CanonicalGCN2' if name == 'GCN' else 'MLP64x2'),
                         output(scorepath.with_name('model.pt')))
                    del result
                scores.append(score)
            write_predictions(m, fold, name, scores, dict(selection=selection, unique_outer_fits=3))
            del ctx
            gc.collect()
            torch.cuda.empty_cache()
    status(worker, 'complete', elapsed_seconds=time.time() - start)


def run_ppr():
    worker = 'ppr'
    m, _, graph, _ = load()
    y = m.primary.to_numpy(bool)
    eligible = y | ~m.excluded_U.to_numpy(bool)
    for fold in range(5):
        allowed = m.outer80.to_numpy(int) != fold
        test = ~allowed
        inner = m['inner80_outer' + str(fold)].to_numpy(int)
        path = ROOT / 'runs/family80' / ('outer' + str(fold)) / 'PPR'
        if (path / 'complete.json').exists():
            continue
        t = transition(graphops.induced(graph, allowed), len(m))
        records = []
        for inner_fold in range(4):
            training = allowed & (inner != inner_fold)
            validation = allowed & (inner == inner_fold) & eligible
            seeds = training & y
            assert not seeds[validation | test].any()
            assert not set(m.loc[training, 'group80']) & set(m.loc[validation | test, 'group80'])
            for index, params in enumerate(CONFIG['PPR_grid']):
                status(worker, 'inner_fit', fold=fold, inner_fold=inner_fold, config_index=index)
                score, audit = ppr(t, seeds, params['restart'])
                record = dict(inner_fold=inner_fold, config_index=index, params=params,
                              AP=float(average_precision_score(y[validation], score[validation])), **audit)
                save_json(path / 'inner' / ('fold%d_config%d.json' % (inner_fold, index)), record)
                records.append(record)
        means = [float(np.mean([r['AP'] for r in records if r['config_index'] == i])) for i in range(len(CONFIG['PPR_grid']))]
        chosen = int(np.argmax(means))
        selection = dict(params=CONFIG['PPR_grid'][chosen], chosen_index=chosen, inner_AP_means=means,
                         selected_by='inner_AP_only', inner_evaluations=len(records), deterministic=True)
        save_json(path / 'selection.json', selection)
        seed_mask = allowed & y
        assert not seed_mask[test].any()
        score, audit = ppr(transition(graph, len(m)), seed_mask, selection['params']['restart'])
        np.save(output(path / 'train_positive_seed_mask.npy'), seed_mask)
        np.save(output(path / 'full_ppr_scores.npy'), score)
        write_predictions(m, fold, 'PPR', [score[test]] * len(CONFIG['seeds']),
             dict(selection=selection, unique_outer_fits=1, deterministic_columns_repeated=True,
                  graph_only=True, outer_test_seed_overlap=0, convergence=audit))
    status(worker, 'complete')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--job', choices=['classical', 'neural', 'ppr'], required=True)
    parser.add_argument('--fold', type=int, choices=range(5))
    args = parser.parse_args()
    verify_training()
    warnings.simplefilter('error', ConvergenceWarning)
    if args.job == 'classical':
        assert args.fold is not None
        run_classical(args.fold)
    elif args.job == 'neural':
        run_neural()
    else:
        run_ppr()
