"""Final registered candidate fit; no sealed labels, no model/site replacement."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import joblib
from common import ROOT, config, load, target, save, verify_freeze, now
from train import directory, context, neural_fit, report_progress
from background import load_full_raw

def main():
    verify_freeze(check_sources=True); c = config(); torch.set_num_threads(c['threads'])
    assert (ROOT / 'reports/decision.json').exists()
    dev, rawdev, _ = load(); old = Path(c['previous_release']); condition = c['registered_main']
    other = pd.read_csv(old / 'sealed/registry.tsv', sep='\t', usecols=['gene_id', 'gene_index'])
    with np.load(old / 'sealed/features.npz') as z:
        rawother = z['main'].astype('float32')
    with np.load(c['graph_source']) as z:
        ids = z['gene_ids'].astype(str)
        graph = {k: z[k].copy() for k in ['ppi_edge_index', 'ppi_base_weight', 'reg_edge_index', 'reg_base_weight']}
    n = len(ids); raw = np.full((n, 92), np.nan, np.float32); di = dev.gene_index.to_numpy(int); oi = other.gene_index.to_numpy(int)
    assert np.array_equal(ids[di], dev.gene_id.to_numpy(str)) and np.array_equal(ids[oi], other.gene_id.to_numpy(str))
    assert len(np.unique(np.r_[di, oi])) == n and not np.intersect1d(di, oi).size
    raw[di] = rawdev; raw[oi] = rawother; bgraw = load_full_raw(); assert np.array_equal(ids, bgraw['ids'])
    allowed = np.zeros(n, bool); allowed[di] = True
    eligible = dev.primary.to_numpy(bool) | ~dev.excluded_U.to_numpy(bool)
    loss_mask = np.zeros(n, bool); loss_mask[di[eligible]] = True
    y = np.zeros(n, bool); y[di] = dev.primary.to_numpy(bool)
    selections = [json.loads((directory(80, fold, condition) / 'selection.json').read_text()) for fold in range(5)]
    chosen = int(np.argmax(np.mean([r['inner_AP_means'] for r in selections], axis=0)))
    params = c['grid'][chosen]; epochs = []
    for fold in range(5):
        for p in (directory(80, fold, condition) / 'inner').glob('*.json'):
            r = json.loads(p.read_text())
            if r['config_index'] == chosen:
                epochs.append(r['best_epoch'])
    epochs = max(c['minimum_epochs'], int(round(np.median(epochs))))
    path = ROOT / 'models/deployment' / condition; scores = []; gates = []; bg_gates = []
    for seed in c['seeds']:
        report_progress(path, 'full_DEV_refit', seed=seed, epochs=epochs)
        ctx = context(raw, graph, allowed, np.ones(n, bool), loss_mask, condition, seed, 'deployment', torch.device('cuda:0'), bgraw)
        fit = neural_fit(ctx, loss_mask, y[loss_mask], params, seed, epochs=epochs)
        p = path / ('seed' + str(seed))
        torch.save(dict(state_dict=fit['state'], input_dim=191, condition=condition, params=params,
                        epochs=epochs, seed=seed, architecture='BackgroundRelGCN64x2'), target(p / 'model.pt'))
        joblib.dump(ctx['preprocessor'], target(p / 'preprocessor.joblib'))
        scores.append(fit['score']); gates.append(fit['gates']); bg_gates.append(fit['background_gates'])
        del ctx, fit; torch.cuda.empty_cache()
    stack = np.stack(scores); mean = stack.mean(0)
    out = pd.DataFrame(dict(gene_id=ids, priority_score=mean,
        priority_percentile=(pd.Series(mean).rank(method='average').to_numpy() - 1) / (n - 1),
        uncertainty_sd=stack.std(0, ddof=1), runs=len(c['seeds']), in_development_training=allowed,
        in_development_positive=y, temporal_modalities_observed=(raw[:, 40:80].reshape(n, 5, 8).sum(2) > 0).sum(1),
        spatial_observed=raw[:, 91] > .5, sequence_observed=bgraw['seq_present'], protein_observed=bgraw['protein_present']))
    for seed, score in zip(c['seeds'], scores):
        out['score_seed' + str(seed)] = score
    out.sort_values('priority_score', ascending=False).to_csv(target(ROOT / 'predictions/genomewide_background_gcn_consensus.tsv.gz'),
        sep='\t', index=False, compression='gzip')
    np.savez_compressed(target(ROOT / 'predictions/genomewide_mean_gates.npz'), gene_id=ids,
                        relation_gates=np.mean(gates, axis=0), background_gates=np.mean(bg_gates, axis=0))
    # Candidate composition is descriptive; no sealed outcomes are evaluated or used for fitting/selection.
    partition = pd.read_csv(old / 'source_data/universe_partition.tsv', sep='\t', usecols=['gene_id', 'ppr', 'excluded_U'])
    from common import boolcol
    partition = partition.set_index('gene_id').loc[ids].reset_index()
    take = np.flatnonzero(~boolcol(partition.excluded_U)); take = take[np.argsort(-mean[take], kind='mergesort')][:100]
    from report import longest_run
    flags = boolcol(partition.ppr)[take]
    top = out.iloc[take].copy(); top['ppr'] = flags; top['candidate_rank'] = np.arange(1, len(take) + 1)
    top.to_csv(target(ROOT / 'source_data/deployment_eligible_U_top100.tsv'), sep='\t', index=False)
    save(ROOT / 'models/deployment/manifest.json', dict(completed=now(), N=n, architecture='BackgroundRelGCN64x2',
        registered_candidate=condition, training_genes=int(allowed.sum()), training_P=int(y.sum()),
        selection='aggregate inner AP only; no outer winner selection', params=params, epochs=epochs,
        priority_not_calibrated_probability=True, PPR_top100=int(flags.sum()), PPR_longest_run=longest_run(flags),
        changed_original_model=False, website_deployed=False, graph=c['graph_source'], independent_validation=False))
    verify_freeze(check_sources=True); report_progress(path, 'complete', genes=n)

if __name__ == '__main__':
    main()
