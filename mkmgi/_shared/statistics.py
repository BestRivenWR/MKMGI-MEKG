"""Tie-aware weighted metrics, paired within-outer-fold family bootstrap."""
import numpy as np
from common import finite_random_ap

def metrics(y, score):
    from sklearn.metrics import average_precision_score, roc_auc_score
    y = np.asarray(y, bool); score = np.asarray(score, float); n = len(y); p = int(y.sum())
    assert n and 0 < p < n and np.isfinite(score).all()
    r = dict(N=n, P=p, AP=float(average_precision_score(y, score)), AUROC=float(roc_auc_score(y, score)),
             prevalence=p/n, finite_random_AP=finite_random_ap(n, p))
    r['AP_prevalence_lift'] = r['AP'] / r['prevalence']
    order = np.argsort(-score, kind='mergesort')
    for k in [20, 50, 100]:
        k = min(k, n); r['precision_at_' + str(k)] = float(y[order[:k]].mean())
        r['known_P_at_' + str(k)] = int(y[order[:k]].sum())
    for pct in [1, 5]:
        k = max(1, int(np.ceil(n * pct / 100)))
        r['recall_top' + str(pct) + 'pct'] = float(y[order[:k]].sum() / p)
    return r

class WeightedScores:
    def __init__(self, y, score, mask=None):
        y = np.asarray(y, bool); s = np.asarray(score, float)
        self.ix = np.flatnonzero(np.ones(len(y), bool) if mask is None else mask)
        self.ix = self.ix[np.argsort(-s[self.ix], kind='mergesort')]
        self.y = y[self.ix].astype(float); ss = s[self.ix]
        self.ends = np.r_[np.flatnonzero(ss[:-1] != ss[1:]), len(ss) - 1] if len(ss) else np.array([], int)

    def calculate(self, weights):
        if not len(self.ix):
            return np.full((len(weights), 2), np.nan)
        w = weights[:, self.ix]; cp = np.cumsum(w * self.y, axis=1)[:, self.ends]
        cn = np.cumsum(w * (1 - self.y), axis=1)[:, self.ends]
        bp = np.diff(cp, prepend=np.zeros((len(w), 1)), axis=1)
        bn = np.diff(cn, prepend=np.zeros((len(w), 1)), axis=1)
        p = cp[:, -1]; n = cn[:, -1]
        with np.errstate(divide='ignore', invalid='ignore'):
            precision = np.divide(cp, cp + cn, out=np.zeros_like(cp), where=cp + cn > 0)
            ap = (bp * precision).sum(1) / p
            auc = (bp * (n[:, None] - cn + .5 * bn)).sum(1) / (p * n)
        ap[p == 0] = np.nan; auc[(p == 0) | (n == 0)] = np.nan
        return np.column_stack([ap, auc])

def interval(v, family=1):
    v = np.asarray(v); v = v[np.isfinite(v)]; a = .05 / family
    return [float(t) for t in np.quantile(v, [a / 2, 1 - a / 2])] if len(v) else [None, None]

def bootstrap(m, scores, grouping):
    from common import config
    c = config(); B = c['bootstrap_replicates']; y = m.primary.to_numpy(bool)
    masks = dict(all=np.ones(len(m), bool), nonPPR=~m.ppr.to_numpy(bool), PPR=m.ppr.to_numpy(bool), TF=m.tf.to_numpy(bool))
    result = {name: {key: np.zeros((B, 2)) for key in masks} for name in scores}
    valid = {name: {key: np.ones((B, 2), bool) for key in masks} for name in scores}
    rng = np.random.default_rng(c['bootstrap_seed'] + int(grouping))
    for fold in range(5):
        ix = np.flatnonzero(m['outer' + str(grouping)].to_numpy(int) == fold)
        _, inverse = np.unique(m['group' + str(grouping)].to_numpy(str)[ix], return_inverse=True)
        ng = int(inverse.max()) + 1
        packs = {name: {key: WeightedScores(y[ix], s[ix], mask[ix]) for key, mask in masks.items()}
                 for name, s in scores.items()}
        for start in range(0, B, 25):
            end = min(B, start + 25)
            multiplicity = rng.multinomial(ng, np.full(ng, 1 / ng), size=end - start)
            weights = multiplicity[:, inverse].astype(float)
            for name, packet in packs.items():
                for key, p in packet.items():
                    vals = p.calculate(weights)
                    valid[name][key][start:end] &= np.isfinite(vals)
                    result[name][key][start:end] += np.nan_to_num(vals) / 5
    for name in scores:
        for key in masks:
            result[name][key][~valid[name][key]] = np.nan
    return result
