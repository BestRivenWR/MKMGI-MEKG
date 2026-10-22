"""Small residual relational GCN, with a parameter-matched no-message control."""
import random
import numpy as np
import torch
from torch import nn
from common import config

def seed_all(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False

class RelationBlock(nn.Module):
    def __init__(self, hidden, dropout, no_message=False):
        super().__init__()
        self.transforms = nn.ModuleList([nn.Linear(hidden, hidden, bias=False) for _ in range(4)])
        self.gate = nn.Linear(hidden, 4)
        self.norm = nn.LayerNorm(hidden)
        self.drop = nn.Dropout(dropout); self.no_message = no_message

    def forward(self, h, matrices, available):
        gate = torch.softmax(self.gate(h).masked_fill(~available, -1e9), dim=1)
        own = self.transforms[0](h)
        message = gate[:, :1] * own
        for i, mat in enumerate(matrices):
            z = self.transforms[i + 1](h)
            if not self.no_message:
                z = torch.sparse.mm(mat, z)
            message = message + gate[:, i + 1:i + 2] * z
        return self.drop(torch.nn.functional.gelu(self.norm(message))), gate

class RelGCN(nn.Module):
    def __init__(self, input_dim, dropout, no_message=False):
        super().__init__(); c = config(); width = c['hidden']
        self.project = nn.Sequential(nn.Linear(input_dim, width), nn.LayerNorm(width), nn.GELU())
        self.blocks = nn.ModuleList([RelationBlock(width, dropout, no_message) for _ in range(c['layers'])])
        self.head = nn.Linear(width, 1); self.alpha = c['residual_alpha']

    def forward(self, x, graph, return_gates=False):
        h0 = self.project(x); h = h0; gates = []
        for block in self.blocks:
            z, gate = block(h, *graph); h = self.alpha * h0 + (1 - self.alpha) * (h + z) / 2
            gates.append(gate)
        score = self.head(h).squeeze(-1)
        return (score, gates) if return_gates else score

def fit_preprocess(raw, train_mask):
    from sklearn.preprocessing import StandardScaler
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        median = np.nanmedian(raw[train_mask], axis=0)
    median = np.nan_to_num(median, nan=0.).astype('float32')
    x = np.where(np.isfinite(raw), raw, median).astype('float32')
    scaler = StandardScaler().fit(x[train_mask])
    return dict(median=median, mean=scaler.mean_.astype('float32'), scale=scaler.scale_.astype('float32'))

def apply_preprocess(raw, fitted):
    x = np.where(np.isfinite(raw), raw, fitted['median']).astype('float32')
    return ((x - fitted['mean']) / fitted['scale']).astype('float32')
