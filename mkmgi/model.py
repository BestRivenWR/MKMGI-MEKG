"""8-DAP Transformer and temporal pooling before unchanged residual relation GCN."""
import torch
from torch import nn
from common import config,core_module
_core=core_module('model')
seed_all=_core.seed_all
fit_preprocess=_core.fit_preprocess
apply_preprocess=_core.apply_preprocess

class TimeEncoder(nn.Module):
    def __init__(self,condition,dropout):
        super().__init__(); c=config(); w=c['time_width']; self.condition=condition
        self.token=nn.Linear(10,w); self.position=nn.Linear(1,w,bias=False)
        self.norm1=nn.LayerNorm(w); self.norm2=nn.LayerNorm(w)
        self.self_attention=nn.MultiheadAttention(w,c['time_heads'],dropout=dropout,batch_first=True)
        self.feedforward=nn.Sequential(nn.Linear(w,c['time_feedforward']),nn.GELU(),nn.Dropout(dropout),nn.Linear(c['time_feedforward'],w))
        self.drop=nn.Dropout(dropout); self.output_norm=nn.LayerNorm(w)
        self.pool=nn.Sequential(nn.Linear(w,16),nn.Tanh(),nn.Linear(16,1))
        dap=torch.tensor(c['DAP'],dtype=torch.float32); self.register_buffer('dap',(dap-dap.min())/(dap.max()-dap.min()))

    def forward(self,values,observed):
        present=observed.any(-1); available=present.clone(); empty=~available.any(1)
        available[empty,0]=True
        values=values.masked_fill(~observed,0.)
        tokens=self.token(torch.cat([values,observed.to(values.dtype)],-1))
        if self.condition!='time_no_pos':
            tokens=tokens+self.position(self.dap[:,None])[None,:,:]
        q=self.norm1(tokens)
        z,_=self.self_attention(q,q,q,key_padding_mask=~available,need_weights=False)
        tokens=tokens+self.drop(z)
        tokens=self.output_norm(tokens+self.drop(self.feedforward(self.norm2(tokens))))
        if self.condition=='time_mean':
            weight=available.to(values.dtype)/available.sum(1,keepdim=True)
        else:
            logits=self.pool(tokens).squeeze(-1).masked_fill(~available,-1e9)
            weight=torch.softmax(logits,1)
        weight=weight.masked_fill(empty[:,None],0.)
        pooled=(tokens*weight[:,:,None]).sum(1)
        return pooled,weight

class BackgroundGCN(nn.Module):
    def __init__(self,condition,dropout):
        super().__init__(); c=config(); w=c['hidden']; self.condition=condition; self.alpha=c['residual_alpha']
        self.time=TimeEncoder(condition,dropout)
        self.project=nn.Sequential(nn.Linear(143,w),nn.LayerNorm(w),nn.GELU())
        self.blocks=nn.ModuleList([_core.RelationBlock(w,dropout,condition=='no_message') for _ in range(c['layers'])])
        self.head=nn.Linear(w,1); self.chunk=c['time_chunk_size']

    def encode(self,x):
        values=x[:,:40].reshape(-1,5,8).transpose(1,2)
        observed=(x[:,40:80].reshape(-1,5,8).transpose(1,2)>.5)
        vectors=[]; weights=[]
        for start in range(0,len(x),self.chunk):
            z,a=self.time(values[start:start+self.chunk],observed[start:start+self.chunk])
            vectors.append(z);weights.append(a)
        pooled=torch.cat(vectors); attention=torch.cat(weights)
        return self.project(torch.cat([pooled,x[:,80:]],1)),attention

    def forward(self,x,graph,return_gates=False):
        h0,time_weights=self.encode(x); h=h0; gates=[]
        for block in self.blocks:
            z,gate=block(h,*graph); h=self.alpha*h0+(1-self.alpha)*(h+z)/2;gates.append(gate)
        score=self.head(h).squeeze(-1)
        return (score,gates,time_weights) if return_gates else score
