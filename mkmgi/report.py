"""Full-grid paired comparison with unchanged candidate pools and baseline OOF."""
import numpy as np
import pandas as pd
from common import ROOT,config,load,target,save,verify_freeze,now,core_module
from train import directory
stats=core_module('statistics')

def longest_run(flags):
    best=current=0
    for f in flags:
        current=current+1 if f else 0;best=max(best,current)
    return int(best)

def main():
    verify_freeze(check_sources=True);c=config();m,_,_=load();eligible=m.primary.to_numpy(bool)|~m.excluded_U.to_numpy(bool)
    scores={};foldrows=[];composition=[];seedrows=[];pooled=[]
    references=['bg_concat','ExtraTrees_background'];conditions=references+c['conditions80']
    for condition in conditions:
        if condition in references:
            z=pd.read_csv(ROOT/'reference'/('family80_'+condition+'_OOF.tsv.gz'),sep='\t')
        else:
            frames=[]
            for fold in range(5):
                p=directory(80,fold,condition);assert (p/'complete.json').exists()
                a=pd.read_csv(p/'outer_predictions.tsv.gz',sep='\t')
                assert np.array_equal(a.gene_id.to_numpy(str),m.loc[m.outer80.eq(fold),'gene_id'].to_numpy(str));frames.append(a)
            z=pd.concat(frames)
        z=z.sort_values('row_index').reset_index(drop=True)
        assert np.array_equal(z.gene_id.to_numpy(str),m.gene_id.to_numpy(str))
        s=z.score.to_numpy(float);assert np.isfinite(s).all();scores[condition]=s
        z.to_csv(target(ROOT/'predictions'/('family80_'+condition+'_OOF.tsv.gz')),sep='\t',index=False,compression='gzip')
        for fold in range(5):
            take=eligible&m.outer80.eq(fold).to_numpy(bool)
            row=dict(condition=condition,fold=fold,**stats.metrics(m.primary.to_numpy(bool)[take],s[take]))
            for subset,mask in [('nonPPR',~m.ppr.to_numpy(bool)),('PPR',m.ppr.to_numpy(bool)),('TF',m.tf.to_numpy(bool))]:
                use=take&mask;yy=m.primary.to_numpy(bool)[use]
                if 0<yy.sum()<len(yy):
                    r=stats.metrics(yy,s[use]);row[subset+'_AP']=r['AP'];row[subset+'_AUROC']=r['AUROC']
            foldrows.append(row)
            for seed in c['seeds']:
                seedrows.append(dict(condition=condition,fold=fold,seed=seed,
                    **stats.metrics(m.primary.to_numpy(bool)[take],z['score_seed'+str(seed)].to_numpy(float)[take])))
        pooled.append(dict(condition=condition,**stats.metrics(m.primary.to_numpy(bool)[eligible],s[eligible])))
        u=np.flatnonzero(~m.excluded_U.to_numpy(bool));order=u[np.argsort(-s[u],kind='mergesort')]
        for k in [10,20,100]:
            ix=order[:k];flags=m.ppr.to_numpy(bool)[ix]
            composition.append(dict(condition=condition,scope='OOF_eligible_U',top_k=k,PPR_count=int(flags.sum()),
                PPR_fraction=float(flags.mean()),PPR_longest_run=longest_run(flags),
                PPR_positions=','.join(str(j+1) for j in np.flatnonzero(flags))))
        top=z.iloc[order[:100]].copy();top['candidate_rank']=np.arange(1,101)
        top.to_csv(target(ROOT/'source_data'/('family80_'+condition+'_OOF_unlabelled_top100.tsv')),sep='\t',index=False)
    fd=pd.DataFrame(foldrows);summary=fd.groupby('condition').agg(AP=('AP','mean'),AP_fold_SD=('AP','std'),
        AUROC=('AUROC','mean'),AUROC_fold_SD=('AUROC','std'),nonPPR_AP=('nonPPR_AP','mean'),
        PPR_AUROC=('PPR_AUROC','mean'),mean_precision_at_20=('precision_at_20','mean')).reset_index()
    boot=stats.bootstrap(m.loc[eligible].reset_index(drop=True),{k:v[eligible] for k,v in scores.items()},80)
    pairs=[('time_attention','bg_concat','time_encoder',2),('time_mean','bg_concat','time_encoder',2),
        ('time_attention','time_mean','time_components',2),('time_attention','time_no_pos','time_components',2),
        ('time_attention','no_message','graph',2),('time_attention','gcn_rewired','graph',2),
        ('time_attention','ExtraTrees_background','strong_baseline_reference',1)]
    comp=[];ci=[]
    for a,b,family,size in pairs:
        for subset in ['all','nonPPR','PPR','TF']:
            for idx,metric in enumerate(['AP','AUROC']):
                col=metric if subset=='all' else subset+'_'+metric
                point=fd.loc[fd.condition.eq(a),col].mean()-fd.loc[fd.condition.eq(b),col].mean()
                delta=boot[a][subset][:,idx]-boot[b][subset][:,idx]
                lo,hi=stats.interval(delta);al,ah=stats.interval(delta,size)
                comp.append(dict(main=a,comparator=b,subset=subset,metric=metric,delta=float(point),
                    conditional95_low=lo,conditional95_high=hi,family_adjusted_low=al,family_adjusted_high=ah,
                    family=family,family_size=size,valid_bootstrap=int(np.isfinite(delta).sum())))
    for condition in conditions:
        for subset in ['all','nonPPR','PPR','TF']:
            for idx,metric in enumerate(['AP','AUROC']):
                lo,hi=stats.interval(boot[condition][subset][:,idx]);ci.append(dict(condition=condition,subset=subset,
                    metric=metric,conditional95_low=lo,conditional95_high=hi))
    cd=pd.DataFrame(comp);cf=pd.DataFrame(composition)
    for name,frame in [('summary',summary),('per_fold_metrics',fd),('seed_metrics',pd.DataFrame(seedrows)),
        ('pooled_OOF_metrics',pd.DataFrame(pooled)),('paired_comparisons',cd),('absolute_intervals',pd.DataFrame(ci)),
        ('candidate_family_concentration',cf)]:
        frame.to_csv(target(ROOT/'source_data'/(name+'.tsv')),sep='\t',index=False)
    pools=pd.read_csv(ROOT/'source_data/matched_candidate_pools.tsv',sep='\t');rows=[]
    for (draw,fold),part in pools.groupby(['draw','outer80']):
        ix=part.row_index.to_numpy(int);yy=m.primary.to_numpy(bool)[ix]
        for condition,s in scores.items():
            rows.append(dict(condition=condition,draw=int(draw),fold=int(fold),**stats.metrics(yy,s[ix])))
    mm=pd.DataFrame(rows);mm.to_csv(target(ROOT/'source_data/matched_candidate_pool_metrics.tsv'),sep='\t',index=False)
    dm=mm.groupby(['condition','draw'])[['AP','AUROC','precision_at_20']].mean().reset_index()
    dm.to_csv(target(ROOT/'source_data/matched_per_draw_means.tsv'),sep='\t',index=False)
    ms=dm.groupby('condition').agg(AP=('AP','mean'),AP_draw_SD=('AP','std'),AUROC=('AUROC','mean'),AUROC_draw_SD=('AUROC','std')).reset_index()
    ms.to_csv(target(ROOT/'source_data/matched_summary.tsv'),sep='\t',index=False)
    def contrast(rival,subset,metric):
        return cd.loc[cd.main.eq('time_attention')&cd.comparator.eq(rival)&cd.subset.eq(subset)&cd.metric.eq(metric)].iloc[0]
    base=cf.loc[cf.condition.eq('bg_concat')&cf.top_k.eq(100)].iloc[0]
    cand=cf.loc[cf.condition.eq('time_attention')&cf.top_k.eq(100)].iloc[0]
    mainrow=summary.loc[summary.condition.eq('time_attention')].iloc[0];gain=contrast('bg_concat','all','AP')
    guards=dict(main_pool_AUROC_target=bool(mainrow.AUROC>.8),paired_AP_gain=bool(gain.delta>=c['minimum_useful_AP_gain'] and gain.family_adjusted_low>0),
        nonPPR_AP_preserved=bool(contrast('bg_concat','nonPPR','AP').delta>=-c['nonPPR_AP_tolerance']),
        within_PPR_AUROC_preserved=bool(contrast('bg_concat','PPR','AUROC').delta>=-c['PPR_AUROC_tolerance']),
        PPR_fraction_not_increased=bool(cand.PPR_fraction<=base.PPR_fraction),PPR_run_reduced=bool(cand.PPR_longest_run<base.PPR_longest_run))
    save(ROOT/'reports/decision.json',dict(completed=now(),registered_main=c['registered_main'],guards=guards,
        meets_all_guards=bool(all(guards.values())),automatic_model_replacement=False,independent_validation=False,
        learned_pooling_increment_supported=bool(contrast('time_mean','all','AP').family_adjusted_low>0),
        explicit_DAP_increment_supported=bool(contrast('time_no_pos','all','AP').family_adjusted_low>0),
        graph_increment_supported=bool(all(contrast(v,'all','AP').family_adjusted_low>0 for v in ['no_message','gcn_rewired']))))
    def table(frame):
        out=['|'+'|'.join(frame.columns)+'|','|'+'|'.join(['---']*len(frame.columns))+'|']
        for row in frame.itertuples(index=False,name=None):
            out.append('|'+'|'.join('%.5f'%v if isinstance(v,(float,np.floating)) else str(v) for v in row)+'|')
        return '\n'.join(out)
    target(ROOT/'reports/README_结果.md').write_text('\n'.join(['# 时间Transformer与时间注意力结果','','本轮是回顾性开发比较，不是新增独立验证。',
        '',table(summary),'','## 匹配1:3候选池（另一个口径）','',table(ms),'','## OOF未标注候选PPR分布','',table(cf),'',
        '时间权重只用于辅助解释，不是逐组学归因或因果机制；原始分数/排名不按家族调分。',
        '协议护栏：'+str(guards),'未自动替换原模型、MEKG或网站。'])+'\n')
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
