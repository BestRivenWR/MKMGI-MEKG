"""Release tests and human-readable handoff. Only run after all imports finalize."""
from pathlib import Path
import json,csv,gzip,sqlite3,random,collections,time
from build_v31 import ROOT,PARENT,old,sha,js
from coordinate_tools import StrictChain,self_test
from query_v31 import query

def main():
    release_db=ROOT/'release/mekg.sqlite';copy_file=ROOT/'audit_local_copy.json'
    local_copy=json.loads(copy_file.read_text()) if copy_file.exists() else None
    db=Path(local_copy['database']) if local_copy else release_db
    c=sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1',uri=True)
    c.execute('PRAGMA cache_size=-262144');c.execute('PRAGMA temp_store=MEMORY')
    checks={};details={}
    def check(name,sql):
        n=c.execute(sql).fetchone()[0];checks[name]=n==0
        if n:details[name]=n
    checks['sqlite_integrity']=c.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    print('QA_INTEGRITY_DONE',flush=True)
    for col in ['subject_node_id','object_node_id']:
        check('edge_'+col+'_exists',f'SELECT COUNT(*) FROM edges e LEFT JOIN nodes n ON e.{col}=n.node_id WHERE n.node_id IS NULL')
    check('evidence_edges_exist','SELECT COUNT(*) FROM evidence v LEFT JOIN edges e ON v.edge_id=e.edge_id WHERE e.edge_id IS NULL')
    check('sources_resolve','SELECT COUNT(*) FROM v3_assertions a LEFT JOIN v3_sources s ON a.source_id=s.source_id WHERE s.source_id IS NULL')
    check('new_evidence_not_training','SELECT COUNT(*) FROM evidence v JOIN v3_assertions a USING(evidence_id) WHERE v.used_for_training!=0')
    check('no_unearned_independent_support',"SELECT COUNT(*) FROM evidence WHERE independence_class='PUBLIC_INDEPENDENT_SUPPORT'")
    check('new_confidence_not_invented','SELECT COUNT(*) FROM evidence v JOIN v3_assertions a USING(evidence_id) WHERE v.confidence_score IS NOT NULL')
    check('peaks_target_genes_exist',"SELECT COUNT(*) FROM v31_peaks p LEFT JOIN nodes n ON n.node_id='gene:'||p.target_gene_id WHERE p.target_gene_id IS NOT NULL AND n.node_id IS NULL")
    check('peaks_native_intervals_valid','SELECT COUNT(*) FROM v31_peaks WHERE native_start0<0 OR native_end0<=native_start0')
    check('peaks_projected_intervals_valid','SELECT COUNT(*) FROM v31_peaks WHERE v5_chrom IS NOT NULL AND (v5_start0<0 OR v5_end0<=v5_start0 OR v5_end0-v5_start0!=native_end0-native_start0)')
    check('unmapped_peaks_not_given_coordinates',"SELECT COUNT(*) FROM v31_peaks WHERE mapping_status='native_only_unmapped_or_gap' AND (v5_chrom IS NOT NULL OR v5_start0 IS NOT NULL OR v5_end0 IS NOT NULL)")
    checks['chain_unit_tests']=self_test()
    chain=StrictChain(ROOT/'raw/v4_to_v5.chain')
    sampled=c.execute('SELECT native_chrom,native_start0,native_end0,v5_chrom,v5_start0,v5_end0 FROM v31_peaks WHERE v5_chrom IS NOT NULL LIMIT 1000').fetchall()
    checks['1000_projected_peaks_recompute']=all(chain.map(*r[:3])[:3]==r[3:] for r in sampled) and len(sampled)==1000
    check('cell_expression_groups_exist','SELECT COUNT(*) FROM v31_cell_expression e LEFT JOIN v31_cell_groups s USING(group_id) WHERE s.group_id IS NULL')
    check('cell_expression_gene_nodes_exist',"SELECT COUNT(*) FROM v31_cell_expression e LEFT JOIN nodes n ON n.node_id='gene:'||e.gene_id WHERE n.node_id IS NULL")
    check('cell_counts_and_detection_valid','SELECT COUNT(*) FROM v31_cell_expression e JOIN v31_cell_groups s USING(group_id) WHERE e.UMI_sum<0 OR e.detected_nuclei<0 OR e.detected_nuclei>s.nuclei OR e.UMI_sum<e.detected_nuclei OR s.UMI_total<=0')
    cellcounts=dict(c.execute("SELECT label,SUM(nuclei) FROM v31_cell_groups WHERE group_level='cell_type' GROUP BY label"))
    expected_cells={'Aleurone':914,'Basal endosperm transfer layer':521,'Endosperm adjacent to scutellum':295,
      'Placenta-chalazal region':303,'Starchy endosperm':29473}
    checks['five_author_cell_types_preserved']=cellcounts==expected_cells
    checks['both_group_levels_total_31506']=dict(c.execute('SELECT group_level,SUM(nuclei) FROM v31_cell_groups GROUP BY group_level'))=={'cell_type':31506,'cluster':31506}
    checks['three_replicates_preserved']=dict(c.execute("SELECT replicate,SUM(nuclei) FROM v31_cell_groups WHERE group_level='cell_type' GROUP BY replicate"))=={'NDO001':4043,'NDO002':14183,'NDO003':13280}
    check('maternal_cells_labeled',"SELECT COUNT(*) FROM v31_cell_groups WHERE label='Placenta-chalazal region' AND json_extract(attributes_json,'$.maternal_tissue')!=1")
    check('contact_anchor_coordinates_valid','SELECT COUNT(*) FROM v31_contact_anchors WHERE start0<0 OR end0<=start0')
    check('GWAS_p_values_valid','SELECT COUNT(*) FROM v31_gwas_loci WHERE pvalue<=0 OR pvalue>1 OR position<1')
    check('KIL_significant_edges_pass_padj_threshold',"SELECT COUNT(*) FROM evidence WHERE source_node_id='source:public:v3:kil_doublemutant_rna_2025' AND (json_extract(attributes_json,'$.padj')>=0.05 OR ABS(effect_size)<1)")
    check('HiChIP_not_mislabeled_endosperm',"SELECT COUNT(*) FROM evidence WHERE source_node_id='source:public:v3:hichip_2025' AND tissue_or_compartment!='seedlings'")
    check('GWAS_no_combined_p_or_causal_gene_claim',"SELECT COUNT(*) FROM evidence WHERE source_node_id='source:public:v3:kernel_gwas_2026' AND (json_extract(attributes_json,'$.no_combined_pvalue')!=1 OR json_extract(attributes_json,'$.causal_gene_proven')!=0)")
    check('TSS_coordinates_and_summits_valid','SELECT COUNT(*) FROM v31_tss_peaks WHERE start0<0 OR end0<=start0 OR summit0<start0 OR summit0>=end0')
    check('TSS_replicate_support_consistent','''SELECT COUNT(*) FROM v31_tss_consensus t JOIN v31_tss_peaks a ON a.peak_id=t.peak1_id JOIN v31_tss_peaks b ON b.peak_id=t.peak2_id
      WHERE a.replicate!=1 OR b.replicate!=2 OR a.chromosome!=b.chromosome OR ABS(a.summit0-b.summit0)>100
      OR (MIN(a.end0,b.end0)-MAX(a.start0,b.start0))*1.0/(a.end0-a.start0)<0.5
      OR (MIN(a.end0,b.end0)-MAX(a.start0,b.start0))*1.0/(b.end0-b.start0)<0.5''')
    with gzip.open(old.PRED,'rt') as f:expected={r['gene_id']:r for r in csv.DictReader(f,delimiter='\t')}
    actual={g:json.loads(s) for g,s in c.execute('SELECT * FROM v3_model_predictions')}
    checks['all_17917_frozen_prediction_rows_identical']=actual==expected and len(actual)==17917
    preserved=json.loads((ROOT/'parent_preservation.json').read_text())
    checks['parent_v3_file_hash_unchanged']=sha(PARENT/'release/mekg.sqlite')==preserved['sha256']
    checks['original_prediction_file_hash_unchanged']=sha(old.PRED)==preserved['model_sha256']
    modules={k:json.loads(v) for k,v in c.execute('SELECT * FROM v31_modules')}
    checks['seven_modules_complete']=set(modules)=={'tf_peaks','myb_peaks','hichip','gwas','cells','kil_rna','tss'}
    mapping_audit=json.loads((ROOT/'release/TF_MAPPING_AUDIT.json').read_text())
    checks['TF_mapping_loss_accounted_for']=sum(mapping_audit['raw_row_counts'].values())==modules['tf_peaks']['input_rows']
    checks['TF_unmapped_loss_matches_builder']=sum(v for k,v in mapping_audit['raw_row_counts'].items() if k!='unique_bijective_mapping')==modules['tf_peaks']['skipped']['tf_peak_target_nonbijective']
    checks['no_source_certified_without_audit']=all(not json.loads(r[0])['certified_independent'] for r in c.execute('SELECT audit_json FROM v31_source_audit'))
    tested=random.Random(20260907).sample(sorted(expected),8)+['Zm00001eb125330','Zm00001eb139600']
    out=ROOT/'examples';out.mkdir(exist_ok=True)
    for gene in tested:
        result=query(gene,limit=8,db=db)
        result['evidence_database']=str(release_db)
        checks['query_'+gene]=result.get('gene_id')==gene and ((result.get('prediction') is not None)==(gene in expected))
        if result.get('prediction') and gene in expected:
            checks['query_score_'+gene]=result['prediction']['priority_score']==float(expected[gene]['priority_score'])
        (out/(gene+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2))
    checks['unknown_gene_not_fabricated']=query('Zm00001eb999999',db=db).get('status')=='not_found'
    alias=query('ereb192',db=db);checks['alias_resolution']=alias.get('gene_id')=='Zm00001eb125330'
    print('QA_QUERIES_DONE',flush=True)
    tables=['nodes','edges','evidence','observations','v31_peaks','v31_contacts','v31_gwas_loci','v31_cell_groups',
      'v31_cell_expression','v31_cluster_markers','v31_tss_peaks','v31_tss_consensus','v31_tss_genes']
    counts={t:c.execute('SELECT COUNT(*) FROM '+t).fetchone()[0] for t in tables}
    parentqa=json.loads((PARENT/'release/QA.json').read_text())
    added={t:counts[t]-parentqa['total_counts'][t] for t in ['nodes','edges','evidence','observations']}
    sourcecounts=dict(c.execute('SELECT source_id,COUNT(*) FROM v3_assertions GROUP BY source_id'))
    source_rows={r[0]:json.loads(r[1]) for r in c.execute('SELECT source_id,source_json FROM v3_sources')}
    parent_sources=set(parentqa['new_source_evidence_counts'])
    for key,src in source_rows.items():
        base=PARENT if key in parent_sources else ROOT
        src['resolved_source_files']={f:str(base/'raw'/f) for f in src['files']}
    details['peak_mapping_by_source']=[{'source_id':r[0],'peak_gene_records':r[1],'projected_v5_records':r[2]} for r in c.execute('SELECT source_id,COUNT(*),SUM(v5_chrom IS NOT NULL) FROM v31_peaks GROUP BY source_id')]
    details['replicated_RAMPAGE_associated_model_genes']=c.execute('SELECT COUNT(DISTINCT g.gene_id) FROM v31_tss_genes g JOIN v3_model_predictions p USING(gene_id)').fetchone()[0]
    details['cell_expression_model_gene_coverage']=c.execute('SELECT COUNT(DISTINCT g.gene_id) FROM v31_cell_expression g JOIN v3_model_predictions p USING(gene_id)').fetchone()[0]
    details['gwas_unique_chromosome_position_loci']=c.execute('SELECT COUNT(*) FROM (SELECT chromosome,position FROM v31_gwas_loci GROUP BY chromosome,position)').fetchone()[0]
    c.close()
    digest=sha(db)
    if local_copy:
        checks['NVMe_QA_copy_matches_source_bytes']=digest==local_copy['source_sha256_during_copy']==sha(release_db)
        checks['source_not_changed_during_QA']=release_db.stat().st_size==local_copy['bytes'] and release_db.stat().st_mtime_ns==local_copy['source_mtime_ns']
    report={'passed':all(checks.values()),'checks':checks,'details':details,'total_counts':counts,'added_to_v3':added,
      'source_evidence_counts':sourcecounts,'modules':modules,'author_cell_type_counts':cellcounts,
      'database_sha256':digest,'live_website_changed':False,'model_scores_modified':False,'read_only_local_validation_copy':local_copy}
    (ROOT/'release/QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    manifest={'release_id':'MEKG-public-v3.1-20260907','database':str(release_db),'database_sha256':digest,
      'parent_preservation':preserved,'query_entrypoint':str(Path(__file__).resolve().parents[2]/'query.py'),
      'prediction_table':'v3_model_predictions','prediction_source':str(old.PRED),
      'source_registry':source_rows,'source_audit':str(ROOT/'release/TRAINING_OVERLAP_AUDIT.json'),
      'new_quantitative_tables':['v31_cell_expression','v31_peaks','v31_contacts','v31_gwas_loci','v31_tss_peaks'],
      'legacy_tables_not_current_query':['gene_cards','public_gene_support'],
      'code_sha256':{p.name:sha(p) for p in Path(__file__).resolve().parent.iterdir() if p.suffix in ['.py','.R','.sh']},
      'query_runtime_files':['release/mekg.sqlite','repository:mekg/query.py','repository:mekg/_builders/v31/query_v31.py','repository:mekg/_builders/v31/query_base_v3.py'],
      'build_dependency':'read-only builder sources at '+str(Path(__file__).resolve().parents[1]/'v3'),
      'parent_code_sha256':{n:sha(Path(__file__).resolve().parents[1]/'v3'/n) for n in ['build_v3.py','table_io.py','query_v3.py']},
      'website_deployed':False,'QA_passed':report['passed']}
    (ROOT/'release/manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    write_report(report,manifest)
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
    assert report['passed'],'Release remains unfrozen until all tests pass'
    (ROOT/'release/COMPLETE.json').write_text(js({'release_id':manifest['release_id'],'database_sha256':digest,
      'qa_sha256':sha(ROOT/'release/QA.json'),'manifest_sha256':sha(ROOT/'release/manifest.json'),'checks_passed':len(checks)}))

def write_report(r,m):
    mods=r['modules'];n=r['total_counts'];added=r['added_to_v3'];tf=mods['tf_peaks'];cell=mods['cells']
    rows='\n'.join(f'| {k} | {v:,} |' for k,v in added.items())
    sources='\n'.join(f'- [{s["title"]}](https://doi.org/{s["doi"]})' for s in m['source_registry'].values())
    text=f'''# MEKG public v3.1 完成记录

本轮补充证据图谱，不重新训练，不改模型分数，不切换 8788 网站。上一版 v3 文件及冻结预测表逐项验证未变。

## 已完成

1. **作者最终单核对象**：B73 19 DAP，31,506 个核、3 个样本、23 个聚类。RNA assay 有 32,815 个特征；默认 SCT assay 的 28,137 个特征不是本次计数分母。
   - 最终对象包含 26,078 个来自 filtered 矩阵的核与 5,428 个来自 raw 矩阵的核。旧 filtered 矩阵共有 32,037 个核，二者不是单纯包含关系。
   - 作者 annotation.predicted：淀粉胚乳 29,473；糊粉层 914；BETL 521；EAS 295；母体胎盘合点区域 303。注释是作者预测，未冒称本项目新验证。
   - 聚类与预测细胞类型不是同一字段：例如 cluster 19 共 290 个核，其中 249 个标为 EAS、41 个标为淀粉胚乳。未把整个聚类强行改成单一类型。
   - 作者保留的 4,725 个非 Singlet 分类核原样保留并标记；未暗中再次筛选或宣称重新完成 QC。
   - 导出 15 个“细胞类型×样本”和 69 个“聚类×样本”汇总，合计 {n['v31_cell_expression']:,} 条基因–组计数。两种分组来自相同核，不是独立重复。
   - 导入 {n['v31_cluster_markers']:,} 条作者聚类差异标记记录。UMAP 是降维坐标，不是空间坐标。
2. **513 个 TF 峰文件**：解析 {tf['input_rows']:,} 条原始行，去重后峰–基因记录及各来源数量见 QA。{tf['resolved_tf_profiles']}/513 个文件的 TF 身份可唯一映射；其余保留原始文件名而不猜测 TF。
   - 仅使用作者主要 geneId 关联，不把 flank_geneIds 全部扩成靶基因。
   - 使用 MaizeGDB 官方 v4→v5 chain，仅接受完整区间落在唯一无缺口块且反向映射一致的转换。未通过的峰保留原生 v4 坐标。
   - 目标基因 ID 的 1,986,023 行可唯一映射；240,796 行为一对多、86,148 行为多对一、417,270 行无可用 v5 对应。未映射行仍可按原始文件和行号追溯，详见 `TF_MAPPING_AUDIT.json`；这与峰坐标是否能转换是两个不同问题。
   - 513 个文件不等于 513 个均已确认的体内实验。逐文件的 ampDAP/ChIP 等实验类型仍需补充来源对应表。
3. **MYB127 CUT&Tag**：两 FLAG 重复峰互惠重叠≥50%，且原始峰均不与 WT 峰重叠；保守坐标转换后关联基因体或上游 2 kb，得到 {mods['myb_peaks']['peak_associated_genes']:,} 个峰关联基因。这是上下文筛选，不是 BAM 层面的差异结合统计或普遍的直接调控证明。
4. **HiChIP 物理接触**：新增 {n['v31_contacts']:,} 条作者接触记录，保留原生双端坐标；通过映射的锚点建立 v5 查询索引。B73 材料来自幼苗，不是胚乳，查询明确标为跨组织线索。作者基因对表与接触表属于同一实验，不重复计算独立支持。
5. **籽粒性状 GWAS**：纳入 187 份材料的关联研究，保留性状、年份、模型及原始 P 值，共 {n['v31_gwas_loci']:,} 条位点–模型–性状–年份记录，对应 {r['details']['gwas_unique_chromosome_position_loci']:,} 个唯一染色体–位置位点。按原文阈值筛选后，在 SNP ±50 kb 范围内关联所有符合条件的基因；不合并相关性状的 P 值，也不把邻近当作因果。
6. **KIL1/KIL2 双突变体 RNA**：补入 A188 背景 EAS 富集样本的作者差异表达清单，区别于 B73 单核数据。双敲结果不拆成某一个 KIL 的单基因效应，差异表达不冒称直接靶标。补充表的 stat 符号与 log2FC 相反，原值保留；上下调方向依据一致的工作表名称与 log2FC，不使用 stat 推断方向。
   - 原清单中 Zm00001d053108 的 padj=0.0503633，未通过本版 padj<0.05 阈值；其连续表达观测保留，但不新增显著差异关系。
7. **MaizeCODE RAMPAGE**：B73 胚乳 15 DAP、两个重复、原生 NAM v5 坐标，{n['v31_tss_peaks']:,} 个原始峰。互惠重叠≥50%且 summit 距离≤100 bp，得到 {n['v31_tss_consensus']:,} 个重复支持的峰对，关联 {mods['tss']['peak_associated_genes']:,} 个基因。保留两个 summit，不擅自替换注释 TSS；发布文件无链方向，不能据此认定转录方向。

## 相比 v3 的增量

| 标准表 | 新增条数 |
|---|---:|
{rows}

细胞计数、原始峰、接触锚点和 GWAS 位点另存专用表；这些记录数不能简单相加当作独立生物学证据数。

## 查询和网站衔接

- 新目录：`{ROOT}`
- 数据库：`{ROOT}/release/mekg.sqlite`
- 查询脚本：`{ROOT}/query_v31.py`；封装入口 `query_gene.sh`。
- 查询部署最小文件：`release/mekg.sqlite`、`query_v31.py`、`query_base_v3.py`。后两者使用 Python 标准库；原始/派生数据用于审查与重建，不需每次查询加载。
- 冻结模型分数仍来自 `{old.PRED}`，数据库 `v3_model_predictions` 为逐行原样副本。
- 旧 `gene_cards` 与 `public_gene_support` 是历史快照，不应作为新网页的当前结果来源。

```bash
{ROOT}/query_gene.sh Zm00001eb003760 --limit 20
```

JSON 保留预测与不确定性、内部关系数、公共证据、实验条件、来源、机制线索；新增细胞类型/聚类表达、TF 峰、RAMPAGE TSS、HiChIP 锚点、标签/家族重叠审查。网页尚未切换到本版。

## 独立性审查与剩余事项

- 已核对正式训练文件清单、源文件清单、标签 DOI、基因家族和旧公共来源注册表。标签文件 SHA 与冻结清单一致；现有 17,917 行标签仅 39 行含 primary_doi。
- 未检出文件哈希/DOI 重合并不能证明原始样本独立；同一论文的多种派生表已归组。`PUBLIC_INDEPENDENT_SUPPORT` 继续保留但不擅自认证来源。
- 还需：513 TF 逐文件的实验/样本对应；旧标签缺失文献依据与样本变换链补全；模型集合外基因的家族净化。此版本可作解释，不是新增独立 benchmark。
- MapMan 官方映射下载需用户登录并接受许可；本轮未绕过。等待用户提供官方文件后再接入。
- 已增加真实 GWAS 位点层；更全面的跨研究 QTL/meta-QTL 更新仍未完成，未以零散位点清单冒充全覆盖。
- 保守坐标转换会舍弃跨 indel 的区间；覆盖率见 QA，不将缺失视为无生物学效应。

## 验收

自动检查 {len(r['checks'])} 项，通过 {sum(r['checks'].values())} 项。详细记录见 `QA.json`、`manifest.json` 和 `TRAINING_OVERLAP_AUDIT.json`。所有 17,917 个冻结预测记录逐行一致。

早期 TF 入库因 SQLite 小缓存产生大量随机读写而中断；未完成数据库保留在 `build_attempt_tf_cache2mb`，仅供过程追踪，不可作为可查询发布版。当前版已从未变的父版本重建。

## 原始来源

{sources}
'''
    (ROOT/'release/UPGRADE_REPORT_ZH.md').write_text(text,encoding='utf8')

if __name__=='__main__':main()
