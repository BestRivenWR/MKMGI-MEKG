# MKMGI × MEKG

本仓库集中提供模型与解释系统代码：首页按“模型训练 → 图谱构建 → 基因查询与解释”阅读即可。

已同步 2026 年 10 月 22 日最终投稿版网页与解释 agent，代码同步日期为 10 月 23 日。
页面默认英文，支持 **EN / 中文** 切换；按 B73 v5 基因 ID 或 Gene symbol 查询，
示例为 **O2**。结果页整合模型结果、研究问答、证据网络和文献与数据，支持桌面和手机使用。

| 目录 | 保留内容 |
|---|---|
| `mkmgi/` | 时间 Transformer、时间注意力、两层残差关系 GCN；预处理、训练、评价和全基因组预测 |
| `evaluation/` | 六种比较方法，以及组学、空间、功能、蛋白、图和结构消融 |
| `mekg/` | MEKG 公共证据导入、v5 坐标处理、图谱组装和查询 |
| `agent/` | DeepSeek 解释 agent、联网文献检索、来源核对、网络展示和网页前端 |

详细测序处理、生信探索、旧模型试验和论文绘图已从这版移出。必要的共享实现收进
`mkmgi/_shared/`，图谱分阶段导入收进 `mekg/_builders/`，不再要求用户找多个历史日期目录。

主模型保持本次论文确定的结构与训练逻辑。重组涉及运行目录和模块导入，未重新训练或计算
论文数值。模型数据、固定划分、权重、预测表和 MEKG 数据库需要单独提供。

- [模型训练、比较与消融](docs/training.md)
- [MEKG 数据组织与构建顺序](docs/mekg.md)
- [解释 agent 和网站部署](docs/agent.md)

网页中的 `priority_score`、`priority_percentile` 和跨种子不确定性来自最终模型的
`predictions/genomewide_time_gcn_consensus.tsv.gz`；时间权重来自
`predictions/genomewide_mean_gates.npz`。MEKG 和联网文献补充证据，不修改模型分数。

代码已发布到 [BestRivenWR/MKMGI-MEKG](https://github.com/BestRivenWR/MKMGI-MEKG)。
API 密钥通过环境变量或部署端密钥文件提供；仓库不包含密钥、运行数据库、对话缓存或
未发表原始数据。许可证暂未确定，软件归档 DOI 和数据存取信息待补充。
代码来源与改动记录保存在 `provenance.json`，投稿版同步核查见 `submission_verification.json`。
