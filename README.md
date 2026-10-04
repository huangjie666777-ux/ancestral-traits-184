# eDNA 固定参考树短序列放置后端

基于 FastAPI + Biopython 的环境 DNA 短序列放置（phylogenetic placement）服务：
在**固定参考拓扑和原始枝长**的前提下，把每条待测序列独立地枚举到参考树的
每一条枝上，用 JC69 模型通过 Felsenstein 递推最大化整树对齐的对数似然，
并以 jplace version 3 格式交付结果。

## 范围

- 参考序列 3–20 条，待测序列 1–5 条；所有序列等长且 ≤ 500 列。
- 允许 ACGT 与 IUPAC 歧义码（R/Y/S/W/K/M/B/D/H/V/N），缺口（`-` `.` `?`）
  视为未知、不提供碱基证据。
- 空序列、非法字符、重名、非正/非有限枝长、树叶与参考 ID 不一一对应等
  问题会被**定位**（记录名、列号或树节点）并**整单拒绝**（HTTP 422）。
- 对每条待测序列：枚举所有枝，插入节点把原枝拆成两段（均非负、总长不变），
  新叶枝长在 [0, 2] 内优化（粗网格 + 模式搜索细化），返回全部枝的最优位置、
  新枝长、对数似然与归一化的相对似然权重（like_weight_ratio）。
- 权重是跨枝的**相对支持度**，不是物种归属正确的概率。
- 完全无有效碱基证据（全缺口/全 N）的待测序列不强制归属，返回无法放置原因。
- 计算采用逐列缩放避免长序列下溢；不使用汉明距离或最近叶启发式。
- DNA 输入只接受 ACGT 与 DNA IUPAC 歧义码；RNA 的 `U` 不做静默 T 转换，
  直接定位为非法字符并整单拒绝（请先在上游转为 DNA）。
- 含逗号、括号、冒号、引号或花括号的叶 ID 在导出 jplace Newick 时按标准
  Newick 规则加单引号（内部引号双写），下游解析可完整还原。

## 群落距离分析（`POST /community/distances`）

环境采样人员可把 2–10 个样本放在一起比较谱系组成。每个样本给出
`sample_id`、已完成的放置任务 `job_id`，以及待测 ID → 非负整数读数的映射
（映射中省略的待测 ID 读数按 0 处理）。

仅当所有任务的**带枝号 Newick 文本完全一致**（同一参考树定位）时才比较；
未知任务、未知待测 ID、非法计数（负数/小数/布尔/非整数）或参考树不一致
都会定位并**整单拒绝**（HTTP 422）。该接口只读取已存储的 jplace 结果，
**不修改任务、不重算放置**。

质量构造：
- 某待测 ID 的读数按该 ID **全部候选枝**的 `like_weight_ratio` 分配到各参考枝
  的插入位置（既不只取最佳枝，也不把质量移到枝端）；插入坐标按 jplace 的
  `distal_length` 解释，即沿枝长方向距子端的距离，`pendant_length` 不参与。
- 无法放置的待测 ID（如全缺口序列）按 `{id, count, reason}` 列入
  `excluded` 并从该样本剔除；每个样本按**可放置读数总和**归一化，
  有效总数（`valid_total`）必须为正，否则整单拒绝。

距离（树上 Kantorovich–Rubstein / Wasserstein-1，p=1）：
- 对每条参考枝，按插入位置分段，积分两样本在该点**远端子树内累计质量差的
  绝对值**，再对全部枝求和；即把质量沿树搬运所需的最小加权枝长工作量。
- 距离为枝长单位（substitutions/site 尺度）的无量纲量：组成完全相同（允许
  总读数不同，归一化后一致）时为 0；同枝内位置不同也会产生小于整枝长的
  贡献，而不是非 0 即 1。

响应保持输入样本顺序：`matrix` 为对称矩阵、对角为 0；`pairs` 对每对样本
列出**全部枝号**及其非负贡献 `edge_contributions`（按枝号升序），
`contributions_sum` 等于该对 `distance`。`samples` 中保留每样本的
`excluded` 排除清单与 `valid_total` 有效总数。

## 模块协作

- `app/validation.py` — FASTA/Newick 解析与定位校验（整单拒绝）
- `app/tree.py` — 无根二叉树图、稳定枝号、带 `{edge_num}` 的 jplace 树序列化
- `app/likelihood.py` — JC69 转移矩阵、Felsenstein 双向消息（rerooting）与逐列缩放
- `app/placement.py` — 逐枝优化（拆分点 + 新枝长）、权重归一化与排序
- `app/jplacefmt.py` — jplace v3 文档组装（五个标准放置字段）
- `app/community.py` — 样本校验、质量构造（全候选权重）与树上 KR(p=1) 距离
- `app/main.py` — FastAPI 入口：`POST /place`、`GET /place/{job_id}/jplace`、
  `POST /community/distances`

## 运行

\`\`\`bash
PYTHONPATH=. .venv/bin/python -m uvicorn app.main:app --port 8169
\`\`\`

## 自测

\`\`\`bash
.venv/bin/python -m compileall -q app scripts
PYTHONPATH=. .venv/bin/python scripts/selftest.py
\`\`\`

自测覆盖：示例放置、权重归一化、似然降序排序、jplace 下载结构、
非法字符与 `U` 的定位拒绝、群落 KR 矩阵（对称性/对角/贡献和等于距离/
相同组成距离为 0）、群落接口的各类整单拒绝、以及逗号叶名 Newick 解析往返。
KR 积分另与 scipy 最小费用流线性规划在小树上交叉验证一致。

## curl 示例

\`\`\`bash
# 组装请求（examples/ 内含 5 参考序列、3 待测序列、示例树）
python3 - <<'EOF' > /tmp/payload.json
import json, pathlib
ex = pathlib.Path('examples')
print(json.dumps({
  "reference_fasta": ex.joinpath('reference.fasta').read_text(),
  "query_fasta": ex.joinpath('queries.fasta').read_text(),
  "newick": ex.joinpath('tree.nwk').read_text()}))
EOF

# 提交放置
curl -s -X POST http://127.0.0.1:8169/place \
  -H 'Content-Type: application/json' --data-binary @/tmp/payload.json

# 下载 jplace（job_id 取自上一步响应）
curl -OJ http://127.0.0.1:8169/place/<job_id>/jplace

# 群落距离：2..10 个样本引用同一参考树上的已有任务
curl -s -X POST http://127.0.0.1:8169/community/distances \
  -H 'Content-Type: application/json' -d '{
    "samples": [
      {"sample_id": "site1", "job_id": "<job_id>",
       "counts": {"q1": 8, "q3": 3}},
      {"sample_id": "site2", "job_id": "<job_id>",
       "counts": {"q1": 2, "q2": 6}}
    ]}'
\`\`\`

## jplace 输出

\`tree\` 字段为带 `{edge_num}` 注释的 Newick 树，枝号与结果中的
`edge_num` 一致；`fields` 为五个标准放置字段：
`edge_num, likelihood, like_weight_ratio, distal_length, pendant_length`。
同一查询的放置按似然降序排列，同值按枝号升序。

## community 响应示例（节选）

```json
{
  "sample_ids": ["site1", "site2"],
  "matrix": [[0.0, 0.2387], [0.2387, 0.0]],
  "pairs": [{
    "i": 0, "j": 1, "sample_a": "site1", "sample_b": "site2",
    "distance": 0.2387, "contributions_sum": 0.2387,
    "edge_contributions": [
      {"edge_num": 0, "contribution": 0.0597},
      {"edge_num": 1, "contribution": 0.0477}
    ]
  }],
  "samples": [
    {"sample_id": "site1", "job_id": "...", "valid_total": 8,
     "excluded": [{"id": "q3", "count": 3, "reason": "no valid base ..."}]},
    {"sample_id": "site2", "job_id": "...", "valid_total": 8, "excluded": []}
  ]
}
```
