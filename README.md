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

## 批次约束 PERMANOVA（`POST /community/permanova`）

实验人员可判断处理组间的群落组成差异是否超出随机分组水平。请求含
4–10 个原格式样本（`sample_id`/`job_id`/`counts`），每样本附非空
`group`；至少 2 组且每组至少 2 个样本。可选 `batch`：任一样本给出则
全体必须给出非空批次。`permutations` 为 1–9999 的整数，`seed` 为整数；
非法设计（样本数、组人数、批次缺失、置换数/种子类型越界）会被定位并
整单拒绝（HTTP 422）。接口复用已有任务、全候选质量构造与 KR 矩阵，
**不重算放置、不修改任务**。

统计量（距离平方和分解）：
- `ss_total` = 全部样本对距离平方和 / 样本数 n
- `ss_within` = 各组内样本对距离平方和 / 组人数，跨组相加
- `ss_between` = 两者之差；pseudo-F = (ss_between/(g-1)) / (ss_within/(n-g))
- `R2` = ss_between / ss_total

置换方案：只交换组别标签，距离矩阵与各组人数固定；有批次时仅在批次内
均匀交换，保留各批次组别数量。若不存在任何能改变分组的合法交换（如
每个批次内组别单一），明确拒绝，**不把跨批次混洗当作有效检验**。合法
分配总数不超过请求的置换数时全部枚举（含原分组），p 值为 F 不小于观测
值的比例；否则按 `seed` 随机抽取 `permutations` 次（允许重复），
p = (极端次数 + 1)/(置换数 + 1)，同请求同种子完全复现。总平方和或组内
平方和为 0 时明确拒绝。

响应含样本顺序（`sample_ids`）、`groups`、`batches`、三个平方和、两个
自由度、`F`、`R2`、`p_value`、`method`（`exact`/`random`）、实际次数
`permutations` 及合法分配总数 `n_assignments`。

## 模块协作

- `app/validation.py` — FASTA/Newick 解析与定位校验（整单拒绝）
- `app/tree.py` — 无根二叉树图、稳定枝号、带 `{edge_num}` 的 jplace 树序列化
- `app/likelihood.py` — JC69 转移矩阵、Felsenstein 双向消息（rerooting）与逐列缩放
- `app/placement.py` — 逐枝优化（拆分点 + 新枝长）、权重归一化与排序
- `app/jplacefmt.py` — jplace v3 文档组装（五个标准放置字段）
- `app/community.py` — 样本校验、质量构造（全候选权重）与树上 KR(p=1) 距离
- `app/permanova.py` — 批次约束的分组设计校验、平方和分解与精确/随机置换检验
- `app/main.py` — FastAPI 入口：`POST /place`、`GET /place/{job_id}/jplace`、
  `POST /community/distances`、`POST /community/permanova`

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
相同组成距离为 0）、深层树远端质量只累计一次、零计数未知 ID 拒绝、
群落接口的各类整单拒绝、PERMANOVA 的精确/随机置换与全部非法设计拒绝、
以及逗号叶名 Newick 解析往返。
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

# 批次约束 PERMANOVA：4 个样本、2 组、仅在批次内交换标签
curl -s -X POST http://127.0.0.1:8169/community/permanova \
  -H 'Content-Type: application/json' -d '{
    "permutations": 999, "seed": 42,
    "samples": [
      {"sample_id": "s1", "job_id": "<job_id>", "group": "treat", "batch": "b1",
       "counts": {"q1": 8, "q2": 1}},
      {"sample_id": "s2", "job_id": "<job_id>", "group": "ctrl", "batch": "b1",
       "counts": {"q1": 1, "q2": 8}},
      {"sample_id": "s3", "job_id": "<job_id>", "group": "treat", "batch": "b2",
       "counts": {"q1": 7, "q2": 2}},
      {"sample_id": "s4", "job_id": "<job_id>", "group": "ctrl", "batch": "b2",
       "counts": {"q1": 2, "q2": 7}}
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
