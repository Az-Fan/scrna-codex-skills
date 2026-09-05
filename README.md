# scRNA-seq Codex Skills

这是一套面向 Codex、Claude Code 和 WispScience 的可审计单细胞 RNA 测序工作流。当前发布版包含 13 个 skill，覆盖输入标准化、QC、人工批准后过滤、整合评估、预处理聚类、marker、人工注释、子集导出、基因程序评分、差异表达、通路富集和样本感知的细胞丰度变化分析。

规范开发仓库位于 `/home/faz_laptop/projects/scrna-codex-skills`；GitHub 仓库 `git@github.com:Az-Fan/scrna-codex-skills.git` 是固定版本的分发来源。科学计算默认复用服务器上已经注册的 pixi 环境，不会自动创建环境、修改环境或安装缺失依赖。

## 一、工作流总览

主流程如下：

```text
01 输入标准化
  ↓
02 计算 QC 指标
  ↓
03 生成 QC 审阅图表（不筛选）
  ↓ 人工确定并批准逐细胞决策
04 应用已批准的 QC 过滤
  ↓
05 整合/批次校正基准评估（按需运行）
  ↓
06 预处理与聚类
  ↓
07 查找 cluster marker
  ↓
08 准备注释审阅 → 人工确认 → 写入正式注释
  ↓
09 导出目标细胞子集（按需运行）
  ├─→ 10 基因集/通路程序评分
  ├─→ 11 样本级差异表达
          ↓
        12 ORA/GSEA 通路富集
  └─→ 13 细胞类型组成与局部邻域丰度变化
```

几个重要边界：

- `03` 只生成 QC 审阅材料，不删除细胞。
- `04` 只执行已经人工批准的逐细胞决策，不自行发明阈值。
- `05` 是方法比较步骤，不是每个项目都必须做；未校正的 `none` 始终作为基线。
- `07` 的 FindAllMarkers 用于注释证据，不能代替样本级疾病差异分析。
- `08` 必须先生成审阅表，再使用人工确认的完整决策表写入标签。
- `09` 只负责导出子集，不隐式完成子集 QC、整合、重聚类或注释。
- `11` 正式推断默认使用样本级 pseudobulk DESeq2，细胞不是生物学重复。
- `12` 读取完整差异结果表进行富集，不重新运行差异分析。
- `13` 比较的是样本间相对组成或 Milo 邻域丰度，必须先声明分母；它不等于绝对组织细胞数量。

## 二、安装固定版本

推荐从固定 release tag 安装，确保不同机器使用同一套代码：

```bash
git clone --branch v3.1.0 --depth 1 git@github.com:Az-Fan/scrna-codex-skills.git
python3 scrna-codex-skills/scripts/install_skills.py --target ~/.codex/skills
```

更新已有安装：

```bash
git -C scrna-codex-skills fetch --tags
git -C scrna-codex-skills checkout v3.1.0
python3 scrna-codex-skills/scripts/install_skills.py --target ~/.codex/skills --force
```

安装到 Claude Code 时，把目标改为 `~/.claude/skills`。安装后重新启动 agent 会话，使 skill discovery 读取新版本。安装脚本只组装 skill 自身的指令和执行器，不修改项目数据与 pixi 环境。

## 三、通用使用方式

### 1. 在项目中保存配置

每个 skill 都提供 `references/config.example.json`。建议将适配后的配置保存在分析项目的 `config/` 目录，并让 `output_dir` 指向新的阶段目录。不要直接修改 skill 内的示例文件。

### 2. 先做 dry-run

以下命令只校验配置并生成执行 manifest，不执行完整分析：

```bash
python3 ~/.codex/skills/<skill-name>/scripts/run.py --config config/example.json
```

检查 manifest 中的输入对象、元数据字段、pixi 环境、参数、输出目录和实际命令，确认无误后再执行：

```bash
python3 ~/.codex/skills/<skill-name>/scripts/run.py \
  --config config/example.json \
  --execute
```

### 3. 长任务使用 tmux supervisor

`02`、`04`–`13` 带有统一的 tmux 监督器。dry-run 仍在前台执行；只有已经确认的 `--execute` 命令才放入 tmux：

```bash
python3 ~/.codex/skills/<skill-name>/scripts/run_in_tmux.py \
  --session unique-session-name \
  --cwd /absolute/path/to/project \
  --log results/<stage>/tmux.log \
  --status results/<stage>/tmux_status.json \
  --execute \
  -- python3 ~/.codex/skills/<skill-name>/scripts/run.py \
     --config /absolute/path/to/config.json \
     --execute
```

监督器拒绝重复的活动 session，并记录独立的日志和状态 JSON。科学完成状态仍以 skill 自己的 `run_manifest.json`、状态表和必要结果文件为准。

### 4. 注册环境

| Skill | 默认复用的 pixi 项目 |
|---|---|
| `01`–`04` | `01-scrna-qc` |
| `05`–`06` | `03-integration` |
| `07`–`09` | `02-annotation` |
| `10` | `05-pathway_program` |
| `11`–`12` | `06-deg-analysis` |
| `13` | `07-cell-abundance`（R 方法为 `default`，scCODA 为 `sccoda`） |

完整兼容性说明见 [toolkit/references/compatibility.md](toolkit/references/compatibility.md)。

## 四、13 个 skill 的输入、用法和输出

### 01-scrna-standardize-input：标准化输入和元数据

**什么时候用**

用于一个数据集首次进入分析流程，或把已有 10x/Seurat 数据迁移到稳定的下游对象契约。它负责确认矩阵方向、基因标识、条形码、样本字段和来源，不应改变原始生物学数值。

**主要输入**

- `input.path`：10x 矩阵目录、10x H5，或 Seurat `.rds/.qs` 对象。
- `input.format`：通常使用 `auto`，也可明确指定格式。
- `metadata.sample`：对象中的样本字段。
- `project.id`、`output_dir`。
- 可选 `output.object_format=auto|qs|rds`。

配置模板：[config.example.json](skills/01-scrna-standardize-input/references/config.example.json)。

**运行**

```bash
python3 ~/.codex/skills/01-scrna-standardize-input/scripts/run.py \
  --config config/01_standardize_input.json
python3 ~/.codex/skills/01-scrna-standardize-input/scripts/run.py \
  --config config/01_standardize_input.json --execute
```

**主要输出**

- `standardized_object.qs` 或 `standardized_object.rds`：标准化 Seurat 对象。
- `cell_metadata.tsv`：完整细胞元数据。
- `samples.tsv`：每个生物学样本一行的样本表。
- `field_mapping.json`：源字段到 sample/condition/batch 等角色的映射。
- `provenance.json`：来源、格式、维度和变换记录。
- `run_manifest.json`：本次运行配置与文件清单。

详细对象约定见 [input-contract.md](skills/01-scrna-standardize-input/references/input-contract.md)。原始输入不会被覆盖。

### 02-scrna-calculate-qc-metrics：计算逐细胞 QC 指标

**什么时候用**

用于过滤前计算可用 QC 指标。支持 STARsolo 目录和含 raw counts 的 Seurat 对象；缺少 GTF、Velocyto 或可选软件包时，相关指标会被明确记录为 skipped，而不是静默消失。

**主要输入**

- 现有 pixi 可执行文件、pixi 项目和 environment。
- `input.type=seurat|starsolo`。
- Seurat 模式：`input.object`；STARsolo 模式：`input.starsolo_dir` 和样本声明。
- `metadata.sample`，可选 `metadata.batch`。
- `parameters.species=mouse|human`、raw-count assay、线粒体基因模式和随机种子。
- 可选 GTF、Velocyto、DecontX 设置和并行样本数。
- 明确的 `output_dir`。

配置模板：[config.example.json](skills/02-scrna-calculate-qc-metrics/references/config.example.json)。

**主要输出**

Seurat 输入时：

- `qc_metrics_object.rds`：添加 QC metadata 的派生对象。
- `metadata.tsv.gz`：完整逐细胞 QC 表。
- `qc_diagnosis.png`：基础 QC 诊断图。
- `_provenance/`：供执行器审计的 `metric_status.tsv` 和 `run_manifest.json`；不作为主要结果展示。

STARsolo 输入时，每个样本目录只保留 `counts.mtx.gz`、`features.tsv.gz`、`metadata.tsv.gz` 和 `qc_diagnosis.png`；所有样本的指标状态合并到输出根目录的 `_provenance/metric_status.tsv`。

可能计算的字段包括 `n_genes`、`n_UMIs`、`mito_frac`、`chrY_frac`、`nuclear_frac`、`ambient_frac_decontx`、`doublet_score`、细胞周期和血红蛋白评分。具体定义见 [input-output.md](skills/02-scrna-calculate-qc-metrics/references/input-output.md)。本 skill 不过滤细胞。

### 03-scrna-review-qc：生成 QC 审阅图表

**什么时候用**

用于过滤前后审阅 QC 分布、样本差异和候选阈值。它自动识别现有 QC 字段，并可结合 condition、batch、cluster、annotation 和 UMAP；缺失的可选字段会跳过。

**主要输入**

- 含 QC metadata 的 Seurat `.rds/.qs`。
- 现有 pixi 项目和 environment。
- 必需 `metadata.sample`；condition、batch、cluster、annotation 可选。
- `output.detail_level=compact|full`。
- 可选候选阈值 `thresholds`，但这些阈值仍需人工审阅。

配置模板：[config.example.json](skills/03-scrna-review-qc/references/config.example.json)。

**主要输出**

compact 模式根目录只保留：

- `qc_atlas.pdf`：样本、分组和 UMAP 感知的 QC 图集。
- `threshold_review.tsv`：候选阈值及待审批字段。
- `qc_summary_by_sample.tsv`：样本级 QC 汇总。
- `run_manifest.json`。

`output.detail_level=full` 时，额外的指标可用性、分位数、候选保留率、绘图状态表和预览 PNG 放在 `details/`。此步骤不会生成过滤对象，也不会自动填写批准状态。

### 04-scrna-apply-qc-filter：应用已批准的 QC 决策

**什么时候用**

仅在逐细胞过滤决策已经完成科学审阅并明确批准后使用。这个 skill 负责可靠执行和审计，不负责推导阈值。

**主要输入**

- `input.object`：含 raw counts 的 Seurat `.rds/.qs`。
- `input.decision_table`：覆盖对象全部细胞且每个 cell ID 唯一的 TSV/TSV.GZ。
- `decision.cell_id_column`。
- `decision.include_all_true`：必须全部为真的列。
- 可选 `decision.exclude_any_true`：任一为真即排除的列。
- `decision.expected_retained_cells`：人工批准的预期保留细胞数。
- `metadata.sample`，可选 `metadata.condition`。
- `approval.status=approved`、批准时间和说明。
- 新的 `output_dir` 与输出对象名。

配置模板：[config.example.json](skills/04-scrna-apply-qc-filter/references/config.example.json)，严格契约见 [filter-contract.md](skills/04-scrna-apply-qc-filter/references/filter-contract.md)。

**主要输出**

- `filtered_object.qs` 或配置的 `.rds`：只含 raw counts 和保留 metadata 的交接对象，不继承 reductions/graphs。
- `cell_filter_decisions.tsv.gz`：完整逐细胞输入决策、最终 retained 和原因。
- `filter_summary_by_sample.tsv`。
- 可选 `filter_summary_by_condition.tsv`。
- `filter_decision_counts.tsv`：各排除/保留原因计数。
- `approved_filter_record.tsv`：批准内容、预期细胞数和决策表 SHA-256。
- `session_info.txt`、`run_manifest.json`。

对象保存采用临时文件、回读验证和原子最终化；实际保留数与批准数不一致时立即停止。

### 05-scrna-benchmark-integration：比较批次校正方案

**什么时候用**

当项目需要判断是否校正批次、选择方法或评估过度校正时使用。支持 `none`、Harmony、Seurat RPCA、scVI、scANVI、BBKNN 和外部预计算 embedding。不是所有项目都需要运行。

**主要输入**

- 过滤后的 Seurat 对象和 assay。
- `metadata.sample`。
- 一个或多个明确的 `metadata.batch_variables`。
- condition 和可信 biological labels，用于混杂与生物学保真审计。
- 用户明确选择的方法及参数网格、dims、邻居数、随机种子。
- 明确选择的 batch-removal 与 biological-conservation 指标。
- 只生成用户请求的图；marker/program 图还需要用户提供基因集。

配置模板：[config.example.json](skills/05-scrna-benchmark-integration/references/config.example.json)。方法和指标语义见 [benchmark-criteria.md](skills/05-scrna-benchmark-integration/references/benchmark-criteria.md)。

**主要输出**

- `method_runs.tsv`/`method_runs_r.tsv`：各场景运行状态。
- `metric_results_long.tsv`：逐场景、batch、label、metric 的长表。
- `method_summary.tsv`、`method_ranking.tsv`、`skipped_metrics.tsv`。
- `design_confounding.tsv`：batch-condition 列联表和混杂指标。
- 用户选择的 UMAP、score heatmap、tradeoff、marker/program retention 图。
- `integration_benchmark_object.qs` 或 embedding bundle。
- `recommendation.md` 与 `recommendation_status.json`。
- `run_manifest.json`。

当 batch 与 condition 完全混杂，或请求的排名指标均不可用时，推荐状态必须为 `unresolved`。不会仅凭 UMAP “混得好”自动选方法。

### 06-scrna-preprocess-and-cluster：预处理和聚类

**什么时候用**

用于 QC 过滤之后的 Normalize/HVG/Scale/PCA、可选回归、可选 Harmony、邻居图、UMAP 和聚类。正常运行一个明确选择的 workflow；只有用户明确要求比较时才运行多个场景。

**主要输入**

- 含 raw counts 的 Seurat `.qs/.rds`。
- `input.assay`、`metadata.sample`、condition/batch 字段。
- `input.qc_status=filtered|unfiltered`；未过滤对象还必须设置 `input.allow_unfiltered=true`，且结果标记为探索性。
- 明确选择标准、细胞周期回归、Harmony、回归后 Harmony，或多方案比较。
- normalization、HVG、PCA dims、邻居参数、UMAP 参数和随机种子。
- fixed resolution，或 guided resolution scan。

固定分辨率配置见 [config.example.json](skills/06-scrna-preprocess-and-cluster/references/config.example.json)；guided 扫描见 [config.guided.example.json](skills/06-scrna-preprocess-and-cluster/references/config.guided.example.json)；确认分辨率后使用 [config.finalize.example.json](skills/06-scrna-preprocess-and-cluster/references/config.finalize.example.json)。

**主要输出**

每次运行：

- `preprocessed_clustered_object.qs`。
- `scenario_summary.tsv`、`scenario_cluster_similarity.tsv`。
- `<scenario>_elbow.png`。
- `workflow_state.json`、`session_info.txt`、追加式 `run.log`。
- `run_manifest_preprocess.json`；finalize 时另写 `run_manifest_finalize.json`。

固定或人工确认聚类还输出：

- `cell_assignments.tsv`。
- `<scenario>_cluster_sizes.tsv`。
- `<scenario>_sample_cluster_counts.tsv`。
- `<scenario>_umap_diagnostics.pdf`。

resolution scan 还输出：

- `<scenario>_umap_clusters_by_resolution.png`。
- `<scenario>_clustree_resolution.png`。
- `<scenario>_resolution_stability.tsv/.png`。

guided 模式只给出稳定性推荐并暂停，必须由用户确认 resolution 后再 finalize；推荐值不等于生物学最优值。完整契约见 [input-output.md](skills/06-scrna-preprocess-and-cluster/references/input-output.md)。

### 07-scrna-find-cluster-markers：查找 cluster marker

**什么时候用**

在聚类完成后、正式注释前运行。它用 FindAllMarkers 比较 cluster，生成无偏 marker 证据；不用于 PAH vs control 等条件差异推断。

**主要输入**

- 已聚类的 Seurat 对象。
- `metadata.cluster` 或活动 identities。
- normalized RNA/SCT assay；缺少 data layer 时可配置在内存中 LogNormalize。
- `test_use`、`only_pos`、`logfc_threshold`、`min_pct`、`min_diff_pct`、`return_thresh`。
- `reporting.top_n` 和 dotplot 每簇基因数。

配置模板：[config.example.json](skills/07-scrna-find-cluster-markers/references/config.example.json)。

**主要输出**

- `cluster_markers.tsv`：Seurat 返回的完整 marker 表，加每簇 rank。
- `top_cluster_markers.tsv`：确定性排序后的展示用 top marker。
- `cluster_marker_summary.tsv`：每簇细胞数和 marker 数，包括零 marker 的簇。
- `top_marker_dotplot.pdf`。
- `run_manifest.json`。

解读时应检查 sample-specific cluster、微小簇、线粒体/核糖体/应激/细胞周期/环境 RNA 等信号，见 [marker-interpretation.md](skills/07-scrna-find-cluster-markers/references/marker-interpretation.md)。

### 08-scrna-annotate-cells：人工审阅并应用细胞注释

这个 skill 有两个必须分开的 action。

#### A. `prepare_review`

**输入**

- 已聚类对象、sample/condition、确认的 cluster 字段和 UMAP reduction。
- 优先提供 `07` 的完整 marker 表 `input.markers`。
- 可选 canonical marker 分组，用于 dot plot。
- 只有需要兼容旧流程时才设置 `clustering.compute_if_missing=true`。

**输出**

- `cluster_markers.tsv`：复用或显式回退计算的完整 marker 表。
- `annotation_review.tsv`：待人工填写的候选 broad/fine/state、证据、冲突、sample bias、QC flag、confidence 和 decision。
- `clustered_object.qs`。
- `cluster_umap.pdf`、`cluster_sample_umap.pdf`。
- 可选 `canonical_marker_dotplot.pdf`。
- `run_manifest.json`。

配置模板：[config.example.json](skills/08-scrna-annotate-cells/references/config.example.json)。

#### B. `apply_confirmed`

**输入**

- 同一个聚类对象和 cluster/reduction 字段。
- 完整人工决策 TSV；每个对象 cluster 必须恰好出现一次。
- broad、fine、可选 state 列，以及 `decision=confirmed`。

**输出**

- `annotated_object.qs`：写入 broad/fine/可选 state 的派生对象。
- `cell_annotations.tsv`：完整逐细胞标签。
- `annotation_summary.tsv`：每簇决定与细胞数。
- `annotated_umap.pdf`。
- `cluster_sample_condition_umap.pdf`。
- `session_info.txt`、`run_manifest.json`。

应用配置见 [config.apply.example.json](skills/08-scrna-annotate-cells/references/config.apply.example.json)，决策表要求见 [annotation-review.md](skills/08-scrna-annotate-cells/references/annotation-review.md)。部分、重复、空白或未确认的决策表会被拒绝；注释过程不删除细胞，也不覆盖 cluster ID。

### 09-scrna-export-subset：导出目标细胞子集

**什么时候用**

从已注释对象中选择一个或多个正式标签，生成一个可交给下游 QC、重聚类、评分或 DE 的派生对象。例如导出 vascular EC，但不在此步骤完成 EC 重分析。

**主要输入**

- `input.object`：通常为 `08` 的 annotated object。
- `metadata.sample`。
- `metadata.cell_type`：用于筛选的正式注释列。
- `subset.include`：精确 inclusion labels。
- 新的 `output_dir`。

配置模板：[config.example.json](skills/09-scrna-export-subset/references/config.example.json)。

**主要输出**

- `subset_object.qs`。
- `subset_counts.mtx`：raw-count 稀疏矩阵。
- `subset_metadata.tsv`。
- `subset_summary.tsv`：按样本等维度的保留汇总。
- `features.tsv`、`barcodes.tsv`。
- `run_manifest.json`。

父对象不会被修改。导出后的 QC、整合、聚类、marker 和注释必须使用相应 skill 单独完成。

### 10-scrna-score-programs：基因集和通路程序评分

**什么时候用**

用于 Hallmark、代谢、信号通路和自定义 signature 的逐细胞评分。支持 VISION、AUCell、UCell、AddModuleScore 和 PROGENy；支持 inline、GMT、MSigDB 和 scMetabolism 资源。

**主要输入**

- Seurat 对象、assay、layer 和 species。
- 一个或多个显式 `tasks`；每个 task 指定 name、method、gene-set source 和 coverage policy。
- 可选 `summarize_by`，建议至少包含 sample，并按需要加入 cell type/condition。
- 随机种子、cores、资源缓存目录和输出格式。

配置模板：[config.example.json](skills/10-scrna-score-programs/references/config.example.json)。方法选择和 layer 要求见 [scoring-methods.md](skills/10-scrna-score-programs/references/scoring-methods.md)。

**主要输出**

- `program_scored_object.qs/.rds`：每个任务作为 namespaced assay 附加到派生对象。
- `scores/<task>_scores.tsv.gz`：完整 cell × program 分数矩阵。
- `scores/<task>_<cache-key>.rds`：可复现缓存。
- `signature_coverage.tsv`：输入、匹配、缺失基因和 coverage status。
- `assay_feature_mapping.tsv`：导出 signature 名和 Seurat assay feature 名映射。
- `score_summary.tsv`：按配置字段汇总的描述性结果。
- `figures/<task>_group_mean_heatmap.png`。
- `task_manifest.json`、`session_info.txt`、`run_manifest.json`。

不同方法的绝对分数不能直接比较。按 condition 的正式推断仍必须以独立样本为统计单位，见 [interpretation-and-inference.md](skills/10-scrna-score-programs/references/interpretation-and-inference.md)。

### 11-scrna-run-differential-analysis：差异表达分析

**什么时候用**

用于全对象、细胞类型、亚型或 cluster 的条件比较。正式推断默认使用 raw counts 按 `sample × population` 聚合的 pseudobulk DESeq2。Seurat Wilcoxon/MAST/LR 只作为明确标注的细胞级探索性分析。

**主要输入**

- 含 raw counts 的 Seurat `.rds/.qs`。
- `metadata.sample`、`metadata.condition`，可选 `metadata.population` 和 covariates。
- `population.mode=all|selected`、include/exclude labels。
- 一个或多个 comparison：`id`、`numerator`、`denominator`。
- `analysis.method=pseudobulk_deseq2`、assay、设计公式和样本/细胞/基因过滤下限。
- `padj_threshold`、`lfc_threshold`。
- 默认 `lfc_shrink=true`，要求 apeglm 成功；正 log2FC 始终表示 numerator 高于 denominator。

配置模板：[config.example.json](skills/11-scrna-run-differential-analysis/references/config.example.json)，设计规则见 [statistical-design.md](skills/11-scrna-run-differential-analysis/references/statistical-design.md)。

**主要输出**

根目录：

- `design_audit.tsv`：重复、设计和混杂审计。
- `task_status.tsv`：每个 population × comparison 的状态。
- `all_comparisons.tsv`：所有任务的完整基因结果。
- `significant_all_comparisons.tsv`：展示/交接用显著结果子表。
- `DEG_count_summary.pdf`。
- `sessionInfo.txt`、`run_manifest.json`。

每个 `comparisons/<population>__<comparison>/`：

- `sample_cell_counts.tsv`、`sample_design.tsv`。
- `effect_size_audit.tsv`：apeglm 是否请求、应用和回退原因。
- `pseudobulk_data.rds`。
- `all_genes.tsv`：保留所有返回基因和 tested/filter_reason。
- `significant_genes.tsv`、`upregulated_genes.tsv`、`downregulated_genes.tsv`。
- `volcano.pdf`、`MA_plot.pdf`、`pseudobulk_PCA.pdf`、`top_DE_heatmap.pdf`。
- 任务失败时保留 `ERROR.txt`，其他成功任务不被丢弃。

apeglm 只替换 log2FC 和 lfcSE；p 值、padj 和 Wald stat 仍来自未收缩 DESeq2 Wald 检验。完整字段与状态契约见 [input-output-contract.md](skills/11-scrna-run-differential-analysis/references/input-output-contract.md)。

### 12-scrna-run-pathway-enrichment：独立 ORA/GSEA 富集

**什么时候用**

用于读取 `11` 的完整差异结果，或外部 CSV/TSV/TXT/XLS/XLSX 差异表，运行 GO-BP/MF/CC、KEGG、Reactome 和 Hallmark ORA/GSEA。它不会重新计算 DE。

**主要输入**

- 一个 `input.differential_table`，或多个 `input.differential_tables`。
- `analysis.stage=enrichment_only`。
- species 和 `gene_id_type`。
- 明确的列映射：gene、log2FC、pvalue、padj、stat、population、comparison。
- ORA 的 padj/LFC 阈值与最小输入基因数。
- GSEA gene-set size 范围。
- 请求的数据库，以及绘图 top-N、标签换行宽度和每页 term 数。

配置模板：[config.example.json](skills/12-scrna-run-pathway-enrichment/references/config.example.json)。

**主要输出**

根目录：

- `task_status.tsv`。
- `enrichment_all_comparisons.tsv`：所有成功任务的完整合并结果。
- `sessionInfo.txt`、`run_manifest.json`。

每个 `comparisons/<task>/`：

- `standardized_input_table.tsv` 和 `input_column_mapping.tsv`：外部输入标准化与列映射审计。
- `enrichment/`：富集结果目录，其中包含：
  - `gene_id_mapping.tsv`：成功映射的输入 ID；可与标准化完整输入表核对未映射 ID。
  - `enrichment_status.tsv`：逐数据库 completed/empty/skipped/failed。
  - 每个数据库、方法和方向的完整 `*_full.tsv`。
  - `enrichment_dotplot_overview.pdf`：每个数据库/方法/方向最多 3 条的固定尺寸总览。
  - `enrichment_dotplot_<database>_<method>[_pageN].pdf`：标签换行、尺寸受限并自动分页的细节图。
  - `gsea_nes_<database>[_pageN].pdf`：各数据库独立 GSEA NES 图。

单个数据库失败不会删除其他数据库的成功结果。绘图只是完整 TSV 的摘要视图；默认每页最多 20 条，可通过 `enrichment.plot_terms_per_page` 调整。设计说明见 [enrichment-design.md](toolkit/references/enrichment-design.md)。

### 13-scrna-test-cell-abundance：细胞组成和局部丰度变化

**什么时候用**

在 QC、聚类和人工确认注释完成后，比较 control、PAH 或其他条件之间的细胞类型/亚型相对比例。也可用 Milo 在连续表达空间中寻找不完全服从 cluster 边界的局部富集或减少状态。正式统计重复始终是样本或供体，不是细胞。

在运行前必须先写清楚分母。例如：

- 全肺对象以全部通过 QC 的细胞为分母：回答“EC 在捕获细胞中的相对比例是否变化”。
- 内皮富集对象以全部确认的 EC 为分母：回答“动脉、静脉、毛细血管亚型在 EC 内部的相对组成是否变化”。
- `selected_cell_types` 只保留指定父群后重新计算分母。不得把这类结果写成绝对组织细胞数、细胞生成或死亡率。

**两种输入方式**

1. `input.object`：带 raw counts、sample、condition 和正式 cell-type 标签的 Seurat `.rds/.qs`。Milo 必须使用这种输入，并要求对象已有明确命名的 reduction。
2. `input.counts_table`：UTF-8 TSV，至少包含 sample、condition、cell_type 和非负整数 `n_cells`。重复的 sample×cell_type 行会求和；缺失组合会补成零。适合复用既有汇总计数，不要求复制大型对象。

每个样本必须唯一对应一个 condition 和一组样本级 covariates。comparison 中正效应始终表示 numerator 相对 denominator 增加，但不同方法的效应尺度不能直接相加或平均。

配置模板：[config.example.json](skills/13-scrna-test-cell-abundance/references/config.example.json)。完整契约见 [input-output-contract.md](skills/13-scrna-test-cell-abundance/references/input-output-contract.md)。

**五种方法如何选择**

| 方法 | 实际调用 | 回答的问题 | 优点 | 主要不足 | 推荐角色 |
|---|---|---|---|---|---|
| propeller | `speckle::propeller.ttest` | 每个已注释群的样本比例是否改变 | 快、透明、支持协变量和 limma 经验贝叶斯；适合作为可复核基线 | 逐群检验；零比例和稀有群受转换影响；展示的 log2 平均比例比是描述性效应 | 常规主分析或敏感性分析 |
| sccomp | `sccomp_estimate` + `sccomp_test` | 在组成约束、过度离散和异常样本下，细胞群组成是否改变 | 稳健多元 Beta-binomial；不要求人为选参考群；可扩展到差异变异性 | 依赖 CmdStan；首次模型编译较慢；posterior FDR 不是 BH padj | 注释级组成的优先正式模型 |
| scCODA | `pertpy.tl.Sccoda` + NUTS | 相对一个参考群，哪些细胞群有可信组成改变 | 联合贝叶斯组成模型；适合稀疏变化；输出后验纳入概率 | 结论依赖参考群；NUTS 慢；必须检查诊断；credible 不是传统 P 值 | 小样本且参考群可解释时补充 |
| DCATS | `DCATS::dcats_GLM` | 在 Beta-binomial 模型下检验组成，可选校正注释误分类 | 可使用细胞类型相似度矩阵；支持协变量；通常较快 | 没有可信相似度矩阵时优势有限；极端零计数和小样本可能收敛不稳 | 注释不确定性/误分类敏感性分析 |
| Milo | `miloR` 邻域计数 + edgeR | 连续低维空间中的局部细胞状态是否富集/减少 | 不依赖硬 cluster；可发现亚型内部激活态或过渡态；样本感知 | 依赖 reduction、批次校正、k、d 和邻域抽样；邻域互相关；不等价于整体细胞类型比例 | 与注释级组成并行的局部状态分析 |

通常使用：

- 一般注释级组成：`sccomp + propeller`。
- 每组只有约 3–5 个独立样本且参考群可审计：再加 `sccoda`，保持参考敏感性分析开启。
- 相近亚型可能误分：再加 `dcats`，并提供有科学依据的 similarity matrix。
- 怀疑一个 cluster 内部存在 PAH 激活态或连续过渡：单独加 `milo`；不要拿 Milo 替代完整组成表。

更细的统计解释和优缺点见 [methods-and-interpretation.md](skills/13-scrna-test-cell-abundance/references/methods-and-interpretation.md)。实现复用各方法的官方包，不重新实现其统计模型。

**配置关键字段**

```json
{
  "input": {"object": "results/GSE/08_annotation/annotated_object.qs"},
  "metadata": {
    "sample": "sample_id",
    "condition": "condition",
    "cell_type": "annotation_fine",
    "covariates": []
  },
  "comparisons": [
    {"id": "pah_vs_control", "numerator": "pah", "denominator": "control"}
  ],
  "analysis": {
    "methods": ["sccomp", "propeller", "sccoda"],
    "denominator": {
      "mode": "selected_cell_types",
      "include": ["Arterial_EC", "Venous_EC", "Capillary_EC_general", "Capillary_EC_aerocyte"],
      "description": "Relative abundance among confirmed vascular endothelial cells"
    },
    "min_samples_per_group": 3,
    "min_cells_per_sample": 20,
    "fdr": 0.05,
    "random_seed": 154959
  }
}
```

`all_input_cells` 不需要 `include`，但仍必须填写 `description`。Milo 还必须填写 `method_options.milo.reduction/k/d/prop`；scCODA 必须明确 reference 策略；DCATS 的 similarity matrix 必须是行列名与细胞类型一致的方阵。

**运行**

```bash
python3 ~/.codex/skills/13-scrna-test-cell-abundance/scripts/check_dependencies.py \
  --skill 13-scrna-test-cell-abundance \
  --pixi-root /home/faz_laptop/projects/scrna_envs

python3 ~/.codex/skills/13-scrna-test-cell-abundance/scripts/run.py \
  --config config/13_cell_abundance.json

python3 ~/.codex/skills/13-scrna-test-cell-abundance/scripts/run.py \
  --config config/13_cell_abundance.json --execute
```

sccomp、scCODA 或多方法运行建议交给前述 tmux supervisor。`07-cell-abundance` 环境的首次部署命令为：

```bash
cd /home/faz_laptop/projects/scrna_envs/07-cell-abundance
pixi install -e default
pixi run -e default install-sccomp-runtime
pixi install -e sccoda
```

环境文件和 `pixi.lock` 固定 R/sccomp/DCATS/Milo、CmdStan 以及 Python/pertpy/scCODA 依赖；skill 自己不会在分析运行时安装或升级包。

**主要输出**

根目录始终保留：

- `sample_cell_counts.tsv`、`sample_cell_proportions.tsv`：完整样本×细胞类型计数和比例。
- `design_audit.tsv`：每个样本的条件、协变量、输入细胞数、分母细胞数和最低细胞数检查。
- `cell_type_eligibility.tsv`：每群总细胞数、存在的样本数和达到最低细胞数的样本数。
- `task_status.tsv`：每个 comparison×method 的 completed、failed、missing_dependency、invalid_design 等状态。
- `all_method_results.tsv`：所有成功方法的完整标准化长表；`significant_method_results.tsv` 只是筛选视图。
- `method_concordance.tsv`：按细胞类型记录 `supported/partial/discordant/not_supported`；scCODA 多个参考敏感性运行只算一个方法证据。
- `sample_composition.pdf`、`cell_type_proportions_by_condition.pdf`、`sample_proportion_heatmap.pdf`：自动分页、样本点带标签的描述图。
- `sessionInfo.txt`、`run.log`、`run_manifest.json`。

每个 `comparisons/<comparison>/<method>/` 还保留官方原始结果、完整标准表、显著子表和分页 `effect_summary.pdf`。scCODA 额外保存逐参考结果与 posterior diagnostics；sccomp 保存抽样/Pathfinder 产物；DCATS 保存系数和似然比；Milo 保存完整 neighborhood 表及 DA 图，可配置是否保存 Milo 对象。

**如何解释跨方法结果**

- `supported`：至少两种不同方法达到各自证据阈值且效应方向一致。
- `partial`：仅一种方法达到阈值。
- `discordant`：达到阈值的方法方向冲突，应优先审计分母、零计数、参考群、异常样本和模型诊断。
- `not_supported`：本次方法均未达到阈值；不等同于证明“绝对无变化”。

propeller/DCATS 的调整 P 值、sccomp posterior FDR、scCODA posterior inclusion/credible、Milo spatial FDR 是不同统计证据。汇总表只比较方向和是否达到各自阈值，不把数值伪装成同一尺度。

## 五、版本 2 到版本 3 的名称迁移

版本 3 保留既有科学计算，同时重新划分职责并增加 QC 过滤和独立富集入口：

| Version 2 | Version 3 |
|---|---|
| `04-scrna-preprocess-and-cluster` | `06-scrna-preprocess-and-cluster` |
| `06-scrna-find-cluster-markers` | `07-scrna-find-cluster-markers` |
| `07-scrna-annotate-cells` | `08-scrna-annotate-cells` |
| `08-scrna-analyze-subset` | `09-scrna-export-subset` |
| `09-scrna-score-programs` | `10-scrna-score-programs` |
| `10-scrna-run-differential-analysis` | `11-scrna-run-differential-analysis` |

新增：

- `04-scrna-apply-qc-filter`：复用项目中已经验证的逐细胞决策表过滤模式。
- `12-scrna-run-pathway-enrichment`：复用原差异分析中的富集实现，提供独立入口。
- `13-scrna-test-cell-abundance`：新增样本感知的注释级组成与 Milo 邻域丰度分析；没有删除或替代既有差异表达、富集和可视化输出。

迁移没有删减既有 DE、marker、subset、program scoring 和 enrichment 的完整计算表。旧版安装目录不应与新版并存，否则相同任务可能被重复或错误路由。

## 六、开发、验证和打包

共享 Python/R 运行时代码只维护在 `toolkit/`。各 skill 的最小运行时由 [release/runtime-manifest.json](release/runtime-manifest.json) 声明，安装和打包时自动注入；不要直接编辑安装后的生成副本。

基础验证：

```bash
python3 scripts/validate_runtime_manifest.py
python3 scripts/smoke_install.py \
  --quick-validate ~/.codex/skills/.system/skill-creator/scripts/quick_validate.py
python3 tests/test_tmux_runner.py
```

完整 fixture E2E：

```bash
python3 tests/e2e/run_fixture_e2e.py
```

生成 `.skill` 包：

```bash
python3 scripts/package_skills.py --output dist
```

确定性 fixture 生成器位于 `tests/fixtures/create_fixture.R`，包含 80 个细胞、4 个样本、2 个 condition、2 个 batch、2 个 cell type、2 个 cluster、整数 counts、QC metadata 和 UMAP。测试必须保证 fixture 运行前后 SHA-256 不变。

仓库同步、发布门禁和 tag 规则见 [AGENTS.md](AGENTS.md)。
