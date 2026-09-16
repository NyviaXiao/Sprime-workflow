# SPrime workflow (Linux)

VCF → 群体 SPrime → 多古人匹配 → 汇总、GMM、个体单倍型 call。
普通 Linux 服务器即可，不要求 SLURM。所有原始代码保存在工程外；本目录是独立实现。

## 安装与运行

```bash
cd Sprime-workflow-main
conda env create -f environment.yaml
conda activate sprime-workflow
cd tools/map_arch
make -B

cd ../..
python -m unittest discover -s tests -v
python run.py run --config config/config.yaml --cores 8 --memory-mb 32000 --dry-run
python run.py run --config config/config.yaml --cores 8 --memory-mb 32000
```

先修改 `config/config.yaml` 的 VCF、遗传图谱、古人 VCF/mask 路径与 `config/samples.tsv`。
配置中的相对路径以配置文件所在目录为基准；命令行路径以当前目录为基准。
VCF 支持单文件或 `chr{chrom}.vcf.gz` 模板；在 shell 中请给模板加引号。
现代人输入通过流式筛选读取，无需预先建立索引；群体子集自动生成 CSI 索引。
古人 VCF 和 mask 必须逐染色体存放，每个古人 VCF 只有一个样本。
染色体命名必须与遗传图谱一致。参考版本由使用者确认，软件不能从 VCF 自动证明组装版本一致。

```bash
python run.py run --config config/config.yaml \
  --vcf '/data/chr{chrom}.vcf.gz' --samples /data/samples.tsv \
  --outdir results/run002 --cores 16 --memory-mb 64000
```

相同命令可断点续跑。`--target matches|gmm|individual` 可只运行对应依赖链。
默认一任务一 CPU；内存预算用于限制并发，不是操作系统硬限制。
map_arch 按染色体坐标分配大量内存，须按实际服务器调整 `resources.map_arch_mem_mb`。
工作目录中保留群体 VCF，供 individual call 和重跑使用；不自动删除大文件。

## 输入与规则

`samples.tsv` 是带表头的 TSV：`sample_id, population, role`；role 为 target 或 outgroup。
任意数量的尼安德特/丹尼索瓦人写入 `archaic_references`；ID/tag 必须唯一。
每个古人从同一个 SPrime score 独立匹配，保留每位古人的计数和 match rate。

汇总分母是 `match + mismatch`，`notcomp` 不计入；默认分母至少 10，低于阈值为 NA。
内部保存完整精度，旧 R 脚本四舍五入到四位，因此展示精度可能不同。
群体 pooled rate 为分子之和/分母之和（按片段记录累计，非全基因组去重估计）。

GMM 预处理从 `SCORE > 150000` 的位点重新计数。每个 reference 只有在
callable 数至少 30 时才参与该片段的阈值判断；每个 Nean/Deni 组至少要有一个
合格 reference。判断使用合格 reference 的最大 rate：max(Nean) < 0.3 且
max(Deni) > 0.3。每个 target Denisovan 都使用自己的 callable 合格片段单独建模，
不对多个 target 求平均或合并。AIC/BIC 输出在 `gmm/model_selection.tsv` 中，仅供参考，
不参与现有 LRT、p-value 或 selected_components 决策。
**不进行片段长度筛选**，`length_bp` 仅为结果描述字段。
基础匹配表、分类表和 GMM 筛选表分别保留，以便核查筛选对结果的影响。

GMM 默认复现现有逐级 1/2/3 成分选择方式：先 1 vs 2，再 2 vs 3 或 1 vs 3。
修正了旧代码将 1 vs 3 的结果标作 2 vs 3 的问题；多次初始化默认 10 次。
Bonferroni 第一阶段按“群体数 × 建模参考数”，第二阶段按其 2 倍校正；
包含数据不足的预定检验。输出原始 p、校正 p、比较名称、均值、权重及标准差。
`legacy_chi_square` 保留旧统计方法；有限混合模型不满足常规卡方 LRT 的一般条件，
可选 `parametric_bootstrap` 校准混合模型比较（并不校准上游筛选造成的偏倚）。
GMM 成分数是分布模型的结果，不能单独当作已证实的历史渗入事件次数。
少于 10 个片段或少于 3 个不同数值时输出 `insufficient_data`，不伪造 p 值。

分类使用配置中的 Nean/Deni reference 列表，并分别取可用 rate 的最大值。
Nean: max(N) > 0.6 且 max(D) < 0.4；Deni: max(D) > 0.3 且 max(N) < 0.3；
所有未命中这两类的片段（包括 match rate 缺失）归为 ambiguous。

个体 call 直接按 `(CHROM, POS, REF, ALT)` 对齐 VCF 与 score，在每个 SPrime 片段内，
分别查找每个目标个体两条单倍型连续携带 SPrime ALLELE 的 marker run，默认至少 2 个。
“连续”指相邻 SPrime marker；没有最大物理距离限制。缺失 GT 中断 run；
非缺失 GT 必须为 phased diploid，且相位应在分析区间内一致。本流程不执行定相。
outgroup 不进入个体 call。样本名与 haplotype 是独立列，避免切割样本名中的点。
坐标表使用 1-based inclusive；BED mask 使用 0-based half-open。个体边界是首末支持 marker，
并非碱基级断点估计。重复运行覆盖目标输出，不进行旧代码的追加写入。

## 输出

```text
results/run001/
  qc/validation.tsv
  resolved_config.json
  run_info.json
  sprime/{population}/chr{chrom}.score
  archaic_match/{population}/{reference}/chr{chrom}.mscore
  tables/segment_match_rates.wide.tsv.gz
  tables/segment_match_rates.long.tsv.gz
  tables/population_match_summary.tsv
  tables/classification_summary.tsv
  classification/{population}.tsv
  gmm/input/{population}/{reference}.tsv
  gmm/model_selection.tsv
  gmm/components.tsv
  gmm/plots/{population}.{reference}.png
  individual_calls/all_individual_calls.tsv.gz
  individual_calls/by_population/{population}/{class}.tsv.gz
  individual_calls/individual_summary.tsv
  plots/landscape/{population}.introgression_landscape.png
  plots/landscape/{population}.neanderthal_affinity_landscape.png
  plots/landscape/{population}.denisovan_affinity_landscape.png
  affinity/{population}/{reference}.tsv.gz
  plots/contour/{population}.{reference1}__{reference2}.png
  adaptive/{variants,core_variants}.tsv.gz
  adaptive/segment_summary.tsv
  adaptive/top2_candidates.tsv
  adaptive/shared_top2_regions.tsv
  provenance/{run_manifest,resolved_config}.json
  provenance/{software_versions,input_manifest,output_manifest}.tsv
  report/report.html
  logs/
  work/
```

wide 表每个古人包含 matched/mismatch/callable/notcomp/match_rate；long 表每行一个片段-古人。
个体表包含 population/sample_id/haplotype/chromosome/start/end/length_bp/marker_count/segment_id/archaic_class。
`adaptive/segment_summary.tsv` 与 `adaptive/top2_candidates.tsv` 的
`candidate_start/candidate_end` 使用 BED-style 0-based half-open 坐标；其余既有
segment summary 表保持原有的 1-based inclusive 坐标定义。
provenance 精简为配置和 run_info：工具版本、小型关键文件 SHA256、输入路径/大小/修改时间。
QC 会扫描古人 VCF/mask，检查单染色体限制；完整验证耗时随参考文件大小增加。

landscape、affinity、adaptive、provenance 和 report 是独立 downstream 分支，
直接消费已有结果；它们不会重新计算 SPrime 或 map_arch。Contour 支持
`all_pairs`（全部古人组合）和 `explicit_pairs`（仅配置组合）。

## 实现与验证范围

```text
Sprime-workflow-main/
  run.py                       # Linux 命令行入口
  environment.yaml             # Conda 运行环境
  config/                      # 路径、参考列表、阈值、样本表
  workflow/
    Snakefile                  # 总入口和最终目标
    rules/
      validate.smk
      prepare_vcf.smk
      sprime.smk
      archaic_match.smk
      summarize.smk
      classify.smk
      gmm.smk
      individual_calls.smk
      report.smk
      new_modules.smk
      tools.smk
    scripts/
      task.py                  # 公共执行入口和日志
      core.py                  # 汇总、筛选、分类等共享函数
      validate_inputs.py
      call_individual_tracts.py
      run_gmm.py
      plot_contour.R
      build_affinity_tables.py
      plot_landscape.R
      plot_affinity_landscape.R
      run_adaptive.py
      collect_adaptive.py
      build_provenance.py
      build_report.py
  tools/                       # JAR、map_arch 源码及编译文件
  tests/
  legacy/                      # 原始分析脚本的只读参考副本
```

`rules/` 按分析阶段拆分依赖，`scripts/` 保存实现；`plot_contour.R` 保留 MASS::kde2d 绘图。
Snakemake 的代码签名按 stage 计算：入口会关闭同一 dispatcher 脚本的整体 code trigger，
改由 stage-local `params.code` 追踪；因此只修改 GMM fitting 不会使上游匹配或 classification
失效，classification 改变会使 individual 分支失效，但不会触发 GMM。
原 `score_summary.r/pre_data.r/selectGene2.py` 的重复处理已合并；原文件副本在 `legacy/`。
`tools/map_arch` 保留提供的 C 程序接口，并修复 VCF 表头读取、未初始化值及缺失 GT 处理；
Linux 上须重新编译。由于修复会影响异常/缺失输入的行为，应与旧结果逐项比较。
本工程未附现代人/古人/遗传图谱真实数据；真实 Linux 全链路和科学结果一致性需用研究数据验收。
