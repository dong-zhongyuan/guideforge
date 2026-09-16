# 预注册登记：细胞杀伤 panel v5（2026-09-09 封存）

> 本登记先于任何 v5 panel 实验结果入库。逐构象预测数值见
> `data/prereg_cellpanel_v5_predictions.json`（生成命令：
> `python scripts/crrna_wetlab_panel_cell.py --prereg`，来源全部为既有
> 干实验产物：Chai 4×8 界面矩阵、16 监督对选型器、同源 trans 终点、
> Sanger-54 耐受邻域、Scholz 标定虚拟细胞、功效预分析）。

## Panel

2 靶 ×（WT + 3 变体）= 8 条（`wetlab_panel_cell8_order.csv`）。
TP53-R248Q（HCT116）与 TP53-R273H（SW480）；变体 = A1C / A1U+A2C / A8C+U15G。

## 预注册假设（H1–H4）

- **H1 非劣性**：A1C 与 A1U+A2C 的杀伤不劣于 WT（界面中性 + 同源顶档；
  预期 |Δkill| < MDD）。
- **H2 茎稳定化机制**：A8C+U15G 若茎稳定化转化为杀伤增强，Δkill>0；
  若蛋白界面主导，Δ≈0 或为负。
- **H3 等位一致性**：两靶上 4 构象排序方向一致（基因层同集假设）；
  不一致 = 等位水平分型信号（探索性，不设显著性口径）。
- **H4 虚拟细胞校验**：实测存活落点 vs Scholz 标定预测
  （HCT116 70.4% / SW480 45.0%，含 EC50 CI 场景带与 20pp 移植边界）。

## 预注册排序假设

杀伤两档：第一档 {WT, A1C, A1U+A2C}（trans 顶档+选型器顶档+界面中性），
第二档 {A8C+U15G}（trans 弱+界面轻降）；**档内不做排序预测**。

## 判定口径（先于数据定死）

- 显著性：靶内成对单侧 t（α=0.05），重复数按实际；**差异 < MDD 一律
  方向性叙述、不下显著性结论**。MDD 与重复数扫描见
  `data/power_analysis_cellpanel.json`（n=3→37.6pp；n=5→26.1pp；n=8→19.7pp；
  预测构象差量级 10–30pp——建议 n≥5）。
- 变体判优：Δkill ≥ MDD 且两靶同向。
- 等位选择性：SI = 突变株杀伤 / WT 株对照杀伤；变体 SI ≥ WT 的 SI 判非劣。
- 分析命令：`crrna_ivt_template.py --anova <结果>`（2×4 交互）；
  `--twin-check <结果> --cell-csv <SI表>`（H4 校验）；
  `crrna_bayesopt.py --ingest-cell <结果>`（回流提案）。

## 边界声明

- 虚拟细胞存活预测为靶基因×细胞系级（模型无骨架项），骨架级差异由
  界面/选型器/trans 三列承载；
- Scholz 标定为 RNP 电转体系，移植到质粒/U6 体系保留 20pp 边界；
- 本登记不预设 A8C+U15G 的方向（H2 双向开放）。
