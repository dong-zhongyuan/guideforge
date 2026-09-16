# 预登记 · 细胞杀伤 panel v6（2026-09-16 封存）

> v6 = v5（2026-09-09 登记，见 `preregistration_cellpanel_v5.md`，全文有效并存档）
> 的面板收窄版：**每靶 WT + 2 变体（A1C、A8C+U15G），共 6 条**；判读逻辑
> H1–H4 与判定口径不变。封存时点先于任何 panel 实验数据（v5 封存后无任何
> 实验数据回填，本登记因此合法）。产物：`data/prereg_cellpanel_v6_predictions.json`
> （`scripts/crrna_wetlab_panel_cell.py --prereg` 生成）。

## 与 v5 的差异（仅面板构成）

- A1U+A2C 退出面板（用户 2026-09-16 指示收窄），其证据与 v5 产物保留于 git 历史；
- H1 表述由"A1C 与 A1U+A2C 不劣于 WT"收窄为"A1C 不劣于 WT"；
- H3 等位一致性的构象数 4 → 3；
- ANOVA 交互因子 2×4 → 2×3；
- 回填模板升级为 n=5 重复列（对齐 `power_analysis_cellpanel.json` 的 rep_scan 建议，
  n=5 时 MDD≈26.1pp，配对+时序 AUC 口径 10–15pp）。

## 2026-09-16 JD12 事件与分型学说核验（先于实验，如实入档）

1. 用户提交 JD12-Tp53-HT29 两条 23nt spacer。核验（`scripts/crrna_jd12_redesign.py`，
   `data/jd12_redesign.json`）：两条均为 **TP53 R273H 等位特异**（c.818G>A），
   即 SW480 突变；HT29 为 R248W，标签不符。crRNA-1 的 PFS（CCUGG）落在 Scholz
   实证景观红旗区（dep 0.628），弃用；crRNA-2（GGGAG，dep 3.90，前 21.4%，
   全转录组 0 完全脱靶）可用，但等位碱基在 spacer 第 17 位（远端），选择性弱于
   canonical R273H 设计（PFS 鉴别 4.68×）。**用户确认细胞系为 SW480**，
   v6 主表维持 canonical spacer；JD12-crRNA-2 × 3 骨架列为**可选附表**
   （`data/wetlab_jd12_sp2_order.csv`），**不在本预登记判读口径内**。

2. 分型学说核验（用户指出：DR×spacer 互作特征承载分型证据，换 spacer 必须重分型；
   `scripts/crrna_panel_retyping.py`，`data/panel_retyping.json`）：
   312 条队列复算与 §A3 记录逐位一致（簇大小 [137,56,119]）。
   - canonical R248Q → **型0**：A1C 3/3 代表进 TOP16（最好第 3）、
     A8C+U15G 3/3（最好第 5）——两变体与分型学说一致；
   - canonical R273H → **型2**：A1C 3/3（最好第 2）一致；
     A8C+U15G 仅 1/3（第 5）——**弱支持，如实标注**（其入板依据为四要素复合分
     与 H2 机制覆盖，非分型 TOP；预登记预测本就为 tier-2 ≤WT）；
   - JD12-crRNA-2 → **型1**（与两个 canonical 分属不同型）：A1C 仅 1/3、
     A8C+U15G 0/3——**分型学说不支持在型1 上下文沿用这两条变体**；
     若正式启用 JD12 附表，需按型1 代表 TOP 谱（WT/U7C/A1G/U13C 等）
     另行登记变体集。

## 判读规则（与 v5 相同，摘要）

靶内成对单侧 t（α=0.05）；Δkill ≥ MDD 且两靶同向判优；SI 等位选择性非劣判；
差异 < MDD 只做方向性叙述。分析命令见预测 JSON 的 decision_rules。
