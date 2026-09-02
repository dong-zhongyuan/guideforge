# GuideForge — Cas12a2 crRNA 骨架优化管线

AI 辅助优化 **Cas12a2 crRNA 的骨架（direct-repeat 茎环）**，用于更高效、更特异的癌细胞杀伤。
设计变量只有骨架 DR 区；spacer 向导序列固定为 TP53-R248Q 靶向序列（24nt,
`GTTCATGCCGCCCATGCAGGAACT`, 来自 Zeng et al. 2026 Nature 补充表, 为该文最优选择性向导)。

## 管线流程

```
注册表(WT 骨架 + 固定 spacer; 双路置信门控 ViennaRNA × RNet-SS)
  → 变体生成: A 突变扫描 | B 模拟退火 | A2 Struct2SeQ 逆折叠 | A3 gRNAde 三维几何逆折叠 | A4 茎区共变枚举
  → 硬过滤(结构保持/spacer 游离/polyT·G/加工位点保护*)
  → 多维打分: ΔΔG · 结构距离 · 蛋白接触(8D4A) · 氢键互补 · 3'保守窗惩罚* · RNet SHAPE 一致性
  → TOP-K + WT 候选库(FASTA 可送合成)
  → 特异性扫描(WT vs 突变转录本, PFS 鉴别)
  → 湿实验闭环(体外旁切初筛 → 细胞杀伤验证)

平台扩展(v1.1): spacer 上下文分型 → 设计智能体
  86 条文献 spacer 按全长折叠特征聚 3 型, 各型最优 DR 位点图谱不同
  → 智能体: 输入 spacer(项目口径)或突变转录本(换靶标) → 骨架分型选择/tilling 设计
  * 加工位点保护与 3'保守窗惩罚的文献依据见下"Dmytrenko 2023"条
```

### 状态竞争模型与交叉特征去混杂(2026-09-01 联审方案落地)

- **活性态-侵占态竞争模型**(`scripts/crrna_state_competition.py`, ViennaRNA 单能量标尺约束配分, hc_add_bp/up): 预注册判别检验 **FAIL**——14 条中 Sp8 的 dG_comp 仅排 6/14; 互斥定义(强制侵占螺旋+禁死茎臂)下所有可枚举侵占态均比茎态贵 5.8+ kcal(含 Sp8), 论文失活折叠在伪结外能量面上不可达。**结论: Sp8 正式定性 `out_of_model_domain`**, 活性态竞争假说的检验需统一假结能量模型(中期, RNAstructure/统一 pk 模型)。
- **交叉特征阶段拆分+组成校正**(`scripts/crrna_cross_decompose.py`, 单碱基组成匹配 N=200 + LOO + GC 偏相关): raw cross_coreDR rho=+0.55(LOO 稳定) 在组成匹配后降至 +0.17~+0.21、GC 偏相关后 +0.08~+0.18; flank 信号 z 后归零(纯成分代理, 且配对对象是被加工丢弃的 U-rich 残留 repeat)。**判据 R1: raw cross_pp 退出主评分获数据支持**; 阶段拆分保留(flank 只进加工模块口径)。
- **v1.7 结构置换口径过滤**: `inv_max_run`(最长连续 DR-spacer 侵占螺旋) <= 同上下文 WT(t1 型 WT 自身 6 对连续侵占故禁用固定阈值) + `partner_switch` 标注(t1 诊断: A8G 类 DDR-alone 稳定化突变在型1 口径强化的恰是侵占螺旋) + `ddG_dr` 语义更名 "DR 单独折叠稳定化(DDR-alone 口径)" + top.json 增加 `model_domain` 声明(伪结外: 仅侵占筛查+同 spacer 相对排序)。
- **补偿突变实验面板**(`scripts/crrna_compensatory_design.py` -> `data/ivt_compensatory_panel.*`): MYCg1 三臂 C_WT(inv_run=6)/A_break(0, GC 不变, 靶向配对保持)/B_break_compensate(6 恢复), 各配同源靶——A vs B 活性差 = 折叠竞争的因果检验, 寡核苷酸级成本。

### NUPACK 3.2.2 含假结配分检验(第四模型, 2026-09-02)

`scripts/crrna_nupack_pk_weight.py`(NUPACK 3.2.2 官方 pfunc/-pseudo, Dirks-Pierce 假结配分; 源码取自 sbi-rostock/nupack-serve vendor 树, 容器 -fcommon 编译; 手册确认 pairs/mfe 不支持 -pseudo)。指标 W_pk=(Z_pk-Z_nopk)/Z_pk(假结构象权重占比)。**Sp8 判别 FAIL(第四模型)**: W_pk 排 5/14, 且与活性 Spearman=-0.327(方向相反)。至此四个独立模型(ViennaRNA 约束配分/ThreshKnot bpp/Knotty DP09-MFE/NUPACK Dirks-Pierce 配分)全部不能复现论文 Sp8 失活态——**out_of_model_domain 为四模型一致的最终口径**; 该机制需蛋白环境或实验手段。平台口径注记: GuideForge WT/茎稳定化 TOP3 的 W_pk 仅 0.9%~1.0%(假结权重可忽略), 与 8D4A/Knotty "Cas12a2 成熟 crRNA 无假结"三方互证, 本体系不受该边界影响。

### DP09 假结能量模型竞争检验(Knotty, 2026-09-01)

`scripts/crrna_pk_competition.py`(容器侧 `crrna_pk_runner.py` + Knotty 官方源码编译, CCJ+DP09): 14 条 Creutzburg 前体/成熟两口径 + GuideForge pdbdr WT/TOP3 的含假结 MFE。**判别检验 FAIL(三模型独立一致)**: DP09 能量面上 Sp8 折为"活性样"(干净 6 对茎[寄存器与 ViennaRNA 滑移 2 位——跨模型模板精确匹配无效的实测教训, 分类已改寄存器无关]+spacer 自发夹), 而高活性 guide(Sp2/Sp4/Sp7/Sp9/Sp12)反呈 DR-spacer 交叉+多重假结(方向与论文机制相反)。至此 ViennaRNA 约束配分 / bpp 阈值(ThreshKnot) / DP09 假结 MFE 三条独立路线均不能复现论文的 Sp8 失活态——**Creutzburg 侧 Sp8 out_of_model_domain 判定升级为三模型证据的最终口径**; 仅剩路径 = NUPACK 3.2.2 pfunc -pseudo 约束配分(需注册, 站点可达)。本体系(Cas12a2)不受影响: Knotty 对 GuideForge WT 预测与 ViennaRNA/8D4A 一致(无假结)。

### 结合态活性模板提取(8D4A, 2026-09-01)

`scripts/crrna_active_template.py` -> `data/cas12a2_active_template.json`: 8D4A 链 B(成熟 crRNA, DR18+spacer23) 三维几何判据(N1-N3<4.0A + C1-C1<12.5A)提取——**结合态活性模板无假结**: DR 区仅 5 对标准茎(1基 4-17/5-16/6-15/7-14/8-13) + 3 个近距摆动接触, DR-spacer 交界零配对。含义: Cas12a2 成熟 crRNA 的活性态在 ViennaRNA 伪结外空间内完整可表示(茎即为活性模板), **假结模型缺口只存在于 Creutzburg/FnCas12a 前体回验侧**(其失活态涉及前体 5prime repeat 侧翼, 论文自述活性态含 canonical pseudoknot)。两体系的模型域边界因此不同, 引用时须分开陈述。

### Chai-1 共折叠: 模板注入打通蛋白-RNA 界面(蛋白预测层正式口径, 2026-09-02)

`scripts/crrna_chai_input.py` + `scripts/crrna_chai_collect.py` + 容器 runner(~/dzy/envs/chai, 权重经 hf-mirror): SuCas12a2(1207aa)+crRNA(DR+固定 spacer)+靶RNA 三元共折叠。能力递进三步:

1. **单序列/ESM 基线(历史首跑)**: 无模板下 aggregate 0.23-0.24、蛋白-crRNA 链间 ipTM 0.02-0.08——1200aa RNA 引导核酸酶的已知短板, 界面不可用。
2. **8D4A 自模板注入(提交 e69fc74, 界面特征正式可用)**: 蛋白链自模板 m8(8D4A 链A, 100% 同一)注入后 **aggregate/ipTM 0.25→0.874, 蛋白-crRNA 链间 ipTM 0.02→0.49**——蛋白-RNA 界面恢复到可作选型特征的置信水平; 模板注入是界面恢复的承载步骤(无模板即回退到基线)。
3. **7 骨架 × 5 模型界面差量表(提交 963745d, `data/chai_cofold_matrix_*.json` / `chai_matrix_iptm.json`)**: WT/6 变体的 aggregate 与 prot-crRNA ipTM 模型间一致(sd≤0.007/0.09), 变体间排序稳定; **口径张力如实声明: crRNA-靶RNA 链对分数低(0.02-0.43)且部分变体模型间方差大, RNA-RNA 双链预测的模型稳定性存疑**——界面差量以 prot-crNA 链对为主口径, crRNA-靶对仅作参考; 与 MD 判据(M1-M3 预注册)的结论各自独立陈述。
4. **四靶 × 8 骨架 32 组合矩阵(运行中, `chai_matrix_4t/`)**: 靶 RNA 统一为 protospacer 窗口 24nt + PFS 5nt(8D4A 排布同构), 供选型分类器的"骨架×靶标界面特征差量"输入(策划案 V3 §4.1/§5.1); 出齐后入库 `data/chai_matrix_4t.json`。

### Han 2025 数据集与选型分类器(同源外部数据, 2026-09-02)

`scripts/crrna_han2025_features.py` + `scripts/crrna_train_selector.py`(数据: `data/han2025_dataset.json`, `data/han2025_sanger_sequences.json`)。

**数据链(全部程序化提取, 无硬编码)**:
- Fig1g 活性全表 **145 值**(F×25/S×23/L×48/FL×48 + Canonical; 修正了旧版漏 48 条 FL 的解析 bug);
- 工具箱 9 条序列(CN/F1/F2/L1-L4/FL1/FL2, MOESM3)与 Fig1g 同尺度锚点 7 条;
- **54 条 Sanger 耐受集**(MOESM1 Sup Fig.2/3, 600dpi 视觉转录 ×2-3 次独立读数共识 + 模板守恒校验 + 与 MOESM3/主文锚点交叉验证, 逐条置信度分级; S8 的 20 条与论文自报数完全吻合);
- cis/trans 切割终点比(Sup Fig.5/6 源数据, 三重复均值)。
- **Fig.1f 编号-序列对照表视觉转录(2026-09-03)**: 主文 Fig.1f 以图像形式给出部分映射, 分块放大两次独立逐字母转录 + Sanger 独立源 loop 六联体交叉校验, 恢复 **9 条新监督对**(S3/S5/S7/S10/S14/S16/S20/S21/FL4, 置信度 high 2 / medium 7, 溯源见 `fig1f_pairs`);
- **查证边界(2026-09-03 更新)**: 145 个编号变体与序列的映射**未随任何数值源数据发表**(MOESM7 全部 25 工作表 + MOESM3/4 已穷尽复查; 库为随机合成, Methods 明示), 公开渠道上限=7 工具箱 + 9 Fig1f 转录 = **16 对**; 完整映射待作者回复(docs/ 邮件草稿)。

**选型分类器(先验排序口径, 2026-09-03 用 16 对重训)**: 决策树(Fig1g 尺度, n=16), 特征重要性 p_fold 0.52 / ddG_dr 0.48; 训练集 Spearman 0.729(样本内读数, 仅作参考); **fig1f 留出验证(训练 7 工具箱 → 检验 9 转录对): Spearman +0.25, MAE 0.217(超预登记可用阈 1/3 值域=0.162, 秩方向为正但水平误差偏大——检验本身仅 7 点树, 功效有限, 两口径并列报告)**; 对保守变体仍无区分力(多数候选落同一叶节点)——定位为文献先验排序, 非活性预测。池化模型扩至 n=46(Han 5 终点×7 + Fig1f×9 + Teng 序数对); LOEO 有升有降(rbs33 0.64→0.96, rbs0 0.39→0.68; trans 0.68→0.25, fig1g 0.29→0.04, cis -0.68→-0.86), 8 员族池化排序**翻转**(旧 U5G+A18C +0.58 居首 → 新 WT/A1C/A1U 0.12 居首, U5G+A18C -1.38 末位; 新数据的 loop 同聚物变体均为中等 ddG_dr/p_fold 而活性不差, 强茎稳定化奖赏被压低)——排序对数据增量敏感, 印证"先验排序"而非稳定预测的定位。多终点版: cis/trans 切割终点各训一树, 8 员族逐终点排序入库(`data/selector_family_ranking.json`); **trans 旁切终点上同源骨架间差 30 倍(WT 39.2 vs L1 1.2)**, 为"DR 可调旁切活性"提供了同家族文献证据; 不同终点给出不同排序, 正是"骨架分型"的同源佐证。Sanger 耐受集体检: 真实耐受变体预测中位 0.609(锚点范围内), 204 候选分布中平均 6.5% 分位(模型不排斥真实功能变体)。

### V3 计算侧(第 8 员骨架 / 四靶扫描 / 统计件 / 智能体, 2026-09-02)

- **8 员跨型候选族**(`scripts/crrna_orientation_library.py`): 四取向代表 6 + WT + compensatory 代表 B_break_compensate(zengDR+G11C/G14C/G17C, 补偿臂因果设计), **IVT 模板 32 行(8 骨架 × 4 靶标: R248Q/G12C/G12D/R273H, `data/ivt_round1_template.csv`)**——正合策划案 V3 §5.1 的 8×4=32 组合矩阵口径。
- **四靶转录组脱靶扫描**(GENCODE v47, 385,659 转录本, 统一口径含 TP53 重扫; `data/*_scan_v47.*`): 四条 spacer 的 0 错配位点全部落在本基因异构体(预期靶点); KRAS 两 spacer 在 **KRASP1 假基因各 1 个 1 错配位点**(错配位置与风险评级见扫描表); 无高危脱靶。参考库 gencode.v47.transcripts.fa 不入 git(`data/raw/`, 源 URL 见 summary)。
- **贝叶斯优化同源先验版**(`scripts/crrna_bayesopt.py --prior-han`): Han 工具箱 7 条作 GP 观测(选择器特征空间), 对 204 候选提 EI 建议; 明确标注"LbCas12a CRISPRi 抑制终点, 非 Cas12a2 杀伤, 待 IVT `--ingest` 替换"。
- **虚拟细胞文献先验场景版**(`scripts/crrna_virtual_cell.py --prior-lit`): 依据 Cas12a2 体外激活浓度量级(Dmytrenko 2023 / 2026 Nature)场景化 EC50(TPM) 网格, 输出杀伤窗口表与主验证细胞系激活率(HCT116/SW480 等), 状态标注"文献先验场景, 非标定预测"。
- **矩阵统计件(数据一到即出)**: `scripts/crrna_ivt_template.py --anova`(骨架×靶标两因素方差分析含交互项——"分型假说是否成立"的统计判定, 合成数据自测通过)与 `--twin-check`(数字孪生预测 vs 细胞实测 Spearman/RMSE; 参数未标定时如实报 UN-CALIBRATED)。
- **智能体模块三接入选型器**(`scripts/crrna_agent_webapp.py`): `/api/design` 与 `/api/panel` 输出 8 员族文献先验活性排序(与 `crrna_train_selector.py` 同一模型同一定义, 随机种子固定可复现)。

### 打分输入端可靠性(系综采样检验)

`scripts/crrna_ensemble_check.py`(输出: `data/ensemble_check.pdbdr.v2.json` / `ensemble_check.zengdr.v6.json`): 对 WT+TOP-12(两口径共 26 构建)做 Boltzmann 系综采样检验(ViennaRNA 2.7.2 pbacktrack 本构建不可用, 改用 subopt 精确枚举 5 kcal 窗口 + 权重采样, 窗口尾部质量均 <2.2%, 采样 100/构建, seed=42)——
- **MFE 口径 ΔΔG 与精确系综自由能差 F_var−F_WT(pf 精确值)全部一致, 偏差 ≤0.03 kcal/mol**, "-2.2 kcal/mol 茎稳定化"等数字不是 MFE 单构象假象;
- 稳定化候选(−2.2/−2.0/−1.4)的 ΔΔG bootstrap 95% CI 均不含 0(PDB 口径 [−2.33,−1.83] 等), 系综口径下主张成立;
- DR 茎完整构象占比 0.81~0.93(WT 0.81~0.83), 无构象二态; 单构象能量涨落 sd≈1.1~1.5 kcal/mol 属系综固有涨落, 在自由能差中抵消。
- 边界: 此检验证明的是 ViennaRNA 最近邻模型内的数值可靠性(指标计算端), 不涉及指标与活性的相关性(见"路径A"条, 后者已证明不构成跨体系活性预测器)。

## 文献依据（对应模块）

| 模块 | 文献依据 |
|---|---|
| 加工位点保护硬过滤 / DR 互换不误杀验证 | Dmytrenko et al. 2023, Nature 613:588-594 (ED Fig.2c Cas12a↔Cas12a2 DR 互换功能保持; ED Fig.3 加工切点在 Cas12a 下游 1nt)。位点敏感性×保守性相关验证**仅 score 口径方向一致(ρ=−0.33, p=0.17, n=19 功效低)**, tolerance 口径为零相关(ρ=−0.01)——引用时只可引 score 口径并说明样本量 |
| 3'保守窗惩罚 | Dmytrenko et al. 2023, Fig.1c 跨家族 DR 3' 端保守/loop 可变(定性结论); **窗口大小 --cons3-window=5nt 为项目设定, 非文献出处**; 功能佐证: Zhang et al. 2025 (PMC11780881, DR 3' 端化学修饰可逆调控 Cas12a 活性——非序列突变, 佐证 3' 端敏感性) |
| spacer-DR 互扰 / spacer 游离度硬过滤 | **Creutzburg et al. 2020, NAR 48(6):3228-3243** (PMC7102956, 同家族直接证据): "base pairing of the direct repeat, other than with itself, was found to be detrimental"—DR 与 spacer 碱基配对损害 Cas12a 活性; Sp8 案例证实 spacer 侵占 DR 致假结破坏; 通过 spacer 尾部回折竞争性保护 DR 挽救活性(与本项目 spacer 游离度约束直接对应)。**Cas12a→Cas12a2 为同家族短迁移(DR 均含保守 5' 假结), 非跨体系假设** |
| spacer 上下文分型的机制前提 | 同上(Creutzburg 2020, 同家族); 补充: Bush et al. 2023, Cell Chem Biol 30:879-892 (SpCas9 体系, sgRNA 80nt 大骨架, 跨体系佐证); Liao et al. 2018 (PMC6546362, FnCas12a DR 下游发卡抑制切割)。分型主张的直接证据为本仓库分型试验计算结果; Cas12a2 本体系暂无直接折叠互扰研究(待验) |
| DR/骨架序列工程先例(Cas12a, 存在性证明) | Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434, DR 突变策略调编辑与检测) |
| guide 层修改跨体系迁移的边界 | **Cas9(II 类)与 Cas12a2(V 类)蛋白无同源性, 不存在也不应主张 Cas9→Cas12a2 序列级迁移**(骨架排布相反: Cas9 骨架在 spacer 3' 侧, Cas12a/Cas12a2 在 5' 侧, Kweon 2017 PMC5700056)。可用证据分三层: ①同家族 DR 互换直接证据 = Dmytrenko 2023 ED Fig.2c (Cas12a↔Cas12a2, 本项目迁移链的唯一同源层); ②同一 guide 策略跨两体系的唯一直接演示 = zCRISPR, Nat Commun 2024 (s41467-024-48012-x, 碱基 Z 修饰同时提升 Cas9 与 Cas12a——化学修饰层, 非骨架序列层); ③类比层 = Cas9 侧 sgRNA 骨架工程 (Chen 2013 eSG; Bush 2023) 与 Cas12a 侧 DR 工程 (Han 2025/Zhang 2025/Teng 2019) 各自独立成立, 仅支持"骨架工程调活性为跨 II/V 类重复出现的原理", 不支持设计移植 |
| "按靶标分型选骨架" | **无先例(空白区)**: 本仓库分型试验提供了反对"单一通用骨架"的计算证据(型间最优位点图谱结构性不同); 功能收益待湿实验 |
| 游离态指标→活性 反向验证(路径A, 阴性校准) | **Creutzburg 2020 Fig 1D 逐条反向打分**(本仓库 `scripts/crrna_creutzburg_reverse.py` / `data/creutzburg_reverse_validation.json`): 14 spacer × 构建内真实 DR(反向引物公共尾 revcomp 推得 AATTTCTACTGTTGTAGA±T, 与注册表 zeng2026 条目一致; 差异位点与正文三重吻合: Sp8/Sp12=0基 10,11,19,20, Sp12/Sp14=5,7,9, 位置19=K(G/U))——**本仓库自行复算未能复现该文自报的正相关**: spacer_up ρ=−0.45 / stem_prob ρ=−0.51, 好坏 AUC 0.23–0.38(方向反转); 仅 5端种子区可及性弱正向(seed5_up ρ=+0.41, AUC 0.68, 置换 p≈0.14); 机制个案复现(Sp8 MFE=spacer 自缠发夹, Sp12 全游离); 好 guide Sp4/Sp9 游离态茎概率≈0.001 但活性 72–83%(蛋白重折叠可挽救茎)。对照: 该文自报特征相关为正(pre-crRNA 口径, n≈25, 茎含手工定义假结 P1 配对, ViennaRNA 不可算假结)。**校准结论: 游离态指标仅用于同 spacer 下 DR 变体相对排序(结构口径), 不主张游离态结构预测活性; 排序主张(TOP1>WT)未经实验标定——Creutzburg 层只支持特征类型与活性相关(论文自报 n约25 显著; 自有复算 n=14 不显著、方向混合, 两口径并列引用), MD 阳性仅为必要条件, 最终确认须体外旁切实验**。Sp8 盲区(复算最有价值发现): 唯一 0% 活性样本未被 spacer_up(全长平均稀释)/stem_prob(DR 自配与 DR-spacer 错配混同)抓住, 已据此新增 **DR-spacer 交叉配对硬过滤** v1.6(变体 MFE 交叉配对 nt 不得多于 WT; Creutzburg 设计守则的保守实现, 定位=无损工程规则, 非活性预测器)。**回验结果为阴性且方向反转**: 三口径 cross 均未抓住 Sp8。机制诊断(2026-09-01): Cas12a 家族 crRNA 的活性态含 5' 假结, 在伪结外折叠空间不可表示, Sp8/Sp12 之争本质是与假结活性态的平衡竞争, 伪结外模型按构造无法打分; 次优空间中侵占型嵌套中间态在好坏 guide 均大量出现(Sp8 93/318 vs Sp12 146/391, Sp12 反而更多), 故任何伪结外交叉指标对该案例无区分力。cross_pp 与活性的显著正相关(mature rho=+0.61, p=0.023)为成分混杂: spacer 5' A 富集与被加工丢弃的 U-rich 残留 repeat 配对, GC vs 活性 rho=-0.32 同向, n=14。Sp8 盲区仍开放。过滤回溯审计(`data/cross_audit.json`): 两主线口径 0 拦截/TOP-12 不受影响; **t1_MYCg1 口径 WT 自身 cross_nt=20**(MYC spacer 与 DR 大量交叉配对, 型1 上下文注脚), 拦截 46/118 条含 TOP-12 中 6 条茎稳定化候选——通用候选 A8G+U15C 在型1 口径不再通过 v1.6, 其通用性主张须降级 |

> 引用层级声明: 机制前提的证据链为 **Cas12a 同家族(Creutzburg 2020)→ Cas12a2(本体系, 暂无直接折叠研究)**, 而非跨体系的 SpCas9→Cas12a2; 硬过滤"结构保持/spacer 游离"两项由此获得同家族文献锚定而非纯启发式。 路径A反向复算为阴性(见上行): 游离态全局指标不构成跨体系活性预测器, 打分体系定位=结构口径相对排序+文献机制锚定。

## 快速开始

**环境前置(重要)**: 主管线需要 **Python ≥3.8 + ViennaRNA ≥2.6**(开发验证环境: `/public/home/mengxl/dzy/envs/guideforge`, py3.10.20 + ViennaRNA 2.7.2)。宿主机裸 `python` 是 2.7, 直接跑会报编码错误——请用上述环境解释器或任何满足前置的 py3 环境。路径类默认值(如 rnet 引擎)全部可用 CLI 参数覆盖; 仓库内少数数据 JSON 的绝对路径为溯源记录, 换机重跑时按需替换。

```bash
PY=/public/home/mengxl/dzy/envs/guideforge/bin/python   # 或你自己的 py3+viennarna 环境

# 注册表契约测试 + 主管线黄金集回归(秒级)
PYTHONPATH=src $PY -m unittest discover -s tests

# 最小运行(纯 ViennaRNA, CPU)
$PY scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 --out-prefix data/run1

# 完整四引擎 + SHAPE 筛选
$PY scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 \
    --use-struct2seq --use-grnade --rnet-screen \
    --s2s-python /path/to/engine-env/bin/python --out-prefix data/run2
```

## 一键复算清单(数值结果第三方复核入口)

```bash
PY=/public/home/mengxl/dzy/envs/guideforge/bin/python

$PY scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 --out-prefix data/repro_v6   # 候选库重生成
$PY scripts/crrna_ensemble_check.py           # ΔΔG 系综可靠性(≤0.03 kcal/mol 主张)
$PY scripts/crrna_creutzburg_reverse.py       # 路径A 阴性校准(方向反转主张)
$PY scripts/crrna_cross_decompose.py          # 交叉特征组成去混杂(rho +0.55→+0.17)
$PY scripts/crrna_state_competition.py        # Sp8 状态竞争 FAIL 判定
$PY scripts/crrna_pk_competition.py           # Knotty DP09 FAIL 判定(容器)
$PY scripts/crrna_nupack_pk_weight.py         # NUPACK 第四模型 FAIL 判定
$PY scripts/crrna_negative_control.py         # 19×上下文负对照
$PY scripts/crrna_han2025_features.py         # Han2025 数据链(Fig1g 145/Sanger 54/cis-trans)
$PY scripts/crrna_train_selector.py           # 选型器训练+204 候选应用+耐受集体检
$PY scripts/crrna_bayesopt.py --cold-start    # 冷启动 EI; --prior-han 同源先验版
$PY scripts/crrna_virtual_cell.py --prior-lit # 文献先验场景表
$PY scripts/crrna_ivt_template.py --anova data/ivt_round1.csv   # 矩阵交互项判定(数据到位后)
```

## 安装

两个环境（清单见仓库根目录）:

- 主环境（跑主管线）: `pip install -r requirements.txt`
- 引擎环境（逆折叠/SHAPE 子进程）: `pip install -r requirements-engine.txt`
  （内含 torch_scatter/torch_cluster 等需匹配 torch 版本的包，安装备注见文件头注释）

更多环境细节与排坑见 `README_移交说明.md`。

## 大文件（模型权重/基准数据）

模型权重与 OpenKnot 基准数据不进 git 历史，存放在本仓库 **Release `data-v1`** 的资产中，
清单（路径/大小/sha256）见 `BIG_FILES.md`；下载后按相同相对路径放回即可。
双源说明：Gitee Release `data-v1` 已附 7 件单附件（RibonanzaNet×4 / gRNAde / OpenKnot M2·M2R）；
其余 3 件超 Gitee 单附件 100MB 上限（Struct2SeQ.pt / Struct2SeQ_SHAPE.pt / OpenKnotBench CSV），
已按 95MB 分片（.part00/.part01）同附于 [Gitee Release data-v1](https://gitee.com/eastern-zhongyuan/guideforge/releases/data-v1)——
下载全部分片后执行 `cat 文件名.part* > 文件名` 拼回，sha256 校验值见 `BIG_FILES.md`（拼回一致性已验证）；
也可从 [GitHub Release data-v1](https://github.com/dong-zhongyuan/guideforge/releases/tag/data-v1) 直接下载全量单件。

## 目录

- `configs/scaffold_registry.json` — 骨架注册表(cas12a2 PDB 结构口径 / cas12a2_zeng2026 细胞实验口径, 均 gate=confident)
- `scripts/` `src/` `tests/` — 管线代码与契约测试; 主管线 `crrna_scaffold_design.py`,
  文献验证 `crrna_dmytrenko_validate.py`, 上下文分型 `crrna_context_typing.py`,
  设计智能体 `crrna_design_agent.py`(--spacer 直输选骨架 / 突变转录本 tilling 设计)
- `data/` — 8D4A 结构、DR-蛋白接触表、TP53 WT/R248Q 转录本、正式候选库产出
  (zengdr.v2 = 加工保护+3'窗+共变口径; context_typing/ = 分型试验; agent/ = 智能体演示)
- `toolbox/` — 结构工具链: rnet-inference (RNet/RibonanzaNet) · Struct2SeQ · geometric-rna-design (gRNAde) · OpenKnotAIDesignData(竞赛基准数据)

## 诚实声明与边界限定

候选库是"候选压缩"结果, 不是活性预测; 打分为透明启发式, 未经实验标定。
最终活性/特异性以体外生化与细胞实验为准。两条必须主动声明的口径边界:

- **假结**: Cas12a2 crRNA 5' 端假结为蛋白诱导的结合态结构(Dmytrenko 2023; Methods in Enzymology 2025),
  游离态预测中不出现——本管线双路结构预测(ViennaRNA × RNet-SS)均为游离态口径,
  ViennaRNA 本身不折叠假结; 注册表条目的 bound_state_note 已如实标注。
- **游离态 ≠ 结合态**: 全部结构指标(含种子区游离 d_seed)为溶液态预测,
  结合态下蛋白会重排 spacer(Bravo 2023 四元复合物结构即证据); 正向项的含义
  是"预测降低溶液态自缠结、加速种子区暴露", RNet SHAPE 一致性为桥接证据,
  不是结合态活性保证。

## 来源

- Cas12a2 结构: PDB 8D49/8D4A/8D4B (Bravo et al. 2023, Nature 613:582)
- Cas12a2 功能/DR 互换/加工位点: Dmytrenko et al. 2023, Nature 613:588-594
- 细胞实验 DR 与 R248Q spacer: Zeng et al. 2026, Nature (10.1038/s41586-026-10738-7) 补充材料
- spacer-DR 互扰(同家族直接证据): **Creutzburg et al. 2020, NAR 48(6):3228-3243** (PMC7102956)
- spacer-骨架误折叠(SpCas9 体系, 跨体系佐证): Bush et al. 2023, Cell Chem Biol 30:879-892 (PMC10529641)
- DR 下游结构敏感性: Liao et al. 2018 (PMC6546362)
- Cas12a DR 序列工程先例: Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434); DR 3' 端敏感性(化学修饰): Zhang et al. 2025 (PMC11780881)
- Cas12a DR 突变提活先例: Teng et al. 2019, Genome Biology (4n96 骨架, PMID 30717767)
- crRNA 预折叠驱动 RNP 组装: Sudhakar et al. 2023 (PMC10200996; 因果方向: RNA 预折叠→诱导蛋白构象变化)
- 结构工具链: Science 2026 (10.1126/science.aeg6829, Das lab, OpenKnot); gRNAde 权重: HuggingFace chaitjo/gRNAde
