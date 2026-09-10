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
  86 条文献 spacer 按全长折叠特征聚 3 型, 各型最优 DR 位点图谱不同(弱证据阳性, 见文献依据表披露)
  → 智能体: 输入 spacer(项目口径)或突变转录本(换靶标) → 骨架分型选择/tilling 设计
    (含阴性退化模式: go/no-go 判据未通过时退化为通用型单骨架推荐, 见 V3 计算侧条)
  * 加工位点保护与 3'保守窗惩罚的文献依据见下"Dmytrenko 2023"条
```

### 状态竞争模型与交叉特征去混杂(2026-09-01 联审方案落地)

- **活性态-侵占态竞争模型**(`scripts/crrna_state_competition.py`, ViennaRNA 单能量标尺约束配分, hc_add_bp/up): 预注册判别检验 **FAIL**——14 条中 Sp8 的 dG_comp 仅排 6/14; 互斥定义(强制侵占螺旋+禁死茎臂)下所有可枚举侵占态均比茎态贵 5.8+ kcal(含 Sp8), 论文失活折叠在伪结外能量面上不可达。**结论: Sp8 正式定性 `out_of_model_domain`**, 活性态竞争假说的检验需统一假结能量模型(中期, RNAstructure/统一 pk 模型)。
- **交叉特征阶段拆分+组成校正**(`scripts/crrna_cross_decompose.py`, 单碱基组成匹配 N=200 + LOO + GC 偏相关): raw cross_coreDR rho=+0.55(LOO 稳定) 在组成匹配后降至 +0.17~+0.21、GC 偏相关后 +0.08~+0.18; flank 信号 z 后归零(纯成分代理, 且配对对象是被加工丢弃的 U-rich 残留 repeat)。**判据 R1: raw cross_pp 退出主评分获数据支持**; 阶段拆分保留(flank 只进加工模块口径)。
- **v1.7 结构置换口径过滤**: `inv_max_run`(最长连续 DR-spacer 侵占螺旋) <= 同上下文 WT(t1 型 WT 自身 6 对连续侵占故禁用固定阈值) + `partner_switch` 标注(t1 诊断: A8G 类 DDR-alone 稳定化突变在型1 口径强化的恰是侵占螺旋) + `ddG_dr` 语义更名 "DR 单独折叠稳定化(DDR-alone 口径)" + top.json 增加 `model_domain` 声明(伪结外: 仅侵占筛查+同 spacer 相对排序)。
- **补偿突变实验面板**(`scripts/crrna_compensatory_design.py` -> `data/ivt_compensatory_panel.*`): MYCg1 三臂 C_WT(inv_run=6)/A_break(0, GC 不变, 靶向配对保持)/B_break_compensate(6 恢复), 各配同源靶——A vs B 活性差 = 折叠竞争的因果检验, 寡核苷酸级成本。**round-3 去混杂重设计(design_mode=deconfounded)**: 旧 B 臂补偿位点落在天然茎 3' 臂(毁掉茎, DR-only p_fold 0.0), A vs B 无法归因单一变量; 新 B 臂(zengDR T10A/T12G/T19G)补偿位点全部避开天然茎, DR-only p_fold 0.93≥WT 0.87、bp_dist 0, A vs B 可归因于折叠竞争(计算口径, 待 IVT 验证)。

### NUPACK 3.2.2 含假结配分检验(第四模型, 2026-09-02)

`scripts/crrna_nupack_pk_weight.py`(NUPACK 3.2.2 官方 pfunc/-pseudo, Dirks-Pierce 假结配分; 源码取自 sbi-rostock/nupack-serve vendor 树, 容器 -fcommon 编译; 手册确认 pairs/mfe 不支持 -pseudo)。指标 W_pk=(Z_pk-Z_nopk)/Z_pk(假结构象权重占比)。**Sp8 判别 FAIL(第四模型)**: W_pk 排 5/14, 且与活性 Spearman=-0.327(方向相反)。至此四个独立模型(ViennaRNA 约束配分/ThreshKnot bpp/Knotty DP09-MFE/NUPACK Dirks-Pierce 配分)全部不能复现论文 Sp8 失活态——**out_of_model_domain 为四模型一致的最终口径**; 该机制需蛋白环境或实验手段。平台口径注记: GuideForge WT/茎稳定化 TOP3 的 W_pk 仅 0.9%~1.0%(假结权重可忽略), 与 8D4A/Knotty "Cas12a2 成熟 crRNA 无假结"三方互证, 本体系不受该边界影响。

### DP09 假结能量模型竞争检验(Knotty, 2026-09-01)

`scripts/crrna_pk_competition.py`(容器侧 `crrna_pk_runner.py` + Knotty 官方源码编译, CCJ+DP09): 14 条 Creutzburg 前体/成熟两口径 + GuideForge pdbdr WT/TOP3 的含假结 MFE。**判别检验 FAIL(三模型独立一致)**: DP09 能量面上 Sp8 折为"活性样"(干净 6 对茎[寄存器与 ViennaRNA 滑移 2 位——跨模型模板精确匹配无效的实测教训, 分类已改寄存器无关]+spacer 自发夹), 而高活性 guide(Sp2/Sp4/Sp7/Sp9/Sp12)反呈 DR-spacer 交叉+多重假结(方向与论文机制相反)。至此 ViennaRNA 约束配分 / bpp 阈值(ThreshKnot) / DP09 假结 MFE 三条独立路线均不能复现论文的 Sp8 失活态——**Creutzburg 侧 Sp8 out_of_model_domain 判定升级为三模型证据的最终口径**; 仅剩路径 = NUPACK 3.2.2 pfunc -pseudo 约束配分(需注册, 站点可达)。本体系(Cas12a2)不受影响: Knotty 对 GuideForge WT 预测与 ViennaRNA/8D4A 一致(无假结)。

### 结合态活性模板提取(8D4A, 2026-09-01)

`scripts/crrna_active_template.py` -> `data/cas12a2_active_template.json`: 8D4A 链 B(成熟 crRNA, DR18+spacer23) 三维几何判据(N1-N3<4.0A + C1-C1<12.5A)提取——**结合态活性模板无假结**: DR 区仅 5 对标准茎(1基 4-17/5-16/6-15/7-14/8-13) + 3 个近距摆动接触, DR-spacer 交界零配对。含义: Cas12a2 成熟 crRNA 的活性态在 ViennaRNA 伪结外空间内完整可表示(茎即为活性模板), **假结模型缺口只存在于 Creutzburg/FnCas12a 前体回验侧**(其失活态涉及前体 5prime repeat 侧翼, 论文自述活性态含 canonical pseudoknot)。两体系的模型域边界因此不同, 引用时须分开陈述。

### Chai-1 共折叠: 模板注入打通蛋白-RNA 界面(蛋白预测层正式口径, 2026-09-02)

`scripts/crrna_chai_input.py` + `scripts/crrna_chai_collect.py` + 容器 runner(~/dzy/envs/chai, 权重经 hf-mirror): SuCas12a2(1207aa)+crRNA(DR+固定 spacer)+靶RNA 三元共折叠。能力递进三步:

1. **单序列/ESM 基线(历史首跑)**: 无模板下 aggregate 0.23-0.24、蛋白-crRNA 链间 ipTM 0.02-0.08——1200aa RNA 引导核酸酶的已知短板, 界面不可用。
2. **8D4A 自模板注入(提交 e69fc74)**: 蛋白链自模板 m8(8D4A 链A, 100% 同一)注入后 aggregate/ipTM 0.25→0.874, 蛋白-crRNA 链间 ipTM 0.02→0.49。**round-3 降级声明: 100% 一致自模板注入的高分是模板复述的必然结果, 只能作为模板合规性 sanity check, 不能作为界面可预测性证据; 且 prot-crRNA ipTM≈0.49 低于 0.5-0.6 可用界面预测区间; template-free / scrambled-template 对照补齐前, 所有"界面恢复"表述均按此口径理解**(声明同步写入各 chai_* JSON)。
3. **AF3 独立对照层(2026-09-04 初版; 2026-09-05 格式勘误重发, 策划案 V3 表1 口径)**: 本机无 AF3 权重(需 DeepMind 审批)与 GPU 机时, 已生成 AlphaFold Server 就绪任务 13 套(WT + TOP-12, `data/af3_inputs/`, 生成器 `scripts/crrna_af3_input.py`, 蛋白/靶 RNA 与 Chai 矩阵同源; 逐任务文件 + 合并单文件 `GF_all13_alphafoldserver.json`)。**勘误两处**: ①旧版误用 AF3 开源代码库方言(dialect=alphafold3, 带链 id), 服务器只认 `alphafoldserver` 方言(顶层任务列表, 实体键 proteinChain/rnaSequence, 官方 server/README.md 口径), 已重发; ②旧版"服务器不接受自定义模板→恰好 template-free"的说法不成立——服务器**默认自动启用 PDB 模板**(8D4A/8D49 在库), 已在蛋白链强制 `useStructureTemplate=false`, 才是真正的 template-free 独立对照。提交操作与判读口径(预登记二分支: 界面恢复→独立佐证 / 界面低迷→主张撤回, 均完整可交付)见 `docs/af3_server_runbook.md`; 结果 zip 放入 `data/af3_results/` 后用 `scripts/crrna_af3_collect.py` 回填汇总(mock 包解析已测)。若 AF3 template-free 下 prot-crRNA 界面同样恢复, "界面可预测"才具备首个独立证据。**AF3 权重不可得的开源替代通道(2026-09-05)**: Protenix-v1(`protenix_base_default_v1.0.0`, ByteDance, Apache-2.0 代码+权重全开放, 368M 参数与 AF3 同规模、训练截止 2021-09-30 与 AF3 对齐, 基准在同等约束下达到/超过 AF3, bioRxiv 2026.02.05.703733)——13 套同一复合物输入已按 Protenix 方言转换(`data/protenix_inputs/`, 生成器 `scripts/crrna_protenix_input.py`, 逐链序列与 AF3 输入 13/13 一致; 不含 templates 字段即 template-free, 与 AF3 口径一致), 云端运行手册 `tools/protenix_cloud_runbook.md`(4090 级 GPU <2 小时跑完), 结果回填 `data/protenix_results/` 后与 Chai/AF3 做跨引擎一致性对照。**Protenix 13 折结果已回收并出对照(2026-09-09)**: 容器 A6000 批跑(2026-09-06, seed_101 × 5 sample)产物 13/13 完整回收, `scripts/crrna_af3_collect.py --protenix` 汇总为 `data/protenix_summary.json`——**预登记分支一触发: template-free 独立引擎下 WT prot-crRNA ipTM 0.835(5 sample sd 0.001), 全部 13 套 0.802–0.832 均 ≥0.6 可用阈, 界面恢复成立**; Chai 自模板 WT 参照 0.489, 两引擎界面恢复方向一致(数值差 0.35 属引擎间离散)。**"蛋白-crRNA 界面可预测"主张自此具备独立引擎佐证**(此前仅有 Chai 自模板单引擎证据); 变体对 WT 差量 -0.003~-0.033 且 WT 居首, 与 Chai 差量表方向一致, 但按预登记口径本层不用于骨架间排序主张。AF3 云端通道(alphafoldserver)**中止(2026-09-09, 用户决定)**: 提交需 Google 登录+逐任务手动上传, 自动化通道全部受阻(computer-use 会话锁/新版 Edge 禁止默认配置远程调试); 已提交的 GF_WT/GF_A1C 两任务在服务器队列自行运行, 结果 zip 若日后回收可作官方 AF3 单点抽查(WT 单点即可复核分支一), 但**不作为本层结论的必要件**——独立引擎佐证已由 Protenix template-free 分支一成立, 预登记判据未绑定官方 AF3, 中止不削弱现有主张。
3. **7 骨架 × 5 模型界面差量表(提交 963745d, `data/chai_cofold_matrix_*.json` / `chai_matrix_iptm.json`)**: WT/6 变体的 aggregate 与 prot-crRNA ipTM 模型间一致(sd≤0.007/0.09), 变体间排序稳定; **口径张力如实声明: crRNA-靶RNA 链对分数低(0.02-0.43)且部分变体模型间方差大, RNA-RNA 双链预测的模型稳定性存疑**——界面差量以 prot-crNA 链对为主口径, crRNA-靶对仅作参考(全部为 100% 自模板口径的模板合规性检查, 非界面可预测性证据——见上条降级声明); 与 MD 判据(M1-M3 预注册)的结论各自独立陈述。
4. **四靶 × 8 骨架 32 组合矩阵(已完成 n=32, `data/chai_matrix_4t.json`)**: 靶 RNA 统一为 protospacer 窗口 24nt + PFS 5nt(8D4A 排布同构)。reading/interpretation 由 `scripts/crrna_chai_matrix_summary.py` 从 scaffold_deltas 数值程序化生成(round-3 修复: 旧版手写 reading 与自身数值方向相反); 结论口径为"与骨架-靶标互作一致的单构建观测(自模板口径, 待对照校准)", 不作分型假说的结构证据。

### Han 2025 数据集与选型分类器(同源外部数据, 2026-09-02)

`scripts/crrna_han2025_features.py` + `scripts/crrna_train_selector.py`(数据: `data/han2025_dataset.json`, `data/han2025_sanger_sequences.json`)。

**数据链(除 Fig.3b RBS0/RBS33 两终点为源数据手动提取——代码内 TOOLBOX_ACTIVITY 如实标注——外, 其余全部程序化提取)**:
- Fig1g 活性全表 **145 值**(F×25/S×23/L×48/FL×48 + Canonical; 修正了旧版漏 48 条 FL 的解析 bug);
- 工具箱 9 条序列(CN/F1/F2/L1-L4/FL1/FL2, MOESM3)与 Fig1g 同尺度锚点 7 条;
- **54 条 Sanger 耐受集**(MOESM1 Sup Fig.2/3, 600dpi 视觉转录 ×2-3 次独立读数共识 + 模板守恒校验 + 与 MOESM3/主文锚点交叉验证, 逐条置信度分级; S8 的 20 条与论文自报数完全吻合);
- cis/trans 切割终点比(Sup Fig.5/6 源数据, 三重复均值)。
- **Fig.1f 编号-序列对照表视觉转录(2026-09-03)**: 主文 Fig.1f 以图像形式给出部分映射, 分块放大两次独立逐字母转录 + Sanger 独立源 loop 六联体交叉校验, 恢复 **9 条新监督对**(S3/S5/S7/S10/S14/S16/S20/S21/FL4, 置信度 high 2 / medium 7, 溯源见 `fig1f_pairs`);
- **查证边界(2026-09-03 更新)**: 145 个编号变体与序列的映射**未随任何数值源数据发表**(MOESM7 全部 25 工作表 + MOESM3/4 已穷尽复查; 库为随机合成, Methods 明示), 公开渠道上限=7 工具箱 + 9 Fig1f 转录 = **16 对**; 完整映射待作者回复(docs/ 邮件草稿)。
- **Tian 2025 RRS 单点扫描入库(2026-09-03, 同源第二个数据集)**: [Tian et al., Nat Commun 16:6694](https://www.nature.com/articles/s41467-025-62082-5) 系统扫描 LbCas12a DR 的 RRS 区(5' 端 4nt)——8 条 crRNA × 12 单点突变 = **96 条带序列+活性配对**(序列=MOESM3, trans-切割 ΔF=MOESM6 Fig.1d 三重复均值, 活性口径=同 crRNA 内 mut/WT 比值), 程序化提取零命名冲突(`parse_tian2025()`); 论文结论: RRS 3/4 位突变近乎灭活(假结承载)。**该集按设计用作选型器的独立外部验证而非训练数据**: 池化模型(训练集不含 Tian)预测 96 条, 整体 Spearman **−0.007**(2026-09-04 tie-aware 全量重算 **−0.040**, 逐 crRNA 中位 −0.092; 两套独立新环境 numpy 2.3.5 与 2.5.2 重算一致——漂移源于统计口径/特征表重生成, 非版本噪声, 两口径并列报告, 详见下行 DeWeirdt 条)、逐 crRNA 中位 −0.094——证实伪结外特征对 RRS/假结区按构造无分辨力(与 Creutzburg Sp8 盲区同源); **选型器排序不适用 DR 5' 端 RRS 位点变体, 候选库中该区突变应按位置规则排除**。入库字段: `han2025_dataset.json` 的 `tian2025_rrs_pairs` / 验证结论 `homolog_training.json` 的 `tian2025_external_check`。
- **DeWeirdt 2020 alt-DR 大库入库(2026-09-03, 同源第三个数据集)**: [DeWeirdt et al., Nat Biotechnol 39:94](https://www.nature.com/articles/s41587-020-0600-6) 的 AsCas12a 替代 DR 深度扫描——**35,883 条 20nt DR 变体**(茎区 ≤3 可变碱基对 + 单链/环区 ≤3 可变核苷酸)× 2 方向构建(BCL2L1×MCL1 合成致死负筛, MELJUSO/2xNLS-AsCas12a, LFC 低=活性强; 序列+读数=MOESM6 alt_DR_reads, spacer 上下文=MOESM1 Supp Table 3 的 MCL1/BCL2L1 guide), 程序化提取 + 逐方向特征(`parse_deweirdt2020()`, 全量特征表 `data/raw/deweirdt2020_dr_scan.json` 不入 git)。**入池实测被拒**: 120 条 LFC 秩分层子样(占池 72%)把其余 5 终点 LOEO 全部拖负(rbs33 +0.96→−0.75、rbs0 +0.68→−0.64、trans +0.25→−0.71、fig1g +0.04→−0.70, 自身 −0.36, 仅 cis 转正)——单终点淹没多终点共识, 与 Tian 同模式改作独立外部验证(拒绝记录见 `homolog_training.json` rejected_endpoints)。**外部验证(训练集不含 DeWeirdt, 全量无剔除)**: pRDA_127(MCL1 spacer 上下文) Spearman **+0.202**, 且不低于同数据集 5 折 CV(+0.182); pRDA_128(BCL2L1 上下文) −0.024(CV +0.109)。**2026-09-04 重算并列报告(特征表 data/raw 重新生成 + tie-aware Spearman 统一口径; `homolog_training.json` 已刷新为新值, 旧值为其保留注记所述"旧 tie-blind 口径"产物)**: pRDA_127 **+0.016**(CV +0.181)、pRDA_128 **−0.102**(CV +0.111)——同数据集 CV 复现一致, 迁移端绝对值漂移; 两套独立新环境(numpy 2.3.5/scipy 1.16.3 与 numpy 2.5.2/scipy 1.18.1)重算数值完全相同, 排除随机/版本噪声。方向与定性结论不变(pRDA_127 弱正/近零、pRDA_128 近零/负、弱信号不足以支撑逐点排序); 但旧值支持的子判读"外部验证不低于同集 CV、瓶颈在特征/终点而非迁移"(+0.202≥+0.182)在新值下**不再成立**(+0.016≪+0.181), 该比较不再作为论据, 迁移损失是否存在留待更多同源终点数据判定。判读: 该筛选动态范围由茎/环随机化细节主导, 5 维伪结外特征仅携带弱信号(秩方差 ~3%), 定性为**方向正确的弱正证据**, 不足以支撑逐点排序。Nguyen 2020(Nat Commun 11:4906)已评估不入库: 其变体为 crRNA 3'/5' 端尾延伸(DNA/RNA/PS), DR 序列不变, 非 DR 序列-活性配对。

**选型分类器(先验排序口径, 2026-09-03 用 16 对重训; 2026-09-04 round-3 全面切 tie-aware Spearman)**: 决策树(Fig1g 尺度, n=16), 特征重要性 p_fold 0.52 / ddG_dr 0.48; 训练集 Spearman 0.501(样本内读数, 仅作参考); **fig1f 留出验证(训练 7 工具箱 → 检验 9 转录对): Spearman +0.388 [bootstrap CI95 0.00, +0.864], MAE 0.217 [CI95 0.161, 0.289](超训练折可用阈 0.0467; 旧"预登记"阈 0.162 用含测试标签的全量 y 计算, 标签泄漏已于 round-3 修复并改名 mae_threshold_trainfold——两口径下 MAE 均超阈, 判定方向不变); MAE 超阈已结构化归因——9/9 条 fig1f 标签高于工具箱训练上界 0.583, 决策树按构造不可外推, 超阈主因为主图/附图两图版间水平偏移而非排序能力失败; 秩方向为正但 n=9 功效有限, 故 fig1f 配对只用于秩级体检、绝对水平不采用, 两口径并列报告**; **线性第二口径(2026-09-08)**: 可外推的岭回归(StandardScaler 仅训练折拟合, Ridge α=1.0)同条件重评 fig1f 留出——MAE=0.2164 与决策树 0.217 一致超阈, 「树不可外推」假设被证伪, 图版间水平偏移坐实为**模型无关的数据因素**; 秩一致性 +0.367 与树口径 +0.388 同向, 转录配对秩级可用的结论获第二模型类加强(产物块 `fig1f_holdout_check.linear_second_caliber`, 契约测试 tests/test_fig1f_linear_caliber.py); 对保守变体仍无区分力(多数候选落同一叶节点)——定位为文献先验排序, 非活性预测。池化模型扩至 n=46(Han 5 终点×7 + Fig1f×9 + Teng 序数对); LOEO 有升有降(rbs33 0.64→0.889, rbs0 0.39→0.667; trans 0.68→0.433, fig1g 0.29→0.216, cis -0.68→-0.692 反预测依旧), 8 员族池化排序**翻转**(旧 U5G+A18C +0.58 居首 → 新 WT/A1C/A1U 0.12 居首, U5G+A18C -1.38 末位; 新数据的 loop 同聚物变体均为中等 ddG_dr/p_fold 而活性不差, 强茎稳定化奖赏被压低)——排序对数据增量敏感, 印证"先验排序"而非稳定预测的定位。多终点版: cis/trans 切割终点各训一树, 8 员族逐终点排序入库(`data/selector_family_ranking.json`); **trans 旁切终点上同源骨架间差 30 倍(WT 39.2 vs L1 1.2)**, 为"DR 可调旁切活性"提供了同家族文献证据; 不同终点给出不同排序, 正是"骨架分型"的同源佐证。Sanger 耐受集体检: 真实耐受变体预测中位 0.609(锚点范围内), 204 候选分布中平均 6.5% 分位(模型不排斥真实功能变体)。**排序权威声明(round-3 R3): 合成候选排序以管线打分为准(`data/selector_model.json` 的 `ranking_authority: "pipeline_score"`); 此前 JSON 中"管线 vs 选型器 ρ=−0.872"为 tie-blind 统计伪影(决策树仅 3 个 distinct 预测值, 统计量随 numpy 版本/输入顺序摆动 +0.587/−0.493/−0.872), tie-aware 重算为 +0.157(p=0.025; 任务⑪保守窗 7nt 重算后, 旧口径 5nt 时为 +0.152), 顶部 TOP-12 重合 10/12, 残余弱分歧源于选型器偏好 DR 单独强稳定化与 LbCas12a→Cas12a2 终点迁移差距**。

- **同源跨上下文 DR 效应交互检验(2026-09-08, 判据登记 docs/preregistration.md §E,
  登记先于评估; `scripts/crrna_homolog_context_interaction.py` → `data/homolog_context_interaction.json`,
  契约测试 `tests/test_homolog_interaction.py`, 独立复算 `scripts/crrna_homolog_interaction_audit.py`)**: 
  分型假说(骨架/DR 效应的最优解随靶标上下文变化)的**原理层**统计判定, 用同源大库替代
  未执行的 IVT 8×4 矩阵——只回答"跨上下文 DR 效应交互是否为该家族普遍现象",
  不替代 Cas12a2 本体系验证(功能收益维持"待湿实验")。**E1 主判据(DeWeirdt 2020,
  AsCas12a alt-DR 负筛): 同一批 35,682 条 test 变体在 pRDA_127(MCL1)/pRDA_128(BCL2L1)
  两个 spacer 上下文下的跨方向 Spearman ρ=0.257 [bootstrap CI95 0.247, 0.267], 
  按 §E 判据(ρ<0.4 且 CI 上界<0.5)判 **PASS**——同批 DR 变体效应仅 ~6.6% 共享方差,
  强烈依赖 spacer 上下文, 分型假说获同源跨上下文大样本统计支持**; n_mut 分层
  ρ 0.33/0.25/0.22(1-3/4-6/7+), 随突变负担增加一致性递减, 方向一致。**E2 佐证
  (区位分层)**: 茎区变体 ρ=0.289(n=4,389) vs 环/单链区 ρ=0.249(n=24,462),
  Fisher z p=0.0096(差 0.04, 大样本小差, 佐证级); 与 §A 本体系计算弱证据阳性
  (公共集 3 条 ≤ 4, 置换 P=0.036)并列, 形成「同源原理统计 + 本体系计算实例」
  双层证据。**E3 佐证(Tian 2025, 12 RRS 单点突变 × 8 crRNA)细化为双层(如实,
  不掩盖)**: E3a 全量 Kendall W=0.769 高, 由**位置梯度**驱动——RRS 3/4 位
  (假结核心)近乎全灭活、位置 1 最耐受, "哪一位敏感"跨 crRNA 稳定(通用层;
  与 Cas12a2 本体系 §A 公共集 A1C/A1G/A1U 同为位置 1 惰性, 跨体系呼应);
  **E3b 位置 1(耐受位)碱基偏好 Kendall W=0.016 < 0.3 → 交互佐证(E1 同向)**——
  同一耐受位点上"哪个碱基最优"随 spacer 上下文翻转(最优碱基分布 1C×4/1G×2/1U×2,
  两两比较 9/20 与多数相反)。位置级敏感图谱通用与位点内碱基最优上下文特异
  是两个层次, 不矛盾, 后者正是分型假说的同源实例; DeWeirdt 仅 2 个上下文
  (是"效应是否普遍一致"的最大样本检验, 非"型"的统计)、负筛 LFC 非动力学、
  AsCas12a 与 Cas12a2 结构域拓扑差异(假结承载 vs 成熟无假结), 迁移引用按本
  README「Cas12a→Cas12a2 同家族短迁移」边界声明执行。
- **IVT 8×4 矩阵统计件(`scripts/crrna_ivt_template.py --anova`)** 保持 data-pending 占位:
  本体系"骨架×靶标交互"的实验统计判定, 回补执行后与 §E 同源层构成完整闭环。

### V3 计算侧(第 8 员骨架 / 四靶扫描 / 统计件 / 智能体, 2026-09-02)

- **8 员跨型候选族**(`scripts/crrna_orientation_library.py`): 四取向代表 6 + WT + compensatory 代表 B_break_compensate(zengDR+T10A/T12G/T19G, round-3 去混杂重设计, 补偿位点避开天然茎; `data/orientation_library.json` 已于 round-3 重跑同步); 构成口径披露——补偿臂为 round-3 评审要求的因果检验臂(非优化取向产物), 计入策划案 V3 "约 8 个跨型候选(含 WT 对照)"的"约 8"口径, 特此说明。**IVT 模板 32 行(8 骨架 × 4 靶标, `data/ivt_round1_template.csv`)2026-09-04 起按策划案 V3 表2 对齐: R248Q/G12D/R273H/APC-Q1328x**(2026-09-08 起经阅读框核查自 Q1312x 修订, CSV 已同步重生成; KRAS-G12C 移出湿实验矩阵, 其 v47 扫描等干实验产物保留作附加证据); 模板与订单表(`data/ivt_round1_order_sheet.csv`)已同步重生成——**同步传播了去混杂新补偿臂**(旧版两文件仍编码 round-3 前退役臂 AATTTCTACTCTTCTACAT, 按旧单合成将买到已撤回分子; 旧订单表备份于 `tmp/ivt_round1_order_sheet.retired_arm.bak.csv`); 面板-库-订单三向一致性由 `tests/test_ivt_panel_sync.py` 契约测试锁定。
- **五靶转录组脱靶扫描**(GENCODE v47, 385,659 转录本; `data/*_scan_v47.*`; 策划案 V3 表2 四靶 R248Q/G12D/R273H/APC-Q1328x + KRAS-G12C 干实验附加): 五条 spacer 的 0 错配位点全部落在本基因异构体(预期靶点); KRAS 两 spacer 在 **KRASP1 假基因各 1 个 1 错配位点**; **APC 无任何 ≤2 错配位点(最干净)**, 仅 6 个 4mm 位点在 HNRNPR; 无高危脱靶。参考库 gencode.v47.transcripts.fa 不入 git(`data/raw/`; 源: EBI GENCODE release_47 FTP, 2026-09-04 重新下载, 385,659 条与旧扫描计数一致)。**口径注记(round-3 R5; v47 表 2026-09-04 已用修复版扫描器重跑为 tiered 口径, APC 同日补扫)**: "R248Q PFS=CAGAG" 与共识 GAAAG 实为 2 个错配(第 1、3 位), 早期文档"1 错配"表述有误; PFS 规则已统一入注册表 `pfs` 字段(scanner/agent/面板均消费同源); **v47 五靶重扫(effector=cas12a2_zeng2026, GAAAG±2)位点计数与旧表逐一一致(四靶), 每站点现含 pfs_mm/pfs_tolerant 分级与双模型计数——全转录组 PFS 匹配脱靶五靶均 exact 0mm=0 / tolerant ≤2mm=0**(参考库为 WT 等位, 突变体转录本不在其中; R248Q 突变体场景由 2 转录本重扫覆盖: exact 0mm=0, tolerant ≤2mm=1, 容忍模型恰好覆盖 CAGAG, `data/tp53_r248q_scan.summary.json`)。注意: 参考库中本基因 0 错配位点的 WT PFS 语境距共识 3~5 错配(CGGAG/GGTGG/GTGGC/CGTGT/CAGAT), KRAS G12C/D、R273H 与 APC Q1328x 突变体的 PFS 语境未入任何扫描——属设计层遗留问题, 不在脱靶扫描口径内。
- **贝叶斯优化同源先验版**(`scripts/crrna_bayesopt.py --prior-han`): Han 工具箱 7 条作 GP 观测(选择器特征空间), 对 204 候选提 EI 建议; 明确标注"LbCas12a CRISPRi 抑制终点, 非 Cas12a2 杀伤; 湿实验=细胞-only 后为正式先验, 细胞实测经 --ingest-cell 回流"。
- **虚拟细胞文献标定场景版**(`scripts/crrna_virtual_cell.py --prior-lit`): Scholz 2026 Fig 1h 剂量曲线标定(n=15 点, R²=0.867)的 EC50(TPM) 网格 + Cas12a2 体外激活浓度量级场景化, 输出杀伤窗口表与主验证细胞系激活率(HCT116/SW480 等)。**round-3 修复与如实标注**: crrna_rel 场景旋钮已接入计算(e50/κ 线性近似, 未标定); 输出 JSON 含 `qa` 字段——15 点为 7 靶标×4 细胞系混合池化(未分层, 仅作包络参考)、含 1 个 >100% 存活点(如实保留)、RPKM≈FPKM 跨源等效与统一 50% 杂合因子为近似假设。
- **矩阵统计件(数据一到即出)**: `scripts/crrna_ivt_template.py --anova`(骨架×靶标两因素方差分析含交互项——"分型假说是否成立"的统计判定, 合成数据自测通过; **ivt_round1.csv 数据待补, 当前为 data-pending 占位**)与 `--twin-check`(数字孪生预测 vs 细胞实测 Spearman/RMSE; 参数未标定时如实报 UN-CALIBRATED)。
- **智能体模块三接入选型器**(`scripts/crrna_agent_webapp.py`): `/api/design` 与 `/api/panel` 输出 8 员族文献先验活性排序。**口径修正(round-3)**: webapp 侧为轻量重训(n=7 工具箱锚点), 与 `crrna_train_selector.py` 的 n=16/n=46 模型**不同规模**——此前"同一模型同一定义"的表述不准确, webapp 排序仅作演示, 正式先验排序以 `data/selector_model.json` 为准; 特征表(`FEATURES`/`TARGET_FEATURES`)2026-09-05 起两端经 import 共用同一定义, 不再各自复制字面量。
- **靶RNA丰度档位接入智能体(策划案 V3 §4.4 模块一特征维, 2026-09-05)**: `scripts/crrna_target_abundance.py` 从 CCLE/DepMap 18q3 RPKM gct 实测提取 panel 基因表达(TP53 HCT116 15.95 / SW480 32.12; KRAS 5.34 / 5.33; APC 7.14 / 2.84 RPKM), 按显式规则换算档位——杂合场景 mut = total×0.5(与虚拟细胞层同口径), 阈值锚定 Scholz 2026 实测 EC50 95%CI [5.2, 17.4]: tier 0 低(<5.2, 低于激活阈)/1 中(阈附近)/2 高(≥17.4), 落盘 `data/target_abundance_tiers.json`(含来源/阈值/实测值/判读, 判读全部由数值按规则生成)。四靶实测档位: R248Q=1、G12D=0、R273H=1、APC-Q1328x=0。模块一: `crrna_design_agent.py --target-key` 把档位写入每条候选的 `target_features`(四个 panel 预计算设计已重生成, 设计内容不变、仅增量字段); 模块三消费方式如实声明: 文献决策树**不**在该维分裂(Han/Fig1f/Teng 训练行无靶标上下文), 档位经显式阈值门控 `activation_gate`(低于 CI 下界→激活预计不足的门控判读等)进入推荐输出与 Web 演示页展示, IVT 8×4 矩阵训练集 schema 将含此列(`selector_model.json` 的 `target_context_features`); 契约测试 `tests/test_target_abundance.py` 锁定 JSON↔gct 一致性、分档边界、agent 特征向量含该维、webapp↔离线 selector 特征列同源。
- **智能体阴性退化模式(策划案 V3 §5.3 末段 / 表5 风险1 对策落地, 2026-09-06)**: 「若矩阵数据显示骨架最优解与靶标特征无显著相关(阴性结果), 智能体退化为单骨架推荐, 科学结论依然完整可交付」。触发判据 = 预登记 go/no-go(docs/preregistration.md §A, 判据本体未改, 本次仅登记下游消费点): `crrna_design_agent.py` 每次运行从分型产物重算各型代表 TOP-8 骨架的全型公共子集, **> 阈 4(topk//2) 即阴性** → 输出从分型推荐退化为**通用型单骨架推荐**(`universal_scaffold()`: WT 0 基线 + 各型代表骨架中取管线分全场最优, 并列 WT 优先, 显式确定性规则); 退化时全部候选统一通用骨架、`scaffold_type=null`, 分型结果降级为 `typing_nearest_type` 诊断(不参与推荐——样本级「置信度 boundary 报两型」的既有行为保留, 与全局模式级退化是两个层级, 见模块 docstring)。输出 JSON 的 `degradation` 块如实携带模式、判据数值(公共集/阈/型间 Jaccard 均值)、事后置换检验证据(post-hoc, 不作触发)与「科学结论依然完整可交付」交付口径, 判读文字全部由数值按显式规则生成; `--degrade on/off` 强制演示/覆盖并以 `forced` 字段如实标注。Web 演示页启动时同源评估, 退化时页面以醒目横幅如实展示(不掩饰)。**当前真实数据判为 typed(弱证据阳性, 公共集 3 条 ≤ 4), 智能体维持分型推荐**; 四个 panel 预计算设计已重生成(设计内容不变、仅增量 `degradation` 块); 契约测试 `tests/test_degradation_mode.py` 锁定判据边界(公共集 4/5)、真实归档判型、触发/不触发两场景端到端行为与 webapp 载荷一致性。
- **突变类型判定接入模块一(最长 ORF 法, 策划案 V3 §4.4 / 表2 APC 入选理由对齐, 2026-09-07, 任务⑩)**: 新模块 `scripts/crrna_mutation_typing.py` 分别在 WT 与突变转录本中定位最长 ORF(正义链 3 框, ATG 起、首个框内终止子止; tie-break 长度→起点→框序; 末端部分密码子丢弃; 规则与边界全部显式声明于模块 docstring), 按突变对参考 ORF 的影响判定: **错义**(ORF 起止/长度不变、单氨基酸替换)/**无义**(ORF 内提前终止子、ORF 变短)/**移码**(|Δlen|%3≠0、阅读框位移), 边界类别(同义/框内 indel/终止子或起始密码子丢失/非编码)如实单列, 复杂多事件拒绝判定; 判读文字由数值按显式规则生成。四靶判定: **R248Q=错义(R248Q, ORF 393→393 aa)、G12D=错义(G12D, 189→189)、R273H=错义(R273H, 393→393)、APC=无义(Q1328*: WT ORF 60..8591=2843 aa → 提前终止于密码子 1328, MUT ORF 60..4043=1327 aa, 截短 53.3%)**。**阅读框核查与口径定稿(2026-09-08 拍板)**: 旧 `apc_q1312x` 序列对的 C>T 编辑(转录本 3946 位)是 panel 脚本以 `seq.find("ATG")` 锚定的——NM_000038.6 的 5'UTR 有上游 ATG(13-15 位), 真实 CDS=最长 ORF(60..8591), 旧编辑在真实阅读框是 A1296V 错义并非无义(真实阅读框 MCR 首个 CAG = 密码子 1328; 真实密码子 1312 为 GGA, "Q1312 无义"在真实阅读框不存在)。经用户拍板, APC 靶点正式整体切换为真实无义 **Q1328\***(c.3982C>T, 转录本 4041 位), panel 键名更名 `apc_q1328x`, 湿实验矩阵 CSV / 演示快照 / Chai 增量输入 / 选型特征全部重生成, 旧键名仅存在于 git 历史。溯源: 该修订源自任务⑩阅读框核查, 策划案旧口径曾短暂锁定, 最终采纳核查结论。接入: `crrna_design_agent.py` 模块一把 `mutation_type` 与紧凑 ORF 证据写入每条候选 `target_features`, 完整 ORF 证据(WT/MUT 最长 ORF 起止、MUT 同起点对应 ORF、提前终止位置、截短比例、判读)入输出 JSON 的 `mutation.typing` 块(四个 panel 预计算设计已重生成); 突变类型维登记进 `crrna_train_selector.py` 的 `TARGET_FEATURES` 契约(`selector_model.json` 已同步, 如实声明文献决策树在该维不能分裂——类目维且文献训练行无突变类型标注), webapp panel 页展示该维(与离线同源于 design.json 证据块, 只读不重构); 契约测试 `tests/test_mutation_typing.py` 锁定 ORF 规则边界、四靶 WT ORF 合理性(经典 CDS 长度)、四靶判定与数值证据、无义/移码等合成边界用例、design.json 字段与 webapp↔离线同源。

### 打分输入端可靠性(系综采样检验)

`scripts/crrna_ensemble_check.py`(输出: `data/ensemble_check.pdbdr.v2.json` / `ensemble_check.zengdr.v6.json`): 对 WT+TOP-12(两口径共 26 构建)做 Boltzmann 系综检验(ViennaRNA 2.7.2 pbacktrack 本构建不可用, 改用 subopt 精确枚举 5 kcal 窗口, 窗口尾部质量均 <2.2%; 能量均值/展宽/序翻转/茎完整率全部为枚举窗口 Boltzmann 加权精确值, 100 采样仅供退役统计量溯源, seed=42)——
- **R1 腿成立(24/24)**: MFE 口径 ΔΔG 与精确系综自由能差 F_var−F_WT(pf 精确值)全部一致, 偏差 ≤0.03 kcal/mol, "-2.2 kcal/mol 茎稳定化"等数字不是 MFE 单构象假象; 稳定化候选(−2.2/−2.0/−1.4)的 ΔΔG bootstrap 95% CI 均不含 0;
- **旧 R2/R3 判据退役 + 重设计 R2′/R3′ 全部通过(round-3, 先登记后评估)**: 旧 R2/R3 统计量按构造病态(独立抽样按采样序号配对, std≈√(sd_v²+sd_w²)≈1.3 kcal, 与观测最大偏差 0.09 kcal——构造性必然), 其普遍触发(19/24、24/24)不构成真实不稳证据, 已退役并留档 JSON `deprecated_ill_posed_R2R3`; 重设计判据(登记于 `docs/preregistration.md` §B, 登记先于评估运行): R2′ 序翻转——6 条实效候选(|ΔΔG_ens|≥1.0) P_contra ≤0.13(阈 0.2, 两枚举窗口卷积精确值, 近中性候选 P_contra≈0.5 属普通热涨落重叠, 按登记不判定、逐候选披露); R3′ 景观粗糙度——展宽比 ≤1.14(阈 2.0, WT 分裂半零分布超阈经验率 1.0%/0.5%, 在 ≤1% 采信线内); R1 0/24、R4 0/24(茎完整率改精确窗口值, 消除 100 采样在 0.8 判界附近的噪声)。**结论强度=`passed_all_rules`——"打分输入端在系综采样下稳定"主张在重设计判据下成立**;
- DR 茎完整构象占比 0.83~0.87(WT 0.87, 枚举窗口精确值), 无构象二态; 单构象能量涨落 sd≈1.0 kcal/mol(WT 窗口精确值 1.00/1.02)属系综固有涨落, 在自由能差中抵消。
- 边界: 此检验证明的是 ViennaRNA 最近邻模型内的数值可靠性(指标计算端), 不涉及指标与活性的相关性(见"路径A"条, 后者已证明不构成跨体系活性预测器)。

### 打分增量信息与排序权威(2026-09 round-3 R3)

- **filter+random 基线**(`scripts/crrna_selection_baseline.py` → `data/selection_baseline.cas12a2_zeng2026.json`): 对 zeng2026 口径重枚举同一突变库(441 变体/203 通过, 与 v6 逐条对拍一致), 管线 TOP-12 均分 −0.400 vs 过滤后随机基线 −1.714±0.218(200 次重复, MC p=0.005 为构造性上界), 效应量 1.70 个通过池 SD; 随机抽取与 TOP-12 重合≈超几何期望(0.72 vs 0.71)——"硬过滤+随机"捞不回同一批候选。独立性质上 TOP-12 显著避开真 loop 突变(0.0 vs 随机中位 0.215, p=0.015); 茎区 GC/选型器先验与随机无显著差异(如实报告)。
- **逐项消融(LOO)与 ddG 不对称**(`scripts/crrna_score_ablation.py` → `data/score_ablation.cas12a2_zeng2026.json`): w_contact 是对 TOP-12 构成影响最大的打分项(置零后重合仅 4/12, TOP-1 易主), w_ens 次之(7/12), 其余项重合 ≥10/12, 严格无效项不存在(w_seed 默认已为 0); ddG 轴 0.1 罚/0.3 赏的 V 形不对称对 TOP-12 构成不敏感(对称化两档重合 ≥10/12, TOP-1 不变), 是否保留随湿实验标定复核。
- **排序权威**: 合成候选排序以管线打分为准; 选型器(n=16 文献先验)仅为先验参考, 不作合成排序依据(分歧细节见上"选型分类器"段末尾声明)。

## 文献依据（对应模块）

| 模块 | 文献依据 |
|---|---|
| 加工位点保护硬过滤 / DR 互换不误杀验证 | Dmytrenko et al. 2023, Nature 613:588-594 (ED Fig.2c Cas12a↔Cas12a2 DR 互换功能保持; ED Fig.3 加工切点在 Cas12a 下游 1nt)。位点敏感性×保守性相关验证(round-3 R2 修复后重算: loop 分区改 flanked-by-paired——旧版把 5' 悬垂 1–4 位错标为 loop; Spearman 改 tie-aware midranks): score-adj 口径 ρ=−0.335 (置换 p=0.16), mean-score 口径 ρ=−0.392 (p=0.10), tolerance 口径 ρ=−0.426 (p=0.08), n=19 功效低、均不显著但三口径方向一致为负——引用时须注明修复口径与样本量 |
| 3'保守窗惩罚 | Dmytrenko et al. 2023, Fig.1c 跨家族 DR 3' 端保守/loop 可变(定性结论); **窗长自任务⑪(2026-09-08)起数据化: 同源 DR 库扩至 n=11(9 基因座, 逐条带来源, `data/dr_homologs.json` v2: SuCas12a2 8D4A + Methods Enzymol 2025 (PMC12975306) array 全长 repeat 复核 + Zeng 2026 + Zetsche 2015 直系同源 Fn/Lb/As/Pm/Mb/Ts/Bs(专利 SEQ ID 多源复核) + Creutzburg 2020 + Teng 2019), 锚定(TTTCTACT/GTAGA)+middle 星形比对(自实现 NW, 显式打分)得 18 列位点保守性谱(`data/dr_conservation.json`; 阈值显式: identity==1.00 完全保守 / >=0.80 保守 / 否则可变; 可变列集中在比对区中段 = Fig.1c loop 可变区, 3' 端锚区完全保守); --cons3-window 默认取先验 cons3_window_prior=7nt(strictly_conserved 3' run 6 列 + trailing3 1nt, trailing 碱基在全部含 trailing 条目均为 T), CLI 可覆盖, 先验缺失回退 5nt; 重算影响: v6 主库 442 变体 pass 集合与 TOP-12 构成不变, 仅含保守窗突变行的 score 微调(详见 dr_conservation.json 与任务⑪报告)**; 功能佐证: Zhang et al. 2025 (PMC11780881, DR 3' 端化学修饰可逆调控 Cas12a 活性——非序列突变, 佐证 3' 端敏感性) |
| spacer-DR 互扰 / spacer 游离度硬过滤 | **Creutzburg et al. 2020, NAR 48(6):3228-3243** (PMC7102956, 同家族直接证据): "base pairing of the direct repeat, other than with itself, was found to be detrimental"—DR 与 spacer 碱基配对损害 Cas12a 活性; Sp8 案例证实 spacer 侵占 DR 致假结破坏; 通过 spacer 尾部回折竞争性保护 DR 挽救活性(与本项目 spacer 游离度约束直接对应)。**Cas12a→Cas12a2 为同家族短迁移(DR 均含保守 5' 假结), 非跨体系假设** |
| spacer 上下文分型的机制前提 | 同上(Creutzburg 2020, 同家族); 补充: Bush et al. 2023, Cell Chem Biol 30:879-892 (SpCas9 体系, sgRNA 80nt 大骨架, 跨体系佐证); Liao et al. 2018 (PMC6546362, FnCas12a DR 下游发卡抑制切割); **综述级背书**: Annu Rev Chem Biomol Eng 2023 (annurev-chembioeng-100522-114706)——"gRNA 二级结构, 无论 spacer 内部还是 **spacer 与 guide 其余部分之间**, 都是 gRNA 活性的关键决定因素"(A–B 耦合前提的综述级陈述); **工程先例**: MIT 学位论文(Cas12a scaffold variants toolbox——"多样化骨架变体使每个 guide 可优化折叠为 Cas12a 可识别结构", 即按上下文选骨架的工程实践); EP4403638(DR 茎环保双链突变不影响切割、破坏双链则灭活——结构保持硬过滤的实验先例)。分型主张的直接证据为本仓库分型试验计算结果; Cas12a2 本体系暂无直接折叠互扰研究(待验) |
| DR/骨架序列工程先例(Cas12a, 存在性证明) | Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434, DR 突变策略调编辑与检测) |
| 共折叠界面差量(结构层选型特征) | 工具链领域地位: Chai-1 (PoseBusters 77%, 与 AF3 报告的 76% 持平), Boltz-1/Protenix 为同族 AF3 级复现; ipTM 阈值惯例(>0.6 中高置信界面几何, >0.8 高置信)适用于**绝对界面断言**; **本仓库仅作差量(deltas)排序使用**——「共折叠初筛 + 第二方法复核」是领域标准工作流(先例: arXiv:2511.14669 用 AF3 ipTM 筛选界面候选 + Rosetta 能量复核), 差量用途不以全称"界面可预测性"为前提; 自模板复述局限与 template-free/scrambled 对照计划见上文 Chai 条 round-3 降级声明, Protenix 跨工具回填后提供一致性检验 |
| guide 层修改跨体系迁移的边界 | **Cas9(II 类)与 Cas12a2(V 类)蛋白无同源性, 不存在也不应主张 Cas9→Cas12a2 序列级迁移**(骨架排布相反: Cas9 骨架在 spacer 3' 侧, Cas12a/Cas12a2 在 5' 侧, Kweon 2017 PMC5700056)。可用证据分三层: ①同家族 DR 互换直接证据 = Dmytrenko 2023 ED Fig.2c (Cas12a↔Cas12a2, 本项目迁移链的唯一同源层); ②同一 guide 策略跨两体系的唯一直接演示 = zCRISPR, Nat Commun 2024 (s41467-024-48012-x, 碱基 Z 修饰同时提升 Cas9 与 Cas12a——化学修饰层, 非骨架序列层); ③类比层 = Cas9 侧 sgRNA 骨架工程 (Chen 2013 eSG; Bush 2023) 与 Cas12a 侧 DR 工程 (Han 2025/Zhang 2025/Teng 2019/DeWeirdt 2020) 各自独立成立, 仅支持"骨架工程调活性为跨 II/V 类重复出现的原理", 不支持设计移植 |
| "按靶标分型选骨架" | **无先例(空白区)**: 本仓库分型试验提供了反对"单一通用骨架"的计算证据(型间最优位点图谱结构性不同); 功能收益待湿实验。**round-3 脆弱性如实披露**: 判据=全型公共 TOP 子集 ≤ topk/2=4 即阳性, 实测公共集 3 条且全为 position-1 惰性突变(A1C/A1G/A1U), 仅以 1 条余量通过; k=3 vs k=2 轮廓系数差仅 0.006; 型 0/2(覆盖 77/86 spacer)TOP 并集 Jaccard 1.00/0.78, 型间差异主要由 junction-pairing 驱动的 n=9 小簇贡献——判定为**弱证据阳性**(判据登记 docs/preregistration.md §A); **round-3 事后敏感性补救**(post-hoc 非预登记, `scripts/crrna_typing_robustness.py` → `data/context_typing/typing.robustness.json`, 判读文字全部由数值按显式阈值规则生成; 档案判据在当前代码下不可原样复现——档案 top-8 中 6 条变体已被晚于档案的 v1.6/v1.7 硬过滤从突变库剔除, 现测 k=8 公共集 4 条非 3 条, 以下均为当前代码口径): ①置换检验——9 代表 TOP 列表随机分组(280 种无序划分精确枚举)零分布中位 8 条, 观测 4 条 P(common≤4)=0.0357, 判据零通过率 3.6%, 判据具判别力; ②topk 敏感性——topk=6/8/10/16 公共集 3/4/5/6 条(阈 3/4/5/8)全部通过, 不随 topk 翻动; ③去 position-1——剔除全部含 pos-1 突变变体后公共集 2 条(U3A/U3C, 同为惰性位点)仍通过, 不依赖 pos-1 单一位点; ④聚类谱——k=3 在手写 KMeans(0.379, 领先仅 0.006)与 Agglomerative-ward(0.435, 领先 0.015)下最优, 但 GaussianMixture 下 k=2 最优(0.392 vs k=3 0.005), 选型依赖方法; 判据面 k=2 不通过(公共集 9 条)、k=3..6 通过——判据随 k 翻动, 阳性是 k≥3 选型下的结果; ⑤n=9 小簇审查——型1 横跨 4 个独立基因(CCNE1/EGFRdel/MYC/R280K), 簇内平均同一性 0.416(背景 0.449, 非同源), 分离由 junction_pairs 主导(标准化差 18.78), 判为"来源分散但 junction-pairing 特征一致的真实亚型"而非批次效应, 可主动扩充候选 M246I gRNA6/R248Q gRNA20; **综合判定: 证据不变**(弱证据阳性评级维持, 触发项: 判据随 k 翻动) |
| 游离态结构特征→guide 活性(多体系方法学先例) | **该规律有反复验证的领域先例, 非本项目独创主张**: Moreno-Mateos et al. 2015, Nat Methods 12:982-988 (CRISPRscan, 被引>1600: 序列+结构特征预测 Cas9 sgRNA 体内活性); Wong et al. 2015, Bioinformatics 31:3105 (WU-CRISPR: 二级结构特征区分功能/非功能 guide); Liu et al. 2020, Bioinformatics 36:1454 (CRISPRpred: 加入二级结构特征后预测更好); Zhang 实验室专利 US8889356/US8932814 (MFE+Boltzmann 配分比较 guide 二级结构——与本仓库同类方法); Creutzburg 2020 论文自报正相关(pre-crRNA 口径 n≈25 显著)。**定位: 多体系先例支持该特征类与活性相关; 本仓库路径A复算(下行)划定其口径边界(pre-crRNA、同体系内), 主张=同上下文相对排序的保守使用** |
| 游离态指标→活性 反向验证(路径A, 阴性校准) | **Creutzburg 2020 Fig 1D 逐条反向打分**(本仓库 `scripts/crrna_creutzburg_reverse.py` / `data/creutzburg_reverse_validation.json`): 14 spacer × 构建内真实 DR(反向引物公共尾 revcomp 推得 AATTTCTACTGTTGTAGA±T, 与注册表 zeng2026 条目一致; 差异位点与正文三重吻合: Sp8/Sp12=0基 10,11,19,20, Sp12/Sp14=5,7,9, 位置19=K(G/U))——**本仓库自行复算未能复现该文自报的正相关**: spacer_up ρ=−0.45 / stem_prob ρ=−0.51(mature19 口径, tie-aware), 好坏 AUC 0.28–0.33(方向反转); 此前"5端种子区可及性弱正向(seed5_up ρ=+0.41)"经 round-3 R1 修复(seed5_unpaired 与注册表同源的 bpp 行-only bug)后**不再成立**——修正后 seed5_up ρ=+0.21/+0.14/−0.10(三口径, 置换 p=0.45/0.64/0.73; 该 JSON 的 stats 块经独立克隆复算发现为"修复 rows + 旧 tie-blind spearman"中间态, 2026-09-04 已用终版脚本重跑刷新), 方向混合不显著, 不作为任何排序依据; 机制个案复现(Sp8 MFE=spacer 自缠发夹, Sp12 全游离); 好 guide Sp4/Sp9 游离态茎概率≈0.001 但活性 72–83%(蛋白重折叠可挽救茎)。对照: 该文自报特征相关为正(pre-crRNA 口径, n≈25, 茎含手工定义假结 P1 配对, ViennaRNA 不可算假结)。**校准结论: 游离态指标仅用于同 spacer 下 DR 变体相对排序(结构口径), 不主张游离态结构预测活性; 排序主张(TOP1>WT)未经实验标定——Creutzburg 层只支持特征类型与活性相关(论文自报 n约25 显著; 自有复算 n=14 不显著、方向混合, 两口径并列引用), MD 阳性仅为必要条件, 最终确认须体外旁切实验**。Sp8 盲区(复算最有价值发现): 唯一 0% 活性样本未被 spacer_up(全长平均稀释)/stem_prob(DR 自配与 DR-spacer 错配混同)抓住, 已据此新增 **DR-spacer 交叉配对硬过滤** v1.6(变体 MFE 交叉配对 nt 不得多于 WT; Creutzburg 设计守则的保守实现, 定位=无损工程规则, 非活性预测器)。**回验结果为阴性且方向反转**: 三口径 cross 均未抓住 Sp8。机制诊断(2026-09-01): Cas12a 家族 crRNA 的活性态含 5' 假结, 在伪结外折叠空间不可表示, Sp8/Sp12 之争本质是与假结活性态的平衡竞争, 伪结外模型按构造无法打分; 次优空间中侵占型嵌套中间态在好坏 guide 均大量出现(Sp8 93/318 vs Sp12 146/391, Sp12 反而更多), 故任何伪结外交叉指标对该案例无区分力。cross_pp 与活性的显著正相关(mature rho=+0.61, p=0.025)为成分混杂: spacer 5' A 富集与被加工丢弃的 U-rich 残留 repeat 配对, GC vs 活性 rho=-0.32 同向, n=14。Sp8 盲区仍开放。过滤回溯审计(`data/cross_audit.json`): 两主线口径 0 拦截/TOP-12 不受影响; **t1_MYCg1 口径 WT 自身 cross_nt=20**(MYC spacer 与 DR 大量交叉配对, 型1 上下文注脚), 拦截 47/118 条含 TOP-12 中 6 条茎稳定化候选——通用候选 A8G+U15C 在型1 口径不再通过 v1.6, 其通用性主张须降级 |

> 引用层级声明: 机制前提的证据链为 **Cas12a 同家族(Creutzburg 2020)→ Cas12a2(本体系, 暂无直接折叠研究)**, 而非跨体系的 SpCas9→Cas12a2; 硬过滤"结构保持/spacer 游离"两项由此获得同家族文献锚定而非纯启发式。 路径A反向复算为阴性(见上行): 游离态全局指标不构成跨体系活性预测器, 打分体系定位=结构口径相对排序+文献机制锚定。

## 快速开始

**环境前置(重要)**: 主管线需要 **Python ≥3.8 + ViennaRNA ≥2.6**(开发验证环境: `/public/home/mengxl/dzy/envs/guideforge`, py3.10.20 + ViennaRNA 2.7.2)。宿主机裸 `python` 是 2.7, 直接跑会报编码错误——请用上述环境解释器或任何满足前置的 py3 环境。路径类默认值(如 rnet 引擎)全部可用 CLI 参数覆盖; 仓库内少数数据 JSON 的绝对路径为溯源记录, 换机重跑时按需替换。

```bash
PY=/public/home/mengxl/dzy/envs/guideforge/bin/python   # 或你自己的 py3+viennarna 环境

# 注册表契约测试 + 主管线黄金集回归(~10 秒级)
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
$PY scripts/crrna_ensemble_check.py           # ΔΔG 系综可靠性(R1+R2'/R3'+R4 全部通过, passed_all_rules; 旧 R2/R3 已退役留档)
$PY scripts/crrna_creutzburg_reverse.py       # 路径A 阴性校准(方向反转主张)
$PY scripts/crrna_cross_decompose.py          # 交叉特征组成去混杂(rho +0.55→+0.17)
$PY scripts/crrna_state_competition.py        # Sp8 状态竞争 FAIL 判定
$PY scripts/crrna_pk_competition.py           # Knotty DP09 FAIL 判定(容器)
$PY scripts/crrna_nupack_pk_weight.py         # NUPACK 第四模型 FAIL 判定
$PY scripts/crrna_negative_control.py         # 19×上下文负对照(pdbdr.v2 口径 V3 FAIL, zengdr.v6 全 PASS; 判据见 docs/preregistration.md §D)
$PY scripts/crrna_han2025_features.py         # Han2025 数据链(Fig1g 145/Sanger 54/cis-trans)
$PY scripts/crrna_train_selector.py           # 选型器训练+204 候选应用+耐受集体检(tie-aware)
$PY scripts/crrna_selection_baseline.py       # R3 filter+random 基线检验
$PY scripts/crrna_score_ablation.py           # R3 逐项消融 + ddG 对称性分析
$PY scripts/crrna_chai_matrix_summary.py      # Chai 4 靶矩阵 summary 程序化重生成
$PY scripts/crrna_af3_input.py                # AF3 服务器就绪检查输入 13 套(V3 表1 口径, template-free 独立对照)
$PY scripts/crrna_protenix_input.py           # Protenix(AF3 开源替代)输入 13 套转换, 云端手册 tools/protenix_cloud_runbook.md
$PY scripts/crrna_compensatory_design.py      # 补偿臂去混杂版面板重生成
$PY scripts/crrna_bayesopt.py --cold-start    # 冷启动 EI; --prior-han 同源先验版
$PY scripts/crrna_virtual_cell.py --prior-lit # 文献标定场景表
$PY scripts/crrna_homolog_context_interaction.py  # 同源跨上下文交互检验(E1 PASS rho=0.257[CI 0.247,0.267]; 判据=项目设定常数(与评估同提交, 2026-09-09 审计后如实改标))
$PY scripts/crrna_homolog_interaction_audit.py    # 上述结果的独立复算审计(不 import 被测脚本)
$PY scripts/crrna_ivt_template.py --anova data/ivt_round1.csv   # 矩阵交互项判定【数据待补: ivt_round1.csv 尚不存在, 为占位标记】
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
**sha256 校验现状(round-3 起)**: BIG_FILES.md 中 RibonanzaNet.pt 为完整 64 位可校验; 其余 9 件目前为截断前缀(† 标注, 仅可前缀比对)——完整 64 位值需从 Release 上传机补录, 见 BIG_FILES.md 注记。
双源说明：Gitee Release `data-v1` 已附 7 件单附件（RibonanzaNet×4 / gRNAde / OpenKnot M2·M2R）；
其余 3 件超 Gitee 单附件 100MB 上限（Struct2SeQ.pt / Struct2SeQ_SHAPE.pt / OpenKnotBench CSV），
已按 95MB 分片（.part00/.part01）同附于 [Gitee Release data-v1](https://gitee.com/eastern-zhongyuan/guideforge/releases/data-v1)——
下载全部分片后执行 `cat 文件名.part* > 文件名` 拼回，sha256 校验值见 `BIG_FILES.md`（拼回一致性已验证）；
也可从 [GitHub Release data-v1](https://github.com/dong-zhongyuan/guideforge/releases/tag/data-v1) 直接下载全量单件。

## 目录

- `configs/scaffold_registry.json` — 骨架注册表(cas12a2 PDB 结构口径 / cas12a2_zeng2026 细胞实验口径, 均 gate=confident; 含可消费 `pfs` 字段(共识+tolerant_mismatches), round-3 R1 已重注册修正基线口径并记录 vienna_version, 漂移测试守护)
- `scripts/` `src/` `tests/` — 管线代码与契约测试; 主管线 `crrna_scaffold_design.py`,
  文献验证 `crrna_dmytrenko_validate.py`, 上下文分型 `crrna_context_typing.py`,
  设计智能体 `crrna_design_agent.py`(--spacer 直输选骨架 / 突变转录本 tilling 设计;
  阴性退化模式: `--degrade auto/on/off`, 判据=预登记 §A go/no-go, 阴性时退化通用型单骨架)
- `data/` — 8D4A 结构、DR-蛋白接触表、TP53 WT/R248Q 转录本、正式候选库产出
  (zengdr.v2 = 加工保护+3'窗+共变口径; context_typing/ = 分型试验; agent/ = 智能体演示;
  **zengdr.v6 已于 round-3 用当前 HEAD 重生成: 排名/分数/library_sha256 与旧版零差异, 补齐 cross_nt/inv_max_run/partner_switch 列与 model_domain 声明**)
- `toolbox/` — 结构工具链: rnet-inference (RNet/RibonanzaNet) · Struct2SeQ · geometric-rna-design (gRNAde) · OpenKnotAIDesignData(竞赛基准数据)

## 诚实声明与边界限定

候选库是"候选压缩"结果, 不是活性预测; 打分为透明启发式, 未经实验标定。
最终活性/特异性以体外生化与细胞实验为准。**判据治理(round-3 起)**: 全部 PASS/FAIL 判据集中登记于 `docs/preregistration.md`(2026-09-04, 区分真预登记/项目设定常数/事后导出三级); 判据变更须先登记后评估; 失败的预登记判据(如旧系综 R2/R3(已退役并重设计为 R2′/R3′)、阴性对照 pdbdr.v2 口径 V3)在头条如实报告而非埋在 JSON。两条必须主动声明的口径边界:

- **假结**: Cas12a2 crRNA 5' 端假结为蛋白诱导的结合态结构(Dmytrenko 2023; Methods in Enzymology 2025),
  游离态预测中不出现——本管线双路结构预测(ViennaRNA × RNet-SS)均为游离态口径,
  ViennaRNA 本身不折叠假结; 注册表条目的 bound_state_note 已如实标注。
- **游离态 ≠ 结合态**: 全部结构指标(含种子区游离 d_seed)为溶液态预测,
  结合态下蛋白会重排 spacer(Bravo 2023 四元复合物结构即证据); 正向项的含义
  是"预测降低溶液态自缠结、加速种子区暴露", RNet SHAPE 一致性为桥接证据,
  不是结合态活性保证。

> 答辩口径: 主结论级质疑的「领域先例 → 我们的实例 → 边界声明」三段式应答
> 见 `docs/defense_prep.md`(2026-09-08, 文献锚定版)。

## 来源

- Cas12a2 结构: PDB 8D49/8D4A/8D4B (Bravo et al. 2023, Nature 613:582)
- Cas12a2 功能/DR 互换/加工位点: Dmytrenko et al. 2023, Nature 613:588-594
- 细胞实验 DR 与 R248Q spacer: Zeng et al. 2026, Nature (10.1038/s41586-026-10738-7) 补充材料
- spacer-DR 互扰(同家族直接证据): **Creutzburg et al. 2020, NAR 48(6):3228-3243** (PMC7102956)
- spacer-骨架误折叠(SpCas9 体系, 跨体系佐证): Bush et al. 2023, Cell Chem Biol 30:879-892 (PMC10529641)
- DR 下游结构敏感性: Liao et al. 2018 (PMC6546362)
- Cas12a DR 序列工程先例: Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434); DR 3' 端敏感性(化学修饰): Zhang et al. 2025 (PMC11780881)
- Cas12a DR 突变提活先例: Teng et al. 2019, Genome Biology (4n96 骨架, PMID 30717767)
- guide 结构特征→活性 先例: Moreno-Mateos et al. 2015, Nat Methods (CRISPRscan); Wong et al. 2015, Bioinformatics (WU-CRISPR); Liu et al. 2020, Bioinformatics (CRISPRpred)
- spacer-guide 耦合为活性决定因素(综述): Annu Rev Chem Biomol Eng 2023 (10.1146/annurev-chembioeng-100522-114706)
- Cas12a 骨架变体工具箱(按上下文选骨架工程先例): MIT 学位论文 (dspace.mit.edu 1142190545); DR 茎环保守性实验先例: EP4403638
- 共折叠工具链: Chai-1 (PoseBusters 77%≈AF3 76%); 差量筛选工作流先例: arXiv:2511.14669 (AF3 ipTM 筛选 + Rosetta 复核)
- crRNA 预折叠驱动 RNP 组装: Sudhakar et al. 2023 (PMC10200996; 因果方向: RNA 预折叠→诱导蛋白构象变化)
- 结构工具链: Science 2026 (10.1126/science.aeg6829, Das lab, OpenKnot); gRNAde 权重: HuggingFace chaitjo/gRNAde
