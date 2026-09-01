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

## 文献依据（对应模块）

| 模块 | 文献依据 |
|---|---|
| 加工位点保护硬过滤 / DR 互换不误杀验证 | Dmytrenko et al. 2023, Nature 613:588-594 (ED Fig.2c Cas12a↔Cas12a2 DR 互换功能保持; ED Fig.3 加工切点在 Cas12a 下游 1nt)。位点敏感性×保守性相关验证**仅 score 口径方向一致(ρ=−0.33, p=0.17, n=19 功效低)**, tolerance 口径为零相关(ρ=−0.01)——引用时只可引 score 口径并说明样本量 |
| 3'保守窗惩罚 | Dmytrenko et al. 2023, Fig.1c 跨家族 DR 3' 端保守/loop 可变(定性结论); **窗口大小 --cons3-window=5nt 为项目设定, 非文献出处**; 功能佐证: Zhang et al. 2025 (PMC11780881, DR 3' 端化学修饰可逆调控 Cas12a 活性——非序列突变, 佐证 3' 端敏感性) |
| spacer-DR 互扰 / spacer 游离度硬过滤 | **Creutzburg et al. 2020, NAR 48(6):3228-3243** (PMC7102956, 同家族直接证据): "base pairing of the direct repeat, other than with itself, was found to be detrimental"—DR 与 spacer 碱基配对损害 Cas12a 活性; Sp8 案例证实 spacer 侵占 DR 致假结破坏; 通过 spacer 尾部回折竞争性保护 DR 挽救活性(与本项目 spacer 游离度约束直接对应)。**Cas12a→Cas12a2 为同家族短迁移(DR 均含保守 5' 假结), 非跨体系假设** |
| spacer 上下文分型的机制前提 | 同上(Creutzburg 2020, 同家族); 补充: Bush et al. 2023, Cell Chem Biol 30:879-892 (SpCas9 体系, sgRNA 80nt 大骨架, 跨体系佐证); Liao et al. 2018 (PMC6546362, FnCas12a DR 下游发卡抑制切割)。分型主张的直接证据为本仓库分型试验计算结果; Cas12a2 本体系暂无直接折叠互扰研究(待验) |
| DR/骨架序列工程先例(Cas12a, 存在性证明) | Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434, DR 突变策略调编辑与检测) |
| "按靶标分型选骨架" | **无先例(空白区)**: 本仓库分型试验提供了反对"单一通用骨架"的计算证据(型间最优位点图谱结构性不同); 功能收益待湿实验 |
| 游离态指标→活性 反向验证(路径A, 阴性校准) | **Creutzburg 2020 Fig 1D 逐条反向打分**(本仓库 `scripts/crrna_creutzburg_reverse.py` / `data/creutzburg_reverse_validation.json`): 14 spacer × 构建内真实 DR(反向引物公共尾 revcomp 推得 AATTTCTACTGTTGTAGA±T, 与注册表 zeng2026 条目一致; 差异位点与正文三重吻合: Sp8/Sp12=0基 10,11,19,20, Sp12/Sp14=5,7,9, 位置19=K(G/U))——**本仓库自行复算未能复现该文自报的正相关**: spacer_up ρ=−0.45 / stem_prob ρ=−0.51, 好坏 AUC 0.23–0.38(方向反转); 仅 5端种子区可及性弱正向(seed5_up ρ=+0.41, AUC 0.68, 置换 p≈0.14); 机制个案复现(Sp8 MFE=spacer 自缠发夹, Sp12 全游离); 好 guide Sp4/Sp9 游离态茎概率≈0.001 但活性 72–83%(蛋白重折叠可挽救茎)。对照: 该文自报特征相关为正(pre-crRNA 口径, n≈25, 茎含手工定义假结 P1 配对, ViennaRNA 不可算假结)。**校准结论: 游离态指标仅用于同 spacer 下 DR 变体相对排序(结构口径), 不主张游离态结构预测活性; 活性主张须待本体系湿实验** |

> 引用层级声明: 机制前提的证据链为 **Cas12a 同家族(Creutzburg 2020)→ Cas12a2(本体系, 暂无直接折叠研究)**, 而非跨体系的 SpCas9→Cas12a2; 硬过滤"结构保持/spacer 游离"两项由此获得同家族文献锚定而非纯启发式。 路径A反向复算为阴性(见上行): 游离态全局指标不构成跨体系活性预测器, 打分体系定位=结构口径相对排序+文献机制锚定。

## 快速开始

```bash
# 注册表契约测试
PYTHONPATH=src python -m unittest discover -s tests

# 最小运行(纯 ViennaRNA, CPU)
python scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 --out-prefix data/run1

# 完整四引擎 + SHAPE 筛选
python scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 \
    --use-struct2seq --use-grnade --rnet-screen \
    --s2s-python /path/to/engine-env/bin/python --out-prefix data/run2
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
双源说明：Gitee Release `data-v1` 已附 7 件（ RibonanzaNet×4 / gRNAde / OpenKnot M2·M2R）；
其余 3 件超 Gitee 单附件 100MB 上限（Struct2SeQ.pt / Struct2SeQ_SHAPE.pt / OpenKnotBench CSV），
请从 [GitHub Release data-v1](https://github.com/dong-zhongyuan/guideforge/releases/tag/data-v1) 下载（全量 10 件都在）。

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
- 结构工具链: Science 2026 (10.1126/science.aeg6829, Das lab, OpenKnot); gRNAde 权重: HuggingFace chaitjo/gRNAde
