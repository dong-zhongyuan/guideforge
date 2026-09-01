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
| 加工位点保护硬过滤 / DR 互换不误杀验证 | Dmytrenko et al. 2023, Nature 613:588-594 (ED Fig.2c Cas12a↔Cas12a2 DR 互换功能保持; ED Fig.3 加工切点在 Cas12a 下游 1nt) |
| 3'保守窗惩罚 | Dmytrenko et al. 2023, Fig.1c 跨家族 DR 3' 端保守/loop 可变(定性结论); **窗口大小 --cons3-window=5nt 为项目设定, 非文献出处**; 功能佐证: Zhang et al. 2025 (PMC11780881, DR 3' 端化学修饰可逆调控 Cas12a 活性——非序列突变, 佐证 3' 端敏感性) |
| spacer 上下文分型的机制前提 | Bush et al. 2023, Cell Chem Biol 30:879-892 (WT 骨架与 spacer 共折叠消除 Stem 1; **SpCas9 体系证据, sgRNA 为 80nt 大骨架而 Cas12a2 DR 仅 19nt 单茎环, 迁移为待验假设**); 综述 De Saeger et al. 2025 (MDPI)。分型主张的直接证据为本仓库分型试验计算结果 |
| DR/骨架序列工程先例(Cas12a, 存在性证明) | Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434, DR 突变策略调编辑与检测) |
| "按靶标分型选骨架" | **无先例(空白区)**: 本仓库分型试验提供了反对"单一通用骨架"的计算证据(型间最优位点图谱结构性不同); 功能收益待湿实验 |

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
- spacer-骨架共折叠互扰(SpCas9 体系机制前提): Bush et al. 2023, Cell Chem Biol 30:879-892 (PMC10529641)
- Cas12a DR 序列工程先例: Lin et al. 2018, Mol Ther; Han et al. 2025 (PMC12508434); DR 3' 端敏感性(化学修饰): Zhang et al. 2025 (PMC11780881)
- 结构工具链: Science 2026 (10.1126/science.aeg6829, Das lab, OpenKnot); gRNAde 权重: HuggingFace chaitjo/gRNAde
