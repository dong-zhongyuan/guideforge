# GuideForge — Cas12a2 crRNA 骨架优化管线

AI 辅助优化 **Cas12a2 crRNA 的骨架（direct-repeat 茎环）**，用于更高效、更特异的癌细胞杀伤。
设计变量只有骨架 DR 区；spacer 向导序列固定为 TP53-R248Q 靶向序列（24nt,
`GTTCATGCCGCCCATGCAGGAACT`, 来自 Zeng et al. 2026 Nature 补充表, 为该文最优选择性向导)。

## 管线流程

```
注册表(WT 骨架 + 固定 spacer)
  → 变体生成: A 突变扫描 | B 模拟退火 | A2 Struct2SeQ 逆折叠 | A3 gRNAde 三维几何逆折叠
  → 硬过滤(结构保持/spacer 游离/polyT·G)
  → 多维打分: ΔΔG · 结构距离 · 蛋白接触(8D4A) · 氢键互补 · RNet SHAPE 一致性
  → TOP-K + WT 候选库(FASTA 可送合成)
  → 特异性扫描(WT vs 突变转录本)
  → 湿实验闭环(体外旁切初筛 → 细胞杀伤验证)
```

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
- `scripts/` `src/` `tests/` — 管线代码与契约测试
- `data/` — 8D4A 结构、DR-蛋白接触表、TP53 WT/R248Q 转录本、正式候选库产出
- `toolbox/` — 结构工具链: rnet-inference (RNet/RibonanzaNet) · Struct2SeQ · geometric-rna-design (gRNAde) · OpenKnotAIDesignData(竞赛基准数据)

## 诚实声明

候选库是"候选压缩"结果, 不是活性预测; 打分为透明启发式, 未经实验标定。
最终活性/特异性以体外生化与细胞实验为准。

## 来源

- Cas12a2 结构: PDB 8D49/8D4A/8D4B (Bravo et al. 2023, Nature 613:582)
- 细胞实验 DR 与 R248Q spacer: Zeng et al. 2026, Nature (10.1038/s41586-026-10738-7) 补充材料
- 结构工具链: Science 2026 (10.1126/science.aeg6829, Das lab, OpenKnot); gRNAde 权重: HuggingFace chaitjo/gRNAde
