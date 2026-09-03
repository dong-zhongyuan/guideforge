# crRNA 管线移交包（Cas12a2 crRNA 骨架优化）2026-08-31

## 这是什么

AI 优化 Cas12a2 crRNA 骨架（DR 茎环）的计算管线：固定 spacer（TP53-R248Q 靶向序列
`GTTCATGCCGCCCATGCAGGAACT`），只改骨架。输入注册表中的野生型骨架，输出 WT + TOP-K
优化骨架变体（FASTA 可直接送合成），外加特异性预评估。

本目录是完整自包含的 crRNA 项目（Cas12a2 crRNA 骨架优化）。

## 目录

- `configs/scaffold_registry.json` — 骨架注册表（cas12a2 PDB 结构口径 18nt、
  cas12a2_zeng2026 细胞实验口径 19nt，均 gate=confident；另有 cas12a 占位与 spcas9 对照条目）
- `scripts/`
  - `crrna_scaffold_design.py` — 主入口：突变库生成 + 打分（ViennaRNA, CPU）
  - `crrna_specificity_scan.py` — 转录组错配扫描（特异性预评估，纯序列）
  - `register_scaffold.py` — 新骨架注册流程（ViennaRNA + RNet-SS 子进程复核）
  - `scaffold_registry.py` — 注册表只读接口
  - `crrna_struct2seq_gen.py` / `crrna_grnade_gen.py` / `crrna_rnet_shape.py` — 三个引擎/筛选的包装器（在 rnet 环境运行，由主入口子进程调用）
- `src/scaffold_registry.py` — 接口副本（项目惯例）
- `tests/test_scaffold_registry.py` — 注册表契约测试 6 项
- `data/` — 8D4A 结构、DR-蛋白接触表（按 DR 序列分文件缓存）、TP53 WT/R248Q 转录本、
  两套正式候选库产出（`tp53_r248q_zengdr.*` / `tp53_r248q_pdbdr.*`）与特异性扫描结果
- `toolbox/` — 结构工具链（见"来源"）

## 环境（两个 conda/venv）

主环境（运行主管线）: python 3.10+, `ViennaRNA`(import RNA)、numpy、biopython。
引擎环境（逆折叠/SHAPE 筛选子进程）: python 3.11, torch(cpu 即可)、einops、pandas、
pyyaml、tqdm、scipy、arnie==0.2.7、matplotlib、polars、torch_geometric、torch_scatter、
torch_cluster、biotite、cpdb-protein、python-dotenv、wandb。

安装备注（服务器实测）:
- pip 建议加 `-i https://pypi.tuna.tsinghua.edu.cn/simple`;装不上的包加
  `--only-binary :all:`；源码编译失败时设 `TMPDIR` 到可写目录。
- torch_scatter/torch_cluster 用 `-f https://data.pyg.org/whl/torch-<torch版本>+cpu.html`。
- `cpdb-protein` 若编译失败：本包 `toolbox/geometric-rna-design/src/data/data_utils.py`
  已含 biotite 兜底补丁，无需处理。

## 快速开始

```bash
# 注册表自检
cd <包根目录> && PYTHONPATH=src python -m unittest discover -s tests

# 生成候选库(最小: 仅 ViennaRNA 突变扫描+模拟退火, 纯主环境)
python scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 --out-prefix data/run1

# 完整四引擎 + SHAPE 筛选(需要引擎环境, 用 --s2s-python 指向引擎环境 python)
python scripts/crrna_scaffold_design.py --effector cas12a2_zeng2026 \
    --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 \
    --use-struct2seq --use-grnade --rnet-screen \
    --s2s-python /path/to/engine-env/bin/python --out-prefix data/run2

# 特异性扫描(对 WT 与突变转录本)
cat data/tp53_mrna_NM000546.fa data/tp53_r248q_mrna_NM000546.fa > data/_wt_mut.fa
python scripts/crrna_specificity_scan.py --spacer GTTCATGCCGCCCATGCAGGAACT \
    --fasta data/_wt_mut.fa --pfs CAGAG --out data/scan1
```

注意:
- `--use-grnade` 仅对 PDB 结构口径 DR(cas12a2 条目, 18nt)有效;cas12a2_zeng2026
  的 DR 与结构差 5' 端 1nt + A10G, 无直接结合态几何, 会自动跳过并打印说明。
- 注册新骨架: `python scripts/register_scaffold.py --effector <名> ... --rnet-python <引擎环境 python>`
  (或设环境变量 RNET_PYTHON / RNET_2D_SCRIPT)。

## 产出物口径(诚实声明)

候选库是"候选压缩"结果, 不是活性预测;打分是透明启发式(ΔΔG、结构距离、spacer 游离度、
8D4A 蛋白接触、氢键互补、RNet SHAPE 一致性), 未经实验标定。最终活性/特异性以体外生化
与细胞实验为准(湿实验方案不随本仓库发布)。

## 来源与血缘

- Cas12a2 骨架/结构: PDB 8D49/8D4A/8D4B (Bravo et al. 2023 Nature 613:582)
- 细胞实验用 DR 与 R248Q spacer: Zeng et al. 2026 Nature (s41586-026-10738-7) 补充表 MOESM3
- toolbox: Science 2026 aeg6829 (Rhiju Das, OpenKnot) 开源生态
  (rnet-inference / Struct2SeQ / geometric-rna-design), gRNAde 权重来自 HF chaitjo/gRNAde

## 大文件

模型权重与基准数据不进 git 历史, 见仓库 Release `data-v1` 与 `BIG_FILES.md` 清单。
