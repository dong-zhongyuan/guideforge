# GuideForge 复现说明（Cas12a2 crRNA 骨架优化计算管线）

本包是 Gitee 仓库 https://gitee.com/eastern-zhongyuan/guideforge 工作树的
自包含快照（代码 + 小数据 + 文档），按 Tier 分级复现。**Tier 0 在全新电脑上
纯 CPU、装依赖后无需网络即可全绿**；Tier 1/2 需要联网下载大数据；Tier 3 需要
外部 GPU/容器/服务器机时，输入文件已备在包内，不要求可跑。

包内全部文件的 sha256 见 `MANIFEST.sha256`。

## Windows 中文环境必看（GBK 坑）

Windows 默认代码页是 GBK，Python 3.13 之前 `open()` 默认按 GBK 读写文本，
本仓库大量 UTF-8 中文注释/数据会直接 `UnicodeDecodeError`。两个办法任选：

- 跑任何脚本前设环境变量（Git Bash: `export PYTHONUTF8=1`；
  CMD: `set PYTHONUTF8=1`；PowerShell: `$env:PYTHONUTF8=1`）——**推荐**；
- 或用 Python >= 3.15（默认 UTF-8 模式）。

包内源码读写文本均已显式 `encoding="utf-8"`，但 ViennaRNA 等 C 扩展内部
打印仍受区域设置影响，设 `PYTHONUTF8=1` 最省事。另外：**用 `python`，不要用
`python3`**（Windows 上 `python3` 常是 WindowsApps 的商店残桩，一跑就弹商店）。

## Tier 0 — 核心复现（纯 CPU、装依赖后无网络）

```bash
python -m venv .venv
# Windows Git Bash:
source .venv/Scripts/activate
# Linux/macOS: source .venv/bin/activate

python -m pip install -r requirements.txt
# 国内网络慢可加: -i https://pypi.tuna.tsinghua.edu.cn/simple

export PYTHONUTF8=1   # Windows 必设
python tools/verify_repro.py
```

`verify_repro.py` 逐项输出 PASS/FAIL 并以退出码汇报：

a. 依赖导入检查（RNA/numpy/Bio/scipy/sklearn/openpyxl/matplotlib/flask）；
b. `python -m pytest tests/ -q` 全量测试（21 项）；
c. 核心管线 smoke（全部输出到系统临时目录，**不污染包内 data/**）：
   - `crrna_scaffold_design.py` 最小运行（仅 ViennaRNA，cas12a2_zeng2026 +
     R248Q spacer，topk 12）成功产出 variants.csv/top.json/top.fasta；
   - `crrna_specificity_scan.py` 对包内 TP53 WT/R248Q 两条小 fasta 扫描跑通；
   - `crrna_ivt_template.py` 重生成 `ivt_round1_template.csv` +
     `ivt_round1_order_sheet.csv`，与包内 `data/` 同名文件**逐字节一致**；
   - `crrna_chai_matrix_summary.py` 重跑，产物与包内
     `data/chai_matrix_4t_summary.json` **逐字节一致**
     （该汇总为确定性程序化生成，无时间戳/随机数）。

已实测环境：Windows 10/11 + Python 3.13.9（Anaconda 底包新建的 venv，
`pip install viennarna` 直接有 2.7.2 wheel，无需 conda）。Linux/macOS
Python >= 3.10 同样适用；ViennaRNA 也可走 `conda install -c bioconda viennarna`。

## Tier 1 — 联网数据下载

```bash
export PYTHONUTF8=1
python tools/fetch_transcriptome.py
```

下载 GENCODE 转录组 FASTA（优先 .gz 再解压，校验文件大小，来源 URL 记入
`data/download_log.json`）：

- `data/raw/gencode.v47.transcripts.fa`（约 613MB，EBI GENCODE release_47）；
- `data/transcriptome/gencode.v49.pc_transcripts.fa`（约 683MB）与
  `gencode.v49.lncRNA_transcripts.fa`（约 224MB）（GENCODE release_49）。

NCBI efetch（联网拉取基因 mRNA，供 `scripts/crrna_panel_targets.py` 的
`--gene` 模式与 `crrna_specificity_scan.py --gene` 使用）无需额外安装，
脚本内置 urllib 调用 https://eutils.ncbi.nlm.nih.gov ；被墙时自备代理或改用
`--fasta` 喂本地序列。

## Tier 2 — 引擎环境（逆折叠 / SHAPE 筛选子进程）

```bash
python -m venv .venv-engine        # 建议 python 3.11
source .venv-engine/Scripts/activate
python -m pip install -r requirements-engine.txt
python tools/fetch_big_files.py    # 10 件权重/基准数据, sha256 校验
```

- `requirements-engine.txt`：torch / einops / arnie==0.2.7 / torch_geometric /
  torch_scatter / torch_cluster / biotite 等，安装备注见文件头注释（CPU 版
  torch 索引、pyg wheel 匹配、cpdb-protein 可跳过——
  `toolbox/geometric-rna-design/src/data/data_utils.py` 已内置 biotite 兜底）。
- `tools/fetch_big_files.py` 按 `BIG_FILES.md` 从 Gitee Release `data-v1`
  下载 10 件资产到对应相对路径；3 件超 100MB 的（OpenKnotBench csv /
  Struct2SeQ.pt / Struct2SeQ_SHAPE.pt）为 part00/part01 分片，自动拼回。
  **校验口径如实声明**：仅 `RibonanzaNet.pt` 有完整 64 位 sha256（✅）做全量
  校验；其余 9 件原记录只有前 16 位截断值（†），只能做前缀比对，脚本输出
  会注明"仅前缀比对，不构成完整完整性校验"。

引擎跑通后按 `README_移交说明.md` 的"完整四引擎 + SHAPE 筛选"命令，用
`--s2s-python` 指向引擎环境 python。

## Tier 3 — 外部 GPU / 容器 / 服务器（输入已备，需外部机时）

以下步骤依赖包外环境，**不要求在新机上可跑**，输入文件已全部备在包内：

- **Chai-1 共折叠**：Linux 容器环境（原机 `~/dzy/envs/chai`）。输入：
  `data/chai_input_WT.fasta`、`data/chai_cofold_*.json`、`data/chai_template_hits.m8`；
  运行脚本：`scripts/crrna_chai_run.sh` / `crrna_chai_run_wt.sh`；
  汇总：`scripts/crrna_chai_matrix_summary.py`（Tier 0 已验证可重算）。
- **Boltz 交叉验证**：输入与清单见 `data/boltz_xval/`。
- **AlphaFold Server (AF3)**：在线服务，输入见 `data/af3_inputs/`，
  生成脚本 `scripts/crrna_af3_input.py`。
- **OpenMM MD**：建议独立 conda 环境
  （`conda install -c conda-forge openmm pdbfixer`，另需 parmed/MDAnalysis）；
  准备/运行/分析脚本：`scripts/crrna_md_prep.py` / `crrna_md_run.py` /
  `crrna_md_analysis.py` / `crrna_md_togmx.py`，清单见
  `data/md_prep_manifests.json`、`data/md_analysis_plan.json`。

## 产出物口径（诚实声明）

候选库是"候选压缩"结果，不是活性预测；打分是透明启发式，未经实验标定。
最终活性/特异性以体外生化与细胞实验为准。详见 `README_移交说明.md`。
