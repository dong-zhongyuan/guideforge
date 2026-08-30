# OpenKnot AI 工具箱 — Science 2026 封面论文复现资源

论文: Townley, Kladwang, Baker … Das R. *De novo design of RNA pseudoknots with deep learning*.
**Science** (2026), DOI: [10.1126/science.aeg6829](https://www.science.org/doi/10.1126/science.aeg6829)
(PMC 全文: https://pmc.ncbi.nlm.nih.gov/articles/PMC13228335/)

本目录收集了该论文涉及的全部开源工具、权重与基准数据（2026-08-30 下载）。

## 目录结构

| 目录 | 内容 | 来源 |
|---|---|---|
| `OpenKnotAI/` | **论文官方代码包**（Zenodo 10.5281/zenodo.20649966，2026-06-11 快照）：`design-methods/`（gRNAde、Struct2SeQ、NA-MPNN、RFDpoly、3DRNA）+ `das-eterna/`（ubr、cmuts、fld、OpenKnotScoreMATLAB、OpenKnotScorePipeline） | Zenodo |
| `rnet-inference/` | **RNet (RibonanzaNet)** 推理封装 + 4 个官方权重（`RibonanzaNet-Weights/`：SHAPE / 二级结构 SS / 降解 Deg / Dropout） | github.com/DasLab/rnet-inference |
| `Struct2SeQ/` | **Struct2SeQ**（深度 Q 学习逆折叠，Round 4 超越人类专家）+ 权重 `Struct2SeQ.pt`、`Struct2SeQ_SHAPE.pt` | github.com/Shujun-He/Struct2SeQ + Kaggle 权重 |
| `NA-MPNN/` | **NA-MPNN**（论文中的 MPNN-fixbb，Baker 实验室），自带权重 `models/` | github.com/baker-laboratory/NA-MPNN |
| `geometric-rna-design/` | **gRNAde**（几何深度学习逆折叠，含官方 `src/openknot_score.py` 评分实现） | github.com/chaitjo/geometric-rna-design |
| `OpenKnotAIDesignData/` | **全部竞赛数据**：`Targets/`（Round 1–4 共 57 个假结靶点）+ `Data/`（约 5 万条序列的 SHAPE / M2 / M2R 实验数据，LFS 已校验） | github.com/eternagame/OpenKnotAIDesignData |
| `OpenKnotScoreMATLAB/` | OpenKnot 评分器 MATLAB 版 | github.com/eternagame/OpenKnotScoreMATLAB |
| `demo/` | 本地复现演示脚本（见下） | 本目录自建 |

## 运行环境

- conda 环境 `rnet`：Python 3.11 + PyTorch (CUDA 12.8，适配 RTX 5060) + einops/matplotlib/pandas/pyyaml/tqdm/arnie
- 激活: `conda activate rnet`

## 快速开始

### 1. RNet 预测任意 RNA 的 SHAPE 图谱（论文的"模拟实验"核心）
```bash
cd rnet-inference
python src/rnet_shape.py GGGGAAAACCCC     # SHAPE (2A3 + DMS)
python src/rnet_2d.py GGGGAAAACCCC       # 二级结构 + 置信度
```

### 2. 复现 P20 演示（实验数据 vs 本地 RNet 预测）
```bash
python demo/demo_openknot_p20.py
```
输出: 各 AI 方法最优设计的 实验OpenKnot分 / RNet预测分 / 结构F1 对照表 + SHAPE 对比图。

### 3. Struct2SeQ 从头设计假结序列
```bash
cd Struct2SeQ
python make_dummy_arnie.py   # 一次性
python generate.py \
    --target_df demo_targets.csv \
    --output_csv results.csv --out_folder output_results \
    --gpu_id 0 --n_structures 100 --up_bias 0.0 \
    --weights_path Struct2SeQ_SHAPE.pt \
    --rnet_weights ../rnet-inference/RibonanzaNet-Weights/RibonanzaNet.pt \
    --rnet_ss_weights ../rnet-inference/RibonanzaNet-Weights/RibonanzaNet-SS.pt
```
`demo_targets.csv` 内含 Round 3 的 P20 "Kissing Multiloops" 靶点（100nt 假结）。
输入 CSV 需含列: `Title`, `Dot-bracket`（支持 `[]{}` 假结记号）, `wild_type_sequence`。

### 4. NA-MPNN (MPNN-fixbb) — 需要 3D 骨架输入
```bash
conda create -n NA-MPNN -y && conda activate NA-MPNN
conda install -c pytorch -c nvidia -c conda-forge openbabel ProDy pyarrow pandas pytorch pytorch-cuda -y
cd NA-MPNN
python ./inference/run.py --model_type "na_mpnn" --mode "design" \
    --pdb_path "./inference/examples/4oqu.pdb" --out_folder "./out/design"
```

### 5. gRNAde — 建议 Linux/WSL2
依赖 PyTorch Geometric 全家桶（torch-scatter/sparse/cluster），官方在 Linux 测试；
模型权重在 HuggingFace（国内镜像: `HF_ENDPOINT=https://hf-mirror.com`）。
OpenKnot 基准示例见 `geometric-rna-design/projects/openknot_benchmark/design.ipynb`。

## 评分定义（OpenKnot score, 0–100）
- 靶结构中**未配对**核苷酸: SHAPE > 0.125 记为正确
- 靶结构中**配对**核苷酸: SHAPE < 0.5 记为正确（交叉对子分数用 < 0.25）
- `openknot_score = 0.5 × eterna_score + 0.5 × crossed_pair_quality_score`
- 实现: `geometric-rna-design/src/openknot_score.py`；分数 > 90 视为设计成功

## 备注
- Struct2SeQ 权重也可从 Kaggle 获取: `kaggle datasets download -d shujun717/struct2seq-weights`
- 冷冻电镜结构: PDB 1OZT (gRNAde) / 1OZU (MPNN-fixbb) / 11EH (Struct2SeQ-SHAPE) / 11AG (MPNN-RFdiff 二聚体)
- SHAPE 原始数据: RNA Mapping Database (rmdb.stanford.edu), accession OK45LIB/OK6LIB/OK7ALIB/OK7BLIB
- Rosetta 因商业许可未再分发（见 `OpenKnotAI/LINKED-TOOLS-rosetta.md`）
