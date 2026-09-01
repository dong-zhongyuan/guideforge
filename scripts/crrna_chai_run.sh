#!/bin/bash
# Chai-1 共折叠运行器(容器, 参数化; 2026-09-02 正式化——此前 ESM 版由 /tmp 临时 sed 生成, 纪律违规)
# 用法: bash crrna_chai_run.sh <fasta路径> <输出目录> <0|1=use_esm> [gpu_id]
# 依赖: ~/dzy/envs/chai; ESM/conformers 资产已软链于 ~/chai_downloads(经 /dawn 只读共享)
set -e
FASTA=${1:?fasta 路径}
OUT=${2:?输出目录}
ESM=${3:-0}
GPU=${4:-0}
export HF_ENDPOINT=https://hf-mirror.com
export CUDA_VISIBLE_DEVICES=$GPU
export CHAI_DOWNLOADS_DIR=/home/node/chai_downloads
mkdir -p "$OUT"
~/dzy/envs/chai/bin/python - "$FASTA" "$OUT" "$ESM" << 'EOF'
import sys
from pathlib import Path
import torch
from chai_lab.chai1 import run_inference
fasta, out, esm = sys.argv[1], Path(sys.argv[2]), bool(int(sys.argv[3]))
r = run_inference(
    fasta_file=Path(fasta), output_dir=out,
    use_msa_server=False, use_esm_embeddings=esm,
    device=torch.device("cuda:0"))
print("scores:", r.scores)
EOF
echo CHAI_RUN_DONE
