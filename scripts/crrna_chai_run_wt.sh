#!/bin/bash
# Chai-1 首次共折叠: SuCas12a2 + crRNA(DR18+R248Q spacer) + 靶RNA (8D4A 口径)
# 单序列模式(无 MSA/ESM, 避免额外下载), GPU0
set -x
export HF_ENDPOINT=https://hf-mirror.com
export CUDA_VISIBLE_DEVICES=0
mkdir -p ~/chai_out/WT2
~/dzy/envs/chai/bin/python - << 'EOF'
from chai_lab.chai1 import run_inference
import torch
from pathlib import Path
r = run_inference(
    fasta_file=Path("/dawn/chai_input/WT.fasta"),
    output_dir=Path("/home/node/chai_out/WT2"),
    use_msa_server=False,
    use_esm_embeddings=False,
    device=torch.device("cuda:0"),
)
print("scores:", r.scores)
EOF
echo CHAI_WT_DONE
