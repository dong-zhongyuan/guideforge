# -*- coding: utf-8 -*-
"""
Struct2SeQ CPU/GPU 兼容启动器
=============================
官方 generate.py 硬编码 .cuda()，本包装器在无 CUDA 时把 .cuda() 降级为 no-op，
使代码在纯 CPU 环境也能运行（有 CUDA 时行为不变，自动走 GPU）。

用法（在 rnet conda 环境中，仓库根目录下）:
  python demo/run_struct2seq.py [--n_structures 1]
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S2S_DIR = os.path.join(ROOT, "Struct2SeQ")
os.chdir(S2S_DIR)                      # generate.py 用相对路径读 test10_configs/
sys.path.insert(0, S2S_DIR)

import torch

if not torch.cuda.is_available():
    print("[run_struct2seq] 未检测到 CUDA，.cuda() 将降级为 CPU 运行（速度较慢）")
    torch.Tensor.cuda = lambda self, *a, **k: self          # type: ignore
    torch.nn.Module.cuda = lambda self, *a, **k: self       # type: ignore
    _orig_load = torch.load
    def _load_cpu(*a, **k):                                  # 权重在 GPU 上保存，需映射到 CPU
        k.setdefault("map_location", "cpu")
        return _orig_load(*a, **k)
    torch.load = _load_cpu                                   # type: ignore
else:
    print(f"[run_struct2seq] 使用 GPU: {torch.cuda.get_device_name(0)}")

RNET_W = os.path.join(ROOT, "rnet-inference", "RibonanzaNet-Weights")

sys.argv = [
    "generate.py",
    "--target_df", "demo_targets.csv",       # P20 'Kissing Multiloops' (100nt 假结)
    "--output_csv", "p20_design_results.csv",
    "--out_folder", "output_results",
    "--gpu_id", "0",
    "--n_structures", "1",
    "--up_bias", "0.0",
    "--weights_path", "Struct2SeQ_SHAPE.pt",  # 论文 Round 4 中超越人类专家的变体
    "--rnet_weights", os.path.join(RNET_W, "RibonanzaNet.pt"),
    "--rnet_ss_weights", os.path.join(RNET_W, "RibonanzaNet-SS.pt"),
]

with open(os.path.join(S2S_DIR, "generate.py"), encoding="utf-8") as f:
    code = compile(f.read(), "generate.py", "exec")
exec(code, {"__name__": "__main__", "__file__": os.path.join(S2S_DIR, "generate.py")})
