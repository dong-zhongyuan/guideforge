# -*- coding: utf-8 -*-
"""
OpenKnot P20 (Kissing Multiloops) 复现演示
==========================================
复现 Science 2026 (aeg6829) 论文的核心 in-silico 管线:
  1) 读取 OpenKnot 竞赛真实实验数据 (Round 3, 靶点 P20 'Kissing Multiloops')
  2) 用 RNet (RibonanzaNet) 在本地预测各方法最优设计的 SHAPE 图谱
  3) 用官方 OpenKnot 评分 (gRNAde src/openknot_score.py) 计算 in-silico 分数
  4) 用 RNet-SS 预测二级结构, 计算与靶结构的碱基对 F1
  5) 与实验 OpenKnot 分数对照

运行 (在 rnet conda 环境中):
  python demo/demo_openknot_p20.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "rnet-inference", "src"))
sys.path.insert(0, os.path.join(ROOT, "rnet-inference", "RibonanzaNet"))
sys.path.insert(0, os.path.join(ROOT, "geometric-rna-design", "src"))

os.environ.setdefault("NUPACKHOME", "/tmp")  # dummy arnie config

import numpy as np
import pandas as pd
import torch
from torch import nn

from utils import FinetunableRibonanzaNet, RNA_Dataset, load_config_from_yaml
from Network import RibonanzaNet
from openknot_score import get_openknot_score
from arnie.pk_predictors import _hungarian
from arnie.utils import convert_dotbracket_to_bp_list

USE_GPU = torch.cuda.is_available()
DEVICE = "cuda" if USE_GPU else "cpu"
RNET_DIR = os.path.join(ROOT, "rnet-inference")
BENCH_CSV = os.path.join(ROOT, "OpenKnotAIDesignData", "Data", "OpenKnotBench_data.v4.5.1.csv")


# ---------------------------------------------------------------- models
def load_shape_model():
    cfg = load_config_from_yaml(os.path.join(RNET_DIR, "RibonanzaNet", "configs", "pairwise.yaml"))
    model = RibonanzaNet(cfg)
    model.load_state_dict(torch.load(os.path.join(RNET_DIR, "RibonanzaNet-Weights", "RibonanzaNet.pt"), map_location="cpu"))
    return model.to(DEVICE).eval()


class FinetunedSS(FinetunableRibonanzaNet):
    """RibonanzaNet-SS: 二级结构(配对矩阵)预测头, 与 rnet-inference/src/rnet_2d.py 一致"""
    def __init__(self, config):
        config.dropout = 0.3
        super().__init__(config)
        self.dropout = nn.Dropout(0.0)
        self.ct_predictor = nn.Linear(64, 1)

    def forward(self, src):
        _, pairwise = super().forward(src, torch.ones_like(src).long().to(src.device))
        pairwise = pairwise + pairwise.permute(0, 2, 1, 3)
        return self.ct_predictor(self.dropout(pairwise)).squeeze(-1)


def load_ss_model():
    cfg = load_config_from_yaml(os.path.join(RNET_DIR, "RibonanzaNet", "configs", "pairwise.yaml"))
    model = FinetunedSS(cfg)
    model.load_state_dict(torch.load(os.path.join(RNET_DIR, "RibonanzaNet-Weights", "RibonanzaNet-SS.pt"), map_location="cpu"))
    return model.to(DEVICE).eval()


@torch.no_grad()
def predict_shape(model, seq):
    tokens = RNA_Dataset(pd.DataFrame([{"sequence": seq}]))[0]["sequence"].unsqueeze(0).to(DEVICE)
    pred = model(tokens, torch.ones_like(tokens)).squeeze().cpu().numpy()
    return pred[:, 0]  # 2A3 (SHAPE) 通道


@torch.no_grad()
def predict_structure(model, seq):
    tokens = RNA_Dataset(pd.DataFrame([{"sequence": seq}]))[0]["sequence"].unsqueeze(0).to(DEVICE)
    mat = model(tokens).sigmoid().cpu().numpy()[0]
    n = len(mat)
    for i in range(n):
        for j in range(n):
            if abs(i - j) < 4:
                mat[i][j] = 0.0
    struct, bps = _hungarian(mat, theta=0.5, min_len_helix=1)
    return struct, bps


def bp_f1(true_struct, pred_struct):
    t = {tuple(sorted(bp)) for bp in convert_dotbracket_to_bp_list(true_struct, allow_pseudoknots=True)}
    p = {tuple(sorted(bp)) for bp in convert_dotbracket_to_bp_list(pred_struct, allow_pseudoknots=True)}
    if not p or not t:
        return 0.0
    tp = len(t & p)
    return 2 * tp / (len(p) + len(t))


# ---------------------------------------------------------------- data
def load_p20_top_designs():
    meta_cols = ["id", "method", "target_openknot_score", "sub_start", "sub_end",
                 "design_length", "design_sequence", "target_structure", "RNet_F1"]
    react_cols = [f"reactivity_{i:04d}" for i in range(1, 308)]
    df = pd.read_csv(BENCH_CSV, usecols=meta_cols + react_cols)
    df = df[(df["round"] == 3) & (df["puzzle"] == "P20")] if "round" in df else df
    # round/puzzle 不在 usecols 里, 重新读一次元数据列做过滤
    df2 = pd.read_csv(BENCH_CSV, usecols=["id", "round", "puzzle"])
    keep = df2[(df2["round"] == 3) & (df2["puzzle"] == "P20")].index
    df = df.loc[keep]
    rows = []
    for method in ["MPNN-fixbb", "gRNAde", "Struct2SeQ-SHAPE", "Eterna", "Starting sequence"]:
        sub = df[df["method"] == method]
        if len(sub) == 0:
            continue
        rows.append(sub.sort_values("target_openknot_score", ascending=False).iloc[0])
    return rows, react_cols


def main():
    print(f"加载模型 (device={DEVICE}) ...")
    shape_model = load_shape_model()
    ss_model = load_ss_model()

    print("读取 P20 (Kissing Multiloops) 各方法最优设计 ...")
    rows, react_cols = load_p20_top_designs()

    results = []
    profiles = {}
    for r in rows:
        seq = r["design_sequence"]
        target = r["target_structure"]
        s0, s1 = int(r["sub_start"]), int(r["sub_end"])
        exp_react = np.array([r[c] for c in react_cols[s0 - 1:s1]], dtype=float)

        pred_react = predict_shape(shape_model, seq)
        _, _, _, ok_score = get_openknot_score(target, pred_react[None, :])
        pred_struct, _ = predict_structure(ss_model, seq)
        f1 = bp_f1(target, pred_struct)

        results.append({
            "method": r["method"],
            "实验OpenKnot分": round(float(r["target_openknot_score"]), 2),
            "RNet预测OpenKnot分": round(float(ok_score[0]), 2),
            "RNet-SS结构F1": round(f1, 3),
            "论文数据RNet_F1": round(float(r["RNet_F1"]), 3),
        })
        profiles[r["method"]] = (exp_react, pred_react, target)

    print("\n===== P20 各方法最优设计: 实验 vs 本地 RNet 复算 =====")
    print(pd.DataFrame(results).to_string(index=False))

    # 图: MPNN-fixbb 最优设计的 实验SHAPE vs RNet预测SHAPE
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    method = "MPNN-fixbb" if "MPNN-fixbb" in profiles else list(profiles)[0]
    exp_react, pred_react, target = profiles[method]
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(exp_react, "o-", ms=3, lw=1, label="Experimental SHAPE (2A3)")
    ax.plot(pred_react, "s-", ms=3, lw=1, alpha=0.8, label="RNet predicted SHAPE")
    paired = np.array([c != "." for c in target])
    ax.fill_between(np.arange(len(target)), 0, max(np.nanmax(exp_react), pred_react.max()),
                    where=paired, color="gray", alpha=0.12, label="Paired in target")
    ax.set_xlabel("Nucleotide position")
    ax.set_ylabel("SHAPE reactivity")
    ax.set_title(f"OpenKnot P20 'Kissing Multiloops' top {method} design — RNet reproduction")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "P20_SHAPE_comparison.png")
    fig.savefig(out, dpi=150)
    print(f"\n图已保存: {out}")


if __name__ == "__main__":
    main()
