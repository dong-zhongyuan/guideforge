"""G6: RRS 位置感知模型(2026-09-11, 判据 docs/preregistration.md §G6, 登记先于运行)。

背景: 池化选型器对 Tian 2025 RRS 96 条外部验证 Spearman≈0(伪结外结构特征对
RRS/假结区按构造无分辨力), 现行处置 = 按位置规则排除 DR 5' 端 4nt 变体。
先例: CRISPRpred(Liu 2020, Bioinformatics 36:1454)等证明位置特异性特征对
guide 活性有分辨力; Tian 2025 自身结论亦以位置为主导(RRS 3/4 位近乎灭活,
位置 1 最耐受, 跨 crRNA 稳定)——即该区活性由「位置×碱基」承载, 可被位置
感知模型捕获。

模型: 特征 = 突变位置 one-hot(4) ⊗ 突变碱基 one-hot(3) = 12 维 + 5 维结构特征
(ddG_dr/bp_dist/cross_nt/p_fold/spacer_up); Ridge(alpha=1.0), 标准化仅在训练折
拟合。验证 = 留一 crRNA 交叉验证(LOCO, 8 折, 组间独立, 杜绝同源 crRNA 泄漏)。
消融: 仅位置12维 / 仅结构5维 两对照。

判读规则(§G6a, 先于数值): LOCO Spearman >= +0.3 -> RRS 区候选从「按位置规则
排除」升级为「位置感知模型打分(仅 Tian 同族 RRS 区、仅相对排序)」;
<+0.3 -> 维持排除规则, 数值入库。

运行: python scripts/crrna_rrs_position_model.py
输出: data/rrs_position_model.json
"""
import json
import os
import re
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

STRUCT_FEATS = ["ddG_dr", "bp_dist", "cross_nt", "p_fold", "spacer_up"]
OUT = os.path.join(ROOT, "data", "rrs_position_model.json")


def parse_tag(tag):
    m = re.match(r"^(\d)([ACGU])_crRNA(\d+)$", tag)
    if not m:
        raise ValueError("tag 非法: %s" % tag)
    return int(m.group(1)), m.group(2), int(m.group(3))


def pos_onehot(pos, base):
    """位置 one-hot ⊗ 碱基 one-hot: 12 维, 顺序 pos1A..pos1U(该位实际出现的
    3 个替代碱基按字母序)。"""
    alt = {1: "CGU", 2: "CGU", 3: "ACG", 4: "ACG"}
    v = [0.0] * 12
    v[(pos - 1) * 3 + alt[pos].index(base)] = 1.0
    return v


def loco_spearman(X, y, groups):
    """留一 crRNA 交叉验证, 返回 (整体 Spearman, 逐折 rho 列表)。"""
    pred = np.zeros(len(y))
    for g in sorted(set(groups)):
        tr = [i for i in range(len(y)) if groups[i] != g]
        te = [i for i in range(len(y)) if groups[i] == g]
        sc = StandardScaler().fit(X[tr])
        m = Ridge(alpha=1.0).fit(sc.transform(X[tr]), y[tr])
        pred[te] = m.predict(sc.transform(X[te]))
    per_fold = []
    for g in sorted(set(groups)):
        te = [i for i in range(len(y)) if groups[i] == g]
        r = float(spearmanr(pred[te], y[te]).statistic)
        per_fold.append({"heldout_crRNA": int(g),
                         "spearman": None if r != r else round(r, 3)})
    rho = float(spearmanr(pred, y).statistic)
    return rho, per_fold, pred


def main():
    d = json.load(open(os.path.join(ROOT, "data", "han2025_dataset.json"),
                       encoding="utf-8"))
    rrs = [r for r in d.get("tian2025_rrs_pairs") or []
           if r.get("rrs_rel") is not None]
    rows, Xpos, Xst, y, groups = [], [], [], [], []
    for r in rrs:
        pos, base, cno = parse_tag(r["tag"])
        rows.append(r)
        Xpos.append(pos_onehot(pos, base))
        Xst.append([float(r[f]) for f in STRUCT_FEATS])
        y.append(float(r["rrs_rel"]))
        groups.append(cno)
    Xpos, Xst, y = np.array(Xpos), np.array(Xst), np.array(y)
    Xfull = np.hstack([Xpos, Xst])

    res = {}
    for name, X in (("position_only", Xpos), ("structure_only", Xst),
                    ("position_plus_structure", Xfull)):
        rho, folds, _ = loco_spearman(X, y, groups)
        res[name] = {"loco_spearman": round(rho, 3), "per_fold": folds,
                     "n_features": X.shape[1]}
        print("%-26s LOCO Spearman = %+.3f (n_feat=%d)" % (
            name, rho, X.shape[1]))
        for f_ in folds:
            print("    留出 crRNA%d: rho=%s" % (f_["heldout_crRNA"],
                                                f_["spearman"]))

    main_rho = res["position_plus_structure"]["loco_spearman"]
    upgrade = main_rho >= 0.3
    reading = ("LOCO Spearman=%+.3f>=+0.3 -> RRS 区候选从「按位置规则排除」升级为"
               "「位置感知模型打分(仅 Tian 同族 RRS 区、仅相对排序)」" % main_rho) \
        if upgrade else \
        ("LOCO Spearman=%+.3f<+0.3 -> 维持「按位置规则排除」, 数值入库" % main_rho)

    payload = {
        "generated_by": "scripts/crrna_rrs_position_model.py (G6, 2026-09-11)",
        "criterion": "§G6a(2026-09-11 登记先于运行): LOCO Spearman>=+0.3 -> 升级; "
                     "否则维持排除",
        "precedent": "CRISPRpred (Liu 2020, Bioinformatics 36:1454) 位置特异性"
                     "特征; Tian 2025 自身结论以位置为主导(RRS 3/4 位近乎灭活)",
        "data": "han2025_dataset.json:tian2025_rrs_pairs (8 crRNA x 12, n=%d)" % len(y),
        "features": {"position_onehot_x_base": 12, "structure": STRUCT_FEATS},
        "model": "Ridge(alpha=1.0), StandardScaler 仅训练折拟合; LOCO 8 折组间独立",
        "endpoint": "rrs_rel = 同 crRNA 内 mut/WT trans 切割比值",
        "results": res,
        "reading": reading,
        "scope_limit": "升级仅在 Tian 同族 RRS 区(DR 5' 端 4nt)且仅相对排序; "
                       "不外推到其他 DR 区或绝对活性",
    }
    json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n判读(§G6a): %s" % reading)
    print("输出 -> %s" % OUT)


if __name__ == "__main__":
    main()
