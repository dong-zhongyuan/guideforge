# -*- coding: utf-8 -*-
"""§I: Cas12a2 同源迁移证据链补强(2026-09-11, 判据 docs/preregistration.md §I,
登记 dc9bfa4 先于运行)。

I1 适用域分析: 选型器 5 特征全部由 Cas12a 文献数据标定, 候选为 Cas12a2 DR
变体。若候选落在训练特征分布包络内, 则迁移=内插而非外推(fig1f 式超阈按构造
不适用)。包络 = 训练集标准化空间内 5-NN 均距, 阈 = 训练集自身留一 5-NN 均距
的 95 分位。

I2 折叠系综同构: 11 条同源 DR(2 Cas12a2 + 9 Cas12a)的 ViennaRNA 配分系综
逐残基配对概率, 按 dr_conservation.json 18 列星形比对对齐, 两两 Spearman。
主口径 = Cas12a2×Cas12a 18 对中位; 对照 = Cas12a 内部 36 对中位。

判读规则(§I1a/§I2a, 先于数值)见 docs/preregistration.md §I。

运行: python scripts/crrna_transfer_domain.py
输出: data/transfer_domain.json
"""
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from crrna_train_selector import (FEATURES, load_our_candidates,  # noqa: E402
                                  load_pooled, rank_family)

HOM = os.path.join(ROOT, "data", "dr_homologs.json")
CONS = os.path.join(ROOT, "data", "dr_conservation.json")
OUT = os.path.join(ROOT, "data", "transfer_domain.json")
KNN_K = 5
ENVELOPE_Q = 95.0


def i1_domain_of_applicability():
    Xp, _yp, _meta = load_pooled()
    sc = StandardScaler().fit(Xp)
    Xs = sc.transform(Xp)
    # 训练集自身留一 5-NN 均距
    loo = []
    for i in range(len(Xs)):
        d = np.linalg.norm(Xs - Xs[i], axis=1)
        d = np.sort(np.delete(d, i))
        loo.append(float(d[:KNN_K].mean()))
    thr = float(np.percentile(loo, ENVELOPE_Q))

    cands = load_our_candidates()
    fam = rank_family({})  # models 置空: 只取特征, 不做预测
    fam_rows = [{"desc": "FAM:" + r["scaffold"],
                 **{f: r["features"][f] for f in FEATURES}} for r in fam]
    all_rows = cands + fam_rows
    Xc = np.array([[float(r[f]) for f in FEATURES] for r in all_rows])
    Xc_s = sc.transform(Xc)
    inside, outside = 0, []
    for r, xc in zip(all_rows, Xc_s):
        d = np.sort(np.linalg.norm(Xs - xc, axis=1))
        knn = float(d[:KNN_K].mean())
        if knn <= thr:
            inside += 1
        else:
            outside.append({"desc": r["desc"], "knn_dist": round(knn, 3)})
    # 逐特征 min-max 越界
    per_feat = {}
    for j, f in enumerate(FEATURES):
        lo, hi = float(Xp[:, j].min()), float(Xp[:, j].max())
        n_out = int(np.sum((Xc[:, j] < lo) | (Xc[:, j] > hi)))
        per_feat[f] = {"train_min": round(lo, 4), "train_max": round(hi, 4),
                       "n_candidates_outside": n_out}
    frac = inside / len(all_rows)
    fam_inside = 8 - sum(1 for o in outside if o["desc"].startswith("FAM:"))
    fam_out_names = [o["desc"][4:] for o in outside
                     if o["desc"].startswith("FAM:")]
    if frac >= 0.95:
        reading = ("%.1f%%(%d/%d) 候选在训练包络内 >=95%% -> 「Cas12a2 迁移在"
                   "特征空间为内插」获支持(§I1a)" % (frac * 100, inside, len(all_rows)))
    else:
        reading = ("%.1f%%(%d/%d) 候选在训练包络内 <95%% -> 如实报告越界子集, "
                   "迁移表述追加适用域边界声明(§I1a); 细分: 8 员族(答辩核心候选)"
                   "%d/8 在包络内, 越界者: %s" % (frac * 100, inside, len(all_rows),
                                                fam_inside,
                                                "/".join(fam_out_names) or "无"))
    return {"design": "训练=池化 46x5(load_pooled); 候选=v6 全 pass 变体+8 员族"
                      "(同一特征代码路径); 包络=标准化空间 %d-NN 均距, 阈=训练留一"
                      " %d-NN 均距的 %g 分位" % (KNN_K, KNN_K, ENVELOPE_Q),
            "n_train": int(len(Xp)), "n_candidates": len(all_rows),
            "n_candidates_pass_v6": len(cands), "n_family": len(fam_rows),
            "envelope_threshold": round(thr, 4),
            "train_loo_knn_median": round(float(np.median(loo)), 4),
            "n_inside": inside, "frac_inside": round(frac, 4),
            "family_inside_count": fam_inside,
            "family_outside": fam_out_names,
            "per_feature_minmax": per_feat,
            "outside_candidates": outside,
            "reading": reading}


def _paired_prob_per_column(name, aln_ent, n_cols):
    """折叠 entry 全长 DR(pf 系综), 返回 18 比对列 -> 残基配对概率(gap=None)。"""
    import RNA
    full_dna = aln_ent["full"]
    seq = full_dna.replace("T", "U")
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    pp = [float(sum(bpp[i][j] for j in range(1, n + 1))) for i in range(1, n + 1)]
    row = aln_ent["aligned_core"]
    start = aln_ent["core_start_in_full"]
    out, ci = [], 0
    for ch in row:
        if ch == "-":
            out.append(None)
        else:
            out.append(pp[start + ci])
            ci += 1
    assert len(out) == n_cols, (name, len(out), n_cols)
    return out


def i2_fold_congruence():
    hom = json.load(open(HOM, encoding="utf-8"))
    cons = json.load(open(CONS, encoding="utf-8"))
    n_cols = len(cons["conservation"])
    aln = {a["name"]: a for a in cons["alignment"]}
    profs = {}
    for e in hom["entries"]:
        name = e["name"]
        if name not in aln:
            continue
        profs[name] = _paired_prob_per_column(name, aln[name], n_cols)
    a2 = sorted(n for n in profs if "Cas12a2" in n)
    a1 = sorted(n for n in profs if "Cas12a2" not in n)

    def pair_rho(x, y):
        xs = [(a, b) for a, b in zip(profs[x], profs[y])
              if a is not None and b is not None]
        if len(xs) < 6:
            return None
        u = [p[0] for p in xs]
        v = [p[1] for p in xs]
        return float(spearmanr(u, v).statistic)

    cross, internal = [], []
    for x in a2:
        for y in a1:
            r = pair_rho(x, y)
            if r is not None:
                cross.append({"pair": [x, y], "rho": round(r, 3)})
    for i in range(len(a1)):
        for j in range(i + 1, len(a1)):
            r = pair_rho(a1[i], a1[j])
            if r is not None:
                internal.append(round(r, 3))
    med_cross = float(np.median([c["rho"] for c in cross]))
    med_internal = float(np.median(internal))
    if med_cross >= 0.7:
        reading = ("Cas12a2×Cas12a 配对谱中位 rho=%.3f>=0.7 -> 「两体系 DR 折叠"
                   "系综同构」获支持(§I2a)" % med_cross)
    elif med_cross >= 0.4 or med_cross >= med_internal - 0.1:
        reading = ("Cas12a2×Cas12a 中位 rho=%.3f(Cas12a 内部对照中位 %.3f) -> "
                   "不低于同源家族内部离散度, 部分支持(§I2a)" % (med_cross,
                                                                  med_internal))
    else:
        reading = ("Cas12a2×Cas12a 中位 rho=%.3f, 低于对照离散度 -> 同构不成立, "
                   "如实报告并触发迁移表述修订(§I2a)" % med_cross)
    return {"design": "11 同源 DR ViennaRNA pf 逐残基配对概率 -> 18 比对列谱; "
                      "两两 tie-aware Spearman; 主口径 Cas12a2xCas12a, "
                      "对照 Cas12a 内部",
            "cas12a2_entries": a2, "cas12a_entries": a1,
            "n_columns": n_cols,
            "cross_pairs": cross, "cross_median_rho": round(med_cross, 3),
            "internal_cas12a_rhos": internal,
            "internal_median_rho": round(med_internal, 3),
            "profiles_paired_prob_per_column": profs,
            "reading": reading}


def main():
    print("=== I1 适用域分析 ===")
    i1 = i1_domain_of_applicability()
    print("包络阈(留一5NN 95分位)=%.4f; 候选 %d 条, 包络内 %d (%.1f%%)" % (
        i1["envelope_threshold"], i1["n_candidates"], i1["n_inside"],
        i1["frac_inside"] * 100))
    for f, d in i1["per_feature_minmax"].items():
        print("  %-10s 训练范围 [%8.3f, %8.3f] 越界候选 %d" % (
            f, d["train_min"], d["train_max"], d["n_candidates_outside"]))
    print("判读(§I1a): %s" % i1["reading"])

    print("\n=== I2 折叠系综同构 ===")
    i2 = i2_fold_congruence()
    for c in i2["cross_pairs"]:
        print("  %-24s x %-24s rho=%+.3f" % (c["pair"][0], c["pair"][1], c["rho"]))
    print("主口径中位 %.3f | Cas12a 内部对照中位 %.3f" % (
        i2["cross_median_rho"], i2["internal_median_rho"]))
    print("判读(§I2a): %s" % i2["reading"])

    payload = {
        "generated_by": "scripts/crrna_transfer_domain.py (§I, 2026-09-11)",
        "criterion": "§I1a/§I2a(2026-09-11 登记 dc9bfa4 先于运行)",
        "context": "迁移链现有 L1 Dmytrenko DR 互换(文献实验) / L2 Protenix 界面"
                   "恢复+cas12a2 活性模板无假结(结构) / L3 11 同源 DR 比对(序列); "
                   "本节补 I1 特征适用域与 I2 折叠系综同构两个可计算环节",
        "i1_domain_of_applicability": i1,
        "i2_fold_congruence": i2,
    }
    json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n输出 -> %s" % OUT)


if __name__ == "__main__":
    main()
