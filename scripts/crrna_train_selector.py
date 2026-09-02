"""选型分类器训练: Han 2025 数据 + 本仓库候选库(2026-09-02)。

两段式:
  第一段(文献预训练): Han 2025 的 7 条工具箱变体 × 结构特征 × 实测活性
    -> 决策树回归(可解释): "什么 DR 结构特征区间活性最优"
  第二段(应用): 把第一段模型应用到本仓库 zengdr 口径全部通过变体
    -> 按文献模型预测活性排序候选库

诚实边界: n=7 训练样本(工具箱级), 预测为先验排序非活性保证;
  LbCas12a CRISPRi ≠ Cas12a2 杀伤, 迁移假设已在 README 声明;
  模型类型选决策树 = 答辩可解释("AI 学到什么"有话可答)。
运行: python scripts/crrna_train_selector.py
输出: data/selector_model.json + data/selector_ranked_candidates.csv
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
from sklearn.tree import DecisionTreeRegressor, export_text

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

DATA = os.path.join(ROOT, "data")

FEATURES = ["ddG_dr", "bp_dist", "cross_nt", "p_fold", "spacer_up"]


def load_han():
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    rows = d["toolbox_features"]
    X, y0, y33, tags = [], [], [], []
    for r in rows:
        if "rbs0" not in r:
            continue
        X.append([float(r[f]) for f in FEATURES])
        y0.append(float(r["rbs0"]))
        y33.append(float(r["rbs33"]))
        tags.append(r["tag"])
    return np.array(X), np.array(y0), np.array(y33), tags


def load_our_candidates():
    import crrna_scaffold_design as core
    import RNA
    path = os.path.join(DATA, "tp53_r248q_zengdr.v6.variants.csv")
    wt_top = json.load(open(
        os.path.join(DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
    wt_dr = core.to_rna(wt_top["top"][0]["dr_seq"])
    spacer = core.to_rna(wt_top["spacer_fixed_dna"])
    wt_full_ss, _ = core.fold(wt_dr + spacer)
    wt_dr_ss, wt_dr_mfe = core.fold(wt_dr)
    wt_stem = core.stem_pairs_of(wt_dr_ss)
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["passed"] not in ("True", "true", "1"):
                continue
            dr = core.to_rna(r["dr_seq"])
            full = dr + spacer
            ss, mfe = core.fold(full)
            dr_ss, dr_mfe = core.fold(dr)
            cross_nt, _ = core.cross_pairs(full, len(dr))
            p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
            _, sp_up, _ = core.pf_stats(full, len(dr), 7)
            bp_d = RNA.bp_distance(wt_full_ss, ss)
            rows.append({
                "desc": r["desc"], "score": float(r["score"]),
                "ddG_dr": round(dr_mfe - wt_dr_mfe, 2),
                "bp_dist": bp_d, "cross_nt": cross_nt,
                "p_fold": round(p_fold, 5),
                "spacer_up": round(sp_up, 3),
                "construct_dna": r["construct_dna"]})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    args = ap.parse_args()

    # === 第一段: 文献预训练 ===
    X, y0, y33, tags = load_han()
    print("=== 训练集 ===")
    print("%-6s" % "var", " ".join("%10s" % f for f in FEATURES),
          "  RBS0(活性)")
    for i, t in enumerate(tags):
        print("%-6s" % t, " ".join("%10.3f" % v for v in X[i]),
              "  %6.1f" % y0[i])

    # 目标 = 两档均值(低=抑制好)
    y = (y0 + y33) / 2

    # 决策树(小样本: max_leaf=4, min_samples_leaf=1)
    dt = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                               random_state=42)
    dt.fit(X, y)

    print("\n=== 决策树规则(可解释) ===")
    print(export_text(dt, feature_names=FEATURES, decimals=3))

    # 训练集拟合度
    pred_train = dt.predict(X)
    rho_train = np.corrcoef(np.argsort(np.argsort(pred_train)),
                            np.argsort(np.argsort(y)))[0, 1]
    print("训练集 Spearman rho = %.3f" % rho_train)

    # 特征重要性
    print("\n=== 特征重要性 ===")
    for f, imp in sorted(zip(FEATURES, dt.feature_importances_),
                         key=lambda x: -x[1]):
        print("  %-12s %.3f" % (f, imp))

    # === 第二段: 应用到候选库 ===
    cands = load_our_candidates()
    X_ours = np.array([[float(r[f]) for f in FEATURES] for r in cands])
    preds = dt.predict(X_ours)

    ranked = sorted(zip([r["desc"] for r in cands], preds,
                        [float(r["ddG_dr"]) for r in cands],
                        [float(r["score"]) for r in cands]),
                    key=lambda x: x[1])  # 低=好

    print("\n=== 文献模型对本仓库 %d 条通过变体的排序(前 20) ===" % len(cands))
    print("%-16s %10s %10s %10s" % ("变体", "文献预测", "DDR-ΔΔG", "管线得分"))
    for name, p, ddg, score in ranked[:20]:
        print("%-16s %10.1f %10.2f %10.4f" % (name, p, float(ddg), float(score)))

    # 与管线综合分的 Spearman
    our_scores = np.array([float(s) for _, _, _, s in ranked])
    lit_preds = np.array([p for _, p, _, _ in ranked])
    rho = np.corrcoef(np.argsort(np.argsort(lit_preds)),
                      np.argsort(np.argsort(our_scores)))[0, 1]
    print("\n文献模型 vs 管线综合分 Spearman = %+.3f" % rho)

    # 输出
    out = {"model": "DecisionTreeRegressor(max_leaf=4)",
           "training_data": "Han 2025 (PMC12508434) LbCas12a 7 DR toolbox",
           "features": FEATURES,
           "feature_importance": {f: float(imp) for f, imp in
                                  zip(FEATURES, dt.feature_importances_)},
           "tree_rules": export_text(dt, feature_names=FEATURES),
           "train_spearman": round(float(rho_train), 3),
           "applied_to": "tp53_r248q_zengdr.v6 (204 passed)",
           "vs_pipeline_score_spearman": round(float(rho), 3),
           "top10_lit_model": [{"desc": n, "lit_pred": round(p, 1),
                                "ddG_dr": float(d), "score": float(s)}
                               for n, p, d, s in ranked[:10]]}
    json.dump(out, open(os.path.join(DATA, "selector_model.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)

    with open(os.path.join(DATA, "selector_ranked_candidates.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank_lit", "desc", "lit_pred_activity",
                    "ddG_dr", "pipeline_score"] +
                   [f"feat_{x}" for x in FEATURES])
        for i, (n, p, d, s) in enumerate(ranked, 1):
            row_idx = [j for j, r in enumerate(cands) if r["desc"] == n]
            feats = [cands[row_idx[0]][f] for f in FEATURES] if row_idx else []
            w.writerow([i, n, round(p, 1), float(d), float(s)] + feats)

    print("\n输出 -> data/selector_model.json / selector_ranked_candidates.csv")


if __name__ == "__main__":
    main()
