"""选型分类器训练: Han 2025 数据 + 本仓库候选库(2026-09-02, Fig1g 尺度版)。

两段式:
  第一段(文献预训练): Han 2025 工具箱变体 × 结构特征 × Fig1g 实测活性
    (归一化 GFP 表达, 低=抑制强; 与筛选图同尺度, 优于旧 RBS0/33 端点)
    -> 决策树回归(可解释): "什么 DR 结构特征区间活性最优"
  第二段(应用): 把第一段模型应用到本仓库 zengdr 口径全部通过变体
    -> 按文献模型预测活性排序候选库
  附加(Sanger 耐受集体检): 54 条真实 Sanger 变体特征经模型预测,
    与 204 候选分布对比 —— 检查模型是否把真实耐受变体排到离谱位置。

诚实边界: n=7 干净监督对(Fig1g 尺度); Han 2025 的 145 个编号变体
  (F/S/L/FL)与序列的映射未随论文发表(2026-09-02 查证 MOESM1-7 全部
  补充材料 + 主文), 96 变体序列级训练集不可得; LbCas12a CRISPRi ≠
  Cas12a2 杀伤, 迁移假设已在 README 声明; 决策树 = 答辩可解释。
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
    X, y, tags = [], [], []
    for r in rows:
        if "fig1g" not in r:
            continue
        X.append([float(r[f]) for f in FEATURES])
        y.append(float(r["fig1g"]))
        tags.append(r["tag"])
    return np.array(X), np.array(y), tags, d


ENDPOINTS = {"fig1g": "fig1g", "cis_end": "cis_end", "trans_end": "trans_end",
             "rbs0": "rbs0", "rbs33": "rbs33"}
# 方向统一: True = 数值越大越优(fig1g/rbs 为抑制读数取负号), 供池化模型共用语义
ENDPOINT_HIGHER_BETTER = {"fig1g": False, "cis_end": True, "trans_end": True,
                          "rbs0": False, "rbs33": False}


def load_pooled():
    """同源池化训练集(等湿实验数据的同源先验版, 2026-09-03):

    全部可用同源活性观测, 逐终点 z 标准化并统一"越大越优"方向:
      - Han 2025 工具箱 5 终点 x 7 条(fig1g 抑制 / cis / trans 切割 / RBS0/RBS33 两档表达)
      - Teng 2019 4n96 序数对(WT=-1, 4n96=+1; 跨物种 Fn/Cas12a 口径, 序数级)
    已尝试并放弃的终点(SPR KD / SupFig8-10 原始荧光): 传感器预处理不可校验 /
    增长混杂方向与 Fig1g 相反——拒绝理由记录于 han2025_dataset.json。
    返回 (X, y, meta), meta 每行 = {dataset, tag, endpoint}。
    """
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    X, y, meta = [], [], []
    for ep, key in ENDPOINTS.items():
        vals = [r for r in d["toolbox_features"] if key in r]
        if len(vals) < 5:
            continue
        ys = np.array([float(r[key]) for r in vals])
        z = (ys - ys.mean()) / max(ys.std(ddof=0), 1e-9)
        if not ENDPOINT_HIGHER_BETTER[ep]:
            z = -z
        for r, zi in zip(vals, z):
            X.append([float(r[f]) for f in FEATURES])
            y.append(float(zi))
            meta.append({"dataset": "han2025", "tag": r["tag"],
                         "endpoint": ep})
    # Teng 2019 4n96 序数对(自身数据集内 z=±1; 特征以 canonical spacer 上下文计算)
    try:
        import crrna_scaffold_design as core
        import RNA
        tj = json.load(open(os.path.join(DATA, "teng4n96_validation.json"),
                            encoding="utf-8"))
        wt_top = json.load(open(os.path.join(
            DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
        spacer = core.to_rna(wt_top["spacer_fixed_dna"])
        # Teng 为 18nt PDB 口径 DR; 以其自身 WT 为基准算特征
        drs = {"Teng_WT": core.to_rna(tj["provenance"]["wt_dr_dna"]),
               "Teng_4n96": core.to_rna(tj["provenance"]["v4n96_dr_dna"])}
        base = drs["Teng_WT"]
        base_ss, _ = core.fold(base + spacer)
        _, base_mfe = core.fold(base)
        base_stem = core.stem_pairs_of(core.fold(base)[0])
        for tag, dr in drs.items():
            full = dr + spacer
            ss, _ = core.fold(full)
            _, dr_mfe = core.fold(dr)
            cross_nt, _ = core.cross_pairs(full, len(dr))
            pf = core.stem_intact_prob(full, base_stem) if base_stem else 1.0
            _, sp_up, _ = core.pf_stats(full, len(dr), 7)
            X.append([round(dr_mfe - base_mfe, 2),
                      RNA.bp_distance(base_ss, ss), cross_nt,
                      round(pf, 5), round(sp_up, 3)])
            y.append(-1.0 if tag == "Teng_WT" else 1.0)
            meta.append({"dataset": "teng2019", "tag": tag,
                         "endpoint": "ordinal_4n96"})
    except Exception as e:
        print("Teng 序数对跳过: %s" % e)
    return np.array(X), np.array(y), meta


def tolerance_neighborhood(row_feats):
    """54 条真实 Sanger 耐受变体的特征分布邻域度(半监督门控)。

    对候选的每个特征取其在耐受集分布中的百分位 p, 邻域度
    = 1 - mean(|p-50|)/50; 1 = 每维都落在真实耐受变体分布中心。
    """
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    tol = d.get("sanger_tolerance") or []
    if len(tol) < 10:
        return None
    devs = []
    for f, v in zip(FEATURES, row_feats):
        col = np.array([float(r[f]) for r in tol])
        p = float((col < v).mean() * 100.0)
        devs.append(abs(p - 50.0) / 50.0)
    return round(1.0 - float(np.mean(devs)), 3)


def load_han_endpoints():
    """同源多终点训练集: fig1g(抑制) + cis_end + trans_end(切割, 与旁切杀伤最近缘)"""
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    out = {}
    for ep, key in ENDPOINTS.items():
        X, y, tags = [], [], []
        for r in d["toolbox_features"]:
            if key not in r:
                continue
            X.append([float(r[f]) for f in FEATURES])
            y.append(float(r[key]))
            tags.append(r["tag"])
        if len(y) >= 5:
            out[ep] = (np.array(X), np.array(y), tags)
    return out


def rank_family(models):
    """对 8 员跨型候选族(orientation_library, 含 compensatory)逐终点预测排序"""
    import crrna_scaffold_design as core
    import RNA
    lib = json.load(open(os.path.join(DATA, "orientation_library.json"),
                         encoding="utf-8"))
    spacer = core.to_rna(lib["spacer"])
    dr_wt = core.to_rna(lib["orientations"][0]["wt_construct"][:19]) if \
        lib["orientations"][0].get("wt_construct") else None
    wt_top = json.load(open(os.path.join(
        DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
    wt_dr = core.to_rna(wt_top["top"][0]["dr_seq"])
    wt_full_ss, _ = core.fold(wt_dr + spacer)
    _, wt_dr_mfe = core.fold(wt_dr)
    wt_stem = core.stem_pairs_of(core.fold(wt_dr)[0])
    fam = [{"desc": "WT", "dr": wt_dr}]
    seen = {"WT"}
    for p in lib["orientations"]:
        b = p["best"]
        if b["desc"] in seen:
            continue
        seen.add(b["desc"])
        fam.append({"desc": b["desc"], "dr": core.to_rna(b["dr_seq"])})
    rows = []
    for m in fam:
        dr = m["dr"]
        full = dr + spacer
        ss, _ = core.fold(full)
        _, dr_mfe = core.fold(dr)
        cross_nt, _ = core.cross_pairs(full, len(dr))
        p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
        _, sp_up, _ = core.pf_stats(full, len(dr), 7)
        feats = [round(dr_mfe - wt_dr_mfe, 2),
                 RNA.bp_distance(wt_full_ss, ss), cross_nt,
                 round(p_fold, 5), round(sp_up, 3)]
        row = {"scaffold": m["desc"],
               "features": dict(zip(FEATURES, feats))}
        for ep, dt in models.items():
            row["pred_" + ep] = round(float(dt.predict([feats])[0]), 4)
        rows.append(row)
    return rows


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

    # === 第一段: 文献预训练(Fig1g 尺度) ===
    X, y, tags, han = load_han()
    print("=== 训练集(Fig1g 归一化 GFP, 低=抑制强) ===")
    print("%-6s" % "var", " ".join("%10s" % f for f in FEATURES),
          "  Fig1g")
    for i, t in enumerate(tags):
        print("%-6s" % t, " ".join("%10.3f" % v for v in X[i]),
              "  %6.4f" % y[i])

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

    # === Sanger 耐受集体检: 真实变体在模型预测分布中的位置 ===
    tol_stats = None
    if han.get("sanger_tolerance"):
        tol = [r for r in han["sanger_tolerance"]
               if r["group"] in ("flank", "loop")]  # 与 Cas12a2 DR 口径近缘
        X_tol = np.array([[float(r[f]) for f in FEATURES] for r in tol])
        tol_pred = dt.predict(X_tol)
        cand_pred = dt.predict(X_ours)
        pct = float(np.mean([(cand_pred < p).mean() for p in tol_pred]))
        print("\n=== Sanger 耐受集体检(n=%d, flank+loop) ===" % len(tol))
        print("真实耐受变体预测活性中位数 = %.4f (工具箱锚点范围 %.4f-%.4f)"
              % (np.median(tol_pred), y.min(), y.max()))
        print("耐受变体在 204 候选预测分布中的平均分位 = %.1f%%" % (pct * 100))
        tol_stats = {
            "n": len(tol),
            "pred_median": round(float(np.median(tol_pred)), 4),
            "pred_min": round(float(tol_pred.min()), 4),
            "pred_max": round(float(tol_pred.max()), 4),
            "avg_percentile_in_candidates": round(pct, 3),
            "by_group": {g: round(float(np.median(
                [p for r, p in zip(tol, tol_pred) if r["group"] == g])), 4)
                for g in sorted({r["group"] for r in tol})}}

    # === 同源多终点: cis/trans 切割终点模型 + 8 员家族逐终点排序 ===
    multi = {ep: (X_e, y_e) for ep, (X_e, y_e, _) in
             load_han_endpoints().items()}
    models = {"fig1g": dt}
    print("\n=== 同源多终点模型(cis/trans = LbCas12a 切割, trans 与旁切杀伤最近缘) ===")
    for ep, (X_e, y_e) in multi.items():
        m = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                                  random_state=42)
        m.fit(X_e, y_e)
        models[ep] = m
        print("  %-9s n=%d  训练范围 %.4f-%.4f" % (ep, len(y_e),
                                                  y_e.min(), y_e.max()))
    fam_rows = rank_family(models)
    # 池化模型 + 耐受集邻域度(同源数据全量训练, 等湿实验回流前的先验版)
    Xp, yp, pmeta = load_pooled()
    pooled = DecisionTreeRegressor(max_leaf_nodes=6, min_samples_leaf=2,
                                   random_state=42)
    pooled.fit(Xp, yp)
    print("\n=== 同源池化模型(n=%d 观测 = Han 5 终点 x 7 + Teng 序数对) ===" % len(yp))
    print(export_text(pooled, feature_names=FEATURES, decimals=3))
    # 留一终点交叉验证(方向一致性)
    loeo = {}
    for ep in sorted({m["endpoint"] for m in pmeta}):
        idx_tr = [i for i, m in enumerate(pmeta) if m["endpoint"] != ep]
        idx_te = [i for i, m in enumerate(pmeta) if m["endpoint"] == ep]
        if len(idx_te) < 5:
            continue
        m2 = DecisionTreeRegressor(max_leaf_nodes=6, min_samples_leaf=2,
                                   random_state=42)
        m2.fit(Xp[idx_tr], yp[idx_tr])
        pr = m2.predict(Xp[idx_te])
        rk = np.argsort(np.argsort(pr))
        ry = np.argsort(np.argsort(yp[idx_te]))
        loeo[ep] = round(float(np.corrcoef(rk, ry)[0, 1]), 3)
        print("  LOEO %-9s (n=%d): Spearman=%+.3f" % (ep, len(idx_te), loeo[ep]))
    for r in fam_rows:
        fv = [r["features"][f] for f in FEATURES]
        r["pred_pooled_z"] = round(float(pooled.predict([fv])[0]), 3)
        r["tolerance_neighborhood"] = tolerance_neighborhood(fv)
    print("\n=== 8 员族: 池化模型 z 预测 + 真实耐受集邻域度 ===")
    print("%-22s %10s %12s" % ("scaffold", "pooled_z", "tol_neighbor"))
    for r in sorted(fam_rows, key=lambda x: -x["pred_pooled_z"]):
        print("%-22s %10.3f %12s" % (r["scaffold"], r["pred_pooled_z"],
                                     r["tolerance_neighborhood"]))
    print("\n=== 8 员跨型候选族逐终点预测(fig1g/cis 低=抑制强; trans 高=切割强) ===")
    print("%-22s %10s %10s %10s" % ("scaffold", "fig1g", "cis_end", "trans_end"))
    for r in fam_rows:
        print("%-22s %10.4f %10.3f %10.3f" % (
            r["scaffold"], r.get("pred_fig1g", float("nan")),
            r.get("pred_cis_end", float("nan")),
            r.get("pred_trans_end", float("nan"))))

    # 输出
    out = {"model": "DecisionTreeRegressor(max_leaf=4)",
           "training_data": "Han 2025 (s41467-025-64010-z) LbCas12a "
                            "toolbox, Fig1g-scale activity (n=%d)" % len(y),
           "features": FEATURES,
           "feature_importance": {f: float(imp) for f, imp in
                                  zip(FEATURES, dt.feature_importances_)},
           "tree_rules": export_text(dt, feature_names=FEATURES),
           "train_spearman": round(float(rho_train), 3),
           "applied_to": "tp53_r248q_zengdr.v6 (%d passed)" % len(cands),
           "vs_pipeline_score_spearman": round(float(rho), 3),
           "sanger_tolerance_check": tol_stats,
           "multi_endpoint": {"endpoints": sorted(models),
                              "note": "cis/trans 终点比 = 同源(LbCas12a)切割先验, "
                                      "trans 与 Cas12a2 旁切杀伤最近缘; 均 n=7"},
           "family_ranking_8": fam_rows,
           "top10_lit_model": [{"desc": n, "lit_pred": round(p, 4),
                                "ddG_dr": float(d), "score": float(s)}
                               for n, p, d, s in ranked[:10]]}
    json.dump(out, open(os.path.join(DATA, "selector_model.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump({"family": fam_rows, "models": sorted(models)},
              open(os.path.join(DATA, "selector_family_ranking.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump({
        "status": "同源先验全量训练(等 IVT 实测替换); 池化目标 = 逐终点 z "
                  "标准化并统一'越大越优'(fig1g/rbs 取负, cis/trans 原向)",
        "n_observations": len(yp),
        "composition": {"han2025_5endpoints_x7": 35, "teng2019_ordinal": 2},
        "pooled_tree_rules": export_text(pooled, feature_names=FEATURES),
        "pooled_feature_importance": {f: float(imp) for f, imp in
                                      zip(FEATURES, pooled.feature_importances_)},
        "loeo_spearman": loeo,
        "rejected_endpoints": {
            "spr_kd": "传感图相位预处理无法对论文定性序(canonical 最强)校验, 拟合不收敛",
            "supfig8_10_raw_fluor": "0mM 诱导下面板间仍 3 倍差、10mM 方向与 "
                                    "Fig1g 相反——增长混杂未归一, 拒绝入库"},
        "family_ranking_pooled": fam_rows},
        open(os.path.join(DATA, "homolog_training.json"), "w",
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
            w.writerow([i, n, round(p, 4), float(d), float(s)] + feats)

    print("\n输出 -> data/selector_model.json / selector_ranked_candidates.csv")


if __name__ == "__main__":
    main()
