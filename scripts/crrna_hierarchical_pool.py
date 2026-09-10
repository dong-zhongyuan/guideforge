"""G7: 层级部分池化(2026-09-11, 判据 docs/preregistration.md §G7, 登记先于运行)。

背景(负面清单 #6/#8): 决策树池化 LOEO 在 cis 端点反预测(-0.692); DeWeirdt
120 条子样平池入池实测把其余 5 终点 LOEO 全部拖负(单终点占池 72% 淹没共识
语义), 被拒入。层级部分池化(组随机截距 + 共享斜率)是多终点整合的标准统计
解法: 各终点自带截距吸收水平差, 斜率跨终点共享。

模型(无外部依赖, 交替最小二乘, 确定性):
  y_i = a[g(i)] + beta'x_i + eps,  g = 终点组
  1) X 全局标准化(整个池化集, 固定);
  2) 初始化 a=0; 迭代 100 轮: 给定 a 拟合 Ridge(beta) 于残差 y-a;
     给定 beta 更新 a[g] = shrink * mean(y - X beta | g),
     shrink = n_g/(n_g + LAM), LAM=10(显式常数, 部分池化收缩);
  3) LOEO: 留出整组训练, 留出组截距收缩至 0(先验), 报留出组 Spearman。

池化集 = crrna_train_selector.load_pooled() 现有观测(Han 5 终点 + Fig1f + Teng
序数对) + DeWeirdt 120 条 LFC 秩分层子样(端点内全量 z 标准化后按秩等距抽样,
确定性; 与 2026-09-03 拒入实验同规模)。

判读规则(§G7a, 先于数值): DeWeirdt 入池后其余 5 终点(fig1g/cis_end/trans_end/
rbs0/rbs33)LOEO 保持正值的终点数 >=3 -> 池化口径从「拒入」修订为「层级部分
池化入池」; 否则维持拒入并报告层级口径数值。cis 终点层级 LOEO 与现值 -0.692
并列报告。层级模型仅作终点整合第二口径, 与决策树池化版并存。

运行: python scripts/crrna_hierarchical_pool.py
输出: data/hierarchical_pool.json
"""
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from crrna_train_selector import FEATURES, load_pooled  # noqa: E402

SCAN = os.path.join(ROOT, "data", "raw", "deweirdt2020_dr_scan.json")
OUT = os.path.join(ROOT, "data", "hierarchical_pool.json")
LAM = 10.0          # 部分池化收缩常数(显式)
N_DW = 120          # DeWeirdt 子样规模(与 2026-09-03 拒入实验同)
MAIN5 = ["fig1g", "cis_end", "trans_end", "rbs0", "rbs33"]


def deweirdt_subsample():
    """DeWeirdt 端点: 全量按 -LFC127 z 标准化, 按 LFC127 秩等距抽 120 条
    (确定性秩分层)。返回 (X, z, meta)。"""
    scan = json.load(open(SCAN, encoding="utf-8"))
    rows = [r for r in scan["rows"]
            if r.get("lfc127") is not None and r.get("features_127")]
    lfc = np.array([float(r["lfc127"]) for r in rows])
    z = -(lfc - lfc.mean()) / max(lfc.std(ddof=0), 1e-9)  # -LFC 大=活性强=优
    order = np.argsort(lfc, kind="stable")     # 按 LFC 秩
    picks = order[np.linspace(0, len(order) - 1, N_DW).round().astype(int)]
    X, zz, meta = [], [], []
    for i in picks:
        r = rows[int(i)]
        X.append([float(r["features_127"][f]) for f in FEATURES])
        zz.append(float(z[int(i)]))
        meta.append({"dataset": "deweirdt2020", "tag": r["tag"],
                     "endpoint": "deweirdt_lfc127"})
    return np.array(X), np.array(zz), meta


def fit_hierarchical(X, y, groups, lam=LAM, n_iter=100):
    """交替最小二乘: 共享 Ridge 斜率 + 组截距(收缩)。返回 (sc, beta 模型, a)。"""
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    a = np.zeros(len(y))
    gid = {g: np.where(groups == g)[0] for g in sorted(set(groups))}
    model = None
    for _ in range(n_iter):
        model = Ridge(alpha=1.0).fit(Xs, y - a)
        resid = y - model.predict(Xs)
        a_new = a.copy()
        for g, idx in gid.items():
            shrink = len(idx) / (len(idx) + lam)
            a_new[idx] = shrink * float(resid[idx].mean())
        if np.max(np.abs(a_new - a)) < 1e-10:
            a = a_new
            break
        a = a_new
    return sc, model, a, gid


def loeo(X, y, groups):
    out = {}
    for g in sorted(set(groups)):
        te = np.where(groups == g)[0]
        tr = np.where(groups != g)[0]
        if len(te) < 5:
            continue
        sc, model, _a, _gid = fit_hierarchical(X[tr], y[tr], groups[tr])
        pr = model.predict(sc.transform(X[te]))  # 留出组截距收缩至 0
        out[g] = round(float(spearmanr(pr, y[te]).statistic), 3)
    return out


def main():
    Xp, yp, pmeta = load_pooled()
    groups = np.array([m["endpoint"] for m in pmeta])
    Xdw, zdw, mdw = deweirdt_subsample()
    Xall = np.vstack([Xp, Xdw])
    yall = np.concatenate([yp, zdw])
    gall = np.concatenate([groups, np.array([m["endpoint"] for m in mdw])])
    print("池化集: 现有 %d 观测 + DeWeirdt %d = %d; 组: %s" % (
        len(yp), len(zdw), len(yall), sorted(set(gall))))

    # 对照1: 不含 DeWeirdt 的层级 LOEO(隔离「层级化」本身的效应)
    lo_wo = loeo(Xp, yp, groups)
    # 对照2: 含 DeWeirdt 的层级 LOEO
    lo_w = loeo(Xall, yall, gall)

    print("\n终点            层级LOEO(无DW)  层级LOEO(含DW)")
    for g in sorted(set(list(lo_wo) + list(lo_w))):
        print("  %-16s %10s %10s" % (
            g, "%+.3f" % lo_wo[g] if g in lo_wo else "n/a(<%d)" % 5,
            "%+.3f" % lo_w[g] if g in lo_w else "n/a"))
    n_pos = sum(1 for g in MAIN5 if lo_w.get(g, -1) > 0)
    revise = n_pos >= 3
    reading = ("DeWeirdt 入池后其余 5 终点 LOEO 正值数=%d>=3 -> 池化口径从「拒入」"
               "修订为「层级部分池化入池」" % n_pos) if revise else \
              ("DeWeirdt 入池后其余 5 终点 LOEO 正值数=%d<3 -> 维持拒入, 层级口径"
               "数值如实报告" % n_pos)
    cis_now = lo_w.get("cis_end")
    cis_note = ("cis 层级 LOEO=%+.3f vs 决策树池化现值 -0.692: %s" % (
        cis_now, "改善" if cis_now is not None and cis_now > -0.692 else "未改善"))

    # 全量拟合的截距与共享斜率(供 README 报告)
    sc_f, model_f, a_f, gid_f = fit_hierarchical(Xall, yall, gall)
    intercepts = {g: round(float(a_f[idx[0]]), 4) for g, idx in gid_f.items()}
    coefs = {f: round(float(c), 4)
             for f, c in zip(FEATURES, model_f.coef_)}

    print("\n组截距:", intercepts)
    print("共享斜率(标准化):", coefs)
    print("判读(§G7a): %s" % reading)
    print("cis: %s" % cis_note)

    payload = {
        "generated_by": "scripts/crrna_hierarchical_pool.py (G7, 2026-09-11)",
        "criterion": "§G7a(2026-09-11 登记先于运行): 入池后其余 5 终点 LOEO 正值"
                     "数>=3 -> 修订为层级部分池化入池; 否则维持拒入",
        "model": "组随机截距+共享 Ridge 斜率, 交替最小二乘, LAM=%g 收缩常数, "
                 "X 全局标准化; LOEO 留出组截距收缩至 0" % LAM,
        "deweirdt_subsample": "全量按 -LFC127 z 标准化后按 LFC127 秩等距抽 %d 条"
                              "(确定性秩分层, 无种子)" % N_DW,
        "n_observations": len(yall),
        "groups": {g: int((gall == g).sum()) for g in sorted(set(gall))},
        "loeo_hierarchical_without_deweirdt": lo_wo,
        "loeo_hierarchical_with_deweirdt": lo_w,
        "loeo_tree_pooled_reference": {"cis_end": -0.692},
        "n_main5_positive_with_dw": n_pos,
        "intercepts_full_fit": intercepts,
        "shared_slopes_standardized": coefs,
        "reading": reading,
        "cis_note": cis_note,
        "scope": "层级模型仅作终点整合第二口径, 与决策树池化版并存, 分歧并列报告",
    }
    json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % OUT)


if __name__ == "__main__":
    main()
