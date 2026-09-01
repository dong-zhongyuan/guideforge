"""贝叶斯优化框架(策划案 V3 §4.2 两轮 DBTL, 2026-09-02)。

状态机(结构就绪, 数据未至——不虚跑):
  --cold-start  以现有计算分对候选库拟合 GP(sklearn GaussianProcessRegressor
                官方实现), 给出"计算代理上的期望改进"排序 = 第二轮加密的
                探索建议(明确声明: 代理-在-代理, 非活性证据);
  --ingest <csv>  读第一轮 IVT 矩阵(desc,Vmax,EC50[,stringency]) 回填观测,
                GP 重拟合并按 EI 提出下一批 --batch N 条(每取向一条);
  状态存 data/bayesopt_state.json(观测累积, 幂等)。

特征向量(与主管线同源): ddG_dr/bp_dist/spacer_up/contact_sum/cons3_frac/seed_up。
运行: python scripts/crrna_bayesopt.py --cold-start
      python scripts/crrna_bayesopt.py --ingest data/ivt_round1.csv --batch 4
"""
import argparse
import csv
import json
import os

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

FEATS = ["ddG_dr", "bp_dist", "spacer_unpaired", "contact_sum",
         "cons3_frac", "seed_up"]
STATE = os.path.join(DATA, "bayesopt_state.json")


def load_library():
    rows = []
    for name in ("tp53_r248q_zengdr.v6", "tp53_r248q_pdbdr.v2"):
        p = os.path.join(DATA, name + ".variants.csv")
        if not os.path.isfile(p):
            continue
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["passed"] in ("True", "true", "1") and r["desc"] != "WT":
                    rows.append({"src": name, **{k: r.get(k) for k in
                                ["desc"] + FEATS + ["score"]}})
    return rows


def X_of(rows):
    return np.array([[float(r[k]) for k in FEATS] for r in rows])


def fit_gp(X, y):
    k = ConstantKernel(1.0) * RBF(length_scale=np.ones(X.shape[1]))
    gp = GaussianProcessRegressor(kernel=k, normalize_y=True, n_restarts_optimizer=3)
    gp.fit(X, y)
    return gp


def expected_improvement(gp, X, xi=0.01):
    mu, sd = gp.predict(X, return_std=True)
    best = mu.max()
    z = (mu - best - xi) / np.maximum(sd, 1e-9)
    from scipy.stats import norm
    return (mu - best) * norm.cdf(z) + sd * norm.pdf(z)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cold-start", action="store_true")
    ap.add_argument("--ingest", default=None, help="IVT 矩阵 CSV(desc,Vmax,EC50)")
    ap.add_argument("--batch", type=int, default=4)
    args = ap.parse_args()

    lib = load_library()
    if not lib:
        raise SystemExit("候选库不可读")
    X = X_of(lib)
    obs = {}
    if os.path.isfile(STATE):
        obs = json.load(open(STATE, encoding="utf-8")).get("observations", {})

    if args.ingest:
        with open(args.ingest, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                obs[r["desc"]] = {"Vmax": float(r["Vmax"]), "EC50": float(r["EC50"])}
        print("回填观测 %d 条(累计 %d)" % (len(obs) - len(json.load(
            open(STATE))["observations"]) if os.path.isfile(STATE) else len(obs),
            len(obs)))
        json.dump({"observations": obs}, open(STATE, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)

    if args.ingest and obs:
        # 已有实测: 目标 = Vmax/EC50 比值(活性代理, 由实测而来)
        obs_desc = list(obs)
        idx = [i for i, r in enumerate(lib) if r["desc"] in obs]
        y = np.array([obs[lib[i]["desc"]]["Vmax"] / max(obs[lib[i]["desc"]]["EC50"], 1e-6)
                      for i in idx])
        gp = fit_gp(X[idx], y)
        ei = expected_improvement(gp, X)
        mode = "实测回填(IVT)"
    elif args.cold_start:
        y = np.array([float(r["score"]) for r in lib])
        gp = fit_gp(X, y)
        ei = expected_improvement(gp, X)
        mode = "冷启动(计算代理-在-代理, 非活性证据)"
    else:
        raise SystemExit("需要 --cold-start 或 --ingest")

    order = np.argsort(-ei)
    seen, prop = set(), []
    for i in order:
        d = lib[i]["desc"]
        if d in seen or d in obs:
            continue
        seen.add(d)
        prop.append({"desc": d, "ei": round(float(ei[i]), 4),
                     "src": lib[i]["src"],
                     **{k: lib[i][k] for k in FEATS}})
        if len(prop) >= args.batch:
            break
    print("[%s] 提出下一批 %d 条(EI 降序):" % (mode, len(prop)))
    for p in prop:
        print("  %-16s EI=%.4f  ddG_dr=%s" % (p["desc"], p["ei"], p["ddG_dr"]))
    json.dump({"mode": mode, "proposals": prop, "n_observations": len(obs),
               "n_library": len(lib), "features": FEATS},
              open(os.path.join(DATA, "bayesopt_proposals.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> data/bayesopt_proposals.json")


if __name__ == "__main__":
    main()
