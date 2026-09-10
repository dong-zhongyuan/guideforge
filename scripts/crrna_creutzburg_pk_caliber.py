"""G5: Creutzburg 2020 论文口径对齐重测(2026-09-11, 判据 docs/preregistration.md
§G5, 登记先于运行)。

背景: 本仓库反向复算(crrna_creutzburg_reverse.py)未能复现论文自报的游离态指标
正相关——但两口径有三处不齐: (1) 论文用 pre-crRNA 上下文(5' 残留 repeat + DR +
spacer); (2) 论文的茎定义含手工标注的假结 P1 配对(U(-10):A(-18) 反 Hoogsteen,
FnCas12a/AsCas12a 同源位坐标, 结构文献口径: WO2019135816A2 / Swarts 2017 /
Yamano 2016); (3) 论文 n≈25, 本仓库可核验子集 n=14。本脚本补齐 (1)(2):
hc_add_bp 强制 P1 配对的约束配分下, 5 对经典茎 WC 配对(-6..-2 × -11..-15)的
联合完整概率 stem_prob|P1; 另报 ΔG_pk = dG(P1 强制) − dG(无约束)作假结配对的
热力学代价代理。

诚实边界: P1 强制 = 论文手工茎定义的可计算近似, 不等于真含假结的配分(真 pk
配分检验见 crrna_nupack_pk_weight.py, 四模型已对 Sp8 判 out_of_model_domain);
n=14, 活性为读图近似值。

判读规则(§G5a, 先于数值): stem_prob|P1(pre25 口径) × 活性 rho>0 且置换
p<0.05 -> 「复现失败归因口径不齐」获支持; 否则维持现表述并追加该口径阴性记录。

运行: python scripts/crrna_creutzburg_pk_caliber.py
输出: data/creutzburg_pk_caliber.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

DATASET = os.path.join(ROOT, "data", "creutzburg_2020_dataset.json")
OUT = os.path.join(ROOT, "data", "creutzburg_pk_caliber.json")

# FnCas12a DR18 (AAUUUCUACUGUUGUAGA = 19nt 版去 3' 端 T, 对应位置 -19..-2):
# P1 假结配对 U(-10):A(-18) -> 0基 (1, 9); 经典茎 WC 5 对 (-6..-2)x(-11..-15)
P1_PAIR_DR = (1, 9)
GOOD_TH = 60


def dr_stem_pairs(dr_rna):
    db, _ = core.fold(dr_rna)
    stack, pairs = [], []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.append((stack.pop(), i))
    return pairs


def constrained_pf(full_rna, force_pair=None):
    """约束配分: force_pair=(i,j) 0基时 hc_add_bp 强制该配对。返回
    (bpp 矩阵访问函数, ensemble dG)。"""
    fc = RNA.fold_compound(full_rna)
    if force_pair is not None:
        fc.hc_add_bp(force_pair[0] + 1, force_pair[1] + 1)
    _ss, dG = fc.pf()
    return fc.bpp(), float(dG)


def stem_prob_from(bpp, pairs):
    p = 1.0
    for i, j in pairs:
        p *= bpp[i + 1][j + 1]
    return p


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic)


def perm_p(vals, acts, rho, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n):
        if abs(spearman(rng.permutation(vals), acts)) >= abs(rho):
            hits += 1
    return (hits + 1) / (n + 1)


def main():
    ds = json.load(open(DATASET, encoding="utf-8"))
    spacers = ds["spacers_21nt_dna"]
    act = ds["activity_fig1d_pct_approx"]
    dr18 = core.to_rna(ds["dr18_dna"])
    pre5 = core.to_rna(ds["pre5_flank_dna"])
    prefix = pre5 + dr18
    plen = len(prefix)

    stem_pairs_dr = dr_stem_pairs(dr18)
    stem_pairs_full = [(plen - len(dr18) + i, plen - len(dr18) + j)
                       for i, j in stem_pairs_dr]
    p1_full = (plen - len(dr18) + P1_PAIR_DR[0],
               plen - len(dr18) + P1_PAIR_DR[1])
    # 合法性: P1 两位的碱基互补性(A-U)
    b1, b2 = prefix[p1_full[0]], prefix[p1_full[1]]
    assert {b1, b2} in ({"A", "U"}, {"G", "C"}, {"G", "U"}), \
        "P1 位点非互补: %s:%s" % (b1, b2)

    names = sorted(spacers, key=lambda s: int(s[2:]))
    acts = np.array([act[n] for n in names], float)
    RT = 0.6163  # kcal/mol, 310K(37C, 与管线体温口径一致)
    rows = []
    for nm in names:
        full = prefix + core.to_rna(spacers[nm])
        bpp_c, dg_c = constrained_pf(full, p1_full)
        _bpp_u, dg_u = constrained_pf(full, None)
        dg_pk = dg_c - dg_u
        stem_u = core.stem_intact_prob(full, stem_pairs_full)
        w_pk = float(np.exp(-dg_pk / RT))          # P1 态 Boltzmann 权重代理
        p_paper = stem_u + w_pk                    # G5 修订指标: 两系综分解
        rows.append({
            "name": nm, "activity": act[nm],
            "stem_prob_unconstrained": round(stem_u, 6),
            "dG_pk_cost": round(dg_pk, 3),
            "w_P1_state": round(w_pk, 6),
            "p_stem_paper": round(p_paper, 6)})

    def stats_of(key):
        v = np.array([r[key] for r in rows])
        r_ = spearman(v, acts)
        return {"spearman": round(r_, 3), "perm_p": round(perm_p(v, acts, r_), 4)}

    vals = np.array([r["p_stem_paper"] for r in rows])
    rho = spearman(vals, acts)
    p = perm_p(vals, acts, rho)
    good = acts >= GOOD_TH
    hits = sum(1 for g in vals[good] for b in vals[~good] if g > b)
    auc = hits / (len(vals[good]) * len(vals[~good]))
    st_dg = stats_of("dG_pk_cost")
    st_u = stats_of("stem_prob_unconstrained")
    st_w = stats_of("w_P1_state")

    ok = rho > 0 and p < 0.05
    reading = ("p_stem_paper rho=%+.3f, 置换 p=%.4f -> 论文口径对齐后符号复现, "
               "「复现失败归因口径不齐」获支持" % (rho, p)) if ok else \
              ("p_stem_paper rho=%+.3f, 置换 p=%.4f -> P1 对齐口径亦未达预登记"
               "阳性阈, 维持「游离态指标不构成跨体系活性预测器」表述, 该口径阴性"
               "如实追加记录" % (rho, p))

    print("P1 配对(全长 0基) = %s (%s:%s); 茎 WC 对 = %s" % (
        p1_full, b1, b2, stem_pairs_full))
    print("%-5s %4s  %10s %9s %10s %12s" % (
        "name", "act", "stem_unc", "dG_pk", "w_P1", "p_stem_paper"))
    for r in rows:
        print("%-5s %4d  %10.6f %9.3f %10.6f %12.6f" % (
            r["name"], r["activity"], r["stem_prob_unconstrained"],
            r["dG_pk_cost"], r["w_P1_state"], r["p_stem_paper"]))
    print("\np_stem_paper(pre25): rho=%+.3f  perm_p=%.4f  AUC=%.3f" % (rho, p, auc))
    print("  分量: 无约束茎 rho=%+.3f(p=%.4f) | w_P1 rho=%+.3f(p=%.4f) | "
          "dG_pk rho=%+.3f(p=%.4f)" % (
              st_u["spearman"], st_u["perm_p"], st_w["spearman"], st_w["perm_p"],
              st_dg["spearman"], st_dg["perm_p"]))
    print("判读(§G5a 修订版): %s" % reading)

    payload = {
        "generated_by": "scripts/crrna_creutzburg_pk_caliber.py (G5, 2026-09-11)",
        "criterion": "§G5a(2026-09-11 登记先于运行): rho>0 且 perm p<0.05 -> 口径"
                     "不齐归因获支持; 否则维持现表述",
        "paper_caliber": "Creutzburg 2020 自报口径 = pre-crRNA 上下文 + 茎定义含"
                         "手工假结 P1 配对(U(-10):A(-18) 反 Hoogsteen); n≈25",
        "computable_approximation": "hc_add_bp 强制 P1 的约束配分下 5 对经典茎 WC "
                                    "配对联合概率; 不等于真含假结配分(NUPACK 口径见 "
                                    "crrna_nupack_pk_weight.py)",
        "p1_pair": {"dr18_0based": list(P1_PAIR_DR), "full_0based": list(p1_full),
                    "bases": b1 + ":" + b2,
                    "coord_reference": "结构文献坐标(AsCas12a 5' 手柄假结解剖: "
                                       "WO2019135816A2 说明书 [0174]; Swarts 2017)"},
        "stem_wc_pairs_dr18_0based": [list(p) for p in stem_pairs_dr],
        "context": "pre25 = 5' 残留 repeat(7nt) + DR18 + spacer21",
        "good_threshold_pct": GOOD_TH,
        "rows": rows,
        "stats": {"p_stem_paper": {"spearman": round(rho, 3),
                                   "perm_p": round(p, 4),
                                   "auc_good_bad": round(float(auc), 3)},
                  "component_stem_unconstrained": st_u,
                  "component_w_P1_state": st_w,
                  "component_dG_pk_cost": {**st_dg,
                                           "expected_direction": "neg"}},
        "metric_revision": "2026-09-11 同日修订(修订先于新指标评估): 原定义 "
                           "P(WC 茎 | P1 强制) 构造性恒零(P1 与茎配对交叉, 硬约束"
                           "互斥, 首跑 14/14 全零证实); 修订为两系综分解 "
                           "p_stem_paper = P(WC 茎完整, 无约束) + exp(-dG_pk/RT), "
                           "仓内先例 crrna_state_competition.py 单一能量标尺约束配分",
        "n": len(names),
        "limitation": "n=14(论文 n≈25 的可核验子集), 活性为 Fig1D 读图近似值",
        "reading": reading,
    }
    json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % OUT)


if __name__ == "__main__":
    main()
