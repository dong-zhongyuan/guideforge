#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""同源跨上下文 DR 效应交互检验(判据登记: docs/preregistration.md §E, 2026-09-08)。

回答: 骨架/DR 变体效应的最优解是否随靶标(spacer)上下文变化——即"分型假说"
在 Cas12a 家族同源大库上的原理层统计判定。只回答原理层, 不替代 Cas12a2 本体系
验证(后者维持"功能收益待湿实验"口径)。

数据:
- DeWeirdt 2020 (AsCas12a alt-DR 负筛): data/raw/deweirdt2020_dr_scan.json
  同一批 35,883 条 DR 变体在两个 spacer 上下文(pRDA_127=MCL1 guide /
  pRDA_128=BCL2L1 guide)各有一个 LFC(低=活性强)。判据 E1/E2。
- Tian 2025 (LbCas12a RRS 单点扫描): data/han2025_dataset.json 的
  tian2025_rrs_pairs: 12 个 DR 5' 端单点突变 × 8 条 crRNA(spacer 上下文),
  活性比 rrs_rel(mut/WT)。判据 E3(佐证)。

判定规则(全部为项目设定常数, 先登记后评估, 见 preregistration.md §E):
- E1 PASS(交互显著): rho < 0.4 且 bootstrap CI95 上界 < 0.5;
  E1 FAIL(高度一致, 支持通用型最优): rho >= 0.7;
  0.4 <= rho < 0.7: 不判(中等一致, 如实报告)。
- E2/E3 为佐证, 不参与 PASS/FAIL, 如实分级报告。

判读文字全部由数值按本文件 READ_RULES 显式规则生成, 禁止手写结论。
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402
from scipy.stats import spearmanr, norm  # noqa: E402

DEWEIRDT_RAW = os.path.join(ROOT, "data", "raw", "deweirdt2020_dr_scan.json")
HAN_JSON = os.path.join(ROOT, "data", "han2025_dataset.json")
OUT_JSON = os.path.join(ROOT, "data", "homolog_context_interaction.json")

RHO_PASS = 0.4   # E1: rho < 0.4 且 CI95 上界 < 0.5 -> PASS(交互显著)
CI_PASS_UB = 0.5
RHO_FAIL = 0.7   # E1: rho >= 0.7 -> FAIL(高度一致, 支持通用型最优)
W_POS1_INTERACTION = 0.3   # E3b: 位置 1 碱基偏好 W < 0.3 -> 交互佐证(E1 同向)
W_POS1_UNIVERSAL = 0.7     # E3b: 位置 1 碱基偏好 W >= 0.7 -> 通用佐证
N_BOOT = 1000
SEED = 20260908

DATE = "2026-09-08"

# ---- 统计工具 ---------------------------------------------------------------

def _rho(x, y):
    """tie-aware Spearman rho(scipy midranks)。常量输入返回 nan。"""
    with np.errstate(all="ignore"):
        rho, p = spearmanr(x, y)
    return float(rho) if np.isfinite(rho) else None, (float(p) if np.isfinite(p) else None)


def _bootstrap_ci(x, y, n_iter=N_BOOT, seed=SEED, alpha=0.05):
    """同批索引重抽两向量, 返回 rho 分布的 [alpha/2, 1-alpha/2] 分位。"""
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = []
    for _ in range(n_iter):
        idx = rng.integers(0, n, n)
        with np.errstate(all="ignore"):
            r, _ = spearmanr(np.asarray(x)[idx], np.asarray(y)[idx])
        if np.isfinite(r):
            vals.append(float(r))
    if not vals:
        return [None, None]
    return [float(np.quantile(vals, alpha / 2)), float(np.quantile(vals, 1 - alpha / 2))]


def _fisher_z_2group(r1, n1, r2, n2):
    """两独立 Spearman rho 差异的双侧检验(Fisher z 转换)。"""
    if None in (r1, r2) or n1 < 4 or n2 < 4:
        return None
    z1 = np.arctanh(r1)
    z2 = np.arctanh(r2)
    se = np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z = (z1 - z2) / se
    p = 2.0 * (1.0 - norm.cdf(abs(z)))
    return {"z": float(z), "p": float(p), "rho_diff": float(r1 - r2)}


def _kendall_w(matrix):
    """Kendall 协和系数 W: 行=评委, 列=对象。每行内秩 1..ncols。
    返回 (W, n_judges, n_objects, col_rank_sums)。"""
    mat = np.asarray(matrix, dtype=float)
    k, n = mat.shape
    if k < 2 or n < 3:
        return None
    # 每行内 rankdata(处理缺失: 该行全 NaN 的行剔除)
    valid_rows = [r for r in mat if np.all(np.isfinite(r))]
    k = len(valid_rows)
    if k < 2:
        return None
    from scipy.stats import rankdata
    ranks = np.vstack([rankdata(r) for r in valid_rows])
    col_sums = ranks.sum(axis=0)
    mean = k * (n + 1) / 2.0
    S = float(np.sum((col_sums - mean) ** 2))
    W = 12.0 * S / (k * k * (n ** 3 - n))
    return {"W": float(W), "n_judges": k, "n_objects": n,
            "col_rank_sums": [float(v) for v in col_sums]}


# ---- E1/E2: DeWeirdt 跨上下文 ----------------------------------------------

def load_deweirdt():
    with open(DEWEIRDT_RAW, encoding="utf-8") as f:
        d = json.load(f)
    return d


def region_of(dr_rna, wt_rna, paired):
    """按突变落区分类: stem_major / loop_major / mixed。区位=WT MFE 折叠(伪结外)。"""
    mut = [i for i, (a, b) in enumerate(zip(wt_rna, dr_rna)) if a != b]
    if not mut:
        return None
    n_pair = sum(1 for i in mut if i in paired)
    n_loop = len(mut) - n_pair
    if n_pair > n_loop:
        return "stem_major"
    if n_loop > n_pair:
        return "loop_major"
    return "mixed"


def analyze_deweirdt(d):
    rows = d["rows"]
    wt_rna = d["wt_dr_dna"].replace("T", "U")
    ss, mfe = core.fold(wt_rna)
    paired = {i for i, c in enumerate(ss) if c in "()"}

    out = {"n_total": len(rows), "wt_dr_dna": d["wt_dr_dna"],
           "wt_dr_fold": {"dotbracket": ss, "mfe_kcal": round(float(mfe), 3),
                          "n_paired_pos": len(paired),
                          "note": "ViennaRNA MFE 伪结外折叠, 区位分层仅描述用"}}

    # E1 主判据: cls=test
    test = [r for r in rows if r["cls"] == "test"]
    x = np.array([r["lfc127"] for r in test], dtype=float)
    y = np.array([r["lfc128"] for r in test], dtype=float)
    rho, p = _rho(x, y)
    ci = _bootstrap_ci(x, y)
    e1_pass = rho is not None and rho < RHO_PASS and ci[1] is not None and ci[1] < CI_PASS_UB
    e1_fail = rho is not None and rho >= RHO_FAIL
    if e1_pass:
        e1_verdict = "PASS"
        e1_reading = ("rho=%.3f < %.1f 且 bootstrap CI95 上界 %.3f < %.1f: "
                      "同批 DR 变体效应强烈依赖 spacer 上下文(交互显著), "
                      "分型假说获同源跨上下文大样本统计支持" % (rho, RHO_PASS, ci[1], CI_PASS_UB))
    elif e1_fail:
        e1_verdict = "FAIL"
        e1_reading = ("rho=%.3f >= %.1f: 同批 DR 变体效应跨两 spacer 上下文高度一致, "
                      "支持通用型最优, 未检出跨上下文交互" % (rho, RHO_FAIL))
    else:
        e1_verdict = "INCONCLUSIVE"
        e1_reading = ("rho=%.3f 落在 [%.1f, %.1f) 中等一致区间: 不判, 如实报告"
                      % (rho, RHO_PASS, RHO_FAIL))
    out["e1"] = {"n": len(test), "cls": "test", "rho_cross": rho,
                 "p": p, "bootstrap_ci95": ci, "thresholds": {"rho_pass": RHO_PASS,
                 "ci_ub_pass": CI_PASS_UB, "rho_fail": RHO_FAIL},
                 "verdict": e1_verdict, "reading": e1_reading}

    # 附加披露: 按 n_mut 桶分层(不参与判定)
    bins = {"1-3": lambda m: 1 <= m <= 3, "4-6": lambda m: 4 <= m <= 6,
            "7+": lambda m: m >= 7}
    by_nmut = {}
    for name, f in bins.items():
        sel = [r for r in test if f(r["n_mut"])]
        if len(sel) >= 10:
            r_, p_ = _rho([r["lfc127"] for r in sel], [r["lfc128"] for r in sel])
            by_nmut[name] = {"n": len(sel), "rho_cross": r_, "p": p_}
    out["e1_by_nmut"] = by_nmut

    # 附加披露: 对照类跨方向 rho(噪声/设计对照, 不参与判定)
    by_cls = {}
    for cls in ("wildtype", "intend-6T", "random"):
        sel = [r for r in rows if r["cls"] == cls]
        if len(sel) >= 2:
            r_, p_ = _rho([r["lfc127"] for r in sel], [r["lfc128"] for r in sel])
            by_cls[cls] = {"n": len(sel), "rho_cross": r_, "p": p_}
    out["e1_by_cls"] = by_cls

    # E2 佐证: 区位分层
    regions = {"stem_major": [], "loop_major": [], "mixed": []}
    for r in test:
        reg = region_of(r["dr_rna"], wt_rna, paired)
        if reg:
            regions[reg].append(r)
    e2 = {}
    for reg, sel in regions.items():
        if len(sel) >= 10:
            r_, p_ = _rho([r["lfc127"] for r in sel], [r["lfc128"] for r in sel])
            e2[reg] = {"n": len(sel), "rho_cross": r_, "p": p_}
    # stem vs loop Fisher z(两侧)
    if "stem_major" in e2 and "loop_major" in e2:
        a = e2["stem_major"]
        b = e2["loop_major"]
        if a["rho_cross"] is not None and b["rho_cross"] is not None:
            e2["stem_vs_loop_fisher_z"] = _fisher_z_2group(
                a["rho_cross"], a["n"], b["rho_cross"], b["n"])
    e2["note"] = "区位=WT DR MFE 折叠配对位(茎)vs 非配对位(环/单链); 佐证不参与 PASS/FAIL"
    out["e2"] = e2
    return out


# ---- E3: Tian 跨 crRNA 一致性 ----------------------------------------------

def analyze_tian():
    with open(HAN_JSON, encoding="utf-8") as f:
        h = json.load(f)
    pairs = h["tian2025_rrs_pairs"]
    # tag 格式 "<mut>_crRNA<k>"(如 1U_crRNA1); 突变名 = _crRNA 前缀, 跨 crRNA 对齐
    by_crna = {}
    for r in pairs:
        tag = r["tag"]
        mut_name = tag.split("_crRNA")[0]
        crna = r.get("crRNA")
        key = "crRNA%d" % crna if crna else tag.split("_crRNA")[-1]
        by_crna.setdefault(key, {})[mut_name] = r["rrs_rel"]
    crna_names = sorted(by_crna.keys(), key=lambda s: int(s.replace("crRNA", "")))
    mut_tags = sorted({t for m in by_crna.values() for t in m})
    matrix = []
    per_crna_n = []
    for c in crna_names:
        row = [by_crna[c].get(t, np.nan) for t in mut_tags]
        matrix.append(row)
        per_crna_n.append((c, int(np.isfinite(row).sum())))
    # 仅保留全 crRNA 均有效的突变列
    mat = np.array(matrix, dtype=float)
    valid_cols = np.all(np.isfinite(mat), axis=0)
    mat = mat[:, valid_cols]
    kept_tags = [t for t, v in zip(mut_tags, valid_cols) if v]
    # E3a: 全量 Kendall W(位置梯度, 通用层)
    kw = _kendall_w(mat)
    # 两两 crRNA Spearman 中位数
    from itertools import combinations
    pair_rhos = []
    for c1, c2 in combinations(range(mat.shape[0]), 2):
        r_, _ = _rho(mat[c1], mat[c2])
        if r_ is not None:
            pair_rhos.append(r_)

    # E3b: 位置分层——位置 1(耐受位)3 突变 × 8 crRNA 的碱基偏好一致性
    def _pos_subset(pos_char):
        col_idx = [i for i, t in enumerate(kept_tags) if t.startswith(pos_char)]
        if len(col_idx) < 3:
            return None
        sub = mat[:, col_idx]
        sub_tags = [kept_tags[i] for i in col_idx]
        w_sub = _kendall_w(sub)
        best = {}
        for r_i in range(sub.shape[0]):
            row = sub[r_i]
            mx = int(np.nanargmax(row))
            best["crRNA%d" % (r_i + 1)] = {"best": sub_tags[mx],
                                           "best_val": float(row[mx])}
        dist = {}
        for v in best.values():
            dist[v["best"]] = dist.get(v["best"], 0) + 1
        return {"cols": sub_tags, "kendall_w": w_sub,
                "best_base_distribution": dist, "per_crna_best": best}

    pos1 = _pos_subset("1")
    pos2 = _pos_subset("2")

    out = {
        "n_crna": len(crna_names),
        "n_mutants_total": len(mut_tags),
        "n_mutants_complete": int(valid_cols.sum()),
        "kept_mut_tags": kept_tags,
        "per_crna_n_valid": per_crna_n,
        # E3a 全量(位置梯度, 通用层)
        "kendall_w": kw,
        "pairwise_spearman_median": float(np.median(pair_rhos)) if pair_rhos else None,
        # E3b 位置分层(碱基级偏好, 交互层)
        "pos1": pos1,
        "pos2": pos2,
        "note": "Tian 2025 RRS 区(DR 5' 端 4nt)单点突变, 8 条 crRNA(spacer)共享同一批 "
                "12 突变。E3a 全量 W 由位置梯度驱动: 位置 3/4(假结核心)近乎全灭活、"
                "位置 1 最耐受, '哪一位敏感'跨 crRNA 稳定(通用层, 与 Cas12a2 本体系 "
                "§A 公共集 A1C/A1G/A1U 同为位置 1 惰性, 跨体系呼应); E3b 看同一耐受位点"
                "内 '哪个碱基最优'是否随 spacer 上下文翻转(交互层, E1 同向佐证)",
    }
    # E3b 判读: pos1(耐受位)碱基偏好一致性
    if pos1 and pos1["kendall_w"] and pos1["kendall_w"]["W"] is not None:
        w1 = pos1["kendall_w"]["W"]
        if w1 < W_POS1_INTERACTION:
            out["reading"] = ("位置 1 碱基偏好 Kendall W=%.3f < %.1f: 同一耐受位点上"
                              "'哪个碱基最优'随 spacer 上下文翻转(最优碱基分布 %s), "
                              "交互佐证(E1 同向); E3a 全量 W=%.3f 仅反映位置梯度通用层, "
                              "两层不矛盾"
                              % (w1, W_POS1_INTERACTION,
                                 pos1["best_base_distribution"], kw["W"]))
        elif w1 >= W_POS1_UNIVERSAL:
            out["reading"] = ("位置 1 碱基偏好 Kendall W=%.3f >= %.1f: 碱基偏好跨 "
                              "spacer 上下文稳定, 通用佐证" % (w1, W_POS1_UNIVERSAL))
        else:
            out["reading"] = ("位置 1 碱基偏好 Kendall W=%.3f 落在 [%.1f, %.1f) 中等区间: "
                              "不判, 如实报告" % (w1, W_POS1_INTERACTION, W_POS1_UNIVERSAL))
    else:
        out["reading"] = "位置 1 子集不足, 无法判定碱基偏好一致性"
    return out


# ---- 主流程 ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT_JSON)
    args = ap.parse_args()

    t0 = time.time()
    deweirdt = analyze_deweirdt(load_deweirdt())
    tian = analyze_tian()

    e1 = deweirdt["e1"]
    verdict = e1["verdict"]
    if verdict == "PASS":
        overall = "homolog_interaction_supported"
        overall_reading = ("E1 PASS: 分型假说获同源跨上下文大样本统计支持(AsCas12a "
                           "35,883 变体 × 2 spacer 上下文)。与 §A 本体系计算弱证据阳性 "
                           "并列, 分型主张形成「同源原理统计 + 本体系计算实例」双层证据。")
    elif verdict == "FAIL":
        overall = "homolog_universal_consistent"
        overall_reading = ("E1 FAIL: 同源大库(AsCas12a)未检出跨上下文交互, 效应高度一致, "
                           "支持通用型最优。如实进 README; 本体系 §A 计算口径与"
                           "「待 IVT 矩阵」标注不变。")
    else:
        overall = "inconclusive_mid_consistency"
        overall_reading = ("E1 中等一致不判。按判据如实报告, 不升级为交互证据, "
                           "也不判为通用型; 本体系维持既有口径。")

    payload = {
        "generated_by": "scripts/crrna_homolog_context_interaction.py",
        "date": DATE,
        "scope": ("同源跨上下文 DR 效应交互检验(分型假说原理层统计判定); "
                  "判据登记 docs/preregistration.md §E, 登记先于本次评估运行; "
                  "判读文字全部由数值按显式阈值规则生成"),
        "prereg_ref": "docs/preregistration.md §E",
        "runtime_sec": round(time.time() - t0, 2),
        "deweirdt": deweirdt,
        "tian": tian,
        "verdict": {"overall": overall, "e1": e1["verdict"],
                    "reading": overall_reading},
        "interpretation": overall_reading,
        "limitation": ("DeWeirdt 仅 2 个 spacer 上下文(n=2 上下文, 回答'效应是否普遍一致'"
                       "而非'型'的统计); 负筛 LFC 非动力学读数; AsCas12a 假结承载结构域 "
                       "与 Cas12a2(成熟 crRNA 无假结)存在拓扑差异, 迁移引用须按 README "
                       "边界声明; Tian 12 变体集中 RRS 假结区。同源层证据不替代 Cas12a2 "
                       "本体系验证(功能收益维持'待湿实验')。"),
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    print(json.dumps({k: payload[k] for k in
                      ("verdict",)}, ensure_ascii=False, indent=1))
    print("E1 rho=%.4f p=%s CI95=%s" % (deweirdt["e1"]["rho_cross"],
                                         deweirdt["e1"]["p"], deweirdt["e1"]["bootstrap_ci95"]))
    if tian["kendall_w"]:
        print("Tian Kendall W=%.4f pairwise_median=%s"
              % (tian["kendall_w"]["W"], tian["pairwise_spearman_median"]))
    print("written:", args.out)


if __name__ == "__main__":
    main()
