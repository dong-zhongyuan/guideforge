"""交叉配对特征按加工阶段拆分 + 组成混杂校正(联审方案 #2, 2026-09-01)。

背景: raw cross_pp 与活性 +0.61(mature)/+0.42(pre) 疑为成分混杂(spacer 5' A
富集 × U 富集 flank/DR 5' 端)。本脚本:
  1. 阶段拆分: cross_pp_flank(仅与被加工丢弃的 5' 残留 repeat) vs
     cross_pp_coreDR(与成熟 DR 核心)——前者属加工模块, 后者才可进成熟活性口径;
  2. 组成匹配零假设: 每条 spacer 做保单碱基组成随机洗牌 N=200, 计算
     z = (obs - E)/SD; (诚实边界: 单碱基组成匹配, 未保二核苷酸组成);
  3. LOO 方向稳定性 + GC 偏相关(rank 残差)。

判据(先定后跑):
  R1 若 |Spearman(z_cross_coreDR)| < 0.3 或 LOO 变号 -> raw 信号归因成分混杂,
     "raw cross_pp 退出主评分" 获数据支持;
  R2 若 cross_pp_flank 与 cross_pp_coreDR 信号方向不一致 -> 阶段拆分必要。
运行: python scripts/crrna_cross_decompose.py
输出: data/cross_decomposition.json
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import numpy as np
import crrna_scaffold_design as core  # noqa: E402

N_SHUFFLE = 200
SEED = 11


def bpp_sum_regions(seq, region_a, region_b):
    """配分函数 bpp 中 region_a × region_b 的配对概率总和(0 基闭开区间列表)。"""
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    ia = set(range(*region_a)) if isinstance(region_a, tuple) else set(region_a)
    ib = set(range(*region_b))
    return float(sum(bpp[i + 1][j + 1] for i in ia for j in ib if i < j))


def spearman(x, y):
    rx = np.argsort(np.argsort(np.asarray(x, float))).astype(float)
    ry = np.argsort(np.argsort(np.asarray(y, float))).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def loo_sign_stability(vals, acts):
    rhos = [spearman(np.delete(vals, k), np.delete(acts, k)) for k in range(len(vals))]
    signs = set(np.sign(rhos))
    return float(np.min(rhos)), float(np.max(rhos)), len(signs) == 1


def partial_spearman_gc(vals, acts, gc):
    """对 GC 做 rank 残差后的偏相关。"""
    def resid(v, g):
        rv = np.argsort(np.argsort(np.asarray(v, float))).astype(float)
        rg = np.argsort(np.argsort(np.asarray(g, float))).astype(float)
        b = np.polyfit(rg, rv, 1)
        return rv - np.polyval(b, rg)
    return float(np.corrcoef(resid(vals, gc), resid(acts, gc))[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dataset", default=os.path.join(ROOT, "data", "creutzburg_2020_dataset.json"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "cross_decomposition.json"))
    args = ap.parse_args()

    ds = json.load(open(args.dataset, encoding="utf-8"))
    spacers = ds["spacers_21nt_dna"]
    act = ds["activity_fig1d_pct_approx"]
    dr18 = core.to_rna(ds["dr18_dna"])
    pre5 = core.to_rna(ds["pre5_flank_dna"])
    rng = np.random.default_rng(SEED)
    names = sorted(spacers, key=lambda s: int(s[2:]))
    acts = np.array([act[n] for n in names], float)
    gcs = np.array([(core.to_rna(spacers[n]).count("G") + core.to_rna(spacers[n]).count("C")) / 21
                    for n in names])

    contexts = {
        "pre25": (pre5 + dr18, (0, 7), (7, 25)),
        "mature18": (dr18, None, (0, 18)),
    }
    report = {"dataset": args.dataset, "n_shuffle": N_SHUFFLE,
              "shuffle_scope": "单碱基组成匹配(二核苷酸未保持, 诚实边界)",
              "contexts": {}}
    for ctx, (prefix, flank_rng, dr_rng) in contexts.items():
        sp_rng = (len(prefix), len(prefix) + 21)
        rows = []
        for nm in names:
            sp = core.to_rna(spacers[nm])
            full = prefix + sp
            f_fl = bpp_sum_regions(full, flank_rng, sp_rng) if flank_rng else 0.0
            f_core = bpp_sum_regions(full, dr_rng, sp_rng)
            null_fl, null_core = [], []
            for _ in range(N_SHUFFLE):
                sh = "".join(rng.permutation(list(sp)))
                sfull = prefix + sh
                null_fl.append(bpp_sum_regions(sfull, flank_rng, sp_rng) if flank_rng else 0.0)
                null_core.append(bpp_sum_regions(sfull, dr_rng, sp_rng))
            z_fl = (f_fl - np.mean(null_fl)) / (np.std(null_fl) or 1)
            z_core = (f_core - np.mean(null_core)) / (np.std(null_core) or 1)
            rows.append({"name": nm, "act": act[nm],
                         "cross_pp_flank": round(f_fl, 4), "z_flank": round(float(z_fl), 3),
                         "cross_pp_coreDR": round(f_core, 4), "z_coreDR": round(float(z_core), 3)})
        print("\n== %s ==" % ctx)
        feats = {}
        for feat, key in (("raw_flank", "cross_pp_flank"), ("z_flank", "z_flank"),
                          ("raw_coreDR", "cross_pp_coreDR"), ("z_coreDR", "z_coreDR")):
            vals = np.array([r[key] for r in rows], float)
            rho = spearman(vals, acts)
            lo, hi, stable = loo_sign_stability(vals, acts)
            p_gc = partial_spearman_gc(vals, acts, gcs)
            feats[feat] = {"spearman": round(rho, 3), "loo_min": round(lo, 3),
                           "loo_max": round(hi, 3), "loo_stable": stable,
                           "partial_given_GC": round(p_gc, 3)}
            print("%-11s rho=%+.3f  LOO[%+.3f,%+.3f] 稳定=%s  GC偏相关=%+.3f" % (
                feat, rho, lo, hi, stable, p_gc))
        report["contexts"][ctx] = {"rows": rows, "features": feats}

    zc = report["contexts"]["mature18"]["features"]["z_coreDR"]
    r1 = (abs(zc["spearman"]) < 0.3) or (not zc["loo_stable"])
    report["verdict"] = {
        "R1_raw_cross_pp_attribution_composition": "支持退出主评分" if r1 else "校正后仍有信号, 需复核",
        "R1_detail": zc,
        "R2_stage_split": ("flank 与 coreDR 方向不一致" if
                           np.sign(report["contexts"]["pre25"]["features"]["raw_flank"]["spearman"]) !=
                           np.sign(report["contexts"]["pre25"]["features"]["raw_coreDR"]["spearman"])
                           else "方向一致, 拆分仍必要(flank 属加工模块)")}
    print("\n判据: R1 %s | R2 %s" % (report["verdict"]["R1_raw_cross_pp_attribution_composition"],
                                    report["verdict"]["R2_stage_split"]))
    json.dump(report, open(args.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % args.out)


if __name__ == "__main__":
    main()
