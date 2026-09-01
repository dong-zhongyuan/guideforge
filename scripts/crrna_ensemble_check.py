"""打分输入端可靠性检验: WT 与 TOP-K 的系综采样稳定性(用户 2026-09-01 指定)。

问题: 打分用的 ΔΔG/stem 指标来自 MFE 单构象点估计。若某候选的 ΔΔG 在系综
不同构象间跳变(如 ±3 kcal/mol), 则 "-2.0 kcal/mol 茎稳定化" 这类数字本身不可靠。

方法(ViennaRNA 2.7.2 pbacktrack 在本构建返回空元组, 改用精确等效方案):
  1. 系综自由能 F = -RT lnZ (pf 精确值); 精确系综 ΔΔG_ens = F_var - F_WT,
     与管线 MFE 口径 ΔΔG_mfe 对比 = MFE 点估计的系统性偏移;
  2. subopt 精确枚举 MFE 以上 5 kcal 窗口内全部构象, Boltzmann 权重有放回
     采样 N=100(窗口尾部质量 >5% 自动加宽到 8 kcal, 仍超则如实报告);
  3. 每采样构象计算: 能量 / DR 茎 5 配对是否完整 / spacer 配对比例;
     变体与 WT 按采样序号配对得 ΔΔG 采样分布(std / max-min);
  4. 判不稳规则(先定后跑):
     R1 |ΔΔG_mfe - ΔΔG_ens| > 1.0  (MFE 点估计偏移过大)
     R2 std(ΔΔG_sample) > |ΔΔG_ens| (采样噪声超过效应量)
     R3 max-min(ΔΔG_sample) > 4.0  (跨构象跳变过大, 用户示例口径 ±3 即 range 6)
     R4 茎完整构象占比在 20%~80% 之间 (构象二态, dp_fold 点估计不稳)

诚实边界: 100 采样的 std 估计自身有约 10% 相对误差; 判不稳阈值为项目设定,
非文献出处; 全部分布原始值随 JSON 输出供复核。
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402

RT = 0.00198717 * 310.15  # kcal/mol, 37C

RULES = {"R1_shift_over_1kcal": 1.0,
         "R2_noise_over_effect": None,
         "R3_sample_range_over_4kcal": 4.0,
         "R4_stem_bimodal_20_80": (0.2, 0.8)}


def parse_pairs(db):
    stack, pairs = [], set()
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.add((stack.pop(), i))
    return pairs


def to_rna(s):
    return str(s).upper().replace("T", "U")


def enumerate_window(fc, deltas=(500, 800)):
    """精确枚举低能窗口构象, 返回 (结构列表, 能量数组, 尾部质量比, delta)。"""
    _x, g_ens = fc.pf()
    z_total = float(np.exp(-g_ens / RT))
    structs, energies, tail, delta_used = [], np.array([]), 1.0, None
    for delta in deltas:
        sub = fc.subopt(delta)
        structs = [s.structure for s in sub if s.structure]
        energies = np.array([s.energy for s in sub if s.structure])
        if len(structs) < 2:
            continue
        w = np.exp(-(energies - energies.min()) / RT)
        # Z_win = e^{-Emin/RT} * sum(w); Z_tot = e^{-F/RT} => 覆盖率 = e^{(F-Emin)/RT} * sum(w)
        cov = float(np.exp((g_ens - energies.min()) / RT) * w.sum())
        tail = max(0.0, 1.0 - cov)
        delta_used = delta
        if tail < 0.05 or delta == deltas[-1]:
            break
    return structs, energies, tail, delta_used


def sample_metrics(structs, energies, n, rng, dr_pairs, dr_len, sp_len):
    """Boltzmann 权重采样 n 个构象, 逐构象算指标。"""
    w = np.exp(-(energies - energies.min()) / RT)
    w = w / w.sum()
    idx = rng.choice(len(structs), size=n, p=w)
    e_s, stem_ok, sp_paired = [], [], []
    for i in idx:
        db = structs[i]
        e_s.append(energies[i])
        pairs = parse_pairs(db)
        stem_ok.append(1.0 if dr_pairs.issubset(pairs) else 0.0)
        sp_paired.append(sum(1 for k in range(dr_len, dr_len + sp_len)
                             if db[k] != ".") / sp_len)
    return np.array(e_s), np.array(stem_ok), np.array(sp_paired)


def p_stem_bpp(seq, dr_pairs, dr_len):
    """全长 bpp 矩阵上 WT DR 茎配对概率乘积(与管线 dp_fold 同口径)。"""
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    p = 1.0
    for i, j in dr_pairs:
        p *= bpp[i + 1][j + 1]
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top-json",
                    default=os.path.join(ROOT, "data", "tp53_r248q_pdbdr.v2.top.json"))
    ap.add_argument("--n-samples", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = json.load(open(args.top_json, encoding="utf-8"))
    entries = d["top"]
    wt_entry = entries[0] if entries[0].get("rank", -1) == 0 or entries[0].get("desc") == "WT" else None
    if wt_entry is None or "construct_dna" not in wt_entry:
        raise SystemExit("top.json 无 WT 构建序列条目, 请检查文件")
    spacer = d["spacer_fixed_dna"]
    variants = [e for e in entries if e.get("desc") not in (None, "WT")][:12]

    wt_seq = to_rna(wt_entry["construct_dna"])
    dr_len = len(variants[0]["dr_seq"])
    sp_len = len(spacer)

    wt_dr_db, _ = RNA.fold(wt_seq[:dr_len])
    dr_pairs = parse_pairs(wt_dr_db)
    print("DR(%dnt) 茎配对(0基): %s | spacer %dnt | 采样 %d/候选, seed=%d\n" % (
        dr_len, sorted(dr_pairs), sp_len, args.n_samples, args.seed))

    rng = np.random.default_rng(args.seed)

    g_w, st_w, en_w, tail_w, delta_w = None, [], np.array([]), 1.0, None
    fc_w = RNA.fold_compound(wt_seq)
    _x, g_w = fc_w.pf()
    st_w, en_w, tail_w, delta_w = enumerate_window(fc_w)
    e_w, stem_w, sp_w = sample_metrics(st_w, en_w, args.n_samples, rng,
                                       dr_pairs, dr_len, sp_len)
    p_stem_w = p_stem_bpp(wt_seq, dr_pairs, dr_len)

    header = "%-14s %7s %7s %6s | %6s %5s %5s | %5s %6s"
    print(header % ("cand", "ddG_mfe", "ddG_ens", "shift", "sd", "range", "tail%",
                    "stemF", "Pstem"))
    rows = [{"name": "WT", "ddG_mfe": 0.0, "ddG_ens": 0.0,
             "tail_mass": round(tail_w, 4), "delta_centigrade": delta_w,
             "stem_intact_frac": round(float(np.mean(stem_w)), 3),
             "p_stem_pf": round(p_stem_w, 5),
             "spacer_paired_mean": round(float(np.mean(sp_w)), 3),
             "spacer_paired_std": round(float(np.std(sp_w)), 3),
             "E_sample_std": round(float(np.std(e_w)), 2)}]
    print("%-14s %7s %7s %6s | %6s %5s %5.1f | %5.2f %6.3f" % (
        "WT", "-", "-", "-", "-", "-", 100 * tail_w,
        rows[0]["stem_intact_frac"], p_stem_w))

    for v in variants:
        seq = to_rna(v["construct_dna"])
        fc_v = RNA.fold_compound(seq)
        _x, g_v = fc_v.pf()
        st_v, en_v, tail_v, delta_v = enumerate_window(fc_v)
        e_v, stem_v, sp_v = sample_metrics(st_v, en_v, args.n_samples, rng,
                                           dr_pairs, dr_len, sp_len)
        p_stem = p_stem_bpp(seq, dr_pairs, dr_len)

        ddg_mfe = v["mfe_kcal"] - d["wt"]["mfe_kcal"]
        ddg_ens = g_v - g_w
        ddg_samp = e_v - e_w  # 按采样序号配对
        stem_frac = float(np.mean(stem_v))
        # bootstrap 95% CI: 有限采样下样本均值的置信区间(点估计本身有精确值 ddG_ens)
        boots = [float(np.mean(ddg_samp[rng.integers(0, len(ddg_samp), len(ddg_samp))]))
                 for _ in range(1000)]
        ci_lo, ci_hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

        flags = []
        if abs(ddg_mfe - ddg_ens) > RULES["R1_shift_over_1kcal"]:
            flags.append("R1")
        if float(np.std(ddg_samp)) > abs(ddg_ens):
            flags.append("R2")
        if float(np.max(ddg_samp) - np.min(ddg_samp)) > RULES["R3_sample_range_over_4kcal"]:
            flags.append("R3")
        lo, hi = RULES["R4_stem_bimodal_20_80"]
        if lo < stem_frac < hi:
            flags.append("R4")

        rows.append({
            "name": v["desc"], "rank": v.get("rank"), "score": v.get("score"),
            "ddG_mfe": round(ddg_mfe, 2), "ddG_ens": round(ddg_ens, 2),
            "ddG_mfe_vs_ens_shift": round(ddg_mfe - ddg_ens, 2),
            "ddG_sample_mean": round(float(np.mean(ddg_samp)), 2),
            "ddG_sample_std": round(float(np.std(ddg_samp)), 2),
            "ddG_sample_min": round(float(np.min(ddg_samp)), 2),
            "ddG_sample_max": round(float(np.max(ddg_samp)), 2),
            "ddG_boot95": [round(ci_lo, 2), round(ci_hi, 2)],
            "tail_mass": round(tail_v, 4), "delta_centigrade": delta_v,
            "stem_intact_frac": round(stem_frac, 3), "p_stem_pf": round(p_stem, 5),
            "spacer_paired_mean": round(float(np.mean(sp_v)), 3),
            "spacer_paired_std": round(float(np.std(sp_v)), 3),
            "E_sample_std": round(float(np.std(e_v)), 2),
            "unstable_flags": flags or ["OK"],
        })
        r = rows[-1]
        rng_ = r["ddG_sample_max"] - r["ddG_sample_min"]
        print("%-14s %7.2f %7.2f %6.2f | %6.2f %5.2f %5.1f | %5.2f %6.3f  [%5.2f,%5.2f] %s" % (
            r["name"], r["ddG_mfe"], r["ddG_ens"], r["ddG_mfe_vs_ens_shift"],
            r["ddG_sample_std"], rng_, 100 * r["tail_mass"],
            r["stem_intact_frac"], r["p_stem_pf"], ci_lo, ci_hi,
            "+".join(r["unstable_flags"])))

    stem_tag = os.path.basename(args.top_json).replace(".top.json", "").replace("tp53_r248q_", "")
    out = args.out or os.path.join(ROOT, "data", "ensemble_check.%s.json" % stem_tag)
    payload = {"source_top_json": args.top_json, "n_samples": args.n_samples,
               "seed": args.seed, "rt_kcal": round(RT, 5),
               "rules": {k: str(v) for k, v in RULES.items()},
               "dr_stem_pairs_0based": sorted(dr_pairs),
               "wt": rows[0], "candidates": rows[1:],
               "note": "ViennaRNA 2.7.2 pbacktrack 本构建返回空元组; 改用 subopt 精确枚举 + Boltzmann 权重采样, 各序列窗口尾部质量已随行报告"}
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n输出 -> %s" % out)
    n_bad = sum(1 for r in rows[1:] if r["unstable_flags"] != ["OK"])
    print("TOP-12 中触发不稳规则: %d 条" % n_bad)


if __name__ == "__main__":
    main()
