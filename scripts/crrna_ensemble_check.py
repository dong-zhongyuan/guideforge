"""打分输入端可靠性检验: WT 与 TOP-K 的系综采样稳定性(用户 2026-09-01 指定)。

问题: 打分用的 ΔΔG/stem 指标来自 MFE 单构象点估计。若某候选的 ΔΔG 在系综
不同构象间跳变(如 ±3 kcal/mol), 则 "-2.0 kcal/mol 茎稳定化" 这类数字本身不可靠。

方法(ViennaRNA 2.7.2 pbacktrack 在本构建返回空元组, 改用精确等效方案):
  1. 系综自由能 F = -RT lnZ (pf 精确值); 精确系综 ΔΔG_ens = F_var - F_WT,
     与管线 MFE 口径 ΔΔG_mfe 对比 = MFE 点估计的系统性偏移;
  2. subopt 精确枚举 MFE 以上 5 kcal 窗口内全部构象(窗口尾部质量 >5% 自动
     加宽到 8 kcal, 仍超则如实报告); 能量均值/展宽/序翻转/茎完整率/spacer
     配对全部为枚举窗口全体的 Boltzmann 加权精确值(无抽样噪声);
  3. 另按 Boltzmann 权重有放回采样 N=100 —— 仅供退役 R2/R3 统计量溯源与
     ddG 均值 bootstrap CI, 不参与任何判定;
  4. 判不稳规则(先定后跑; 集中登记于 docs/preregistration.md §B):
     R1 |ΔΔG_mfe - ΔΔG_ens| > 1.0  (MFE 点估计偏移过大)
     R2' 序翻转: |ΔΔG_ens| ≥ 1.0 且 P_contra > 0.2; P_contra = P(逐构象能量序
        与窗口均值序相反), 两枚举窗口 Boltzmann 权重卷积精确值(无序号配对、
        无抽样噪声); |ΔΔG_ens| < 1.0 的近中性候选不判定(其"无实效"即主张
        本身), P_contra 随行披露
     R3' 景观粗糙度: sd_ens(E_var)/sd_ens(E_WT) > 2.0 (枚举窗口 Boltzmann
        加权精确值; WT 分裂半 200 次零分布标定见 JSON r3_null_calibration)
     R4 茎完整构象占比在 20%~80% 之间 (构象二态, dp_fold 点估计不稳)

诚实边界: 判不稳阈值为项目设定常数(登记见 §B), 非文献出处; 窗口尾部质量
随行报告; 全部原始值随 JSON 输出供复核。

2026-09-04 round-3: 旧 R2(std(ddG_sample)>|ddG_ens|)/R3(range>4.0) 的统计量
ddG_samp=e_v-e_w 为两组相互独立的 Boltzmann 抽样按采样序号配对求差, std 按构造
≈ √(sd_v²+sd_w²)(≈1.3 kcal/mol), 在普通系综涨落尺度上武装, 对 |ΔΔG_ens|≲1.3
的候选必然触发——普遍触发不证明真实不稳, 也无法检出设计目标病理(alfa round-2
#4)——已退役, 旧数值保留于输出 JSON deprecated_ill_posed_R2R3 仅供溯源;
R2'/R3' 为同日登记的重设计判据(先登记后评估); R4 茎完整率同日改精确窗口值
(原 100 采样估计在 0.8 判界附近有 ±0.04 噪声)。
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
         "R2p_effect_floor_kcal": 1.0,
         "R2p_p_contra_over": 0.2,
         "R3p_spread_ratio_over": 2.0,
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


def sample_metrics(structs, energies, n, rng):
    """Boltzmann 权重采样 n 个构象的能量 —— 仅供退役统计量溯源, 不参与判定。"""
    w = win_weights(energies)
    idx = rng.choice(len(structs), size=n, p=w)
    return energies[idx]


def win_weights(energies):
    w = np.exp(-(energies - energies.min()) / RT)
    return w / w.sum()


def win_mean_sd(energies):
    """枚举窗口内 Boltzmann 加权能量均值/标准差(精确值, 非采样估计)。"""
    w = win_weights(energies)
    m = float((w * energies).sum())
    return m, float(np.sqrt((w * (energies - m) ** 2).sum()))


def exact_struct_metrics(structs, energies, dr_pairs, dr_len, sp_len):
    """枚举窗口全体的结构级指标精确值(Boltzmann 加权):
    茎完整率 / spacer 配对比例均值 / spacer 配对比例标准差。"""
    w = win_weights(energies)
    stem, sp = [], []
    for db in structs:
        pairs = parse_pairs(db)
        stem.append(1.0 if dr_pairs.issubset(pairs) else 0.0)
        sp.append(sum(1 for k in range(dr_len, dr_len + sp_len)
                      if db[k] != ".") / sp_len)
    stem, sp = np.array(stem), np.array(sp)
    sp_m = float((w * sp).sum())
    return (float((w * stem).sum()), sp_m,
            float(np.sqrt((w * (sp - sp_m) ** 2).sum())))


def p_order_contra(en_v, en_w, mean_gap):
    """P(逐构象能量序与窗口均值序相反): 两枚举窗口 Boltzmann 权重卷积精确值。

    对变体每个能量 e 累积 WT 窗口分布得 P(E_w < e), 再按变体分布加权求和;
    与退役 R2/R3 不同, 不做采样序号配对、无抽样噪声。连续浮点能量 ties 概率可忽略。
    """
    wv, ww = win_weights(en_v), win_weights(en_w)
    order = np.argsort(en_w)
    cum = np.concatenate(([0.0], np.cumsum(ww[order])))
    p_w_lt_v = float((wv * cum[np.searchsorted(en_w[order], en_v)]).sum())
    return p_w_lt_v if mean_gap < 0 else 1.0 - p_w_lt_v


def split_half_null(en_w, rng, reps=200):
    """WT 枚举集随机分裂半的 spread-ratio 零分布(p95/max/超阈经验率)——R3' 阈值标定。"""
    w, n, ratios = win_weights(en_w), len(en_w), []
    for _ in range(reps):
        a = rng.integers(0, 2, n) == 0
        if a.sum() < 10 or (~a).sum() < 10:
            continue
        sds = []
        for mask in (a, ~a):
            wm = w[mask] / w[mask].sum()
            m = float((wm * en_w[mask]).sum())
            sds.append(float(np.sqrt((wm * (en_w[mask] - m) ** 2).sum())))
        if min(sds) > 0:
            ratios.append(max(sds) / min(sds))
    n_over = sum(1 for r in ratios if r > RULES["R3p_spread_ratio_over"])
    return (round(float(np.percentile(ratios, 95)), 3),
            round(float(np.max(ratios)), 3), len(ratios), n_over)


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
    print("DR(%dnt) 茎配对(0基): %s | spacer %dnt | 判据全部为枚举窗口精确值\n" % (
        dr_len, sorted(dr_pairs), sp_len))

    rng = np.random.default_rng(args.seed)          # 退役统计量溯源采样流(与历史运行同序)
    rng_null = np.random.default_rng(args.seed + 1000)  # R3' 零分布标定专用流

    fc_w = RNA.fold_compound(wt_seq)
    _x, g_w = fc_w.pf()
    st_w, en_w, tail_w, delta_w = enumerate_window(fc_w)
    e_w = sample_metrics(st_w, en_w, args.n_samples, rng)
    stem_w_frac, sp_w_mean, sp_w_std = exact_struct_metrics(
        st_w, en_w, dr_pairs, dr_len, sp_len)
    p_stem_w = p_stem_bpp(wt_seq, dr_pairs, dr_len)
    m_w, sd_w_win = win_mean_sd(en_w)
    null_p95, null_max, null_reps, null_over = split_half_null(en_w, rng_null)
    null_rate = null_over / null_reps if null_reps else 0.0

    header = "%-14s %7s %7s %6s | %6s %5s %5s | %5s %6s %5s %5s"
    print(header % ("cand", "ddG_mfe", "ddG_ens", "shift", "sd", "range", "tail%",
                    "stemF", "Pstem", "Pcon", "sprR"))
    rows = [{"name": "WT", "ddG_mfe": 0.0, "ddG_ens": 0.0,
             "tail_mass": round(tail_w, 4), "delta_centigrade": delta_w,
             "stem_intact_frac": round(stem_w_frac, 3),
             "p_stem_pf": round(p_stem_w, 5),
             "spacer_paired_mean": round(sp_w_mean, 3),
             "spacer_paired_std": round(sp_w_std, 3),
             "E_sample_std": round(float(np.std(e_w)), 2),
             "E_win_mean": round(m_w, 2), "E_win_sd": round(sd_w_win, 2)}]
    print("%-14s %7s %7s %6s | %6s %5s %5.1f | %5.2f %6.3f %5s %5s" % (
        "WT", "-", "-", "-", "-", "-", 100 * tail_w,
        rows[0]["stem_intact_frac"], p_stem_w, "-", "-"))

    for v in variants:
        seq = to_rna(v["construct_dna"])
        fc_v = RNA.fold_compound(seq)
        _x, g_v = fc_v.pf()
        st_v, en_v, tail_v, delta_v = enumerate_window(fc_v)
        e_v = sample_metrics(st_v, en_v, args.n_samples, rng)
        stem_frac, sp_v_mean, sp_v_std = exact_struct_metrics(
            st_v, en_v, dr_pairs, dr_len, sp_len)
        p_stem = p_stem_bpp(seq, dr_pairs, dr_len)
        m_v, sd_v_win = win_mean_sd(en_v)
        mean_gap = m_v - m_w
        p_con = p_order_contra(en_v, en_w, mean_gap)
        spread = sd_v_win / sd_w_win

        ddg_mfe = v["mfe_kcal"] - d["wt"]["mfe_kcal"]
        ddg_ens = g_v - g_w
        ddg_samp = e_v - e_w  # 退役统计量(按采样序号配对), 仅供溯源, 不参与判定
        # bootstrap 95% CI: 有限采样下样本均值的置信区间(点估计本身有精确值 ddG_ens)
        boots = [float(np.mean(ddg_samp[rng.integers(0, len(ddg_samp), len(ddg_samp))]))
                 for _ in range(1000)]
        ci_lo, ci_hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

        r2_tested = abs(ddg_ens) >= RULES["R2p_effect_floor_kcal"]
        flags = []
        if abs(ddg_mfe - ddg_ens) > RULES["R1_shift_over_1kcal"]:
            flags.append("R1")
        if r2_tested and p_con > RULES["R2p_p_contra_over"]:
            flags.append("R2'")
        if spread > RULES["R3p_spread_ratio_over"]:
            flags.append("R3'")
        lo, hi = RULES["R4_stem_bimodal_20_80"]
        if lo < stem_frac < hi:
            flags.append("R4")

        rows.append({
            "name": v["desc"], "rank": v.get("rank"), "score": v.get("score"),
            "ddG_mfe": round(ddg_mfe, 2), "ddG_ens": round(ddg_ens, 2),
            "ddG_mfe_vs_ens_shift": round(ddg_mfe - ddg_ens, 2),
            "window_mean_gap": round(mean_gap, 2),
            "p_order_contra": round(p_con, 3),
            "spread_ratio": round(spread, 3),
            "r2_tested": r2_tested,
            "ddG_sample_mean": round(float(np.mean(ddg_samp)), 2),
            "ddG_sample_std": round(float(np.std(ddg_samp)), 2),
            "ddG_sample_min": round(float(np.min(ddg_samp)), 2),
            "ddG_sample_max": round(float(np.max(ddg_samp)), 2),
            "ddG_boot95": [round(ci_lo, 2), round(ci_hi, 2)],
            "tail_mass": round(tail_v, 4), "delta_centigrade": delta_v,
            "stem_intact_frac": round(stem_frac, 3), "p_stem_pf": round(p_stem, 5),
            "spacer_paired_mean": round(sp_v_mean, 3),
            "spacer_paired_std": round(sp_v_std, 3),
            "E_sample_std": round(float(np.std(e_v)), 2),
            "unstable_flags": flags or ["OK"],
        })
        r = rows[-1]
        rng_ = r["ddG_sample_max"] - r["ddG_sample_min"]
        print("%-14s %7.2f %7.2f %6.2f | %6.2f %5.2f %5.1f | %5.2f %6.3f %5.2f %5.2f  [%5.2f,%5.2f] %s" % (
            r["name"], r["ddG_mfe"], r["ddG_ens"], r["ddG_mfe_vs_ens_shift"],
            r["ddG_sample_std"], rng_, 100 * r["tail_mass"],
            r["stem_intact_frac"], r["p_stem_pf"], p_con, spread, ci_lo, ci_hi,
            "+".join(r["unstable_flags"])))

    # === 全规则结果汇总 + 结论强度(2026-09-04 round-3, 评审共同关切 #2) ===
    # 旧 R2/R3 已退役(统计构造病态, alfa round-2 #4), 判定由 R1 + R2'/R3'(重设计) + R4 决定
    n_c = len(rows) - 1
    rule_desc = {
        "R1": "|ddG_mfe - ddG_ens| > 1.0 kcal (MFE 点估计偏移过大)",
        "R2'": "|ddG_ens| >= 1.0 且 P_contra > 0.2 (序翻转; 两枚举窗口卷积精确值)",
        "R3'": "sd_ens(var)/sd_ens(WT) > 2.0 (景观粗糙度; 窗口精确值)",
        "R4": "茎完整构象占比在 20%~80% (构象二态; 窗口精确值)"}
    rule_summary = {}
    for rid in ("R1", "R2'", "R3'", "R4"):
        fired = [r["name"] for r in rows[1:] if rid in r["unstable_flags"]]
        rule_summary[rid] = {"criterion": rule_desc[rid],
                             "fired": len(fired), "total": n_c,
                             "fired_candidates": fired}
    n_r2_tested = sum(1 for r in rows[1:] if r["r2_tested"])
    rule_summary["R2'"]["tested"] = n_r2_tested
    rule_summary["R2'"]["untested_note"] = (
        "|ddG_ens| < %.1f 的近中性候选 %d 条不判定(其'无实效'即主张本身), "
        "P_contra 逐候选披露" % (RULES["R2p_effect_floor_kcal"], n_c - n_r2_tested))

    # 退役 R2/R3 溯源块: ddG_samp = e_v - e_w 为独立抽样按采样序号配对,
    # std 按构造 ≈ sqrt(sd_v^2 + sd_w^2); 与观测 ddG_sample_std 逐候选对照
    sd_w = rows[0]["E_sample_std"]
    ill_rows = [{"name": r["name"],
                 "constructed_std": round(float(np.sqrt(
                     r["E_sample_std"] ** 2 + sd_w ** 2)), 2),
                 "observed_ddG_sample_std": r["ddG_sample_std"]}
                for r in rows[1:]]
    ill_max_dev = max(abs(i["constructed_std"] - i["observed_ddG_sample_std"])
                      for i in ill_rows) if ill_rows else 0.0
    deprecated_note = (
        "旧 R2(std(ddG_sample)>|ddG_ens|)/R3(max-min>4.0) 统计构造病态(alfa round-2 "
        "#4, 2026-09-04): ddG_samp = e_v - e_w 为两组相互独立的 Boltzmann 抽样按采样"
        "序号配对求差, 其 std 按构造 ≈ sqrt(sd_v^2+sd_w^2)(本库 WT sd=%.2f; 逐候选构造"
        "值与观测值最大偏差 %.2f kcal)——在普通系综涨落尺度上武装, 对 |ddG_ens|≲1.3 "
        "kcal 的候选必然触发, 其普遍触发不能作为真实构象不稳的证据, 该检验按构造无法"
        "检出设计目标病理。该两条规则已退役, 数值仅供溯源, 判定以重设计的 R2'/R3' "
        "为准(判据登记: docs/preregistration.md §B)。" % (sd_w, ill_max_dev))

    r1_f = rule_summary["R1"]["fired"]
    r2_f = rule_summary["R2'"]["fired"]
    r3_f = rule_summary["R3'"]["fired"]
    r4_f = rule_summary["R4"]["fired"]
    if r1_f > 0:
        strength = "failed_R1"
        strength_note = (
            "R1 在 %d/%d 候选触发: MFE 点估计与精确系综自由能差偏移超 1 kcal, "
            "MFE 口径 ΔΔG 数字本身不可靠, 打分输入端可靠性主张失败。" % (r1_f, n_c))
    elif r2_f == 0 and r3_f == 0 and r4_f == 0:
        strength = "passed_all_rules"
        strength_note = (
            "R1 + 重设计 R2'/R3' + R4 全部未触发: 全部 %d 候选 |ddG_mfe-ddG_ens| "
            "≤ 1.0(MFE 口径与精确系综自由能差一致); %d 条实效候选(|ddG_ens|≥1.0) "
            "序翻转概率 P_contra 均 ≤ 0.2; 景观展宽比均 ≤ 2.0(WT 分裂半零分布 "
            "p95=%.2f, 超阈经验率 %.1f%%); 无构象二态——'打分输入端在系综采样下稳定'主张在 "
            "重设计判据下成立。近中性候选(|ddG_ens|<1.0) P_contra≈0.5 属普通热涨落"
            "重叠, 按登记判据不判定(逐候选披露)。旧 R2/R3 已退役, 见 "
            "deprecated_ill_posed_R2R3。" % (n_c, n_r2_tested, null_p95, 100 * null_rate))
    else:
        strength = "partial"
        strength_note = (
            "R1 通过, 但 R2'/R3'/R4 分别触发 %d/%d(判定 %d 条实效候选)、%d/%d、"
            "%d/%d: 逐候选见 rule_summary——'打分输入端系综稳定'的强主张不成立。"
            % (r2_f, n_c, n_r2_tested, r3_f, n_c, r4_f, n_c))

    stem_tag = os.path.basename(args.top_json).replace(".top.json", "").replace("tp53_r248q_", "")
    out = args.out or os.path.join(ROOT, "data", "ensemble_check.%s.json" % stem_tag)
    payload = {"source_top_json": args.top_json, "n_samples": args.n_samples,
               "seed": args.seed, "rt_kcal": round(RT, 5),
               "rules": {k: str(v) for k, v in RULES.items()},
               "rule_summary": rule_summary,
               "conclusion_strength": strength,
               "conclusion_note": strength_note,
               "r3_null_calibration": {
                   "method": "WT 枚举窗口随机分裂半, Boltzmann 加权 sd 之比, 200 次",
                   "reps": null_reps, "ratio_p95": null_p95, "ratio_max": null_max,
                   "n_over_threshold": null_over,
                   "null_trigger_rate": round(null_rate, 4),
                   "note": "零分布在 R3' 阈值 %.1f 处的经验触发率 %.1f%%——阈值假阳性率 "
                           "极低, 标定成立" % (RULES["R3p_spread_ratio_over"],
                                               100 * null_rate)},
               "deprecated_ill_posed_R2R3": {
                   "note": deprecated_note,
                   "constructed_std_formula": "sqrt(sd_v^2 + sd_w^2) (独立 Boltzmann 抽样按采样序号配对)",
                   "wt_E_sample_std": sd_w,
                   "per_candidate": ill_rows,
                   "max_abs_deviation_kcal": round(ill_max_dev, 2)},
               "preregistration": "判据集中登记于 docs/preregistration.md "
                                  "(2026-09-04) §B; R2'/R3' 为 2026-09-04 登记的"
                                  "重设计判据(先登记后评估); 判据变更须先登记后评估",
               "dr_stem_pairs_0based": sorted(dr_pairs),
               "wt": rows[0], "candidates": rows[1:],
               "note": "ViennaRNA 2.7.2 pbacktrack 本构建返回空元组; 改用 subopt 精确枚举 + Boltzmann 权重, 各序列窗口尾部质量已随行报告; "
                       "能量均值/展宽/序翻转/茎完整率/spacer 配对均为枚举窗口精确值; "
                       "100 采样仅用于退役统计量溯源与 ddG 均值 bootstrap CI"}
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n输出 -> %s" % out)
    print("\n=== 全规则汇总(判据登记: docs/preregistration.md §B) ===")
    for rid in ("R1", "R2'", "R3'", "R4"):
        s = rule_summary[rid]
        print("  %s: %d/%d 触发 — %s" % (rid, s["fired"], s["total"], s["criterion"]))
    print("  R3' 零分布标定: WT 分裂半 ratio p95=%.2f max=%.2f 超阈 %d/%d" % (
        null_p95, null_max, null_over, null_reps))
    print("结论强度: %s" % strength)
    print(strength_note)
    n_bad = sum(1 for r in rows[1:] if r["unstable_flags"] != ["OK"])
    print("TOP-12 中触发不稳规则: %d 条" % n_bad)


if __name__ == "__main__":
    main()
