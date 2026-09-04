"""Chai-1 4靶x8骨架界面矩阵汇总(程序化生成, round-3 评审整改 common concern 3)。

读取 data/chai_matrix_4t.json(32 行原始均值), 重算 WT 参照 / 差量 / 跨靶 spread,
reading 文字全部由 scaffold_deltas 数值按下列阈值规则程序化生成(方向与量级不手写,
消除手写叙述与数据脱节的可能):
  NEUTRAL_ABS   max|Δ| <= 0.03  -> 界面中性且靶标无关
  COLLAPSE_ABS  max|Δ| >= 0.25  -> 界面塌陷
  SPREAD_FLAG   跨靶 spread >= 0.10 -> 靶标相关(比较 KRAS/TP53 两族 |Δ| 均值的方向)
  NULL_SD_K     |Δ| <= 2 * iptm_sd -> 判为模型间噪声内(null); iptm_sd 为 aggregate
                ipTM 的 5 模型 sd, 作 prot-crRNA 链对噪声代理(后者 sd 未存盘)

显式口径声明(写入产物): 本矩阵全部分数来自 100% 一致自模板(8D4A链A)注入,
高 ipTM 为自模板复述的平凡结果, 只能作为模板合规性/健全性检查(sanity check),
不能作为界面可预测性证据; 在获得 template-free 或 scrambled-template 对照之前,
所有相关结论按此口径陈述。WT prot-crRNA ipTM≈0.49 低于可用界面预测区间(~0.5-0.6)。

运行: python scripts/crrna_chai_matrix_summary.py
"""
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

NEUTRAL_ABS = 0.03
COLLAPSE_ABS = 0.25
SPREAD_FLAG = 0.10
NULL_SD_K = 2.0

CLAIM_DOWNGRADE = (
    "Chai-1 '界面恢复' 为 100% 一致自模板(8D4A链A, 全长1207aa)复述: 注入与查询完全相同的"
    "蛋白坐标后, 高 aggregate/ipTM 是模板合规性的平凡结果, 只能作为模板合规性/健全性检查"
    "(sanity check), 不能作为界面可预测性证据。在获得 template-free 或 scrambled-template "
    "对照(以相同流程打分)之前, 本文件所有 reading 仅描述自模板口径下的数值表型, "
    "不得引为'分型假说的结构证据'或'骨架-靶标互作的直接观测'。"
)
PROT_CRRNA_BAND = (0.5, 0.6)  # 可用界面预测置信区间下沿的常用口径


def family_of(target):
    return target.split("_")[0]


def scaffold_reading(deltas, sds):
    """由单骨架的四靶差量(dict)与对应 iptm_sd(dict) 生成解读文字。"""
    mx = max(abs(v) for v in deltas.values())
    spread = round(max(deltas.values()) - min(deltas.values()), 4)
    nulls = sorted(t for t, v in deltas.items()
                   if abs(v) <= NULL_SD_K * sds[t])
    delta_str = ", ".join("%s %+.4f" % (t, v) for t, v in sorted(deltas.items()))
    parts = ["Δprot-crRNA ipTM vs 同靶WT: %s" % delta_str]
    if nulls:
        parts.append("%s 在 %dsd 噪声内(判 null)" % ("/".join(nulls), NULL_SD_K))
    if mx <= NEUTRAL_ABS:
        parts.append("max|Δ|=%.4f <= %.2f → 界面中性且靶标无关(通用型候选)" % (mx, NEUTRAL_ABS))
    else:
        direction = "下降" if all(v < 0 for v in deltas.values()) else \
                    ("上升" if all(v > 0 for v in deltas.values()) else "方向不一")
        lo, hi = min(abs(v) for v in deltas.values()), mx
        if mx >= COLLAPSE_ABS:
            parts.append("四靶界面均塌陷(|Δ| %.4f-%.4f, 方向: %s)" % (lo, hi, direction))
        else:
            grade = {"下降": "中等下降", "上升": "中等上升"}.get(direction, "中等变化且方向不一")
            parts.append("四靶界面%s(|Δ| %.4f-%.4f)" % (grade, lo, hi))
    if spread >= SPREAD_FLAG:
        fams = {}
        for t, v in deltas.items():
            fams.setdefault(family_of(t), []).append(abs(v))
        ranked = sorted(fams.items(), key=lambda kv: -sum(kv[1]) / len(kv[1]))
        a, b = ranked[0], ranked[1]
        ma, mb = sum(a[1]) / len(a[1]), sum(b[1]) / len(b[1])
        deepest = min(deltas, key=lambda t: deltas[t])
        if ma - mb >= SPREAD_FLAG / 2:
            parts.append("跨靶 spread %.4f ≥ %.2f → 靶标相关: %s |Δ|均值 %.4f > %s %.4f; "
                         "与骨架-靶标互作一致的单构建观测(自模板口径, 待对照校准)"
                         % (spread, SPREAD_FLAG, a[0], ma, b[0], mb))
        else:
            parts.append("跨靶 spread %.4f ≥ %.2f, 但 %s/%.4f 与 %s/%.4f 族均值几乎相同, "
                         "spread 由靶点个体差异(最深 %s)驱动, 非族间互作信号"
                         % (spread, SPREAD_FLAG, a[0], ma, b[0], mb, deepest))
    else:
        parts.append("跨靶 spread %.4f < %.2f → 未见靶标相关" % (spread, SPREAD_FLAG))
    return "; ".join(parts)


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "data",
                                                  "chai_matrix_4t_summary.json"),
                    help="汇总 JSON 输出路径(默认 data/chai_matrix_4t_summary.json)")
    args = ap.parse_args()
    src = os.path.join(ROOT, "data", "chai_matrix_4t.json")
    rows = json.load(open(src, encoding="utf-8"))["rows"]
    targets = sorted({r["target"] for r in rows})
    scaffolds = sorted({r["scaffold"] for r in rows if r["scaffold"] != "WT"})

    by = {(r["scaffold"], r["target"]): r for r in rows}
    wt_reference = {t: by[("WT", t)]["prot_crRNA"] for t in targets}
    wt_pc_mean = sum(wt_reference.values()) / len(wt_reference)

    scaffold_deltas, reading = [], {}
    for s in scaffolds:
        deltas, sds = {}, {}
        for t in targets:
            r = by[(s, t)]
            deltas[t] = round(r["prot_crRNA"] - wt_reference[t], 4)
            sds[t] = r["iptm_sd"]
        spread = round(max(deltas.values()) - min(deltas.values()), 4)
        scaffold_deltas.append({
            "scaffold": s,
            "delta_vs_wt": deltas,
            "cross_target_spread": spread,
            "iptm_sd_max": round(max(sds.values()), 4),
            "null_within_%dsd" % NULL_SD_K: [t for t in targets
                                             if abs(deltas[t]) <= NULL_SD_K * sds[t]]})
        reading[s] = scaffold_reading(deltas, sds)

    ct = [r["crRNA_target"] for r in rows]
    reading["note"] = ("crRNA_target 链对本矩阵全部 %.2f-%.2f(合成靶完全互补), "
                       "与旧矩阵原生靶(0.03-0.2)的差异=靶构造口径差, 两表分开引用"
                       % (min(ct), max(ct)))

    payload = {
        "engine": "chai-1 + 8D4A自模板m8 + ESM",
        "generated_by": "scripts/crrna_chai_matrix_summary.py (reading 由 scaffold_deltas 程序化生成)",
        "design": "4靶 x 8骨架(含WT与B_break_compensate); 靶RNA=protospacer窗口24nt+PFS 5nt(8D4A同构)",
        "n": len(rows),
        "claim_downgrade": CLAIM_DOWNGRADE,
        "prot_crRNA_usability": {
            "wt_mean": round(wt_pc_mean, 4),
            "usable_band": "~%.1f-%.1f" % PROT_CRRNA_BAND,
            "usable": wt_pc_mean >= PROT_CRRNA_BAND[0],
            "note": "WT prot-crRNA ipTM 四靶均值 %.3f 低于可用界面预测区间(~%.1f-%.1f), "
                    "即使撇开自模板问题, 该链对分数本身也不足以支撑界面水平结论"
                    % (wt_pc_mean, PROT_CRRNA_BAND[0], PROT_CRRNA_BAND[1])},
        "staleness": {
            "checked": "2026-09-04",
            "B_break_compensate": (
                "chai_matrix_4t.json 该行折叠对象为 round-3 前的旧补偿臂(DR "
                "AATTTCTACTCTTCTACAT, 茎破坏混杂臂, p_fold=0.0); 2026-09-04 起 "
                "orientation_library / ivt_round1_template / ivt_round1_order_sheet "
                "均已切换为去混杂新臂(zengDR+T10A/T12G/T19G, DR AATTTCTACAGGTGTAGAG, "
                "p_fold=0.89)。本矩阵 B 臂 reading 是**已退役分子**的数值表型(其"
                "'界面塌陷'正是该混杂的可预测表型), 不得引为当前分子的界面结论; "
                "重折叠待 Chai 环境机时。")},
        "reading_thresholds": {"neutral_abs": NEUTRAL_ABS, "collapse_abs": COLLAPSE_ABS,
                               "spread_flag": SPREAD_FLAG, "null_sd_k": NULL_SD_K,
                               "sd_note": "iptm_sd 为 aggregate ipTM 的 5 模型 sd, 作 "
                                          "prot-crRNA 噪声代理(后者 model sd 未存盘)"},
        "controls_pending": ["template-free 运行(同流程打分)",
                             "scrambled-template 运行(同流程打分)",
                             "prot-crRNA 链对的 per-model sd 存盘"],
        "wt_reference": wt_reference,
        "scaffold_deltas": scaffold_deltas,
        "reading": reading,
    }
    out = args.out
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)
    for s, txt in reading.items():
        print("[%s] %s" % (s, txt))


if __name__ == "__main__":
    main()
