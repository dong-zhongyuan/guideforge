# -*- coding: utf-8 -*-
"""MD 预登记判据汇总(判据 data/md_analysis_plan.json, 2026-09-01 先于轨迹固定)。

读取 4 体系 analysis.json(crrna_md_analysis.py 产出), 应用预登记判据:
  M1 主判据: 候选茎 C1' RMSD 中位数 <= WT + 0.5 Å -> 判「不劣于」
  M2 主判据: 候选 DR 接触存活率 >= WT - 10 个百分点 -> 判「不劣于」
  M3 次要(氢键计数), 不作判据, 如实并列。
判读文字全部由数值按显式规则生成; 阴性结果同样报告。

运行: python scripts/crrna_md_verdict.py
输出: data/md_analysis_summary.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

TAGS = ("wt", "A7C_U14G", "A7G_U14C", "U6C_A15G")
M1_DELTA_A = 0.5    # 预注册: 候选 <= WT + 0.5 Å
M2_DELTA_PP = 10.0  # 预注册: 候选 >= WT - 10 个百分点


def main():
    ana = {}
    for tag in TAGS:
        p = os.path.join(ROOT, "data", "md", tag, "analysis.json")
        ana[tag] = json.load(open(p, encoding="utf-8"))
    wt = ana["wt"]
    rows, n_pass = [], 0
    for tag in TAGS[1:]:
        a = ana[tag]
        m1 = a["M1_stem_c1p_rmsd_A_median"]
        m2 = a["M2_contact_survival_frac"]
        m1_pass = m1 <= wt["M1_stem_c1p_rmsd_A_median"] + M1_DELTA_A
        m2_pass = m2 >= wt["M2_contact_survival_frac"] - M2_DELTA_PP / 100.0
        verdict = "不劣于WT" if (m1_pass and m2_pass) else "劣于WT"
        n_pass += m1_pass and m2_pass
        rows.append({"system": tag,
                     "M1_stem_rmsd_median_A": m1,
                     "M1_vs_WT_delta_A": round(
                         m1 - wt["M1_stem_c1p_rmsd_A_median"], 3),
                     "M1_pass": m1_pass,
                     "M2_contact_survival": m2,
                     "M2_vs_WT_delta_pp": round(
                         (m2 - wt["M2_contact_survival_frac"]) * 100, 1),
                     "M2_pass": m2_pass,
                     "M3_hbond_mean": a["M3_dr_donor_hbond_mean"],
                     "verdict": verdict})
        print("%-10s M1 %.3f(%+.3f) %s | M2 %.3f(%+.1fpp) %s | M3 %.1f -> %s" % (
            tag, m1, rows[-1]["M1_vs_WT_delta_A"],
            "过" if m1_pass else "不过",
            m2, rows[-1]["M2_vs_WT_delta_pp"],
            "过" if m2_pass else "不过",
            rows[-1]["M3_hbond_mean"], verdict))
    if n_pass == 3:
        reading = ("3/3 候选双主判据(M1/M2)不劣于 WT -> MD 层支持「双突变骨架在"
                   "结合态保持天然构象与界面接触」, 为结构保持过滤提供动力学"
                   "必要条件证据")
    elif n_pass > 0:
        reading = ("%d/3 候选双主判据不劣于 WT -> 部分支持, 未过体系如实列名" % n_pass)
    else:
        reading = "0/3 候选通过 -> MD 层阴性, 如实报告"
    payload = {
        "generated_by": "scripts/crrna_md_verdict.py (2026-09-11)",
        "criteria": "data/md_analysis_plan.json(2026-09-01 先于轨迹固定): "
                    "M1 候选茎RMSD中位<=WT+0.5A; M2 接触存活>=WT-10pp; M3 次要",
        "engine": "OpenMM(amber14/tip3p) 10ns x4 体系, 2-10ns 窗口; "
                  "轨迹分析 MDAnalysis 2.10.0 官方原语; "
                  "crRNA 按链 ID(B)定位(system.pdb 各链 resid 独立编号, "
                  "裸 resid 选择会串链, 2026-09-11 审计修正)",
        "wt_reference": {"M1_stem_rmsd_median_A": wt["M1_stem_c1p_rmsd_A_median"],
                         "M1_stem_rmsd_p90_A": wt["M1_stem_c1p_rmsd_A_p90"],
                         "M2_contact_survival": wt["M2_contact_survival_frac"],
                         "M3_hbond_mean": wt["M3_dr_donor_hbond_mean"]},
        "candidates": rows,
        "n_pass": n_pass,
        "reading": reading,
        "scope": "MD 阳性仅为必要条件(结合态结构保持可行), 非活性预测; "
                 "最终确认须体外/细胞实验",
    }
    out = os.path.join(ROOT, "data", "md_analysis_summary.json")
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n判读: %s" % reading)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
