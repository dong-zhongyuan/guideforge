#!/usr/bin/env python3
"""细胞杀伤 panel v6(2 靶 × 3 条 = 6 条)——湿实验唯一执行实验的下单用产物。

v6 口径变更(2026-09-16 用户指示): 变体收窄为 A1C 与 A8C+U15G 两条
(每靶 WT + 2 变体, 原 v5 为 WT + 3 变体); 细胞系经用户确认
TP53-R248Q -> HCT116(内源杂合)   TP53-R273H -> SW480(内源)。
同日另呈用户提供的 JD12 两条 23nt spacer: 核验均为 R273H 等位特异
(见 scripts/crrna_jd12_redesign.py 与 data/jd12_redesign.json),
其中 crRNA-2 可用但等位碱基在 spacer 第 17 位(远端), 选择性弱于
本表沿用的 canonical R273H 设计(PFS 鉴别 4.68x), 故仅作可选附表:
  data/wetlab_jd12_sp2_order.csv —— JD12-crRNA-2 x 3 骨架(可选, 不在 v6 预登记内)

变体选型依据(v5 四证据复合分沿用, 2026-09-09 全栈复核):
  A1C       保守通用型: 界面 Δ 0.000/-0.005(中性), 选型器顶档, trans 顶档,
            Protenix Δ-0.003 交叉一致, 耐受邻域 0.385
  A8C+U15G  茎稳定化取向代表: 界面 Δ -0.019/-0.006 轻降, ddG_dr -2.2,
            机制覆盖(茎稳定化)
v5 的 A1U+A2C(保守双突变) 因用户收窄面板退出, 证据与产物保留于 git 历史;
v5 8 条产物(wetlab_panel_cell8_*)与 v4 4 条保留为历史。

序列来源: data/ivt_round1_order_sheet.csv(32 条全量), 本脚本只做子集抽取
+ T7 模板拼接(与 crrna_dna_template_sheet.py 同口径), 不重算任何设计;
JD12 附表来源 data/jd12_redesign.json(脚本生成产物)。

输出:
  data/wetlab_panel_cell6_order.csv   —— 6 条 DNA 合成模板下单表(T7+DR+spacer)
  data/wetlab_panel_cell6.fasta       —— 同 6 条合成模板 FASTA
  data/wetlab_results_template_cell6.csv —— 结果回填模板(n=5 重复, 对齐功效建议)
  data/wetlab_jd12_sp2_order.csv      —— JD12-crRNA-2 可选附表(3 条)

Run: python scripts/crrna_wetlab_panel_cell.py [--prereg]
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "data", "ivt_round1_order_sheet.csv")
DST_CSV = os.path.join(ROOT, "data", "wetlab_panel_cell6_order.csv")
DST_FA = os.path.join(ROOT, "data", "wetlab_panel_cell6.fasta")
DST_TPL = os.path.join(ROOT, "data", "wetlab_results_template_cell6.csv")
DST_JD12 = os.path.join(ROOT, "data", "wetlab_jd12_sp2_order.csv")

T7 = "TAATACGACTCACTATAGG"  # 19nt, wetlab_plan.md 统一口径

# ---- panel 定义(v6: 2026-09-16 每靶 WT+2 变体, 用户收窄)----
PANEL = {  # 靶标 -> (细胞株, [骨架], 声明序即下单序)
    "TP53-R248Q": ("HCT116(内源 R248Q 杂合)",
                   ["WT", "A1C", "A8C+U15G"]),
    "TP53-R273H": ("SW480(内源 R273H)",
                   ["WT", "A1C", "A8C+U15G"]),
}
SCAFFOLD_ROLE = {
    "WT": "参考株(阳性基线, 野生型骨架)",
    "A1C": "保守通用型(界面中性+选型器/trans 顶档+Protenix 一致)",
    "A8C+U15G": "茎稳定化取向代表(ddG_dr -2.2, 机制覆盖)",
}


def main():
    with open(SRC, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    panel = []
    for r in rows:
        tgt = r["target"].strip()
        scaf = r["scaffold_desc"].strip()
        if tgt not in PANEL or scaf not in PANEL[tgt][1]:
            continue
        rna = r["sequence_5to3_RNA"].strip().upper()
        if "U" not in rna:
            raise SystemExit(f"not RNA: {r['oligo_name']}")
        dna_core = rna.replace("U", "T")
        panel.append({
            "oligo_name": r["oligo_name"],
            "target": tgt,
            "cell_line": PANEL[tgt][0],
            "scaffold_desc": scaf,
            "scaffold_role": SCAFFOLD_ROLE[scaf],
            "spacer_dna_24nt": r["spacer_dna_ref"].strip().upper(),
            "rna_crRNA_43nt": rna,
            "dna_core_43nt": dna_core,
            "synthesis_template_62nt": T7 + dna_core,
            "note": "模板=T7(19)+DR(19)+spacer(24); 若IVT试剂盒要求转录起点G, "
                    "用T7(17nt,去尾G)+G+core方案并核验5'端",
        })

    # 排序: 靶标按 PANEL 声明序, 骨架 WT 在前
    def sort_key(x):
        tgts = list(PANEL)
        return (tgts.index(x["target"]), PANEL[x["target"]][1].index(x["scaffold_desc"]))
    panel.sort(key=sort_key)

    # ---- self-checks ----
    assert len(panel) == 6, f"expect 6 rows, got {len(panel)}"
    per_target = {}
    for x in panel:
        per_target[x["target"]] = per_target.get(x["target"], 0) + 1
        assert len(x["dna_core_43nt"]) == 43, f"core!=43nt: {x['oligo_name']}"
        assert len(x["synthesis_template_62nt"]) == 62, f"tpl!=62nt: {x['oligo_name']}"
        assert len(x["spacer_dna_24nt"]) == 24, f"spacer!=24nt: {x['oligo_name']}"
        assert x["synthesis_template_62nt"].startswith(T7)
        assert x["dna_core_43nt"][19:] == x["spacer_dna_24nt"], \
            f"spacer mismatch: {x['oligo_name']}"
    assert per_target == {"TP53-R248Q": 3, "TP53-R273H": 3}, per_target
    spacers = {x["target"]: x["spacer_dna_24nt"] for x in panel}
    assert spacers["TP53-R248Q"] == "GTTCATGCCGCCCATGCAGGAACT"
    assert spacers["TP53-R273H"] == "CACCTCAAAGCTGTTCCGTCCCAG"

    fields = list(panel[0].keys())
    with open(DST_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(panel)
    with open(DST_FA, "w", encoding="utf-8") as f:
        for x in panel:
            f.write(f">{x['oligo_name']} {x['target']}|{x['scaffold_desc']}|"
                    f"{x['cell_line']}|len62|T7+DR19+spacer24\n"
                    f"{x['synthesis_template_62nt']}\n")

    # ---- 结果回填模板(n=5, 对齐 power_analysis_cellpanel 的 rep_scan 建议) ----
    with open(DST_TPL, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["target", "scaffold_desc"] +
                   ["survival_%d" % i for i in range(1, 6)] + ["note"])
        for x in panel:
            w.writerow([x["target"], x["scaffold_desc"], "", "", "", "", "",
                        "%s | 存活率%%(或相对蛋白量%%, 100=无杀伤); 填满后: "
                        "--ingest-cell / --anova / --twin-check" % x["cell_line"]])

    # ---- JD12-crRNA-2 可选附表(不在 v6 预登记内, 头对头备选) ----
    import json
    jd = json.load(open(os.path.join(ROOT, "data", "jd12_redesign.json"),
                        encoding="utf-8"))
    jd_rows = []
    for c in jd["designs"]["JD12_sp2_variants"]["constructs"]:
        dna = c["dna_core"]
        tpl = T7 + dna
        assert len(c["spacer_dna"]) == 23 and len(dna) == 42 and len(tpl) == 61
        assert dna[19:] == c["spacer_dna"]
        jd_rows.append({
            "oligo_name": "GF_SW480_JD12sp2_" + c["scaffold"],
            "target": "TP53-R273H", "cell_line": PANEL["TP53-R273H"][0],
            "scaffold_desc": c["scaffold"],
            "scaffold_role": SCAFFOLD_ROLE[c["scaffold"]] + "(JD12 附表)",
            "spacer_dna_23nt": c["spacer_dna"],
            "rna_crRNA_42nt": c["rna_crRNA"], "dna_core_42nt": dna,
            "synthesis_template_61nt": tpl,
            "note": "可选附表: JD12-crRNA-2(用户序列, R273H 等位第17位, "
                    "PFS GGGAG 良好, 全转录组 0 完全脱靶); 不在 v6 预登记内; "
                    "模板=T7(19)+DR(19)+spacer(23), IVT 后核验 5' 端",
        })
    with open(DST_JD12, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(jd_rows[0].keys()))
        w.writeheader()
        w.writerows(jd_rows)

    print(f"OK: 6 条细胞 panel v6 | {per_target}")
    for x in panel:
        print(f"  {x['oligo_name']:<32} {x['cell_line']}")
    print(f"wrote: {DST_CSV}")
    print(f"wrote: {DST_FA}")
    print(f"wrote: {DST_TPL}")
    print(f"wrote: {DST_JD12} (可选附表 3 条)")


def prereg_mode():
    """v6 6 条的实验前预测封存(2026-09-16, 湿实验开跑前登记)。

    汇集既有干实验产物(Chai 界面矩阵/16 对选型器/同源 trans/耐受邻域/
    Scholz 标定虚拟细胞/功效预分析)为逐构象预测 + 预设判定口径,
    写 data/prereg_cellpanel_v6_predictions.json; 配套登记文档
    docs/preregistration_cellpanel_v6.md(v5 存档保留)。先于任何实验结果入库;
    v5 封存(2026-09-09)后无实验数据回填, v6 仅收窄面板(A1U+A2C 退出),
    H1-H4 判读逻辑不变。
    """
    import json

    def load(fn):
        return json.load(open(os.path.join(ROOT, "data", fn), encoding="utf-8"))

    matrix = {(r["target"], r["scaffold"]): r
              for r in load("chai_matrix_4t.json")["rows"]}
    summary = {r["scaffold"]: r for r in load("chai_matrix_4t_summary.json")
               ["scaffold_deltas"]}
    fam = {r["scaffold"]: r for r in load("selector_family_ranking.json")["family"]}
    vc = load("virtual_cell_prior.json")
    vc_line = {r["cell_line"].split("_")[0]: r for r in vc["named_cell_lines"]}
    power = load("power_analysis_cellpanel.json")

    cell_of = {"TP53-R248Q": "HCT116", "TP53-R273H": "SW480"}
    entries = []
    for tgt, (cell, scafs) in PANEL.items():
        cl = cell_of[tgt]
        vcl = vc_line.get(cl, {})
        tag = "tp53"
        base_surv = vcl.get("surv_%s_mid" % tag)
        for s in scafs:
            m = matrix.get((tgt, s), {})
            f = fam.get(s, {})
            entries.append({
                "target": tgt, "cell_line": cell, "scaffold": s,
                "predictions": {
                    "chai_iptm": m.get("iptm_mean"),
                    "chai_prot_crRNA": m.get("prot_crRNA"),
                    "interface_delta_vs_wt": summary.get(s, {}).get(
                        "delta_vs_wt", {}).get(tgt),
                    "selector_fig1g_pred": f.get("pred_fig1g"),
                    "homolog_trans_end": f.get("pred_trans_end"),
                    "tolerance_neighborhood": f.get("tolerance_neighborhood"),
                    "virtual_cell_survival_mid_pct": base_surv,
                    "note": "虚拟细胞存活为靶基因x细胞系级预测(模型无骨架项), "
                            "骨架级差异由界面/选型器/trans 三列承载"},
                "preregistered_order_hypothesis": {
                    "kill_tier_1": ["WT", "A1C"],
                    "kill_tier_2": ["A8C+U15G"],
                    "reading": "同源证据预测 WT/A1C 杀伤同档(trans 顶档"
                               "+选型器顶档+界面中性), A8C+U15G 预测<=WT 档"
                               "(trans 弱+界面轻降); 档内排序不做预测"}})
    out = {
        "sealed": "2026-09-16, 先于任何 v6 panel 实验结果(v5 口径 2026-09-09 "
                  "登记后无任何实验数据回填, 收窄为 6 条重新封存)",
        "panel": "2 靶 x (WT+2 变体) = 6 条(wetlab_panel_cell6_*)",
        "hypotheses": {
            "H1_非劣性": "A1C 的杀伤不劣于 WT(界面中性+同源顶档;"
                        "预期 |Δkill| < MDD)",
            "H2_茎稳定化机制": "A8C+U15G 若茎稳定化转化为杀伤增强, "
                              "Δkill(=WT存活-变体存活)>0; 若蛋白界面主导, "
                              "Δ≈0 或为负",
            "H3_等位一致性": "两靶上 3 构象排序方向一致(基因层同集假设); "
                           "不一致 = 等位水平分型信号(探索性, 不设显著性口径)",
            "H4_虚拟细胞校验": "HCT116/R248Q 与 SW480/R273H 实测存活落点 vs "
                             "Scholz 标定预测(HCT116 70.4%/SW480 45.0%, "
                             "EC50 CI 场景带 + 20pp RNP→质粒移植边界)",
            "边界": "JD12-crRNA-2 可选附表(wetlab_jd12_sp2_order.csv)不在本"
                   "预登记判读口径内; 如纳入实验, 按探索性头对头叙述, 不做"
                   "显著性主张"},
        "decision_rules": {
            "显著性": "靶内成对单侧 t(α=0.05), 重复数按实际; 差异 <MDD 一律"
                     "方向性叙述不下显著性结论",
            "MDD_pp": power["pairwise_mdd_pp"],
            "rep_scan": power["rep_scan_mdd_pp"],
            "变体判优": "Δkill ≥ MDD 且两靶同向",
            "等位选择性": "SI=突变株杀伤/WT株对照杀伤, 变体 SI ≥ WT 的 SI 判非劣",
            "分析命令": "crrna_ivt_template.py --anova <结果CSV>(2x3 交互) 与 "
                       "--twin-check <结果CSV> --cell-csv <SI表>(虚拟细胞校验); "
                       "crrna_bayesopt.py --ingest-cell <结果CSV>(回流)"},
        "entries": entries}
    dst = os.path.join(ROOT, "data", "prereg_cellpanel_v6_predictions.json")
    json.dump(out, open(dst, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("预注册预测封存 -> %s(6 条逐构象预测 + H1-H4 + 判定口径)" % dst)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--prereg":
        prereg_mode()
    else:
        main()
