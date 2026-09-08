#!/usr/bin/env python3
"""细胞杀伤 panel(2 靶 × 2 骨架 = 4 条)——湿实验唯一执行实验的下单用产物。

导师口径(2026-09-07): 湿实验只做实验三(细胞杀伤验证), 靶点收窄为 2 个。
选靶硬约束 = 必须有内源携带该突变的现成细胞株:
  TP53-R248Q -> HCT116(内源杂合)   TP53-R273H -> SW480(内源)
每靶 2 条: WT 骨架(参考基线) + 综合取向代表 A1G+U3A+A8G+U15C(干实验排序
最优推荐; 无实验一体外裁决, 直接按干实验排序定)。KRAS-G12D / APC-Q1328x
本轮退出湿实验(内源匹配细胞株不可得/需工程化), 干实验产物全部保留。

序列来源: data/ivt_round1_order_sheet.csv(32 条全量), 本脚本只做子集抽取
+ T7 模板拼接(与 crrna_dna_template_sheet.py 同口径), 不重算任何设计。
10 条体外 panel(wetlab_panel10_*)保留为回补备用, 本脚本不改动它。

输出:
  data/wetlab_panel_cell4_order.csv  —— 4 条 DNA 合成模板下单表(T7+DR+spacer)
  data/wetlab_panel_cell4.fasta      —— 同 4 条 62nt 合成模板 FASTA

Run: python scripts/crrna_wetlab_panel_cell.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "data", "ivt_round1_order_sheet.csv")
DST_CSV = os.path.join(ROOT, "data", "wetlab_panel_cell4_order.csv")
DST_FA = os.path.join(ROOT, "data", "wetlab_panel_cell4.fasta")

T7 = "TAATACGACTCACTATAGG"  # 19nt, wetlab_plan.md 统一口径

# ---- panel 定义(导师 2026-09-07 口径)----
PANEL = {  # 靶标 -> (细胞株, [骨架])
    "TP53-R248Q": ("HCT116(内源 R248Q 杂合)",
                   ["WT", "A1G+U3A+A8G+U15C"]),
    "TP53-R273H": ("SW480(内源 R273H)",
                   ["WT", "A1G+U3A+A8G+U15C"]),
}
SCAFFOLD_ROLE = {
    "WT": "参考株(阳性基线, 野生型骨架)",
    "A1G+U3A+A8G+U15C": "综合取向代表(干实验排序最优推荐)",
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
    assert len(panel) == 4, f"expect 4 rows, got {len(panel)}"
    per_target = {}
    for x in panel:
        per_target[x["target"]] = per_target.get(x["target"], 0) + 1
        assert len(x["dna_core_43nt"]) == 43, f"core!=43nt: {x['oligo_name']}"
        assert len(x["synthesis_template_62nt"]) == 62, f"tpl!=62nt: {x['oligo_name']}"
        assert len(x["spacer_dna_24nt"]) == 24, f"spacer!=24nt: {x['oligo_name']}"
        assert x["synthesis_template_62nt"].startswith(T7)
        assert x["dna_core_43nt"][19:] == x["spacer_dna_24nt"], \
            f"spacer mismatch: {x['oligo_name']}"
    assert per_target == {"TP53-R248Q": 2, "TP53-R273H": 2}, per_target
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

    print(f"OK: 4 条细胞 panel | {per_target}")
    for x in panel:
        print(f"  {x['oligo_name']:<32} {x['cell_line']}")
    print(f"wrote: {DST_CSV}")
    print(f"wrote: {DST_FA}")


if __name__ == "__main__":
    main()
