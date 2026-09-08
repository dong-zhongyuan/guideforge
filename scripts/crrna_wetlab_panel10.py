#!/usr/bin/env python3
"""压缩版湿实验 panel(≤10 条)——体外切割(旁切/cis)验证下单用。

导师审定口径(docs/wetlab_plan.md §2, v3):体外旁切初筛候选控制在 10 条
以内——主靶 TP53-R248Q 上 4 骨架全测(裁决主战场 + 因果检验臂),其余 3 靶
各测 WT + 1 条综合取向代表(A1G+U3A+A8G+U15C),合计 4 + 3×2 = 10 条。
细胞实验只带体外裁决出的最优 2-3 条 + WT 对照(不在本 panel)。

选取按 (靶标, 骨架) 逻辑而非硬编码 oligo 名——对 APC 键名口径变更
(Q1312x->Q1328x)与骨架命名稳健。序列来源: data/ivt_round1_order_sheet.csv
(32 条全量, 由 crrna_ivt_template.py 生成), 本脚本只做子集抽取 + T7 模板
拼接(与 crrna_dna_template_sheet.py 同口径), 不重算任何设计。

输出:
  data/wetlab_panel10_order.csv   —— 10 条 DNA 合成模板下单表(T7+DR+spacer)
  data/wetlab_panel10.fasta       —— 同 10 条 62nt 合成模板 FASTA(可直接贴给合成商)

Run: python scripts/crrna_wetlab_panel10.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "data", "ivt_round1_order_sheet.csv")
DST_CSV = os.path.join(ROOT, "data", "wetlab_panel10_order.csv")
DST_FA = os.path.join(ROOT, "data", "wetlab_panel10.fasta")

T7 = "TAATACGACTCACTATAGG"  # 19nt, wetlab_plan.md 统一口径

# ---- panel 选取逻辑(wetlab_plan.md §2)----
MAIN_TARGET = "TP53-R248Q"
# 主靶 4 骨架: WT 参考 + 茎稳定化代表 + 综合取向代表 + 补偿因果臂
MAIN_SCAFFOLDS = ["WT", "A8C+U15G", "A1G+U3A+A8G+U15C", "B_break_compensate"]
# 其余 3 靶各 2 条: WT 参考 + 综合取向代表
OTHER_SCAFFOLDS = ["WT", "A1G+U3A+A8G+U15C"]
# 靶标角色标注(与 wetlab_plan §2 表一致)
TARGET_ROLE = {
    "TP53-R248Q": "主验证·HCT116·4骨架全测(裁决主战场+因果臂)",
    "KRAS-G12D": "泛化·SW480·WT+综合取向代表",
    "TP53-R273H": "同基因第二位点·WT+综合取向代表",
    "APC-Q1328x": "无义·最干净鉴别·WT+综合取向代表",
}
SCAFFOLD_ROLE = {
    "WT": "参考株(阳性基线)",
    "A8C+U15G": "茎稳定化取向代表",
    "A1G+U3A+A8G+U15C": "综合取向代表",
    "B_break_compensate": "补偿臂·因果检验(round-3 去混杂重设计)",
}


def selected(target, scaffold):
    if target == MAIN_TARGET:
        return scaffold in MAIN_SCAFFOLDS
    return scaffold in OTHER_SCAFFOLDS


def main():
    with open(SRC, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    panel = []
    for r in rows:
        tgt = r["target"].strip()
        scaf = r["scaffold_desc"].strip()
        if not selected(tgt, scaf):
            continue
        rna = r["sequence_5to3_RNA"].strip().upper()
        if "U" not in rna:
            raise SystemExit(f"not RNA: {r['oligo_name']}")
        dna_core = rna.replace("U", "T")
        template = T7 + dna_core
        panel.append({
            "oligo_name": r["oligo_name"],
            "target": tgt,
            "target_role": TARGET_ROLE.get(tgt, ""),
            "scaffold_desc": scaf,
            "scaffold_role": SCAFFOLD_ROLE.get(scaf, ""),
            "spacer_dna_24nt": r["spacer_dna_ref"].strip().upper(),
            "rna_crRNA_43nt": rna,
            "dna_core_43nt": dna_core,
            "synthesis_template_62nt": template,
            "note": "模板=T7(19)+DR(19)+spacer(24); 若IVT试剂盒要求转录起点G, "
                    "用T7(17nt,去尾G)+G+core方案并核验5'端",
        })

    # ---- 排序: 主靶 4 条(按 MAIN_SCAFFOLDS 序)在前, 其余靶按 OTHER_SCAFFOLDS 序 ----
    def sort_key(x):
        torder = 0 if x["target"] == MAIN_TARGET else 1
        slist = MAIN_SCAFFOLDS if x["target"] == MAIN_TARGET else OTHER_SCAFFOLDS
        return (torder, x["target"], slist.index(x["scaffold_desc"]))
    panel.sort(key=sort_key)

    # ---- self-checks ----
    assert len(panel) == 10, f"expect 10 rows, got {len(panel)}"
    tgt_counts = {}
    for x in panel:
        tgt_counts[x["target"]] = tgt_counts.get(x["target"], 0) + 1
        assert len(x["dna_core_43nt"]) == 43, f"core!=43nt: {x['oligo_name']}"
        assert len(x["synthesis_template_62nt"]) == 62, f"tpl!=62nt: {x['oligo_name']}"
        assert len(x["spacer_dna_24nt"]) == 24, f"spacer!=24nt: {x['oligo_name']}"
        assert x["synthesis_template_62nt"].startswith(T7)
        assert set(x["dna_core_43nt"]) <= {"A", "C", "G", "T"}
        # crRNA 43nt = DR19 + spacer24; 核对 spacer 落在 3' 端
        assert x["dna_core_43nt"][19:] == x["spacer_dna_24nt"], \
            f"spacer mismatch: {x['oligo_name']}"
    assert tgt_counts.get(MAIN_TARGET) == 4, f"main target != 4: {tgt_counts}"
    assert len(tgt_counts) == 4, f"expect 4 targets, got {tgt_counts}"
    for t, n in tgt_counts.items():
        if t != MAIN_TARGET:
            assert n == 2, f"{t} != 2 rows: {n}"
    # APC 必须是新口径 Q1328x + 新 spacer(防旧 Q1312x 回流)
    apc = [x for x in panel if x["target"].startswith("APC")]
    assert apc and all("Q1328x" in x["target"] for x in apc), \
        f"APC 口径错误(应为 Q1328x): {[x['target'] for x in apc]}"
    assert all(x["spacer_dna_24nt"] == "TGACACTGCTGGAACTTCGCTCAC" for x in apc), \
        "APC spacer 非 Q1328* 口径"

    fields = list(panel[0].keys())
    with open(DST_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(panel)

    with open(DST_FA, "w", encoding="utf-8") as f:
        for x in panel:
            f.write(f">{x['oligo_name']} {x['target']}|{x['scaffold_desc']}|"
                    f"len62|T7+DR19+spacer24\n{x['synthesis_template_62nt']}\n")

    print(f"OK: 10 条 panel | 靶标分布={tgt_counts}")
    print(f"主靶 {MAIN_TARGET} 骨架: "
          f"{[x['scaffold_desc'] for x in panel if x['target'] == MAIN_TARGET]}")
    print(f"wrote: {DST_CSV}")
    print(f"wrote: {DST_FA}")


if __name__ == "__main__":
    main()
