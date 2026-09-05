#!/usr/bin/env python3
"""Convert ivt_round1 RNA order sheet (43nt crRNA) -> DNA synthesis template sheet.

Adds T7 promoter prefix for IVT-based synthesis (wetlab_plan.md convention:
TAATACGACTCACTATAGG). Output columns: oligo name / target / scaffold /
43nt DNA core (U->T) / full synthesis template (61nt, T7+core) / note.

Run: python3 scripts/crrna_dna_template_sheet.py
"""
import csv

SRC = "data/ivt_round1_order_sheet.csv"
DST = "data/ivt_round1_dna_template_order.csv"
T7 = "TAATACGACTCACTATAGG"  # 19nt, wetlab_plan.md 统一口径


def main():
    rows = []
    with open(SRC, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rna = r["sequence_5to3_RNA"].strip().upper()
            if "U" not in rna:
                raise SystemExit(f"not RNA: {r['oligo_name']}")
            dna_core = rna.replace("U", "T")
            template = T7 + dna_core
            rows.append({
                "oligo_name": r["oligo_name"],
                "target": r["target"],
                "scaffold_desc": r["scaffold_desc"],
                "rna_43nt_ref": rna,
                "dna_core_43nt": dna_core,
                "synthesis_template_62nt": template,
                "note": "模板核心=DR19+spacer24; 若IVT试剂盒要求转录起点G, 用T7(17nt,去尾G)+G+core方案并核验5'端",
            })

    # ---- self-checks ----
    assert len(rows) == 32, f"expect 32 rows, got {len(rows)}"
    for x in rows:
        assert len(x["dna_core_43nt"]) == 43, f"core != 43nt: {x['oligo_name']}"
        assert len(x["synthesis_template_62nt"]) == 62, \
            f"template != 62nt: {x['oligo_name']}"
        assert x["synthesis_template_62nt"].startswith(T7)
        assert set(x["dna_core_43nt"]) <= {"A", "C", "G", "T"}
    targets = sorted({x["target"] for x in rows})
    scaffolds = sorted({x["scaffold_desc"] for x in rows})
    print(f"OK: {len(rows)} rows | targets={targets}")
    print(f"scaffolds({len(scaffolds)})={scaffolds}")

    with open(DST, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(list(rows[0].keys()))
        for x in rows:
            w.writerow(list(x.values()))
    print(f"wrote: {DST}")


if __name__ == "__main__":
    main()
