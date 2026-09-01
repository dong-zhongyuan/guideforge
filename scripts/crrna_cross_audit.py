"""旧候选库交叉配对审计(v1.6 硬过滤回溯适用性检查)。

v1.6 新增 DR-spacer 交叉配对硬过滤(cross_nt <= WT)。本脚本对既有主线运行的
通过变体重算 cross_nt/cross_pp, 回答: (1) 新过滤会拦截多少已通过变体;
(2) TOP-12 是否受影响。若 TOP-12 全部保留, 既有候选库结论不变, 仅需注明口径。

用法: python scripts/crrna_cross_audit.py --top-json data/tp53_r248q_pdbdr.v2.top.json
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402


def audit(top_json):
    d = json.load(open(top_json, encoding="utf-8"))
    wt_construct = d["top"][0]["construct_dna"]
    dr_len = len(d["top"][0]["dr_seq"])
    wt_full = core.to_rna(wt_construct)
    wt_cross_nt, wt_cross_pp = core.cross_pairs(wt_full, dr_len)

    rejected, kept, rows = [], [], []
    csv_path = top_json.replace(".top.json", ".variants.csv")
    with open(csv_path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["passed"] not in ("True", "true", "1"):
                continue
            full = core.to_rna(r["construct_dna"])
            cnt, pp = core.cross_pairs(full, dr_len)
            rows.append({"desc": r["desc"], "cross_nt": cnt, "cross_pp": round(pp, 3)})
            (rejected if cnt > wt_cross_nt else kept).append(r["desc"])

    top12 = [e["desc"] for e in d["top"] if e.get("desc") not in (None, "WT")][:12]
    top_hit = [n for n in top12 if n in set(rejected)]
    return {"source": os.path.basename(top_json), "wt_cross_nt": wt_cross_nt,
            "wt_cross_pp": round(wt_cross_pp, 3), "n_passed": len(rows),
            "n_rejected_by_v16": len(rejected), "rejected": rejected[:20],
            "top12_affected": top_hit, "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top-json", nargs="+", default=[
        os.path.join(ROOT, "data", "tp53_r248q_pdbdr.v2.top.json"),
        os.path.join(ROOT, "data", "tp53_r248q_zengdr.v6.top.json"),
        os.path.join(ROOT, "data", "cross_spacer.t0_R248Qg9.top.json"),
        os.path.join(ROOT, "data", "cross_spacer.t1_MYCg1.top.json")])
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "cross_audit.json"))
    args = ap.parse_args()

    report = []
    for tj in args.top_json:
        r = audit(tj)
        report.append(r)
        print("%s: WT cross_nt=%d | 通过 %d 条, v1.6 拦截 %d 条, TOP-12 受影响: %s" % (
            r["source"], r["wt_cross_nt"], r["n_passed"], r["n_rejected_by_v16"],
            r["top12_affected"] or "无"))
        if r["rejected"]:
            print("   拦截例: %s" % ", ".join(r["rejected"][:10]))
    json.dump(report, open(args.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % args.out)


if __name__ == "__main__":
    main()
