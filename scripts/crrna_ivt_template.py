"""IVT 首轮矩阵模板与回填校验(V3 §5.1 8x4 矩阵的回流前置件, 2026-09-02)。

模板: 8 骨架(WT + 四取向候选族去重 7 条) x 4 靶标(R248Q/G12C/G12D/R273H;
APC 为无义截短型, 骨架-向导互扰口径与错义一致, 可选 --with-apc 扩 10x8)
= 32 行, 每行含: 组合 id / 靶标 / 骨架 desc(与 variants.csv 及 BO ingest 一致) /
完整构建序列 / 空白测量列(Vmax_1-3, EC50_1-3, mismatch_stringency)。

校验: --check 对已填 CSV 查 schema、desc 白名单、数值范围(Vmax>0, EC50>0),
输出可直接喂 scripts/crrna_bayesopt.py --ingest 的聚合表(desc,Vmax,EC50 均值)。

用法:
  python scripts/crrna_ivt_template.py                 # 生成模板
  python scripts/crrna_ivt_template.py --check ivt_round1.csv --out-round1 bo_round1.csv
"""
import argparse
import csv
import json
import os
from statistics import mean

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

TARGETS = [("TP53-R248Q", "tp53_r248q"), ("KRAS-G12C", "kras_g12c"),
           ("KRAS-G12D", "kras_g12d"), ("TP53-R273H", "tp53_r273h")]

# 主靶标锚定项目正典 spacer(Zeng 2026 gRNA3), 保证与既有全部工作可比;
# 其余靶标用智能体 tilling 最优设计(源: data/agent/*.design.json)。
CANONICAL_SPACER = {"TP53-R248Q": "GTTCATGCCGCCCATGCAGGAACT"}


def load_scaffolds():
    lib = json.load(open(os.path.join(DATA, "orientation_library.json"),
                         encoding="utf-8"))
    seen, rows = {}, []
    for p in lib["orientations"]:
        b = p["best"]
        if b["desc"] in seen:
            continue
        seen[b["desc"]] = p["orientation"]
        rows.append({"desc": b["desc"], "dr_dna": b["construct_dna"][:19],
                     "orientation": p["orientation"]})
    wt = json.load(open(os.path.join(DATA, "tp53_r248q_zengdr.v6.top.json"),
                        encoding="utf-8"))["top"][0]
    return [{"desc": "WT", "dr_dna": wt["dr_seq"], "orientation": "reference"}] + rows


def load_spacers():
    out = []
    for label, key in TARGETS:
        d = json.load(open(os.path.join(DATA, "agent", key + ".design.json"),
                           encoding="utf-8"))
        top = d["designs"][0]
        sp = CANONICAL_SPACER.get(label) or top.get("spacer_dna") or top.get("spacer", "")
        out.append((label, sp))
    return out


def gen_template(path):
    scaf = load_scaffolds()
    sp = load_spacers()
    rows = []
    for label, spacer in sp:
        for s in scaf:
            rows.append({
                "combo_id": "%s|%s" % (label.replace("-", "_"), s["desc"].replace("+", "")),
                "target": label, "scaffold_desc": s["desc"],
                "orientation": s["orientation"],
                "spacer_dna": spacer,
                "construct_dna": s["dr_dna"] + spacer,
                "Vmax_1": "", "Vmax_2": "", "Vmax_3": "",
                "EC50_1": "", "EC50_2": "", "EC50_3": "",
                "mismatch_stringency": "",
                "note": ""})
    cols = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("模板 %d 行(骨架 %d x 靶标 %d) -> %s" % (len(rows), len(scaf), len(sp), path))


def check_filled(path, out):
    scaf = {s["desc"] for s in load_scaffolds()}
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    errs, agg = [], {}
    for i, r in enumerate(rows, 2):
        if r["scaffold_desc"] not in scaf:
            errs.append("行%d: 未知骨架 %s" % (i, r["scaffold_desc"]))
            continue
        vs = [float(r[k]) for k in ("Vmax_1", "Vmax_2", "Vmax_3") if r[k]]
        es = [float(r[k]) for k in ("EC50_1", "EC50_2", "EC50_3") if r[k]]
        if not vs or not es:
            continue
        if min(vs) <= 0 or min(es) <= 0:
            errs.append("行%d: Vmax/EC50 须为正" % i)
            continue
        agg[r["scaffold_desc"]] = {"Vmax": round(mean(vs), 4),
                                   "EC50": round(mean(es), 4)}
    if errs:
        print("校验失败:")
        for e in errs[:10]:
            print(" ", e)
        raise SystemExit(1)
    if out:
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["desc", "Vmax", "EC50"])
            for d, v in agg.items():
                w.writerow([d, v["Vmax"], v["EC50"]])
        print("聚合表(%d 骨架) -> %s (可喂 crrna_bayesopt.py --ingest)" % (len(agg), out))
    print("校验通过: %d 行, %d 骨架有数据" % (len(rows), len(agg)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", default=None, help="已填 CSV 路径")
    ap.add_argument("--out-round1", default=None, help="校验后聚合表输出")
    ap.add_argument("--out", default=os.path.join(DATA, "ivt_round1_template.csv"))
    args = ap.parse_args()
    if args.check:
        check_filled(args.check, args.out_round1)
    else:
        gen_template(args.out)


if __name__ == "__main__":
    main()
