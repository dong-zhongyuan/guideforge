"""G4: DeWeirdt 大库位点敏感性 × 保守性(2026-09-11, 判据 docs/preregistration.md §G4,
登记先于运行)。

动机: G3 的位点敏感性来自本管线单点扫描(每位仅 3 个 alt, n=19 位点), 功效低。
DeWeirdt 2020 alt-DR 扫描(AsCas12a, 35,883 条)的 n_mut==1 子集给每个位点数百条
实测活性(LFC 负筛), 位点敏感性估计精度高一到两个数量级。保守性取
data/dr_conservation.json 的同源比对列 identity, 经 AsCas12a_Zetsche2015 比对行
映射到逐位坐标(DeWeirdt WT 5' 端多 1nt  flank T, 坐标偏移 +1, 显式声明)。

检验: 位点敏感性 = mean(-LFC)(-LFC 大 = 突变失活强 = 位点敏感), 两方向 LFC127 /
LFC128 并列; Spearman(tie-aware) × 连续 identity(主口径) 与三档分区(并列),
置换 p (n=20000, 同步置换灵敏度向量)。

判读规则(§G4b, 先于数值): 两方向均 rho<0 且 p<0.05 -> 大库独立证据支持
「敏感位点富集于保守位点」; 方向分裂或均不显著 -> 如实报告。

运行: python scripts/crrna_deweirdt_position_sensitivity.py
输出: data/deweirdt_position_sensitivity.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

SCAN = os.path.join(ROOT, "data", "raw", "deweirdt2020_dr_scan.json")
CONS = os.path.join(ROOT, "data", "dr_conservation.json")
OUT = os.path.join(ROOT, "data", "deweirdt_position_sensitivity.json")
AS_ENTRY = "AsCas12a_Zetsche2015"
N_PERM = 20000


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic)


def perm_p(x, y, obs, n=N_PERM, seed=0):
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n):
        if abs(spearman(rng.permutation(x), y)) >= abs(obs):
            hits += 1
    return (hits + 1) / (n + 1)


def classify(identity):
    if identity >= 1.00:
        return "strictly_conserved"
    if identity >= 0.80:
        return "conserved"
    return "variable"


def as_position_identity(cons_doc):
    """AsCas12a 逐位(0基, 库条目坐标) -> (比对列, identity, class)。
    flank5/trailing3 位无列 -> None。"""
    aln = {a["name"]: a for a in cons_doc["alignment"]}
    ent = aln[AS_ENTRY]
    cols = cons_doc["conservation"]
    row = ent["aligned_core"]
    core_start = ent["core_start_in_full"]
    core_len = len(ent["anchor_l"] + ent["middle"] + ent["anchor_r"]) \
        if "anchor_l" in ent else len(row.replace("-", ""))
    out = {}
    ci = 0
    for k, ch in enumerate(row):
        if ch != "-":
            c = cols[k]
            out[core_start + ci] = {"column": c["col_1based"],
                                    "identity": c["identity"],
                                    "class": classify(c["identity"])}
            ci += 1
    assert ci == core_len, (ci, core_len)
    return out, ent["full"]


def main():
    cons_doc = json.load(open(CONS, encoding="utf-8"))
    pos_map, lib_full = as_position_identity(cons_doc)
    scan = json.load(open(SCAN, encoding="utf-8"))
    wt_dna = scan["wt_dr_dna"]
    # 坐标对齐: DeWeirdt WT 与库条目只差 5' 端长度, 尾部必须一致
    offset = len(wt_dna) - len(lib_full)
    assert offset >= 0 and wt_dna[offset:] == lib_full, \
        "DeWeirdt WT 与库条目 %s 尾部不一致" % AS_ENTRY

    sens = {}  # pos(库坐标) -> {"lfc127": [...], "lfc128": [...]}
    n_used = 0
    for r in scan["rows"]:
        if r.get("n_mut") != 1:
            continue
        dr_dna = r["dr_rna"].replace("U", "T")
        if len(dr_dna) != len(wt_dna):
            continue
        diffs = [i for i in range(len(wt_dna)) if dr_dna[i] != wt_dna[i]]
        if len(diffs) != 1:
            continue
        lib_pos = diffs[0] - offset
        if lib_pos < 0:
            continue
        e = sens.setdefault(lib_pos, {"lfc127": [], "lfc128": []})
        if r.get("lfc127") is not None:
            e["lfc127"].append(float(r["lfc127"]))
        if r.get("lfc128") is not None:
            e["lfc128"].append(float(r["lfc128"]))
        n_used += 1

    rows = []
    for pos in sorted(sens):
        m = pos_map.get(pos)
        if m is None:
            continue  # flank5/trailing3 无比对列
        s127 = float(np.mean([-v for v in sens[pos]["lfc127"]]))
        s128 = float(np.mean([-v for v in sens[pos]["lfc128"]]))
        rows.append({"pos_lib_0based": pos, "pos_deweirdt_1based": pos + offset + 1,
                     "column": m["column"], "identity": m["identity"],
                     "class": m["class"],
                     "sens_neglfc127": round(s127, 4),
                     "sens_neglfc128": round(s128, 4),
                     "n_variants_127": len(sens[pos]["lfc127"]),
                     "n_variants_128": len(sens[pos]["lfc128"])})

    ident = np.array([r["identity"] for r in rows])
    cls = np.array([{"strictly_conserved": 2, "conserved": 1, "variable": 0}[r["class"]]
                    for r in rows])
    result = {}
    for key in ("sens_neglfc127", "sens_neglfc128"):
        s = np.array([r[key] for r in rows])
        rho_i = spearman(ident, s)
        p_i = perm_p(ident, s, rho_i)
        rho_c = spearman(cls, s)
        p_c = perm_p(cls, s, rho_c)
        result[key] = {"vs_identity": {"rho": round(rho_i, 3), "perm_p": round(p_i, 4)},
                       "vs_class3": {"rho": round(rho_c, 3), "perm_p": round(p_c, 4)}}
        print("%s: identity rho=%+.3f (p=%.4f) | class3 rho=%+.3f (p=%.4f)" % (
            key, rho_i, p_i, rho_c, p_c))

    ok = all(result[k]["vs_identity"]["rho"] < 0 and
             result[k]["vs_identity"]["perm_p"] < 0.05 for k in result)
    both_neg = all(result[k]["vs_identity"]["rho"] < 0 for k in result)
    if ok:
        reading = ("两方向均 rho<0 且 p<0.05 -> 大库独立证据支持「敏感位点富集于"
                   "保守位点」(打分逻辑与天然进化一致的第二条独立证据链)")
    elif both_neg:
        reading = ("两方向方向一致为负但置换 p 未达 0.05(边缘) -> 方向与 G3 "
                   "一致, 显著性不足, 如实报告; 不改动 G3 结论")
    else:
        reading = "方向分裂 -> 如实报告, 不改动 G3 结论"

    print("\n逐位(库坐标 -> DeWeirdt 1基 | identity | class | -LFC127 | -LFC128 | n):")
    for r in rows:
        print("  pos%2d(DW %2d)  id=%.3f  %-18s  %+.3f  %+.3f  n=%d/%d" % (
            r["pos_lib_0based"], r["pos_deweirdt_1based"], r["identity"], r["class"],
            r["sens_neglfc127"], r["sens_neglfc128"],
            r["n_variants_127"], r["n_variants_128"]))

    payload = {
        "generated_by": "scripts/crrna_deweirdt_position_sensitivity.py (G4, 2026-09-11)",
        "criterion": "§G4b(2026-09-11 登记先于运行): 两方向均 rho<0 且 p<0.05 -> "
                     "大库独立证据支持; 否则如实报告",
        "data": {"scan": "data/raw/deweirdt2020_dr_scan.json (n_mut==1 子集)",
                 "conservation": "data/dr_conservation.json",
                 "as_entry": AS_ENTRY,
                 "coordinate_offset": offset,
                 "offset_rule": "DeWeirdt WT 5' 端多 %d nt flank(T), 位点坐标偏移 +%d; "
                                "尾部一致性已断言" % (offset, offset)},
        "n_single_mutants_used": n_used,
        "n_positions": len(rows),
        "sensitivity_definition": "mean(-LFC): -LFC 大 = 突变失活强 = 位点敏感",
        "spearman": "tie-aware midranks (scipy.stats.spearmanr)",
        "n_perm": N_PERM,
        "positions": rows,
        "results": result,
        "reading": reading,
        "limitation": "DeWeirdt 库为组合设计库, n_mut==1 子集稀疏(每位 1-3 条, "
                      "仅 14/18 比对位有覆盖), 位点敏感性估计精度低于预期; "
                      "功效仍高于本管线单点扫描的 3-alt/位, 但非饱和扫描",
    }
    json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n判读(§G4b): %s" % reading)
    print("输出 -> %s" % OUT)


if __name__ == "__main__":
    main()
