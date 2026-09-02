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

import numpy as np

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


def two_way_anova(path):
    """骨架×靶标 两因素方差分析(含交互项)——V3 §5.3 答辩口径"交互项显著性"。

    平衡设计(每格 n=技术重复), 经典 SS 分解; 输入 = --check 通过的已填 CSV。
    """
    from scipy.stats import f as fdist

    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    cells = {}  # (scaffold, target) -> [Vmax reps], [EC50 reps]
    for r in rows:
        vs = [float(r[k]) for k in ("Vmax_1", "Vmax_2", "Vmax_3") if r[k]]
        es = [float(r[k]) for k in ("EC50_1", "EC50_2", "EC50_3") if r[k]]
        if vs and es:
            cells.setdefault((r["scaffold_desc"], r["target"]), [[], []])
            cells[(r["scaffold_desc"], r["target"])][0] += vs
            cells[(r["scaffold_desc"], r["target"])][1] += es
    if len(cells) < 8:
        raise SystemExit("数据不足: 需要至少 8 个(骨架,靶标)组合有重复测量, "
                         "当前 %d" % len(cells))
    scafs = sorted({k[0] for k in cells})
    tgts = sorted({k[1] for k in cells})
    a, b = len(scafs), len(tgts)
    out = {}
    for mi, metric in enumerate(("Vmax", "EC50")):
        ns = [len(cells.get((s, t), [[], []])[mi]) for s in scafs for t in tgts]
        n = min(ns)
        if n == 0 or set(ns) != {n}:
            print("[%s] 跳过: 组合重复数不平衡 %s" % (metric, ns))
            continue
        y = np.array([[cells[(s, t)][mi] for t in tgts] for s in scafs])  # a×b×n
        gm = y.mean()
        ss_t = ((y - gm) ** 2).sum()
        ss_a = b * n * ((y.mean(axis=(1, 2)) - gm) ** 2).sum()
        ss_b = a * n * ((y.mean(axis=(0, 2)) - gm) ** 2).sum()
        ss_ab = (n * ((y.mean(axis=2) - y.mean(axis=(1, 2))[:, None]
                       - y.mean(axis=(0, 2))[None, :] + gm) ** 2)).sum()
        ss_e = ss_t - ss_a - ss_b - ss_ab
        df_a, df_b, df_ab, df_e = a - 1, b - 1, (a - 1) * (b - 1), a * b * (n - 1)
        res = {}
        for name, ss, df in (("scaffold", ss_a, df_a), ("target", ss_b, df_b),
                             ("interaction", ss_ab, df_ab)):
            ms = ss / df
            f_stat = ms / (ss_e / df_e)
            res[name] = {"SS": round(ss, 4), "df": df, "F": round(f_stat, 3),
                         "p": float(fdist.sf(f_stat, df, df_e))}
        res["error"] = {"SS": round(ss_e, 4), "df": df_e}
        out[metric] = res
        print("\n=== %s 两因素方差分析(a=%d 骨架, b=%d 靶标, n=%d 重复) ===" % (
            metric, a, b, n))
        for name in ("scaffold", "target", "interaction"):
            r = res[name]
            print("  %-12s SS=%10.2f df=%2d F=%8.3f p=%.4g %s" % (
                name, r["SS"], r["df"], r["F"], r["p"],
                "**" if r["p"] < 0.01 else ("*" if r["p"] < 0.05 else "")))
    if not out:
        raise SystemExit("无可分析指标")
    dst = os.path.join(DATA, "ivt_matrix_anova.json")
    json.dump({"input": path, "anova": out,
               "reading": "interaction p<0.05 => 骨架最优解随靶标变化(分型假说成立)"
                          "; 否则交互弱, 结论退化为通用型最优"},
              open(dst, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n结果 -> %s" % dst)


def twin_check(path, cell_csv=None):
    """数字孪生预测 vs 实测校验(V3 §4.3: 实验前预测/实验后校验)。

    用已填 IVT 的 (Vmax, EC50) + CCLE 丰度 + Hill 阈值模型出各骨架杀伤预测;
    若给 --cell-csv(scaffold_desc,target,SI_measured) 则报 Spearman/RMSE。
    """
    prm = json.load(open(os.path.join(DATA, "virtual_cell_params.json"),
                         encoding="utf-8"))
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    meas = {}
    for r in rows:
        vs = [float(r[k]) for k in ("Vmax_1", "Vmax_2", "Vmax_3") if r[k]]
        es = [float(r[k]) for k in ("EC50_1", "EC50_2", "EC50_3") if r[k]]
        if vs and es:
            meas[(r["scaffold_desc"], r["target"])] = (mean(vs), mean(es))
    if not meas:
        raise SystemExit("已填 CSV 无可用 (Vmax, EC50) 行")
    if prm.get("ec50_tpm") is None or prm.get("vmax") is None:
        print("虚拟细胞参数 UN-CALIBRATED(ec50_tpm/vmax 为 null):")
        print("  回填 data/virtual_cell_params.json 后本命令输出杀伤预测;")
        print("  当前 IVT 已测 %d 组合, 参数标定材料齐备时即可运行。" % len(meas))
        return
    ab = json.load(open(os.path.join(DATA, "virtual_cell_abundance.json"),
                        encoding="utf-8"))
    n_hill = float(prm.get("n_hill", 2.0))

    def hill(rna, ec50):
        return rna ** n_hill / (rna ** n_hill + ec50 ** n_hill)

    pred = {}
    for (sc, tg), (vmax, ec50) in meas.items():
        gene = "TP53" if tg.startswith("TP53") else "KRAS"
        mut_tpm = ab.get(gene, {}).get("mut_tpm_p50", 0) or 0
        wt_tpm = ab.get(gene, {}).get("wt_tpm_p50", 0) or 0
        pred[(sc, tg)] = {"SI_pred": round(
            (vmax * hill(mut_tpm, ec50) + 1e-9) /
            (vmax * hill(wt_tpm, ec50) + 1e-9), 3)}
    print("杀伤预测(SI=突变型/野生型激活比):")
    for k, v in sorted(pred.items()):
        print("  %-22s %-12s SI_pred=%.3f" % (k[0], k[1], v["SI_pred"]))
    if cell_csv:
        from scipy.stats import spearmanr

        got = {}
        with open(cell_csv, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                got[(r["scaffold_desc"], r["target"])] = float(r["SI_measured"])
        common = sorted(set(got) & set(pred))
        if len(common) < 3:
            raise SystemExit("重叠样本 %d < 3, 无法校验" % len(common))
        p = [pred[k]["SI_pred"] for k in common]
        g = [got[k] for k in common]
        rho = spearmanr(p, g).statistic
        rmse = (sum((x - y) ** 2 for x, y in zip(p, g)) / len(p)) ** 0.5
        print("\n对比 %d 样本: Spearman=%.3f RMSE=%.3f" % (len(common), rho, rmse))
        json.dump({"pred": {"/".join(k): v for k, v in pred.items()},
                   "measured": {"/".join(k): v for k, v in got.items()},
                   "spearman": float(rho), "rmse": round(rmse, 3)},
                  open(os.path.join(DATA, "virtual_cell_check.json"), "w",
                       encoding="utf-8"), indent=1)
    else:
        print("\n(未给 --cell-csv: 仅输出预测; 细胞 SI 数据到手后加 "
              "--cell-csv 出 Spearman/RMSE)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", default=None, help="已填 CSV 路径")
    ap.add_argument("--out-round1", default=None, help="校验后聚合表输出")
    ap.add_argument("--anova", default=None,
                    help="已填 CSV -> 骨架×靶标两因素方差分析(交互项)")
    ap.add_argument("--twin-check", default=None,
                    help="已填 CSV -> 数字孪生预测(+可选 --cell-csi 对比)")
    ap.add_argument("--cell-csv", default=None,
                    help="细胞 SI 实测(scaffold_desc,target,SI_measured)")
    ap.add_argument("--out", default=os.path.join(DATA, "ivt_round1_template.csv"))
    args = ap.parse_args()
    if args.check:
        check_filled(args.check, args.out_round1)
    elif args.anova:
        two_way_anova(args.anova)
    elif args.twin_check:
        twin_check(args.twin_check, args.cell_csv)
    else:
        gen_template(args.out)


if __name__ == "__main__":
    main()
