"""Han 2025 (Nat Commun, PMC12508434) LbCas12a DR 变体数据提取与特征计算。

数据来源(全部已下载入库):
  MOESM3 XLSX: 7 条工具箱变体序列(CN/F1/F2/L1/L2/FL1/FL2)
  MOESM7 XLSX Fig.1 列63-64: 96 变体筛选活性(归一化 GFP, 低=抑制强)
  MOESM7 XLSX Fig.3b: 7 条工具箱两档表达(RBS0/RBS33)活性(n=3)
本脚本:
  1. 解析 MOESM3 提取 7 条 DR 序列;
  2. 对每条跑本管线特征(DDR-alone 折叠/ΔΔG/bp_dist/cross_pairs/stem_intact_prob/系综);
  3. 与 Fig.3b 活性做 Spearman 相关;
  4. 全部 97 条筛选活性值存 data/han2025_screen.json。
诚实边界: LbCas12a CRISPRi(GFP 抑制)体系, 非 Cas12a2 RNA 触发杀伤;
  n=7 有序列为工具箱级(非全 96 变体), 全量序列需补充 PDF(MOESM1)。
运行: python scripts/crrna_han2025_features.py
"""
import argparse
import json
import os
import sys

import numpy as np
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402

DATA = os.path.join(ROOT, "data")

# Fig.3b 工具箱活性(均值, 低=抑制强, 来自源数据手动提取)
TOOLBOX_ACTIVITY = {
    "CN": {"rbs0": 36.5, "rbs33": 50.9},
    "F1": {"rbs0": 29.6, "rbs33": 51.9},
    "F2": {"rbs0": 36.0, "rbs33": 59.5},
    "FL1": {"rbs0": 82.2, "rbs33": 92.3},
    "FL2": {"rbs0": 79.3, "rbs33": 87.1},
    "L1": {"rbs0": 82.8, "rbs33": 98.6},
    "L2": {"rbs0": 38.5, "rbs33": 67.8},
}


def parse_sequences():
    wb = openpyxl.load_workbook(
        os.path.join(DATA, "41467_2025_64010_MOESM3_ESM.xlsx"), data_only=True)
    ws = wb["Sheet1"]
    seqs = {}
    for r in range(3, ws.max_row + 1):
        name = ws.cell(r, 1).value
        seq = ws.cell(r, 2).value
        if not name or not seq:
            continue
        name = str(name).strip()
        seq = str(seq).replace(" ", "").replace("T", "U").upper()
        # 只取 gfp-1 系列的(同一 spacer, 不同 DR)——保证单变量
        if name.startswith("gfp-1-"):
            tag = name.replace("gfp-1-", "")
            seqs[tag] = seq
    return seqs


def extract_screen():
    """Fig.1 列63-64 的 96 变体筛选活性"""
    wb = openpyxl.load_workbook(
        os.path.join(DATA, "41467_2025_64010_MOESM7_ESM.xlsx"),
        data_only=True, read_only=True)
    ws = wb["Fig. 1"]
    seen = {}
    for r in range(1, ws.max_row + 1):
        for c in range(55, min(65, ws.max_column + 1)):
            v = ws.cell(r, c).value
            if v is not None:
                s = str(v).strip()
                if s == "Canonical" or (len(s) <= 5 and s[0] in "SLF"
                                         and s[1:].isdigit() and s[1:] != ""):
                    for c2 in range(c + 1, min(c + 3, ws.max_column + 1)):
                        v2 = ws.cell(r, c2).value
                        if v2 is not None and isinstance(v2, (int, float)):
                            if s not in seen:
                                seen[s] = float(v2)
                            break
                    break
    return seen


def compute_features(dr_rna, spacer_rna):
    full = dr_rna + spacer_rna
    ss, mfe = core.fold(full)
    dr_ss, dr_mfe = core.fold(dr_rna)
    cross_nt, cross_pp = core.cross_pairs(full, len(dr_rna))
    inv_run = core.inv_max_run(full, len(dr_rna))
    _, spacer_up, seed_up = core.pf_stats(full, len(dr_rna), 7)
    stem_pairs = core.stem_pairs_of(dr_ss)
    p_fold = core.stem_intact_prob(full, stem_pairs) if stem_pairs else 1.0
    return {
        "mfe": round(mfe, 2), "ddG_dr": round(dr_mfe, 2),
        "bp_dist": None,  # 需 WT 基准, 主循环里算
        "cross_nt": cross_nt, "cross_pp": round(cross_pp, 4),
        "inv_max_run": inv_run,
        "spacer_up": round(spacer_up, 3), "seed_up": round(seed_up, 3),
        "p_fold": round(p_fold, 5),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spacer-len", type=int, default=20,
                    help="spacer 用于折叠的截取长度(默认 20)")
    args = ap.parse_args()

    seqs = parse_sequences()
    if "CN" not in seqs:
        raise SystemExit("MOESM3 解析失败: 未找到 gfp-1-CN")
    spacer = seqs["CN"][21:21 + args.spacer_len]  # DR 21nt 后取 spacer
    wt_dr = seqs["CN"][:21]
    print("WT DR(21nt): %s" % wt_dr)
    print("Spacer(前%dnt): %s" % (len(spacer), spacer))

    wt_dr_ss, wt_dr_mfe = core.fold(wt_dr)
    wt_full_ss, wt_mfe = core.fold(wt_dr + spacer)
    print("WT DR MFE=%.2f, 全长 MFE=%.2f" % (wt_dr_mfe, wt_mfe))

    rows = []
    for tag, full_seq in sorted(seqs.items()):
        dr = full_seq[:21]
        feats = compute_features(dr, spacer)
        feats["bp_dist"] = None  # 后填
        act = TOOLBOX_ACTIVITY.get(tag, {})
        rows.append({"tag": tag, "dr_rna": dr, "n_mut": sum(
            1 for a, b in zip(wt_dr, dr) if a != b), **feats, **act})

    # bp_dist vs WT
    for r in rows:
        v_ss, _ = core.fold(r["dr_rna"] + spacer)
        r["bp_dist"] = RNA_bp_distance(wt_full_ss, v_ss)

    print("\n%-6s %4s %8s %8s %6s %6s %7s %7s %8s %8s" % (
        "var", "nmut", "ddG_dr", "MFE", "bp_d", "x_nt", "P_fold", "sp_up",
        "RBS0", "RBS33"))
    for r in rows:
        print("%-6s %4d %8.2f %8.2f %6d %6d %7.4f %7.3f %8.1f %8.1f" % (
            r["tag"], r["n_mut"], r["ddG_dr"], r["mfe"], r["bp_dist"],
            r["cross_nt"], r["p_fold"], r["spacer_up"],
            r.get("rbs0", 0), r.get("rbs33", 0)))

    # Spearman 相关(仅 7 点, 方向参考)
    acts0 = np.array([r.get("rbs0", np.nan) for r in rows])
    acts33 = np.array([r.get("rbs33", np.nan) for r in rows])
    for feat in ("ddG_dr", "bp_dist", "cross_nt", "p_fold", "spacer_up"):
        vals = np.array([float(r[feat]) for r in rows])
        mask = ~np.isnan(acts0)
        rho0 = spearman(vals[mask], acts0[mask])
        rho33 = spearman(vals[mask], acts33[mask])
        print("%-10s  rho(RBS0)=%+.3f  rho(RBS33)=%+.3f" % (feat, rho0, rho33))

    screen = extract_screen()
    json.dump({"toolbox_sequences": {t: s for t, s in seqs.items()},
               "toolbox_features": rows,
               "toolbox_activity": TOOLBOX_ACTIVITY,
               "screen_96": screen,
               "source": "Han et al. 2025 Nat Commun PMC12508434 (LbCas12a CRISPRi)"},
              open(os.path.join(DATA, "han2025_dataset.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n数据集 -> data/han2025_dataset.json (工具箱%d + 筛选%d)" % (
        len(rows), len(screen)))


def RNA_bp_distance(s1, s2):
    import RNA
    return int(RNA.bp_distance(s1, s2))


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


if __name__ == "__main__":
    main()
