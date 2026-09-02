"""Han 2025 (Nat Commun, s41467-025-64010-z) LbCas12a DR 变体数据提取与特征计算。

数据来源(全部已下载入库):
  MOESM3 XLSX: 工具箱 DR 序列(gfp-1 系 7 条 + P1 系 L3/L4)
  MOESM7 XLSX Fig.1 列63-64: Fig1g 活性表(归一化 GFP, 低=抑制强,
    2026-09-02 修正: 含 FL 组共 145 条, 旧版漏 48 条 FL)
  MOESM7 XLSX Fig.3b: 7 条工具箱两档表达(RBS0/RBS33)活性(n=3)
  data/han2025_sanger_sequences.json: MOESM1 Sup Fig.2/3 的 54 条
    Sanger 序列(600dpi 视觉转录+模板校验, 含置信度标注)
本脚本:
  1. 解析 MOESM3 提取 9 条工具箱 DR 序列(含 L3/L4);
  2. 对每条跑本管线特征(DDR-alone 折叠/ΔΔG/bp_dist/cross_pairs/stem_intact_prob/系综);
  3. Fig1g 145 条活性全表 + 工具箱 7 条 Fig1g 同尺度锚点值;
  4. 54 条 Sanger 耐受集特征计算(耐受先验, 无逐条活性标签);
  5. 与 Fig.3b 活性做 Spearman 相关(参考)。
诚实边界: LbCas12a CRISPRi(GFP 抑制)体系, 非 Cas12a2 RNA 触发杀伤;
  Fig1g 的 145 个编号(F/S/L/FL)与序列的映射未随论文发表,
  干净监督对 = 工具箱 7 条(Fig1g 同尺度); 96 变体序列不可得(2026-09-02 查证)。
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
        # gfp-1 系列(同一 spacer, 不同 DR)保证单变量;
        # L3/L4 只在 P1 系列出(取 DR 21nt, 同一 DR 跨系列一致)
        for prefix in ("gfp-1-", "P1-"):
            if name.startswith(prefix):
                tag = name.replace(prefix, "")
                if tag not in seqs:
                    seqs[tag] = seq
    return seqs


def extract_screen():
    """Fig.1 列63-64 的 Fig1g 活性表(含 FL 组, 共 145 条 + Canonical)"""
    import re
    wb = openpyxl.load_workbook(
        os.path.join(DATA, "41467_2025_64010_MOESM7_ESM.xlsx"),
        data_only=True, read_only=True)
    ws = wb["Fig. 1"]
    pat = re.compile(r"^(FL|F|S|L)\d+$")
    seen = {}
    for row in ws.iter_rows(min_row=1, min_col=63, max_col=64,
                            values_only=True):
        s, v = row
        if s is None or v is None:
            continue
        s = str(s).strip()
        if s == "Canonical" or pat.match(s):
            if s not in seen and isinstance(v, (int, float)):
                seen[s] = float(v)
    return seen


def compute_features(dr_rna, spacer_rna, wt_dr_mfe=0.0):
    full = dr_rna + spacer_rna
    ss, mfe = core.fold(full)
    dr_ss, dr_mfe = core.fold(dr_rna)
    cross_nt, cross_pp = core.cross_pairs(full, len(dr_rna))
    inv_run = core.inv_max_run(full, len(dr_rna))
    _, spacer_up, seed_up = core.pf_stats(full, len(dr_rna), 7)
    stem_pairs = core.stem_pairs_of(dr_ss)
    p_fold = core.stem_intact_prob(full, stem_pairs) if stem_pairs else 1.0
    return {
        "mfe": round(mfe, 2), "ddG_dr": round(dr_mfe - wt_dr_mfe, 2),
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
        feats = compute_features(dr, spacer, wt_dr_mfe)
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

    # === Fig1g 活性全表 + 工具箱同尺度锚点(程序化提取, 不硬编码) ===
    fig1g = extract_screen()
    tag_map = {"CN": "Canonical"}
    anchors_fig1g = {}
    for r in rows:
        key = tag_map.get(r["tag"], r["tag"])
        if key in fig1g:
            r["fig1g"] = fig1g[key]
            anchors_fig1g[r["tag"]] = fig1g[key]
    print("\nFig1g 全表 %d 条; 工具箱同尺度锚点 %d 条: %s" % (
        len(fig1g), len(anchors_fig1g), anchors_fig1g))

    # === Sanger 耐受集(MOESM1 Sup Fig.2/3, 54 条)特征计算 ===
    sanger_path = os.path.join(DATA, "han2025_sanger_sequences.json")
    sanger_rows = []
    if os.path.exists(sanger_path):
        sanger = json.load(open(sanger_path, encoding="utf-8"))
        for e in sanger["entries"]:
            dr = e["seq"]
            feats = compute_features(dr, spacer, wt_dr_mfe)
            v_ss, _ = core.fold(dr + spacer)
            feats["bp_dist"] = RNA_bp_distance(wt_full_ss, v_ss)
            sanger_rows.append({
                "label": e["label"], "group": e["group"],
                "dr_rna": dr, "confidence": e["confidence"],
                "reads": e.get("reads", 1),
                "matches_toolbox": e.get("matches_toolbox"),
                "n_mut": sum(1 for a, b in zip(wt_dr, dr) if a != b),
                **feats})
        print("Sanger 耐受集 %d 条特征计算完成" % len(sanger_rows))

    json.dump({"toolbox_sequences": {t: s for t, s in seqs.items()},
               "toolbox_features": rows,
               "toolbox_activity": TOOLBOX_ACTIVITY,
               "toolbox_activity_fig1g": anchors_fig1g,
               "fig1g_activity": fig1g,
               "sanger_tolerance": sanger_rows,
               "sanger_source": (sanger["source"] + "; " + sanger["extraction"]
                                 if sanger_rows else None),
               "source": "Han et al. 2025 Nat Commun s41467-025-64010-z "
                         "(LbCas12a CRISPRi)"},
              open(os.path.join(DATA, "han2025_dataset.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n数据集 -> data/han2025_dataset.json (工具箱%d + Fig1g %d + "
          "Sanger耐受%d)" % (len(rows), len(fig1g), len(sanger_rows)))


def RNA_bp_distance(s1, s2):
    import RNA
    return int(RNA.bp_distance(s1, s2))


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


if __name__ == "__main__":
    main()
