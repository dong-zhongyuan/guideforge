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
  5. 与 Fig.3b 活性做 Spearman 相关(参考);
  6. (2026-09-03) 主文 Fig.1f 序列表视觉转录: 从编号-序列桥中恢复
     9 条新监督对(见 FIG1F_TRANSCRIPTION 溯源), 干净监督对 7 -> 16。
诚实边界: LbCas12a CRISPRi(GFP 抑制)体系, 非 Cas12a2 RNA 触发杀伤;
  Fig1g 的 145 个编号(F/S/L/FL)与序列的映射未随论文数值源数据发表
  (2026-09-02 查证 MOESM1-7; 2026-09-03 复查 MOESM7 全部 25 个工作表 +
  MOESM3/MOESM4 全文, 数值侧确无映射)。主文 Fig.1f 以图像形式给出部分
  编号-序列对照(典型突变体表), 经视觉转录+双重独立读取+Sanger 交叉校验
  恢复 9 条新监督对(置信度分级见 fig1f_pairs); 其余 129 个编号仍不可得,
  完整映射需向作者索取。
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

# 主文 Fig.1f 序列表视觉转录(2026-09-03)。
# 溯源: nature.com 文章页 Fig1 全图 PNG(2017x2201) -> panel f 右侧编号表
#   区域裁切(scripts/tmp/fig1f/, gitignore) -> 分块 2x-3x 放大 -> 视觉模型
#   两次独立逐字母转录 + 与 han2025_sanger_sequences.json(MOESM1 Sup Fig.2/3
#   独立转录源)的 loop 六联体(11-16 位)交叉校验。
# 置信度: high = 两次完整独立读取一致 + loop 六联体与 Sanger 独立源吻合;
#   medium = 单次完整读取 + loop 六联体与 Sanger 独立源吻合。
# 边界: 其余编号(含 FL3/FL5, 不在 Fig1g 活性表内)未收录; 视觉转录存在
#   单字母误差风险(如 S16 七连 C), 训练使用时建议按置信度加权或敏感性剔除。
FIG1F_TRANSCRIPTION = {
    "S3":  {"dr_rna": "UAAUUUCUACUCAUAGUAGAU", "confidence": "high",
            "cross_check": "loop UCAUAG = Sanger B19"},
    "S5":  {"dr_rna": "UAAUUUCUACUGUUUGUAGAU", "confidence": "high",
            "cross_check": "loop UGUUUG = Sanger A36/B19b"},
    "S7":  {"dr_rna": "UAAUUUCUACCCUAAGUAGAU", "confidence": "medium",
            "cross_check": "loop CCUAAG = Sanger B24(二读截边, 单完整读)"},
    "S10": {"dr_rna": "UAAUUUCUACCCGAGGUAGAU", "confidence": "medium",
            "cross_check": "loop CCGAGG = Sanger C48/C37"},
    "S14": {"dr_rna": "UAAUUUCUACGGGCGGUAGAU", "confidence": "medium",
            "cross_check": "loop GGGCGG = Sanger B58/D47"},
    "S16": {"dr_rna": "UAAUUUCUACCCCCCGUAGAU", "confidence": "medium",
            "cross_check": "loop CCCCCG = Sanger F33(完全一致, 同聚段易误读已复核)"},
    "S20": {"dr_rna": "UAAUUUCUACCAAUGGUAGAU", "confidence": "medium",
            "cross_check": "loop CAAUGG = Sanger B29"},
    "S21": {"dr_rna": "UAAUUUCUACUGGAUGUAGAU", "confidence": "medium",
            "cross_check": "loop UGGAUG = Sanger A40"},
    "FL4": {"dr_rna": "CAAUGUCUACCUUCUGUAGAU", "confidence": "medium",
            "cross_check": "flank CAAUGU+loop CUUCUG = Sanger B38"},
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


def extract_cleavage(sheet="Sup Fig. 5", col0=1):
    """Sup Fig.5(cis)/6(trans, b 面板) 实时切割曲线终点/起点比(三重复均值)。

    返回 {标签: 终点比}; 标签 Canonical -> CN。b 面板列布局: 7 个 crRNA x 3 重复。
    """
    wb = openpyxl.load_workbook(
        os.path.join(DATA, "41467_2025_64010_MOESM7_ESM.xlsx"),
        data_only=True, read_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c) if c is not None else "" for c in rows[1]]
    out = {}
    for i, h in enumerate(hdr):
        if h in ("Canonical", "L1", "L2", "FL1", "FL2", "F1", "F2"):
            tri = [j for j in range(i, min(i + 3, len(hdr)))
                   if rows[2][j] is not None]
            first = np.mean([float(rows[2][j]) for j in tri])
            last = np.mean([float(rows[-1][j]) for j in tri])
            tag = "CN" if h == "Canonical" else h
            out[tag] = round(float(last / max(first, 1e-9)), 3)
    return out


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

    # === Fig.1f 视觉转录配对(2026-09-03): 编号->序列->Fig1g 活性 ===
    fig1f_pairs = []
    for tag, e in sorted(FIG1F_TRANSCRIPTION.items()):
        dr = e["dr_rna"]
        assert len(dr) == 21 and set(dr) <= set("AUCG"), tag
        feats = compute_features(dr, spacer, wt_dr_mfe)
        v_ss, _ = core.fold(dr + spacer)
        feats["bp_dist"] = RNA_bp_distance(wt_full_ss, v_ss)
        fig1f_pairs.append({
            "tag": tag, "dr_rna": dr, "confidence": e["confidence"],
            "cross_check": e["cross_check"],
            "n_mut": sum(1 for a, b in zip(wt_dr, dr) if a != b),
            "fig1g": fig1g.get(tag), **feats})
    new_tags = [p["tag"] for p in fig1f_pairs if p["fig1g"] is not None]
    print("Fig1f 转录配对 %d 条, 其中带 Fig1g 活性 %d 条: %s" % (
        len(fig1f_pairs), len(new_tags), new_tags))

    # === 同源多终点: cis(Sup Fig.5b) / trans(Sup Fig.6b) 切割终点比 ===
    cis = extract_cleavage("Sup Fig. 5")
    trans = extract_cleavage("Sup Fig. 6")
    for r in rows:
        if r["tag"] in cis:
            r["cis_end"] = cis[r["tag"]]
        if r["tag"] in trans:
            r["trans_end"] = trans[r["tag"]]
    print("cis 终点比: %s" % cis)
    print("trans 终点比: %s" % trans)

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

    out = {"toolbox_sequences": {t: s for t, s in seqs.items()},
           "toolbox_features": rows,
           "toolbox_activity": TOOLBOX_ACTIVITY,
           "toolbox_activity_fig1g": anchors_fig1g,
           "toolbox_cleavage": {"cis_end": cis, "trans_end": trans},
           "fig1g_activity": fig1g,
           "fig1f_pairs": fig1f_pairs,
           "fig1f_source": (
               "主文 Fig.1f 序列表视觉转录(2026-09-03): nature.com Fig1 全图 "
               "panel f 右侧编号表, 分块放大两次独立逐字母读取 + Sanger "
               "独立转录源 loop 六联体交叉校验; 置信度分级见各条 confidence; "
               "其余 129 个编号的序列映射仍不可得(数值源数据无, 需向作者索取)"),
           "sanger_tolerance": sanger_rows,
           "sanger_source": (sanger["source"] + "; " + sanger["extraction"]
                             if sanger_rows else None),
           "source": "Han et al. 2025 Nat Commun s41467-025-64010-z "
                     "(LbCas12a CRISPRi)"}
    # 合并写: 保留 rejected_endpoints 等后加键, 避免整体重写丢数据
    ds_path = os.path.join(DATA, "han2025_dataset.json")
    if os.path.exists(ds_path):
        try:
            old = json.load(open(ds_path, encoding="utf-8"))
            if isinstance(old, dict):
                old.update(out)
                out = old
        except (ValueError, OSError):
            pass
    json.dump(out, open(ds_path, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\n数据集 -> data/han2025_dataset.json (工具箱%d + Fig1g %d + "
          "Fig1f配对%d + Sanger耐受%d)" % (
              len(rows), len(fig1g), len(fig1f_pairs), len(sanger_rows)))


def RNA_bp_distance(s1, s2):
    import RNA
    return int(RNA.bp_distance(s1, s2))


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


if __name__ == "__main__":
    main()
