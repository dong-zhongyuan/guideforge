"""靶RNA丰度档位特征(策划案 V3 §4.4 模块一特征向量 / 表2 panel 选择标准, 2026-09-05)。

从 CCLE/DepMap 18q3 RPKM gct 实测提取 panel 基因(TP53/KRAS/APC)在对应
细胞系(HCT116/SW480)的表达量, 按显式阈值规则换算丰度档位, 落盘
data/target_abundance_tiers.json; 模块一(靶标解析)特征向量与模块三
(选型/推荐链路)消费同一文件, 阈值唯一定义于本模块, 不允许下游各自换算。

分档规则(显式声明, 可复现):
  quantity  = mut_rpkm_het50 = total_rpkm x 0.5
              (杂合突变场景等位因子 0.5, 与 scripts/crrna_virtual_cell.py 同口径)
  tier 0 低  : mut_rpkm_het50 <  5.2   (低于标定激活阈 95%CI 下界)
  tier 1 中  : 5.2 <= mut_rpkm_het50 < 17.4  (激活阈附近)
  tier 2 高  : mut_rpkm_het50 >= 17.4  (高于标定激活阈 95%CI 上界)
  阈值锚点: Scholz et al. 2026 Nature Fig 1h 实测 EC50=11.31 FPKM,
  95%CI [5.2, 17.4](data/scholz2026_fig1h_calibration.json);
  RPKM≈FPKM 跨源等效为近似(口径同 virtual_cell qa.cross_source_caveat)。
  换表达数据源或校准更新时须重标定阈值并同步本文件与 JSON。

运行: python scripts/crrna_target_abundance.py
输出: data/target_abundance_tiers.json
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

CCLE_GCT = os.path.join(DATA, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")
CALIB_JSON = os.path.join(DATA, "scholz2026_fig1h_calibration.json")
OUT_JSON = os.path.join(DATA, "target_abundance_tiers.json")

GENES = ("TP53", "KRAS", "APC")
CELL_LINES = ("HCT116", "SW480")

ALLELE_FACTOR_HET = 0.5   # 杂合突变场景, 与 crrna_virtual_cell.py 统一 50% 同口径
TIER_LO, TIER_HI = 5.2, 17.4   # Scholz 2026 Fig1h EC50 95%CI(见模块 docstring)
TIER_LABELS = {0: "低(低于激活阈)", 1: "中(近激活阈)", 2: "高(超激活阈)"}

# panel 靶标 -> 代表细胞系(V3 表2 来源癌种列 + 仓库既有验证设计);
# 细胞系突变状态以 DepMap mutation calls 复核(仓库统一告诫, 见 virtual_cell note)
PANEL_TARGETS = {
    "tp53_r248q": {"gene": "TP53", "cell_line": "HCT116",
                   "note": "V3 表2 主验证靶标, 来源癌种列=结直肠癌(HCT116)"},
    "kras_g12d": {"gene": "KRAS", "cell_line": "SW480",
                  "note": "V3 表2 泛化验证靶标, 来源癌种列=胰腺/结直肠(SW480等)"},
    "tp53_r273h": {"gene": "TP53", "cell_line": "SW480",
                   "note": "V3 表2 同基因不同位点; SW480 为 TP53 R273H 突变结直肠系"
                           "(突变状态以 DepMap mutation calls 复核)"},
    "apc_q1328x": {"gene": "APC", "cell_line": "HCT116",
                   "note": "设计无义突变(非特定 COSMIC 条目, V3 表2 来源=结直肠癌); "
                           "代表系取主验证系 HCT116, SW480 值见 genes 块"},
}


def tier_of(mut_rpkm_het50):
    """显式分档规则(唯一实现): 阈值 TIER_LO/TIER_HI, 见模块 docstring。"""
    if mut_rpkm_het50 < TIER_LO:
        return 0
    if mut_rpkm_het50 < TIER_HI:
        return 1
    return 2


def parse_gct(gct_path=CCLE_GCT, genes=GENES, cell_lines=CELL_LINES):
    """解析 CCLE gct, 返回 {gene: {"ensembl":..., "lines": {line: (ach, rpkm)}}}。
    只取目标基因行与目标细胞系列; 数值为文件原始浮点, 不舍入。"""
    out = {}
    with open(gct_path, "rt", encoding="utf-8", errors="ignore") as f:
        header = None
        for line in f:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if header is None:
                if parts[0] == "Name":
                    header = parts
                    col = {}
                    for i, h in enumerate(header[2:]):
                        # 列头形如 HCT116_LARGE_INTESTINE (ACH-000971):
                        # 取系名(首个 "_" 前段)精确匹配, 避免 SW48/SW480 混淆
                        name = h.split(" (")[0].split("_")[0]
                        if name in cell_lines and name not in col:
                            ach = h.split(" (")[1].rstrip(")") if " (" in h else ""
                            col[name] = (i, ach)
                    missing = set(cell_lines) - set(col)
                    if missing:
                        raise SystemExit("gct 缺目标细胞系列: %s" % sorted(missing))
                continue
            sym = parts[1] if len(parts) > 1 else ""
            if sym in genes and sym not in out:
                out[sym] = {"ensembl": parts[0], "lines": {}}
                for cl in cell_lines:
                    i, ach = col[cl]
                    out[sym]["lines"][cl] = (ach, float(parts[i + 2]))
    missing = set(genes) - set(out)
    if missing:
        raise SystemExit("gct 缺目标基因行: %s" % sorted(missing))
    return out


def _line_record(total_rpkm):
    mut = total_rpkm * ALLELE_FACTOR_HET
    t = tier_of(mut)
    return {"total_rpkm": total_rpkm,
            "mut_rpkm_het50": round(mut, 5),
            "abundance_tier": t, "tier_label": TIER_LABELS[t]}


def build(gct_path=CCLE_GCT):
    """构建丰度档位表(全部数值来自 gct 实测, 档位由 tier_of 换算)。"""
    raw = parse_gct(gct_path)
    genes = {}
    for g in GENES:
        genes[g] = {"ensembl": raw[g]["ensembl"],
                    "lines": {cl: dict(zip(("ach",), (raw[g]["lines"][cl][0],)),
                                       **_line_record(raw[g]["lines"][cl][1]))
                              for cl in CELL_LINES},
                    "panel_targets": [k for k, v in PANEL_TARGETS.items()
                                      if v["gene"] == g]}
    panel = {}
    for key, spec in PANEL_TARGETS.items():
        rec = genes[spec["gene"]]["lines"][spec["cell_line"]]
        panel[key] = {"gene": spec["gene"], "cell_line": spec["cell_line"],
                      "ach": rec["ach"], "total_rpkm": rec["total_rpkm"],
                      "mut_rpkm_het50": rec["mut_rpkm_het50"],
                      "abundance_tier": rec["abundance_tier"],
                      "tier_label": rec["tier_label"], "note": spec["note"]}
    # 判读文字: 由数值按显式规则生成(排序+阈值比较), 不手写结论
    order = sorted(panel, key=lambda k: -panel[k]["mut_rpkm_het50"])
    reading = ["panel 靶RNA丰度(杂合场景 mut RPKM)降序: " + " > ".join(
        "%s %.2f(%s)" % (k, panel[k]["mut_rpkm_het50"], panel[k]["tier_label"])
        for k in order)]
    for k in order:
        p = panel[k]
        if p["abundance_tier"] == 0:
            reading.append("%s: %.2f < %.1f -> 低于标定激活阈 CI 下界" % (
                k, p["mut_rpkm_het50"], TIER_LO))
        elif p["abundance_tier"] == 1:
            reading.append("%s: %.1f <= %.2f < %.1f -> 激活阈附近" % (
                k, TIER_LO, p["mut_rpkm_het50"], TIER_HI))
        else:
            reading.append("%s: %.2f >= %.1f -> 高于标定激活阈 CI 上界" % (
                k, p["mut_rpkm_het50"], TIER_HI))
    return {
        "source": {
            "file": os.path.basename(gct_path),
            "release": "CCLE/DepMap 18q3 RNAseq RPKM (2018-07-18)",
            "units": "RPKM",
            "columns": {"HCT116": "HCT116_LARGE_INTESTINE (ACH-000971)",
                        "SW480": "SW480_LARGE_INTESTINE (ACH-000842)"},
            "note": "2018q3 RPKM 旧版口径; 换源(如 DepMap 新版 TPM)须重标定分档阈值"},
        "tier_rule": {
            "quantity": "mut_rpkm_het50 = total_rpkm x allele_factor",
            "allele_factor": ALLELE_FACTOR_HET,
            "allele_factor_meaning": "杂合突变场景(突变等位基因=总丰度 50%), "
                                     "与 scripts/crrna_virtual_cell.py 同口径",
            "edges_mut_rpkm": [TIER_LO, TIER_HI],
            "tiers": {"0": "mut < %.1f: %s" % (TIER_LO, TIER_LABELS[0]),
                      "1": "%.1f <= mut < %.1f: %s" % (TIER_LO, TIER_HI,
                                                       TIER_LABELS[1]),
                      "2": "mut >= %.1f: %s" % (TIER_HI, TIER_LABELS[2])},
            "anchor": "Scholz et al. 2026 Nature Fig 1h 实测 EC50=11.31 FPKM, "
                      "95%%CI [%.1f, %.1f](data/scholz2026_fig1h_calibration.json); "
                      "RPKM≈FPKM 跨源等效为近似" % (TIER_LO, TIER_HI),
            "implementation": "scripts/crrna_target_abundance.py:tier_of(唯一换算点)"},
        "genes": genes,
        "panel_targets": panel,
        "panel_reading": reading,
    }


def load_abundance(path=OUT_JSON):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def target_feature(abn, target_key):
    """模块一靶标解析特征向量的丰度维: 返回 None(未知靶标)或特征 dict。"""
    p = (abn or {}).get("panel_targets", {}).get(target_key)
    if not p:
        return None
    return {"target_key": target_key, "gene": p["gene"],
            "cell_line": p["cell_line"], "total_rpkm": p["total_rpkm"],
            "mut_rpkm_het50": p["mut_rpkm_het50"],
            "target_abundance_tier": p["abundance_tier"],
            "target_abundance_tier_label": p["tier_label"]}


def activation_gate(feat):
    """模块三推荐链路对丰度档位的消费: 显式阈值规则生成判读, 不手写结论。

    返回 {"verdict":..., "note":...}; verdict 仅由 tier 决定:
      0 below_ci_lo / 1 within_ci / 2 above_ci_hi。"""
    mut = feat["mut_rpkm_het50"]
    t = feat["target_abundance_tier"]
    ctx = "%s@%s" % (feat["gene"], feat["cell_line"])
    if t == 0:
        return {"verdict": "below_ci_lo",
                "note": "靶RNA丰度 %s 杂合场景 %.2f RPKM < 标定激活阈 CI 下界 %.1f: "
                        "胞内激活预计不足, 骨架选型仅结构口径, 细胞层建议提高 "
                        "crRNA 剂量或改选高丰度系" % (ctx, mut, TIER_LO)}
    if t == 1:
        return {"verdict": "within_ci",
                "note": "靶RNA丰度 %s 杂合场景 %.2f RPKM 落在标定 EC50 95%%CI "
                        "[%.1f, %.1f] 内: 激活阈附近, 杀伤窗口对剂量与骨架增益敏感"
                        % (ctx, mut, TIER_LO, TIER_HI)}
    return {"verdict": "above_ci_hi",
            "note": "靶RNA丰度 %s 杂合场景 %.2f RPKM > 标定激活阈 CI 上界 %.1f: "
                    "激活条件充裕" % (ctx, mut, TIER_HI)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gct", default=CCLE_GCT)
    ap.add_argument("--out", default=OUT_JSON)
    args = ap.parse_args()

    out = build(args.gct)
    # 阈值锚点与标定文件交叉核验(防阈值手写漂移)
    if os.path.isfile(CALIB_JSON):
        cal = json.load(open(CALIB_JSON, encoding="utf-8"))
        ci = cal.get("ec50_ci95")
        assert ci and abs(ci[0] - TIER_LO) < 1e-9 and abs(ci[1] - TIER_HI) < 1e-9, \
            "分档阈值 %s 与标定文件 CI %s 不一致, 须重标定" % ((TIER_LO, TIER_HI), ci)
    json.dump(out, open(args.out, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("丰度档位表 -> %s" % args.out)
    for key in PANEL_TARGETS:
        p = out["panel_targets"][key]
        print("  %-12s %s@%-7s total=%.5f mut(het50)=%.2f -> tier %d %s" % (
            key, p["gene"], p["cell_line"], p["total_rpkm"],
            p["mut_rpkm_het50"], p["abundance_tier"], p["tier_label"]))
    for line in out["panel_reading"]:
        print("  [判读] %s" % line)


if __name__ == "__main__":
    main()
