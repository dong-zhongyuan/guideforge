"""虚拟细胞层: 杀伤过程数字孪生框架(策划案 V3 §4.3, 2026-09-02)。

数据源(运行时下载+缓存, 不硬编码数值):
  Human Protein Atlas RNA 细胞系 TPM(proteinatlas.org 可达已探明);
  DepMap 门户可达但批量矩阵下载(figshare)403, 留作校准期替换源。

模型(Hill 激活阈值框架, 参数留接口待 IVT 实测回填——未标定前只输出
参数扫描包络, 不给点预测, 不写活性结论):
  activation = [RNA]^n / (EC50^n + [RNA]^n) × crRNA_rel
  杀伤窗口 = 靶基因高表达系激活率 - 低表达系(如同基因 WT 背景/低丰度系)激活率

运行: python scripts/crrna_virtual_cell.py [--genes TP53,KRAS,APC,MYC]
输出: data/virtual_cell_abundance.json(丰度分档表)
      data/virtual_cell_envelope.json(EC50 网格上的杀伤窗口包络, 预校准)
      data/virtual_cell_params.json(参数接口, 待 IVT 回填: ec50_tpm/vmax/crrna_rel)
"""
import argparse
import csv
import io
import json
import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

HPA_TSV = ("https://www.proteinatlas.org/api/search_download.php"
           "?search=%s&format=tsv&columns=g,rnacs"
           "&compress=no")

PARAM_INTERFACE = {
    "ec50_tpm": None,      # 待 IVT: 靶 RNA 激活阈值(TPM), 未知
    "vmax": None,          # 待 IVT: 旁切饱和速率
    "n_hill": 2.0,         # Hill 系数(默认 2, 待标定)
    "crrna_rel": None,     # 待测: 细胞内 crRNA 相对丰度(qPCR 归一)
    "status": "UN-CALIBRATED: 仅参数包络输出, 不构成活性预测",
}
GENE2TARGETS = {"TP53": ["R248Q", "R273H"], "KRAS": ["G12C", "G12D"],
                "APC": ["MCR-nonsense"], "MYC": ["g1"]}


def fetch_hpa(gene, cache):
    if os.path.isfile(cache):
        return open(cache, encoding="utf-8").read()
    req = urllib.request.Request(HPA_TSV % gene,
                                 headers={"User-Agent": "guideforge/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        txt = r.read().decode("utf-8", "ignore")
    open(cache, "w", encoding="utf-8").write(txt)
    return txt


def load_local_ccle(genes):
    """若 data/ 下有完整 CCLE 18q3 gct(手动下载放入, 带宽受限自动抓取不通),
    就地解析四基因行; 返回 {gene: {cell_line: RPKM}}。"""
    path = os.path.join(DATA, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")
    if not os.path.isfile(path):
        return None
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt", encoding="utf-8", errors="ignore") as f:
        header = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if header is None:
                if parts[0] == "Name":
                    header = parts
                continue
            sym = parts[1] if len(parts) > 1 else ""
            if sym in genes:
                vals = [float(x) if x not in ("", "NA") else 0.0 for x in parts[2:]]
                out[sym] = {header[i + 2].split(" (")[0]: v
                            for i, v in enumerate(vals)}
    return out or None


def prior_lit_mode(ec50_grid):
    """文献实测标定版(2026-09-03 升级: 从数量级场景改为 Scholz 2026 实测剂量曲线)。

    标定源: Scholz et al. 2026 Nature (s41586-026-10466-y) Fig 1h 源数据
    (data/scholz2026_fig1h_calibration.json): 4 细胞系 x 7 靶转录本的
    (FPKM, 存活率) 15 点, Hill 拟合 EC50=11.3 FPKM (95%CI 5.2-17.4),
    有效 Hill h=1.78, R2=0.867——胞内杀伤半数阈值的实测锚点。
    分子层注记: Kunwar 2026 (bioRxiv v2) SEC 纯化后结合无双曲线协同性
    (binary-靶 RNA Kd=18±4 nM), 细胞层 h=1.78 为有效剂量斜率非分子协同。
    移植边界: 标定体系为 RNP 电转; 本项目质粒/U6 表达体系的有效 RNP 剂量
    不同, crrna_rel 保留场景因子 {0.3, 1, 3}。
    """
    cal = json.load(open(os.path.join(
        DATA, "scholz2026_fig1h_calibration.json"), encoding="utf-8"))
    ec50, h = cal["ec50_fpkm"], cal["hill"]
    ec50_lo, ec50_hi = cal["ec50_ci95"]
    scenarios = [ec50_lo, ec50, ec50_hi]

    abundance = json.load(open(os.path.join(DATA, "virtual_cell_abundance.json"),
                               encoding="utf-8"))["genes"]
    gct = os.path.join(DATA, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")
    named = {}
    if os.path.isfile(gct):
        import gzip
        op = gzip.open if gct.endswith(".gz") else open
        want = ("HCT116", "SW480", "SW620", "AsPC", "BxPC", "PANC1", "PANC-1")
        with op(gct, "rt", encoding="utf-8", errors="ignore") as f:
            header = None
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    if parts[0] == "Name":
                        header = parts
                    continue
                sym = parts[1] if len(parts) > 1 else ""
                if sym not in abundance:
                    continue
                for i, v in enumerate(parts[2:]):
                    col = header[i + 2]
                    if any(w in col for w in want):
                        try:
                            named.setdefault(col.split(" (")[0], {})[sym] = float(v)
                        except ValueError:
                            pass

    def survival(rpkm, e50):
        # Scholz 拟合: 存活% = 平台 + 跨度/(1+(F/EC50)^h); F≈RPKM(同类长度归一)
        return cal["surv_min_pct"] + cal["surv_span_pct"] / (
            1.0 + (max(rpkm, 1e-9) / e50) ** h)

    lines_tbl = []
    for cl, genes in sorted(named.items()):
        row = {"cell_line": cl, "rpkm": genes}
        for lab, e50 in zip(("lo", "mid", "hi"), scenarios):
            row["surv_tp53_%s" % lab] = round(survival(
                genes.get("TP53", 0.0) * 0.5, e50), 1)  # 突变本按杂合 50%
            row["surv_kras_%s" % lab] = round(survival(
                genes.get("KRAS", 0.0) * 0.5, e50), 1)
        lines_tbl.append(row)
    out = {
        "status": "文献实测标定(Scholz 2026 Fig 1h, n=15 剂量点, R2=0.867); "
                  "移植边界=RNP电转体系, 质粒体系 crrna_rel 场景因子未定; "
                  "IVT/细胞实测回填 virtual_cell_params.json 后 --twin-check 校验",
        "calibration": {k: cal[k] for k in (
            "ec50_fpkm", "hill", "surv_min_pct", "surv_span_pct", "r2",
            "ec50_ci95", "n_points", "source")},
        "scenario_ec50_fpkm": {"lo": ec50_lo, "mid": ec50, "hi": ec50_hi},
        "crrna_rel_scenarios": [0.3, 1.0, 3.0],
        "named_cell_lines": lines_tbl}
    p = os.path.join(DATA, "virtual_cell_prior.json")
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("文献实测标定场景 -> %s" % p)
    print("  EC50=%.1f FPKM (CI %.1f-%.1f), h=%.2f" % (ec50, ec50_lo, ec50_hi, h))
    for r in lines_tbl:
        print("  %-16s TP53_RPKM=%-8s 预测存活(mid)=%5s%%  KRAS_RPKM=%-8s 存活=%5s%%" % (
            r["cell_line"], r["rpkm"].get("TP53"), r["surv_tp53_mid"],
            r["rpkm"].get("KRAS"), r["surv_kras_mid"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--genes", default="TP53,KRAS,APC,MYC")
    ap.add_argument("--ec50-grid", default="0.5,1,2,5,10,20,50",
                    help="EC50(TPM) 扫描网格, 预校准包络用")
    ap.add_argument("--prior-lit", action="store_true",
                    help="文献先验场景表(替代 UN-CALIBRATED 空转: 数量级口径)")
    args = ap.parse_args()

    if args.prior_lit:
        prior_lit_mode(args.ec50_grid)
        return

    abundance = {}
    genes = [g.strip() for g in args.genes.split(",")]
    local = load_local_ccle(set(genes))
    for gene in genes:
        if local and gene in local:
            vals = sorted(local[gene].values())
            abundance[gene] = {
                "n_cell_lines": len(vals),
                "median_tpm": round(vals[len(vals) // 2], 2),
                "p10_tpm": round(vals[int(0.1 * len(vals))], 2),
                "p90_tpm": round(vals[int(0.9 * len(vals))], 2),
                "targets": GENE2TARGETS.get(gene, []),
                "source": "CCLE 18q3 RPKM(本地 data/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct); 口径注记: 2018q3 RPKM 旧版, DepMap 新版 TPM 被验证墙挡, 换源时需重标定丰度分档"}
            print("[%s] CCLE 本地解析: n=%d median=%s" % (
                gene, len(vals), abundance[gene]["median_tpm"]))
            continue
        cache = os.path.join(DATA, "hpa_rna_celline_%s.tsv" % gene)
        try:
            txt = fetch_hpa(gene, cache)
        except Exception as e:
            print("[%s] HPA 拉取失败: %s" % (gene, e))
            continue
        rows = list(csv.DictReader(io.StringIO(txt), delimiter="\t"))
        cells = {}
        key = None
        for row in rows:
            for k, v in row.items():
                if k and "RNA cell line" in k and v not in (None, "", "NA"):
                    key = k
                    m = re.match(r"([^:]+)\s*:\s*([0-9.]+)", v)
                    if m:
                        cells[m.group(1)] = float(m.group(2))
        vals = sorted(cells.values())
        if not vals:
            print("[%s] HPA 返回无细胞系 RNA 数据(格式或列名变动), 保留缓存待查" % gene)
            continue
        abundance[gene] = {
            "n_cell_lines": len(cells),
            "median_tpm": round(vals[len(vals) // 2], 2),
            "p10_tpm": round(vals[int(0.1 * len(vals))], 2),
            "p90_tpm": round(vals[int(0.9 * len(vals))], 2),
            "targets": GENE2TARGETS.get(gene, []),
            "source": "Human Protein Atlas RNA cell line (缓存 %s)" %
                      os.path.basename(cache)}

    out = os.path.join(DATA, "virtual_cell_abundance.json")
    json.dump({"genes": abundance,
               "note": "丰度分档供杀伤窗口计算; DepMap 全矩阵待校准期接入"},
              open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("丰度表 -> %s" % out)
    for g, a in abundance.items():
        print("  %-5s n=%d median=%s TPM (p10=%s p90=%s)" % (
            g, a["n_cell_lines"], a["median_tpm"], a["p10_tpm"], a["p90_tpm"]))

    # 预校准包络: EC50 网格 × 基因, 高低分位激活差
    def hill(rna, ec50, n=2.0):
        return rna ** n / (ec50 ** n + rna ** n)

    env = []
    for ec50 in [float(x) for x in args.ec50_grid.split(",")]:
        for g, a in abundance.items():
            hi = hill(a["p90_tpm"], ec50, PARAM_INTERFACE["n_hill"])
            lo = hill(a["p10_tpm"], ec50, PARAM_INTERFACE["n_hill"])
            env.append({"gene": g, "ec50_tpm": ec50,
                        "act_p90": round(hi, 3), "act_p10": round(lo, 3),
                        "window": round(hi - lo, 3)})
    json.dump({"model": "Hill 阈值框架; 参数未标定, 此为包络非预测",
               "params_interface": PARAM_INTERFACE, "envelope": env},
              open(os.path.join(DATA, "virtual_cell_envelope.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump(PARAM_INTERFACE,
              open(os.path.join(DATA, "virtual_cell_params.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    print("包络 -> data/virtual_cell_envelope.json; 参数接口 -> data/virtual_cell_params.json")


if __name__ == "__main__":
    main()
