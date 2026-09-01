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


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--genes", default="TP53,KRAS,APC,MYC")
    ap.add_argument("--ec50-grid", default="0.5,1,2,5,10,20,50",
                    help="EC50(TPM) 扫描网格, 预校准包络用")
    args = ap.parse_args()

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
                "source": "CCLE 18q3 RPKM(本地 data/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct)"}
            print("[%s] CCLE 本地解析: n=%d median=%s" % (
                gene, len(vals), abundance[gene]["median"]))
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
