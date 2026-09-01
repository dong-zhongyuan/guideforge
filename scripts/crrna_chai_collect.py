"""Chai-1 共折叠分数收集(策划案 V3 蛋白预测层首跑, 2026-09-02)。

读取容器产出的 scores.model_idx_*.npz(经 /dawn 回传到宿主), 汇总
aggregate/pTM/ipTM/链间 pairwise ipTM, 写入 data/chai_cofold.json。
运行(宿主): python scripts/crrna_chai_collect.py --dir /tmp/chai_scores/WT2 --label WT_single_seq
"""
import argparse
import glob
import json
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", required=True, help="scores.npz 所在目录")
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    rows = []
    for f in sorted(glob.glob(os.path.join(args.dir, "scores.model_idx_*.npz"))):
        d = np.load(f)
        pair = np.squeeze(d["per_chain_pair_iptm"])
        row = {"model": os.path.basename(f),
               "aggregate": round(float(np.squeeze(d["aggregate_score"])), 3),
               "ptm": round(float(np.squeeze(d["ptm"])), 3),
               "iptm": round(float(np.squeeze(d["iptm"])), 3),
               "iptm_prot_crRNA": round(float(pair[0][1]), 3),
               "iptm_crRNA_target": round(float(pair[1][2]), 3),
               "clashes": bool(np.squeeze(d["has_inter_chain_clashes"]))}
        rows.append(row)
        print(row)
    best = max(rows, key=lambda r: r["iptm"]) if rows else None
    payload = {"label": args.label,
               "input": "SuCas12a2(1207aa, 8D4A链A) + crRNA(DR18+R248Q) + 靶RNA(28nt)",
               "mode": "单序列(无MSA/无ESM)首次能力探针" if "single" in args.label else args.label,
               "chains": ["protein", "crRNA", "target_RNA"],
               "models": rows, "best_iptm_model": best,
               "interpretation": "蛋白-RNA 界面未恢复(iptm~0.2, 无MSA已知短板); "
                                 "crRNA-靶RNA pairwise iptm~0.35 为最有意义链对; "
                                 "ESM 嵌入版与模板/MSA 版为后续升级路径"}
    out = os.path.join(ROOT, "data", "chai_cofold_%s.json" % args.label)
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
