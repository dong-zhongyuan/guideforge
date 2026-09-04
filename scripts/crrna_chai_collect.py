"""Chai-1 共折叠分数收集(策划案 V3 蛋白预测层首跑, 2026-09-02)。

读取容器产出的 scores.model_idx_*.npz(经 /dawn 回传到宿主), 汇总
aggregate/pTM/ipTM/链间 pairwise ipTM, 写入 data/chai_cofold.json。
interpretation 文字由模型分数程序化生成(方向与量级阈值规则见 build_interpretation),
不再手写——消除"文字与数据脱节"及残留旧解读串的可能(round-3 评审整改)。
运行(宿主): python scripts/crrna_chai_collect.py --dir /tmp/chai_scores/WT2 --label WT_single_seq
重生成既有 JSON 的 interpretation: python scripts/crrna_chai_collect.py --reinterpret data/chai_cofold_X.json
"""
import argparse
import glob
import json
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

IPTM_BAND = (0.5, 0.6)  # 界面预测可用置信区间的常用口径


def _mean(rows, key):
    return sum(r[key] for r in rows) / len(rows) if rows else 0.0


def uses_self_template(label):
    """该 label 对应的运行是否注入了 8D4A 自模板。

    依据提交史: template* 与 matrix_* 两批均带 m8 自模板(commit e69fc74/963745d),
    single_seq / esm_local 为无模板裸序列探针。"""
    return "template" in label or label.startswith("matrix")


def build_interpretation(rows, label):
    """由 5 模型分数均值程序化生成解读文字。

    阈值规则: prot-crRNA ipTM < 0.5 判"低于可用区间"; 0.5-0.6 判"下沿, 弱参考";
    >= 0.6 判"达可用区间"。自模板运行追加降级声明(100% 一致自模板的高分为
    模板合规性平凡结果, 仅 sanity check, 不作界面可预测性证据)。"""
    if not rows:
        return "无模型分数可读"
    iptm = _mean(rows, "iptm")
    pc = _mean(rows, "iptm_prot_crRNA")
    ct = _mean(rows, "iptm_crRNA_target")
    parts = ["%d模型均值: ipTM %.3f, prot-crRNA ipTM %.3f, crRNA-靶RNA ipTM %.3f"
             % (len(rows), iptm, pc, ct)]
    if pc < IPTM_BAND[0]:
        parts.append("prot-crRNA ipTM %.2f 低于可用界面预测区间(~%.1f-%.1f)"
                     % (pc, IPTM_BAND[0], IPTM_BAND[1]))
    elif pc < IPTM_BAND[1]:
        parts.append("prot-crRNA ipTM %.2f 处于可用区间(~%.1f-%.1f)下沿, 仅可作弱参考"
                     % (pc, IPTM_BAND[0], IPTM_BAND[1]))
    else:
        parts.append("prot-crRNA ipTM %.2f 达可用区间(>=%.1f)" % (pc, IPTM_BAND[1]))
    if uses_self_template(label):
        parts.append("本运行注入 100% 一致自模板(8D4A链A): 上述分数为模板合规性/"
                     "健全性检查(sanity check), 高 ipTM 是自模板复述的平凡结果, "
                     "不能作为界面可预测性证据; 相关结论待 template-free 或 "
                     "scrambled-template 对照校准")
    else:
        parts.append("单序列口径(无MSA/无模板): 分数为裸序列预测能力下限")
    return "; ".join(parts)


def reinterpret(path):
    payload = json.load(open(path, encoding="utf-8"))
    payload["self_template"] = uses_self_template(payload.get("label", ""))
    payload["interpretation"] = build_interpretation(payload.get("models", []),
                                                     payload.get("label", ""))
    json.dump(payload, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("重生成 interpretation -> %s" % path)
    print(payload["interpretation"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", help="scores.npz 所在目录")
    ap.add_argument("--label")
    ap.add_argument("--reinterpret", metavar="JSON",
                    help="重生成既有收集 JSON 的 interpretation 字段(由 models 数值推导)")
    args = ap.parse_args()

    if args.reinterpret:
        reinterpret(args.reinterpret)
        return
    if not args.dir or not args.label:
        ap.error("--dir 与 --label 为收集模式必填")

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
               "self_template": uses_self_template(args.label),
               "models": rows, "best_iptm_model": best,
               "interpretation": build_interpretation(rows, args.label)}
    out = os.path.join(ROOT, "data", "chai_cofold_%s.json" % args.label)
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
