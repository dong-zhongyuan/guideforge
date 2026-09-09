# -*- coding: utf-8 -*-
"""Chai-1 增量折叠结果回收合并器(2026-09-06, V3 对齐矩阵落地)。

输入: 服务器跑完 data/chai_increment_inputs/ 11 折后回收的 scores npz 目录树
(每折一个子目录, 目录名 = manifest.json 里 job 的 fasta 文件名去 .fasta,
内含 scores.model_idx_*.npz × 5)。

合并规则(与 scripts/crrna_chai_increment_input.py 的增量策略互为逆过程):
  旧矩阵 data/chai_matrix_4t.json 32 行中
    - 复用 21 行: 3 个保留靶(TP53_R248Q/KRAS_G12D/TP53_R273H) × 7 个非补偿骨架
      (协议一致, 不重跑);
    - 丢弃 8 行: KRAS_G12C 整列(V3 对齐后 G12C 退出湿实验矩阵, 仅留干实验证据);
    - 替换 4 行: 旧 B_break_compensate 行折叠对象是 round-3 前已退役的茎破坏
      混杂臂(DR AATTTCTACTCTTCTACAT), 由新跑的去混杂新臂(zengDR+T10A/T12G/T19G,
      DR AATTTCTACAGGTGTAGAG) × 3 保留靶替换;
    - 新增 8 行: APC_Q1328x × 8 骨架(含 WT)。
  输出 4 靶 × 8 骨架 = 32 行, 行结构与旧矩阵完全一致
  (target/scaffold/n_models/iptm_mean/iptm_sd/prot_crRNA/crRNA_target/clash_frac),
  写入新文件(默认 data/chai_matrix_v3aligned.json), 不覆盖旧矩阵。

口径声明(诚实并列, 与 crrna_chai_matrix_summary.py 的 claim_downgrade 一致):
  1. 矩阵全部分数来自 100% 一致自模板(8D4A链A, data/chai_template_hits.m8)注入
     + ESM, 高 ipTM 为模板合规性的平凡结果, 仅 sanity check, 不作界面可预测性
     证据; 结论待 template-free / scrambled-template 对照校准;
  2. 复用的 21 行与新增/替换的 11 行都是真实 Chai-1 运行的 5 模型聚合值——
     "程序化生成"的是 summary 文件里的 reading 文字(由 scaffold_deltas 按阈值
     规则生成), 不是矩阵数值本身;
  3. iptm_sd 为 5 模型 aggregate ipTM 的总体标准差(np.std, ddof=0), 作
     prot-crRNA 噪声代理(后者 per-model sd 未存盘)。旧 21 行的 sd 生成脚本
     未入库(commit 722d8ea 仅提交数据), 原始口径不可考, 新旧行混用时对
     sd 口径差异保持知情。

运行:
  合并: python scripts/crrna_chai_increment_merge.py --runs-dir <回收目录>
  自检: python scripts/crrna_chai_increment_merge.py --selftest
        (合成 11 折 npz 验证合并逻辑, 不需要真实 Chai 输出)
"""
import argparse
import glob
import json
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

DROP_TARGET = "KRAS_G12C"
REPLACE_SCAFFOLD = "B_break_compensate"
NEW_TARGET = "APC_Q1328x"

TEMPLATE_CAVEAT = (
    "全部 32 行为 100% 一致自模板(8D4A链A)注入 + ESM 结果: 仅模板合规性 "
    "sanity check, 非界面可预测性证据; 结论待 template-free/scrambled-template "
    "对照校准")

PROVENANCE_NOTE = (
    "行数值均为真实 Chai-1 运行的 5 模型聚合(复用行取自 data/chai_matrix_4t.json "
    "commit 722d8ea; 新增/替换 11 行取自本次增量回收 npz); 矩阵数值本身不存在 "
    "程序推算成分——'程序化生成'仅指 chai_matrix_4t_summary.json 的 reading 文字 "
    "(由 scaffold_deltas 按阈值规则生成)。iptm_sd = aggregate ipTM 5 模型总体标准差 "
    "(np.std ddof=0); 旧复用行 sd 原始口径未入库不可考, 混用时保持知情。")


def _nrm(s):
    """命名归一(2026-09-09): 兼容 +/-/_ 三种历史形式的键比较与去重。"""
    return s.replace("+", "").replace("-", "_").replace("_", "") if False else         s.replace("+", "").replace("-", "_").replace("_", "")


def canonical_scaffold(name):
    """增量 fasta 的下划线骨架名 -> 旧矩阵的紧凑名(B_break_compensate 除外)。"""
    return name if name == REPLACE_SCAFFOLD else name.replace("_", "")


def parse_fold(fold_dir):
    """一折目录的 scores.model_idx_*.npz -> 矩阵行数值(链序: 0=protein,1=crRNA,2=target)。"""
    files = sorted(glob.glob(os.path.join(fold_dir, "scores.model_idx_*.npz")))
    if not files:
        raise FileNotFoundError("缺 scores npz: %s" % fold_dir)
    iptms, pc, ct, clashes = [], [], [], 0
    for f in files:
        d = np.load(f)
        pair = np.squeeze(d["per_chain_pair_iptm"])
        iptms.append(float(np.squeeze(d["iptm"])))
        pc.append(float(pair[0][1]))
        ct.append(float(pair[1][2]))
        clashes += int(bool(np.squeeze(d["has_inter_chain_clashes"])))
    n = len(files)
    return {"n_models": n,
            "iptm_mean": round(float(np.mean(iptms)), 4),
            "iptm_sd": round(float(np.std(iptms)), 4),
            "prot_crRNA": round(float(np.mean(pc)), 4),
            "crRNA_target": round(float(np.mean(ct)), 4),
            "clash_frac": round(clashes / n, 4)}


def merge(matrix_path, runs_dir, manifest_path):
    old = json.load(open(matrix_path, encoding="utf-8"))["rows"]
    reused = [dict(r) for r in old
              if _nrm(r["target"]) != _nrm(DROP_TARGET) and _nrm(r["scaffold"]) != _nrm(REPLACE_SCAFFOLD)]
    if len(reused) != 21:
        raise ValueError("旧矩阵复用行应为 21(3 保留靶 x 7 非补偿骨架; 命名归一化匹配), 实得 %d" % len(reused))

    jobs = json.load(open(manifest_path, encoding="utf-8"))["jobs"]
    new_rows = []
    for job in jobs:
        stem = os.path.splitext(job["file"])[0]
        vals = parse_fold(os.path.join(runs_dir, stem))
        new_rows.append({"target": job["target"],
                         "scaffold": canonical_scaffold(job["scaffold"]),
                         **vals})

    rows = reused + new_rows
    targets = sorted({r["target"] for r in rows})
    scaffolds = sorted({r["scaffold"] for r in rows})
    key = lambda r: (_nrm(r["target"]), _nrm(r["scaffold"]))
    if len(rows) != 32 or len({key(r) for r in rows}) != 32:
        raise ValueError("合并后应为 4靶x8骨架=32 不重复行, 实得 %d 行" % len(rows))
    rows.sort(key=lambda r: (scaffolds.index(r["scaffold"]),
                             targets.index(r["target"])))
    payload = {
        "engine": "chai-1 + 8D4A自模板m8 + ESM",
        "template_caveat": TEMPLATE_CAVEAT,
        "design": "4靶 x 8骨架(含WT与B_break_compensate); 靶RNA=protospacer窗口24nt+PFS 5nt(8D4A同构)",
        "v3_alignment": {
            "dropped": "%s 整列 8 行(V3 对齐: 仅留干实验证据, 不入矩阵)" % DROP_TARGET,
            "replaced": ("B_break_compensate 4 行: 旧行为已退役茎破坏混杂臂 "
                         "(DR AATTTCTACTCTTCTACAT), 新行为去混杂新臂 "
                         "(zengDR+T10A/T12G/T19G, DR AATTTCTACAGGTGTAGAG)"),
            "added": "%s x 8 骨架(新靶列)" % NEW_TARGET,
            "reused": "3 保留靶 x 7 非补偿骨架 = 21 行(协议一致不重跑)"},
        "provenance": PROVENANCE_NOTE,
        "rows": rows}
    return payload, reused, new_rows


def print_summary(payload, old_matrix_path):
    rows = payload["rows"]
    old = {(r["target"], r["scaffold"]): r
           for r in json.load(open(old_matrix_path, encoding="utf-8"))["rows"]}
    by = {(r["target"], r["scaffold"]): r for r in rows}
    targets = sorted({r["target"] for r in rows})
    print("=" * 60)
    print("合并后矩阵: %d 行, 靶标 %s" % (len(rows), targets))
    print("[WT 参照 prot_crRNA] 旧 -> 新")
    for t in targets:
        o = old.get((t, "WT"))
        n = by[(t, "WT")]
        print("  %-12s %s -> %.4f" % (
            t, ("%.4f" % o["prot_crRNA"]) if o else "(新列) --", n["prot_crRNA"]))
    print("[B_break_compensate 换臂] 旧(退役混杂臂) -> 新(去混杂新臂)")
    for t in targets:
        o = old.get((t, REPLACE_SCAFFOLD))
        n = by[(t, REPLACE_SCAFFOLD)]
        if o is None:
            print("  %-12s (新列) -> %.4f" % (t, n["prot_crRNA"]))
        else:
            print("  %-12s %.4f -> %.4f (Δ%+.4f)" % (
                t, o["prot_crRNA"], n["prot_crRNA"],
                n["prot_crRNA"] - o["prot_crRNA"]))
    fields = ("n_models", "iptm_mean", "iptm_sd", "prot_crRNA",
              "crRNA_target", "clash_frac")
    print("[复用行] 21 行数值与旧矩阵一致: %s" % all(
        by[k][f] == o[f] for k, o in old.items() if o["target"] != DROP_TARGET
        and o["scaffold"] != REPLACE_SCAFFOLD for f in fields))


def selftest():
    """合成 11 折 npz, 对真实旧矩阵跑合并, 校验结构与数值。"""
    manifest_path = os.path.join(DATA, "chai_increment_inputs", "manifest.json")
    matrix_path = os.path.join(DATA, "chai_matrix_4t.json")
    jobs = json.load(open(manifest_path, encoding="utf-8"))["jobs"]
    tmp = tempfile.mkdtemp(prefix="chai_merge_selftest_")
    # 每折 5 模型: iptm = 0.80 + j*0.001, 链对分数取可识别合成值
    for j, job in enumerate(jobs):
        stem = os.path.splitext(job["file"])[0]
        d = os.path.join(tmp, stem)
        os.makedirs(d)
        for i in range(5):
            pair = np.zeros((3, 3))
            pair[0][1] = 0.50 + j * 0.01
            pair[1][2] = 0.60 + j * 0.01
            np.savez(os.path.join(d, "scores.model_idx_%d.npz" % i),
                     aggregate_score=np.array(0.85), ptm=np.array(0.9),
                     iptm=np.array(0.80 + i * 0.001),
                     per_chain_pair_iptm=pair,
                     has_inter_chain_clashes=np.array(i == 0))
    payload, reused, new_rows = merge(matrix_path, tmp, manifest_path)
    rows = payload["rows"]
    assert len(rows) == 32, len(rows)
    assert sorted({_nrm(r["target"]) for r in rows}) == sorted(map(_nrm, [
        "APC_Q1328x", "KRAS_G12D", "TP53_R248Q", "TP53_R273H"]))
    assert all(len([r for r in rows if _nrm(r["target"]) == _nrm(t)]) == 8
               for t in {r["target"] for r in rows})
    # 新行数值与合成输入一致(以 job0 = TP53_R248Q x B_break_compensate 为例)
    r0 = next(r for r in rows
              if r["target"] == "TP53_R248Q" and r["scaffold"] == REPLACE_SCAFFOLD)
    assert r0["n_models"] == 5 and abs(r0["iptm_mean"] - 0.802) < 1e-4, r0
    assert abs(r0["prot_crRNA"] - 0.50) < 1e-4 and abs(r0["crRNA_target"] - 0.60) < 1e-4
    assert abs(r0["clash_frac"] - 0.2) < 1e-4
    # 骨架名映射: 下划线名 -> 紧凑名
    apc = {r["scaffold"] for r in rows if _nrm(r["target"]) == _nrm(NEW_TARGET)}
    assert "A1GU3AA8GU15C" in apc and "A1G_U3A_A8G_U15C" not in apc, apc
    # 复用 21 行与旧矩阵逐值一致; 旧 B 臂行与 G12C 列不带入
    old = {(_nrm(r["target"]), _nrm(r["scaffold"])): r
           for r in json.load(open(matrix_path, encoding="utf-8"))["rows"]}
    for r in reused:
        o = old[(_nrm(r["target"]), _nrm(r["scaffold"]))]
        assert all(r[f] == o[f] for f in
                   ("n_models", "iptm_mean", "iptm_sd", "prot_crRNA",
                    "crRNA_target", "clash_frac")), (r, o)
    assert not any(_nrm(r["target"]) == _nrm(DROP_TARGET) for r in rows)
    old_b = old[(_nrm("TP53_R248Q"), _nrm(REPLACE_SCAFFOLD))]["prot_crRNA"]
    assert abs(r0["prot_crRNA"] - old_b) > 1e-9, "替换行不应等于退役臂旧值"
    print("[selftest] 合成 11 折合并校验全部通过 (32 行, 4靶x8骨架)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs-dir",
                    help="回收的 Chai 输出根目录(每折一个子目录, 名=fasta 去后缀)")
    ap.add_argument("--matrix", default=os.path.join(DATA, "chai_matrix_4t.json"),
                    help="旧矩阵路径(默认 data/chai_matrix_4t.json)")
    ap.add_argument("--manifest",
                    default=os.path.join(DATA, "chai_increment_inputs",
                                         "manifest.json"),
                    help="增量 manifest(默认 data/chai_increment_inputs/manifest.json)")
    ap.add_argument("--out", default=os.path.join(DATA, "chai_matrix_v3aligned.json"),
                    help="输出路径(默认 data/chai_matrix_v3aligned.json, 不覆盖旧矩阵)")
    ap.add_argument("--selftest", action="store_true",
                    help="用合成 npz 验证合并逻辑(无需真实 Chai 输出)")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(selftest())
    if not args.runs_dir:
        ap.error("--runs-dir 为合并模式必填(或改用 --selftest)")
    payload, _, _ = merge(args.matrix, args.runs_dir, args.manifest)
    json.dump(payload, open(args.out, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("输出 -> %s" % args.out)
    print_summary(payload, args.matrix)


if __name__ == "__main__":
    main()
