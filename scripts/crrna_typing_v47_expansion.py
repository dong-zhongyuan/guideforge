# -*- coding: utf-8 -*-
"""上下文分型 v47 扩库(§A2a 配套, 2026-09-10): spacer 上下文 85 -> 数百条重聚类。

动机: §A/§A2 的分型证据由 n=86(实际 85 条文献 spacer)承载, n=9 小簇主导。
本脚本把上下文库扩到数百条: 从 GENCODE v47 转录本库(data/raw/, 不入 git)
按确定性规则抽取 23nt spacer 上下文——癌基因热点 30 基因 × 至多 6 窗 +
背景转录本 40 条 × 3 窗; 与文献 85 条合并成扩库队列, 同一 DR(zeng2026)
下算同一套 5 维特征, 按 §A2a 预登记 k=3 重聚类(轮廓谱仅诊断), 每型 3 代表
重跑主管线(topk=16), 按 §A2b/A2c 口径做置换检验(含/去 position-1 两层)。

判读文字全部由数值按显式阈值规则生成, 程序内不手写结论方向。
全部抽样规则确定(固定起点/步长, 无随机数), 同一 v47 文件可逐条复现。

用法:
  PYTHONUTF8=1 python scripts/crrna_typing_v47_expansion.py
输入: data/raw/gencode.v47.transcripts.fa(本地, 380MB, 不入库)
输出: data/context_typing/spacers_v47_expansion.txt(新增 spacer 全列)
      data/context_typing/typing_v47_expansion.json(聚类+判据+判读)
"""
import itertools
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_context_typing as typing_mod  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402
import crrna_typing_robustness as rb  # noqa: E402
from scaffold_registry import get_scaffold  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'context_typing')
V47 = os.path.join(ROOT, 'data', 'raw', 'gencode.v47.transcripts.fa')
LIT_SPACERS = os.path.join(DATA, 'spacers.txt')
OUT_TXT = os.path.join(DATA, 'spacers_v47_expansion.txt')
OUT_JSON = os.path.join(DATA, 'typing_v47_expansion.json')

K_FIXED = 3            # §A2a: k 预登记固定
REPS_PER_TYPE = 3
SPACER_LEN = 23        # 与文献队列同长
SEEDS = tuple(range(50))
P_THRESH = 0.05        # §A2b/A2c 项目设定常数

CANCER_GENES = ["TP53", "KRAS", "APC", "EGFR", "MYC", "BRAF", "PIK3CA",
                "PTEN", "RB1", "NRAS", "IDH1", "IDH2", "CTNNB1", "SMAD4",
                "FBXW7", "ERBB2", "AKT1", "CDKN2A", "VHL", "MLH1", "MSH2",
                "STK11", "KEAP1", "NF1", "AR", "GATA3", "PIK3R1", "HRAS",
                "ALK", "RET"]
CANCER_WINDOWS = 6     # 每基因至多 6 窗
CANCER_START, CANCER_STRIDE = 137, 389     # 质数步长, 确定性
BG_EVERY, BG_MAX, BG_WINDOWS = 997, 40, 3  # 每 997 条转录本取 1 条背景
BG_START, BG_STRIDE = 151, 271

HOMOPOLYMER = re.compile(r"(.)\1{4,}")


def parse_v47(path):
    """流式解析 GENCODE fasta, 产出 (transcript_id, gene, seq)。"""
    tid = gene = None
    buf = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if tid is not None:
                    yield tid, gene, "".join(buf)
                parts = line[1:].strip().split("|")
                tid = parts[0]
                gene = parts[5] if len(parts) > 5 else ""
                buf = []
            else:
                buf.append(line.strip())
    if tid is not None:
        yield tid, gene, "".join(buf)


def window_ok(w):
    if len(w) != SPACER_LEN or set(w) - set("ACGT"):
        return False
    gc = (w.count("G") + w.count("C")) / len(w)
    if not 0.25 <= gc <= 0.75:
        return False
    return not HOMOPOLYMER.search(w)


def take_windows(seq, start, stride, n_max):
    out = []
    lo, hi = 100, max(100, len(seq) - 100 - SPACER_LEN)
    pos = start
    while pos <= hi and len(out) < n_max:
        w = seq[pos:pos + SPACER_LEN]
        if window_ok(w):
            out.append((pos, w))
        pos += stride
    return out


def sample_spacers():
    """确定性抽样: 返回 [(name, dna, provenance)]。"""
    picked = []
    seen_genes = set()
    bg_count = 0
    for i, (tid, gene, seq) in enumerate(parse_v47(V47)):
        if gene in CANCER_GENES and gene not in seen_genes:
            seen_genes.add(gene)
            for pos, w in take_windows(seq, CANCER_START, CANCER_STRIDE,
                                       CANCER_WINDOWS):
                picked.append(("v47-%s-%s-%d" % (gene, tid.split(".")[0], pos),
                               w, {"kind": "cancer_gene", "gene": gene,
                                   "transcript": tid, "pos": pos}))
        elif i % BG_EVERY == 0 and bg_count < BG_MAX and len(seq) > 400:
            bg_count += 1
            for pos, w in take_windows(seq, BG_START, BG_STRIDE, BG_WINDOWS):
                picked.append(("v47-bg-%s-%d" % (tid.split(".")[0], pos),
                               w, {"kind": "background", "gene": gene,
                                   "transcript": tid, "pos": pos}))
    return picked


def perm_null(lists_flat, idx, topk):
    """§A2 置换零分布: 9 代表 3/3/3 无序划分(280 种)的公共集大小。"""
    null = []
    for rest_a in itertools.combinations(idx[1:], 2):
        ga = (0,) + rest_a
        rem = [i for i in idx if i not in ga]
        for gb in itertools.combinations(rem, 3):
            gc = tuple(i for i in rem if i not in gb)
            if gb[0] > gc[0]:
                continue
            unions = [set().union(*(lists_flat[i] for i in g))
                      for g in (ga, gb, gc)]
            null.append(len(set.intersection(*unions)))
    return np.array(null)


def main():
    t0 = time.time()
    dr = core.to_rna(get_scaffold(rb.EFFECTOR))

    lit = [l.strip().split("\t") for l in open(LIT_SPACERS, encoding="utf-8")
           if l.strip()]
    lit_names = [n for n, _ in lit]
    print("[v47] 文献队列 %d 条; 解析 v47 转录本库抽样 ..." % len(lit_names))
    new = sample_spacers()
    new_names = [n for n, _, _ in new]
    print("[v47] 新增 %d 条(癌基因 %d / 背景 %d)"
          % (len(new),
             sum(1 for _, _, p in new if p["kind"] == "cancer_gene"),
             sum(1 for _, _, p in new if p["kind"] == "background")))

    with open(OUT_TXT, "w", encoding="utf-8") as fh:
        for n, dna, _ in new:
            fh.write("%s\t%s\n" % (n, dna))

    names = lit_names + new_names
    seqs = [core.to_rna(s) for _, s in lit] + [core.to_rna(d) for _, d, _ in new]
    print("[v47] 计算 %d 条特征矩阵 ..." % len(names))
    F = np.array([typing_mod.spacer_features(dr, n, s)
                  for n, s in zip(names, seqs)], dtype=float)
    Z = (F - F.mean(0)) / F.std(0)

    lab, centers, inertia = typing_mod.kmeans(Z, K_FIXED, SEEDS)
    sizes = [int((lab == c).sum()) for c in range(K_FIXED)]
    sil = float(typing_mod.silhouette(Z, lab, K_FIXED))
    print("[v47] k=3(预登记固定) 簇大小 %s, 轮廓系数 %.3f(仅诊断)" %
          (sizes, sil))

    reps = rb.reps_from_clusters(lab, Z, names, seqs, K_FIXED, REPS_PER_TYPE)
    prefixes, n_fresh = {}, 0
    for c, rl in reps.items():
        for nm, dna in rl:
            prefixes[nm], fresh = rb.run_pipeline(nm, dna)
            n_fresh += fresh
    print("[v47] 9 代表管线重跑: 新跑 %d, 缓存 %d" % (n_fresh, 9 - n_fresh))

    types = sorted(reps.keys())
    lists9 = {t: [rb.top_descs(prefixes[nm], rb.CRIT_TOPK) for nm, _ in reps[t]]
              for t in types}
    lists9_np = {t: [rb.top_descs_no_pos1(prefixes[nm], rb.CRIT_TOPK)
                     for nm, _ in reps[t]] for t in types}
    obs = rb.criterion(lists9, rb.CRIT_TOPK)
    obs_np = rb.criterion(lists9_np, rb.CRIT_TOPK)
    idx = list(range(9))
    flat = [rb.top_descs(prefixes[nm], rb.CRIT_TOPK)
            for t in types for nm, _ in reps[t]]
    flat_np = [rb.top_descs_no_pos1(prefixes[nm], rb.CRIT_TOPK)
               for t in types for nm, _ in reps[t]]
    null, null_np_ = perm_null(flat, idx, rb.CRIT_TOPK), perm_null(flat_np, idx, rb.CRIT_TOPK)
    p_b = float((null <= obs["n_common"]).mean())
    p_c = float((null_np_ <= obs_np["n_common"]).mean())
    pos_b, pos_c = bool(p_b < P_THRESH), bool(p_c < P_THRESH)

    reading = ("扩库队列 n=%d(k=3 预登记固定): A2b(含 pos-1) 公共集 %d 条 "
               "P=%.4f %s %.2f;A2c(去 pos-1, 主结论口径) 公共集 %d 条 P=%.4f "
               "%s %.2f → %s"
               % (len(names), obs["n_common"], p_b, "<" if pos_b else ">=",
                  P_THRESH, obs_np["n_common"], p_c, "<" if pos_c else ">=",
                  P_THRESH,
                  "扩库后预登记口径下阳性" if (pos_b and pos_c)
                  else ("主结论口径 A2c 阳性" if pos_c
                        else "主结论口径 A2c 阴性")))
    print("[§A2-v47]", reading)

    # ---- §A3 特征增强版(2026-09-10 预登记, docs/preregistration.md §A3) ----
    print("[v47] §A3: 计算 10 维增强特征矩阵 ...")
    F2 = np.array([typing_mod.spacer_features_v2(dr, n, s)
                   for n, s in zip(names, seqs)], dtype=float)
    Z2 = (F2 - F2.mean(0)) / F2.std(0)
    lab2, _, inertia2 = typing_mod.kmeans(Z2, K_FIXED, SEEDS)
    sizes2 = [int((lab2 == c).sum()) for c in range(K_FIXED)]
    sil2 = float(typing_mod.silhouette(Z2, lab2, K_FIXED))
    reps2 = rb.reps_from_clusters(lab2, Z2, names, seqs, K_FIXED, REPS_PER_TYPE)
    n_fresh2 = 0
    for c, rl in reps2.items():
        for nm, dna in rl:
            prefixes[nm], fresh = rb.run_pipeline(nm, dna)
            n_fresh2 += fresh
    print("[v47] §A3 九代表管线重跑: 新跑 %d" % n_fresh2)
    types2 = sorted(reps2.keys())
    obs2 = rb.criterion({t: [rb.top_descs(prefixes[nm], rb.CRIT_TOPK)
                             for nm, _ in reps2[t]] for t in types2},
                        rb.CRIT_TOPK)
    obs2_np = rb.criterion({t: [rb.top_descs_no_pos1(prefixes[nm], rb.CRIT_TOPK)
                                for nm, _ in reps2[t]] for t in types2},
                           rb.CRIT_TOPK)
    flat2 = [rb.top_descs(prefixes[nm], rb.CRIT_TOPK)
             for t in types2 for nm, _ in reps2[t]]
    flat2np = [rb.top_descs_no_pos1(prefixes[nm], rb.CRIT_TOPK)
               for t in types2 for nm, _ in reps2[t]]
    null2 = perm_null(flat2, idx, rb.CRIT_TOPK)
    null2np = perm_null(flat2np, idx, rb.CRIT_TOPK)
    p2_b = float((null2 <= obs2["n_common"]).mean())
    p2_c = float((null2np <= obs2_np["n_common"]).mean())
    pos2_b, pos2_c = bool(p2_b < P_THRESH), bool(p2_c < P_THRESH)
    # §A3c 判读规则(登记先于数值)
    if pos2_c and not pos_c:
        a3_read = ("A3 主口径阳性(P=%.4f)而 A2(v47) 阴性(P=%.4f) → 「特征分辨率"
                   "不足」归因成立, 计算层分型证据升级" % (p2_c, p_c))
    elif not pos2_c and not pos_c:
        a3_read = ("A3 主口径亦阴性(P=%.4f) → 计算层不支持分型, 主张转交 §E "
                   "同源层(E1 已 PASS)与 §F 湿实验交互轴裁决" % p2_c)
    elif pos2_c and pos_c:
        a3_read = ("A3 与 A2(v47) 均阳性(P=%.4f/%.4f) → 分型信号对特征增强稳健"
                   % (p2_c, p_c))
    else:
        a3_read = ("A3 阴性(P=%.4f)而 A2(v47) 阳性(P=%.4f) → 特征增强未带来"
                   "一致信号, 如实并列" % (p2_c, p_c))
    print("[§A3]", a3_read)
    a3 = {
        "registered": "docs/preregistration.md §A3(2026-09-10, 登记先于运行)",
        "feature_names": rb.FEAT_NAMES + typing_mod.FEAT_V2_EXTRA,
        "clustering": {"k_fixed": K_FIXED, "cluster_sizes": sizes2,
                       "inertia": round(float(inertia2), 1),
                       "silhouette_diagnostic_only": round(sil2, 4)},
        "representatives": {t: [nm for nm, _ in reps2[t]] for t in types2},
        "a3b_with_pos1": {"obs_common": obs2["n_common"],
                          "common": obs2["common"], "p_obs": round(p2_b, 4),
                          "threshold_p": P_THRESH, "positive": pos2_b,
                          "null_median": float(np.median(null2))},
        "a3c_no_pos1": {"obs_common": obs2_np["n_common"],
                        "common": obs2_np["common"], "p_obs": round(p2_c, 4),
                        "threshold_p": P_THRESH, "positive": pos2_c,
                        "null_median": float(np.median(null2np))},
        "reading": a3_read,
    }

    report = {
        "generated_by": "scripts/crrna_typing_v47_expansion.py",
        "date": time.strftime("%Y-%m-%d"),
        "registered": "docs/preregistration.md §A2(2026-09-10 登记; k=3 固定, "
                      "置换检验主判据, 去 pos-1 主结论口径)",
        "cohort": {"n_literature": len(lit_names), "n_v47_new": len(new),
                   "n_total": len(names), "spacer_len": SPACER_LEN,
                   "sampling_rule": "确定性: 癌基因 %d 基因×≤%d 窗(起点%d/步长%d) "
                       "+ 背景每 %d 条转录本×%d 窗(起点%d/步长%d); 过滤: ACGT, "
                       "GC∈[0.25,0.75], 无≥5 同聚物"
                       % (len(CANCER_GENES), CANCER_WINDOWS, CANCER_START,
                          CANCER_STRIDE, BG_EVERY, BG_WINDOWS, BG_START,
                          BG_STRIDE)},
        "provenance": {n: p for n, _, p in new},
        "clustering": {"k_fixed": K_FIXED, "kmeans_seeds": list(SEEDS),
                       "cluster_sizes": sizes, "inertia": round(float(inertia), 1),
                       "silhouette_diagnostic_only": round(sil, 4)},
        "representatives": {t: [nm for nm, _ in reps[t]] for t in types},
        "a2b_permutation": {"obs_common": obs["n_common"], "common": obs["common"],
                            "p_obs": round(p_b, 4), "threshold_p": P_THRESH,
                            "positive": pos_b,
                            "null_median": float(np.median(null)),
                            "n_partitions": int(len(null))},
        "a2c_no_pos1": {"obs_common": obs_np["n_common"],
                        "common": obs_np["common"], "p_obs": round(p_c, 4),
                        "threshold_p": P_THRESH, "positive": pos_c,
                        "null_median": float(np.median(null_np_)),
                        "n_partitions": int(len(null_np_))},
        "reading": reading,
        "a3_enriched": a3,
        "comparison_literature_cohort": {
            "n": 85,
            "a2b_p": 0.0357, "a2c_common": 2, "a2c_p": 0.0357,
            "source": "data/context_typing/typing.robustness.json a2_prereg"},
        "runtime_sec": round(time.time() - t0, 1),
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print("[v47] -> %s (%.1fs)" % (OUT_JSON, time.time() - t0))


if __name__ == "__main__":
    main()
