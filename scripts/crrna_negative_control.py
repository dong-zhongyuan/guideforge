"""阴性对照方向性检验(证据链补强 #5, 用户 2026-09-01 指定)。

问题: 打分是否只是"随机排列"? 若构造明知会被破坏的 DR(茎区全打散/随机替换/
组成洗牌), 打分体系应给它们最低排名且硬过滤拦截——方向性的一次性证明。

阴性对照集(每口径):
  1. stem_break_5arm_A   : 5' 茎臂(0基 3-7) 全改 A → 5 对全部失配
  2. stem_arms_swap      : 3' 茎臂序列按配对位置拷到 5' 臂(同向) → 全失配
  3. stem_break_3arm     : 3' 茎臂换为 5' 臂反排 → 全失配
  4. revcomp_dr          : 整条 DR 换为其反向互补
  5. shuffle×N           : WT DR 组成内部洗牌(组成不变, 排列破坏)
  6. random×N            : 等长均匀随机 DR(无 TTTT/GGGG)

判定(先定后跑, 全部阴性对照):
  V1 工程茎破坏 4 条全部被硬过滤拦截(passed=False);
  V2 每条阴性对照分数低于通过库中位数;
  V3 ≥80% 阴性对照分数低于通过库第 10 百分位;
  方向性: 阴性对照 bp_dist 均值 > TOP-12 均值, p_fold 均值 < TOP-12 均值。

诚实边界: 随机/洗牌 DR 与 WT 差异位点数远多于常规突变库(1-2 点),
"打分会给垫底"部分来自突变数量本身; 因此判定以硬过滤拦截与方向性为主,
分数百分位为辅。wt-stab 口径参数与所对照的主线运行完全一致(取自 top.json params)。
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402

# 打分参数一律取自所对照主线运行的 top.json params(单一事实源), 不在本脚本重复默认值;
# 缺键即报错终止。


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def build_negatives(wt_dr_dna, n_shuffle, n_random, rng):
    dr = wt_dr_dna
    arm1, arm2 = dr[3:8], dr[12:17]  # 0基 茎臂(5bp 茎: (3,16)(4,15)(5,14)(6,13)(7,12))
    negs = {
        "stem_break_5arm_A": dr[:3] + "AAAAA" + dr[8:],
        "stem_arms_swap": dr[:3] + arm2[::-1] + dr[8:],
        "stem_break_3arm": dr[:12] + arm1[::-1] + dr[17:],
        "revcomp_dr": revcomp(dr),
    }
    for i in range(n_shuffle):
        while True:
            s = list(dr)
            rng.shuffle(s)
            s = "".join(s)
            if s != dr and "TTTT" not in s and "GGGG" not in s:
                break
        negs["shuffle_%d" % (i + 1)] = s
    for i in range(n_random):
        while True:
            s = "".join(rng.choice("ACGT") for _ in dr)
            if "TTTT" not in s and "GGGG" not in s:
                break
        negs["random_%d" % (i + 1)] = s
    return negs


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top-json", default=os.path.join(ROOT, "data", "tp53_r248q_pdbdr.v2.top.json"))
    ap.add_argument("--n-shuffle", type=int, default=5)
    ap.add_argument("--n-random", type=int, default=10)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = json.load(open(args.top_json, encoding="utf-8"))
    spacer_dna = d["spacer_fixed_dna"]
    wt_entry = d["top"][0]
    wt_dr_dna = wt_entry["dr_seq"]
    dr_rna, spacer_rna = core.to_rna(wt_dr_dna), core.to_rna(spacer_dna)

    # 与主线运行完全一致的参数: 全部取自 top.json params, 缺键报错
    needed = ("max_bp_dist", "spacer_unpaired_margin", "proc_window", "seed_len",
              "cons3_window", "w_bp", "w_ddg", "w_contact", "w_ens", "w_hbond",
              "w_stab", "w_cons3", "w_fold", "w_seed", "stab_dr_only",
              "no_protect_processing")
    params = d.get("params", {})
    missing = [k for k in needed if k not in params]
    if missing:
        raise SystemExit("top.json params 缺键 %s(该运行早于对应参数引入), "
                         "请改用含全部键的运行输出" % missing)
    ns = {k: params[k] for k in needed}
    # 早于 v1.6 的运行无交叉配对过滤, 等价 allow_cross_pairing=True
    ns["allow_cross_pairing"] = params.get("allow_cross_pairing", True)
    score_args = argparse.Namespace(**ns)
    seed_len = ns["seed_len"]

    wt = core.wt_reference(dr_rna, spacer_rna, seed_len)
    contact = core.load_contacts(dr_rna)
    stem_pos = core.stem_positions(wt["dr_only_struct"], len(dr_rna))

    rng = __import__("random").Random(args.seed)
    negs = build_negatives(wt_dr_dna, args.n_shuffle, args.n_random, rng)

    # 参照库: 该主线运行的通过变体分数分布
    csv_path = args.top_json.replace(".top.json", ".variants.csv")
    passed_scores, all_scores = [], []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            s = float(row["score"])
            all_scores.append(s)
            if row["passed"] in ("True", "true", "1"):
                passed_scores.append(s)
    passed_scores.sort()
    import numpy as np
    ps = np.array(passed_scores)

    rows = []
    for name, nd in negs.items():
        r = core.score_variant(dr_rna, core.to_rna(nd), wt, spacer_rna,
                               score_args, contact, stem_pos)
        pct = float((ps < r["score"]).mean() * 100)  # 通过库中分数低于它的比例
        rows.append({"name": name, "n_mut": r["n_mut"], "passed": r["passed"],
                     "score": r["score"], "bp_dist": r["bp_dist"],
                     "p_fold": r["p_fold"], "spacer_unpaired": r["spacer_unpaired"],
                     "ddG_dr": r["ddG_dr"], "pct_of_passed": round(pct, 1),
                     "dr_seq": nd})
    rows.sort(key=lambda r: -r["score"])

    top12 = [e for e in d["top"] if e.get("desc") not in (None, "WT")][:12]
    v1 = all(not r["passed"] for r in rows if r["name"].startswith(("stem", "revcomp")))
    v2 = all(r["score"] < float(np.median(ps)) for r in rows)
    v3 = sum(r["pct_of_passed"] < 10 for r in rows) >= 0.8 * len(rows)
    neg_bp = float(np.mean([r["bp_dist"] for r in rows]))
    top_bp = float(np.mean([e["bp_dist"] for e in top12]))
    neg_pf = float(np.mean([r["p_fold"] for r in rows]))
    top_pf = float(np.mean([e["p_fold"] for e in top12]))
    direction_ok = neg_bp > top_bp and neg_pf < top_pf

    print("口径: %s (DR %s, %dnt) | 通过库 n=%d, 中位分=%.3f, 10百分位=%.3f" % (
        os.path.basename(args.top_json), wt_dr_dna, len(wt_dr_dna),
        len(ps), float(np.median(ps)), float(np.percentile(ps, 10))))
    print("%-20s %5s %5s %8s %6s %7s %8s %8s" % (
        "negative", "nmut", "pass", "score", "bp_d", "p_fold", "sp_up", "pct<10"))
    for r in rows:
        print("%-20s %5d %5s %8.3f %6d %7.4f %8.3f %7s" % (
            r["name"], r["n_mut"], "Y" if r["passed"] else "N", r["score"],
            r["bp_dist"], r["p_fold"], r["spacer_unpaired"],
            "YES" if r["pct_of_passed"] < 10 else "no"))
    print("\nV1 茎破坏全部被硬过滤拦截: %s" % ("PASS" if v1 else "FAIL"))
    print("V2 全部低于通过库中位分: %s" % ("PASS" if v2 else "FAIL"))
    print("V3 >=80%% 低于通过库第10百分位: %s" % ("PASS" if v3 else "FAIL"))
    print("方向性 bp_dist 阴性均值 %.1f > TOP12 均值 %.1f; p_fold %.4f < %.4f: %s" % (
        neg_bp, top_bp, neg_pf, top_pf, "PASS" if direction_ok else "FAIL"))

    out = args.out or os.path.join(
        ROOT, "data", "negative_control.%s.json" % os.path.basename(args.top_json).replace(".top.json", ""))
    json.dump({"source_top_json": args.top_json, "n_negative": len(rows),
               "n_passed_library": len(ps), "rows": rows,
               "verdict": {"V1_stem_break_all_filtered": v1,
                           "V2_all_below_median": v2,
                           "V3_80pct_below_p10": v3,
                           "directionality": direction_ok},
               "negative_bp_mean": neg_bp, "top12_bp_mean": top_bp,
               "negative_pfold_mean": neg_pf, "top12_pfold_mean": top_pf},
              open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
