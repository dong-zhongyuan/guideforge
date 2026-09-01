"""多 spacer 合成优先级清单(证据链补强 #9, 用户 2026-09-01 指定)。

输入: 同一 DR 口径下、三个 spacer 上下文各自的全变体库与 TOP 榜——
  主靶标 TP53-R248Q + 跨 spacer 稳健性两代表(型0/型1)。
通用候选 = 在全部上下文都通过硬过滤的 DR 变体;
优先级 = 最差名次(跨上下文 max rank)升序 → 平均分降序。
输出: JSON + 主靶标合成 FASTA(通用候选)。

诚实边界: "通用元件"主张的计算边界——三上下文均同 DR 口径、同打分参数;
型空间覆盖仅 3/3 型采样, 未穷尽 spacer 多样性; 湿实验前的优先级建议,
非活性保证(见路径A校准)。
"""
import argparse
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def load_context(top_json):
    tag = os.path.basename(top_json).replace(".top.json", "")
    d = json.load(open(top_json, encoding="utf-8"))
    wt_dr = d["top"][0]["dr_seq"]
    rows = []
    with open(top_json.replace(".top.json", ".variants.csv"), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["passed"] in ("True", "true", "1"):
                rows.append({"desc": r["desc"], "score": float(r["score"]),
                             "construct_dna": r["construct_dna"], "ddG_dr": float(r["ddG_dr"])})
    rows.sort(key=lambda r: -r["score"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return tag, wt_dr, d["spacer_fixed_dna"], {r["desc"]: r for r in rows}, len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--contexts", nargs="+", default=[
        os.path.join(ROOT, "data", "tp53_r248q_zengdr.v6.top.json"),
        os.path.join(ROOT, "data", "cross_spacer.t0_R248Qg9.top.json"),
        os.path.join(ROOT, "data", "cross_spacer.t1_MYCg1.top.json")])
    ap.add_argument("--max-worst-rank", type=int, default=40)
    ap.add_argument("--cross-audit", default=os.path.join(ROOT, "data", "cross_audit.json"),
                    help="v1.6 交叉配对回溯审计(存在则自动注记各上下文拦截情况; 传空字符串跳过)")
    ap.add_argument("--out-prefix", default=None)
    args = ap.parse_args()

    ctxs = []
    wt_drs = set()
    for tj in args.contexts:
        tag, wt_dr, spacer, by_desc, n_passed = load_context(tj)
        print("上下文 %-28s spacer=%s WT_DR=%s 通过 n=%d" % (tag, spacer[:24], wt_dr, n_passed))
        wt_drs.add(wt_dr)
        ctxs.append({"tag": tag, "wt_dr": wt_dr, "spacer": spacer,
                     "by_desc": by_desc, "n_passed": n_passed})
    if len(wt_drs) != 1:
        raise SystemExit("各上下文 WT DR 口径不一致: %s" % wt_drs)

    common = set(ctxs[0]["by_desc"])
    for c in ctxs[1:]:
        common &= set(c["by_desc"])
    print("\n跨上下文全部通过硬过滤的通用候选: %d / %d(主上下文通过数)" % (
        len(common), ctxs[0]["n_passed"]))

    out_rows = []
    for desc in common:
        ranks = [c["by_desc"][desc]["rank"] for c in ctxs]
        scores = [c["by_desc"][desc]["score"] for c in ctxs]
        main_r = ctxs[0]["by_desc"][desc]
        out_rows.append({
            "desc": desc, "worst_rank": max(ranks),
            "ranks_per_context": {c["tag"]: c["by_desc"][desc]["rank"] for c in ctxs},
            "mean_score": round(sum(scores) / len(scores), 4),
            "min_score": round(min(scores), 4),
            "ddG_dr": main_r["ddG_dr"],
            "construct_dna_main_target": main_r["construct_dna"],
        })
    out_rows.sort(key=lambda r: (r["worst_rank"], -r["mean_score"]))

    shown = [r for r in out_rows if r["worst_rank"] <= args.max_worst_rank]

    # v1.6 交叉配对审计自动注记(数据来自 cross_audit.json, 非人工写入)
    cross_note = None
    if args.cross_audit and os.path.isfile(args.cross_audit):
        audits = {a["source"].replace(".top.json", ""): a for a in
                  json.load(open(args.cross_audit, encoding="utf-8"))}
        affected = {}
        for a in audits.values():
            tag = a["source"].replace(".top.json", "")
            for row in a.get("rows", []):
                if row["cross_nt"] > a["wt_cross_nt"]:
                    affected.setdefault(row["desc"], []).append(tag)
        for r in out_rows:
            tag_hits = affected.get(r["desc"])
            if tag_hits:
                r["v16_cross_filtered_in"] = sorted(set(tag_hits))
        if affected:
            n_uni = sum(1 for r in out_rows if r.get("v16_cross_filtered_in"))
            cross_note = (
                "v1.6 交叉配对硬过滤回溯审计(cross_audit.json): %d 条通用候选在部分上下文 "
                "cross_nt>WT 不再通过 v1.6(含 TOP-12: %s); 这些候选的跨型通用性主张降级为"
                "\"主靶标/型0 口径通过, 型1 口径结构保持存疑\", 优先级仅供合成排序" % (
                    n_uni, ", ".join(sorted(set(k for k in affected))[:6])))
    print("\n优先级清单(最差名次 <= %d, 共 %d 条):" % (args.max_worst_rank, len(shown)))
    print("%-16s %6s %10s %8s  %s" % ("DR_variant", "worst", "mean_score", "ddG_dr", "ranks"))
    for r in shown[:20]:
        rk = "/".join(str(r["ranks_per_context"][c["tag"]]) for c in ctxs)
        print("%-16s %6d %10.3f %8.1f  %s" % (r["desc"], r["worst_rank"],
                                              r["mean_score"], r["ddG_dr"], rk))

    prefix = args.out_prefix or os.path.join(ROOT, "data", "synthesis_priority")
    json.dump({"wt_dr": ctxs[0]["wt_dr"], "contexts": [
                  {"tag": c["tag"], "spacer": c["spacer"], "n_passed": c["n_passed"]}
                  for c in ctxs],
               "n_universal": len(out_rows),
               "rules": "worst_rank asc, mean_score desc; synthesis = main-target constructs",
               "v16_cross_caveat": cross_note,
               "priority": out_rows},
              open(prefix + ".json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    with open(prefix + ".fasta", "w", encoding="utf-8") as f:
        for r in shown:
            f.write(">%s|worst_rank=%d|%s\n%s\n" % (
                r["desc"], r["worst_rank"],
                "|".join("%s#%d" % (t, r["ranks_per_context"][t]) for t in
                         [c["tag"] for c in ctxs]),
                r["construct_dna_main_target"]))
    print("\n输出 -> %s.json / .fasta (FASTA 含最差名次<=%d 的通用候选, 主靶标构建)" % (
        prefix, args.max_worst_rank))


if __name__ == "__main__":
    main()
