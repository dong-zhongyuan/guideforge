"""Pareto 前沿可视化(证据链补强 #6, 答辩 PPT 用, 用户 2026-09-01 指定)。

四维编码到一张散点图:
  x = ddG_dr (DR 单独折叠 ΔΔG, 茎稳定化轴, 越负越稳)
  y = spacer_unpaired (spacer 游离度, 越高越好)
  颜色 = contact_sum + hbond_loss (蛋白接触+氢键损失惩罚, 越浅越好)
  点大小 = bp_dist (与 WT 结构距离, 越小越好)
叠加元素:
  灰点 = 未通过硬过滤的变体(展示过滤器效果); 彩点 = 通过变体;
  金边 = Pareto 前沿; 红星描边 = 茎稳定化候选(ddG_dr <= -1);
  黑菱形 = WT; 标注 = TOP-12。
英文标签(服务器无中文字体)。输出 SVG + PNG。
"""
import argparse
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def load_rows(csv_path):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append({
                "desc": r["desc"], "passed": r["passed"] in ("True", "true", "1"),
                "pareto": r.get("pareto") in ("True", "true", "1"),
                "ddG_dr": float(r["ddG_dr"]), "sp_up": float(r["spacer_unpaired"]),
                "pen": float(r["contact_sum"]) + float(r["hbond_loss"]),
                "bp_dist": int(r["bp_dist"]), "score": float(r["score"]),
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top-json", default=os.path.join(ROOT, "data", "tp53_r248q_pdbdr.v2.top.json"))
    ap.add_argument("--stab-threshold", type=float, default=-1.0)
    ap.add_argument("--label-top", type=int, default=12)
    args = ap.parse_args()

    tag = os.path.basename(args.top_json).replace(".top.json", "")
    d = json.load(open(args.top_json, encoding="utf-8"))
    rows = load_rows(args.top_json.replace(".top.json", ".variants.csv"))
    wt_sp = d["wt"]["spacer_mean_unpaired"]
    tops = [e for e in d["top"] if e.get("desc") not in (None, "WT")][:args.label_top]
    top_names = {e["desc"] for e in tops}
    by_desc = {r["desc"]: r for r in rows}

    failed = [r for r in rows if not r["passed"]]
    passed = [r for r in rows if r["passed"]]

    fig, ax = plt.subplots(figsize=(9.5, 6.8), dpi=200)
    if failed:
        ax.scatter([r["ddG_dr"] for r in failed], [r["sp_up"] for r in failed],
                   s=8, c="#d8d8d8", alpha=0.55, linewidths=0,
                   label="filtered out (n=%d)" % len(failed))
    sc = ax.scatter([r["ddG_dr"] for r in passed], [r["sp_up"] for r in passed],
                    c=[r["pen"] for r in passed], cmap="viridis_r",
                    s=[16 + 26 * r["bp_dist"] for r in passed],
                    alpha=0.85, edgecolors="none",
                    label="passed (n=%d)" % len(passed))
    front = [r for r in passed if r["pareto"]]
    ax.scatter([r["ddG_dr"] for r in front], [r["sp_up"] for r in front],
               s=[16 + 26 * r["bp_dist"] for r in front], facecolors="none",
               edgecolors="#e6a700", linewidths=1.4,
               label="Pareto front (n=%d)" % len(front))
    stab = [r for r in passed if r["ddG_dr"] <= args.stab_threshold]
    ax.scatter([r["ddG_dr"] for r in stab], [r["sp_up"] for r in stab],
               marker="*", s=230, facecolors="none", edgecolors="#d62728",
               linewidths=1.6, label="stem-stabilized ddG_dr<=%.0f (n=%d)" % (
                   args.stab_threshold, len(stab)))
    ax.scatter([0.0], [wt_sp], marker="D", s=90, c="black", zorder=5,
               label="WT (spacer_up=%.2f)" % wt_sp)

    for e in tops:
        r = by_desc.get(e["desc"])
        if r:
            ax.annotate(e["desc"], (r["ddG_dr"], r["sp_up"]),
                        textcoords="offset points", xytext=(7, 4),
                        fontsize=7.5, color="#222222")

    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("protein contact + H-bond penalty (lower = better)", fontsize=9)
    ax.set_xlabel("DR-alone $\\Delta\\Delta$G (kcal/mol)  <-  stem stabilization", fontsize=10.5)
    ax.set_ylabel("spacer unpaired fraction (higher = better)", fontsize=10.5)
    ax.set_title("GuideForge variant landscape — %s\nsize = bp_dist to WT (smaller = better); gold edge = Pareto front; red star = stabilized" % tag,
                 fontsize=11)
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
    fig.tight_layout()

    for ext in ("svg", "png"):
        p = os.path.join(ROOT, "data", "pareto.%s.%s" % (tag, ext))
        fig.savefig(p, dpi=200)
        print("输出 -> %s" % p)


if __name__ == "__main__":
    main()
