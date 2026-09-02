"""跨型候选骨架族生成(策划案 V3 §4.2: 四取向引导, 2026-09-01)。

第一轮生成不以单一最优为目标, 按四种优化取向各引导一轮主管线, 每取向取
最优通过变体, 组成含 WT 的跨型候选族(目标 8 条, 不足如实记录):

  low_threshold  低激活阈值取向: 种子区暴露+茎稳定化奖赏放大 (w_seed/w_stab ↑)
  high_stringency 高错配严格取向: 保守窗+蛋白接触惩罚放大 (w_cons3/w_contact ↑)
  robust         结构稳健取向: 结构距离锁零 (max-bp-dist 0) + 系综多样性惩罚
  balanced       平衡取向: 主管线默认权重

运行: python scripts/crrna_orientation_library.py
输出: data/orientation_library.json / .fasta (WT + 各取向代表; 供体外矩阵)
"""
import argparse
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")
SPACER = "GTTCATGCCGCCCATGCAGGAACT"

ORIENTATIONS = {
    "low_threshold": ["--w-seed", "0.5", "--w-stab", "0.5", "--stab-dr-only"],
    "high_stringency": ["--w-cons3", "0.8", "--w-contact", "1.5"],
    "robust": ["--max-bp-dist", "0", "--w-ens", "0.3"],
    "balanced": [],
}

# 各取向的选择键(在通过硬过滤的变体内按取向目标挑代表, 而非同一综合分)
def pick_key(orientation):
    if orientation == "low_threshold":
        return lambda r: (-(float(r["ddG_dr"])), float(r["seed_up"]))
    if orientation == "high_stringency":
        return lambda r: (-int(r["n_mut"]), -(float(r["contact_sum"]) +
                                              float(r["cons3_frac"])))
    if orientation == "robust":
        return lambda r: (-(int(r["bp_dist"])), -float(r["d_ens"]))
    return lambda r: float(r["score"])  # balanced


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--effector", default="cas12a2_zeng2026")
    ap.add_argument("--spacer", default=SPACER)
    args = ap.parse_args()

    picked = []
    for name, extra in ORIENTATIONS.items():
        prefix = os.path.join(DATA, "orientation.%s" % name)
        cmd = [sys.executable, os.path.join(HERE, "crrna_scaffold_design.py"),
               "--effector", args.effector, "--spacer", args.spacer,
               "--topk", "12", "--use-covariation"] + extra + ["--out-prefix", prefix]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print("[%s] 管线失败: %s" % (name, r.stderr[-300:]))
            continue
        with open(prefix + ".variants.csv", newline="", encoding="utf-8") as f:
            rows = [r for r in csv.DictReader(f)
                    if r["passed"] in ("True", "true", "1")
                    and r["desc"] != "WT" and int(r["n_mut"]) > 0]
        if not rows:
            print("[%s] 无通过变体" % name)
            continue
        rows.sort(key=pick_key(name), reverse=True)
        seen, chosen = set(), []
        for r in rows:
            if r["desc"] in seen:
                continue
            seen.add(r["desc"])
            chosen.append(r)
            if len(chosen) == 2:
                break
        for best in chosen:
            print("[%s] 代表: %s (ddG_dr=%s, n_mut=%s, bp_dist=%s, score=%s)" % (
                name, best["desc"], best["ddG_dr"], best["n_mut"],
                best["bp_dist"], best["score"]))
            picked.append({"orientation": name, "extra_args": extra,
                           "n_passed": len(rows), "best": best,
                           "wt_construct": None})

    library = []
    seen = set()
    for p in picked:
        b = p["best"]
        if b["desc"] not in seen:
            seen.add(b["desc"])
            library.append(p)
    # V3 口径补第 8 员: compensatory 代表(B_break_compensate, zengDR+G11C/G14C/G17C,
    # 源: data/ivt_compensatory_panel.json, 机制故事=侵占配对破坏-补偿因果臂)
    comp = json.load(open(os.path.join(DATA, "ivt_compensatory_panel.json"),
                          encoding="utf-8"))
    for arm in comp["arms"]:
        if arm["arm"] == "B_break_compensate":
            dr19 = arm["dr_dna"][:19]
            library.append({
                "orientation": "compensatory",
                "extra_args": ["(外部设计: 补偿臂, 见 ivt_compensatory_panel.json)"],
                "n_passed": 1, "wt_construct": None,
                "best": {"desc": "B_break_compensate", "dr_seq": dr19,
                         "construct_dna": dr19 + args.spacer,
                         "score": "0", "ddG_dr": "NA",
                         "note": "DR-spacer 侵占配对补偿设计, 非四取向管线产物"}})
            seen.add("B_break_compensate")
            print("[compensatory] 代表: B_break_compensate (DR=%s)" % dr19)
            break
    wt_construct = None
    for name in ORIENTATIONS:
        tj = os.path.join(DATA, "orientation.%s.top.json" % name)
        if os.path.isfile(tj):
            wt_construct = json.load(open(tj, encoding="utf-8"))["top"][0]["construct_dna"]
            break
    out = os.path.join(DATA, "orientation_library.json")
    json.dump({"effector": args.effector, "spacer": args.spacer,
               "rule": "四取向各取前2去重 + WT + compensatory 代表 = 8 员跨型"
                       "候选族(V3 §4.2 口径: 8 骨架 x 4 靶标 = 32 组合)",
               "n_library_with_wt": len(library) + 1,
               "orientations": picked + [library[-1]]},
              open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    with open(os.path.join(DATA, "orientation_library.fasta"), "w") as f:
        f.write(">WT|orientation_reference\n%s\n" % wt_construct)
        for p in library:
            b = p["best"]
            f.write(">%s|%s|score=%s|ddG_dr=%s\n%s\n" % (
                p["orientation"], b["desc"], b["score"], b["ddG_dr"],
                b["construct_dna"]))
    print("\n候选族(含 WT): %d 条 -> data/orientation_library.json/.fasta" % (len(library) + 1))


if __name__ == "__main__":
    main()
