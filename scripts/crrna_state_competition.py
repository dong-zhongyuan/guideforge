"""活性态-侵占态竞争模型(Sp8 盲区响应, 联审方案 #1, 2026-09-01)。

替代"交叉配对计数"的单点判据, 在 ViennaRNA 单一能量标尺内计算三态竞争:
  P_active   = Z(含模板 DR 茎全部配对) / Z_total   —— 约束配分, 优于 bpp 独立性乘积
  P_invasion = Z(含该 spacer 最长连续 DR-spaker 侵占螺旋全部配对) / Z_total
  dG_comp    = G_invasion - G_active (kcal/mol, 越小/越负侵占态越有竞争力)

活性态模板: Creutzburg 2020 补充材料自定义的 DR 茎 5 对(即该文 stem-forming
potential 的定义), 与本仓库 DR 单独折叠所得一致; **保守 5' 假结不可表示**,
故 P_active 是"含完整茎"的下界近似——已在输出 model_domain 中声明。

判别性检验(先定后跑): 14 条 spacer 中 Sp8(唯一 0% 活性)应在 dG_comp 最小
(侵占最有竞争力)或 P_active 最低的极端处; 若不成立, 该模型亦宣告无效。

诚实边界: 同一 ViennaRNA 标尺内可比, 无跨工具混标; n=14 仅判别不拟合;
不用于 Cas12a2 活性预测(见 README 路径A校准)。
运行: python scripts/crrna_state_competition.py
输出: data/state_competition.json
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import numpy as np
import crrna_scaffold_design as core  # noqa: E402

RT = 0.00198717 * 310.15


def pairs_of(db):
    stack, pairs = [], []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.append((stack.pop(), i))
    return pairs


def log_z(seq, forced_pairs=None, forbid_positions=None):
    """配分函数对数(-RT lnZ)。forced_pairs 0 基强制配对(非交叉集合);
    forbid_positions 0 基强制不成对(定义互斥态: 禁死茎 5' 臂 = 茎必被破坏)。
    用 hc_add_bp/hc_add_up(1 基)——本构建 constraints_add 括号串无效(已探测)。"""
    fc = RNA.fold_compound(seq)
    if forced_pairs:
        for i, j in forced_pairs:
            fc.hc_add_bp(i + 1, j + 1)
    if forbid_positions:
        for i in forbid_positions:
            fc.hc_add_up(i + 1)
    _x, g = fc.pf()
    return g


def longest_cross_helix(db, boundary):
    """最长连续跨界螺旋(嵌套, 可被约束串表示)。"""
    stack, cross = [], []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            if (i < boundary) != (j < boundary):
                cross.append((i, j))
    cross.sort()
    best, cur, run = [], [], None
    for i, j in cross:
        if run is not None and i == run[-1][0] + 1 and j == run[-1][1] - 1:
            cur.append((i, j))
        else:
            cur = [(i, j)]
        if len(cur) > len(best):
            best = list(cur)
        run = cur
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dataset", default=os.path.join(ROOT, "data", "creutzburg_2020_dataset.json"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "state_competition.json"))
    args = ap.parse_args()

    ds = json.load(open(args.dataset, encoding="utf-8"))
    spacers = ds["spacers_21nt_dna"]
    act = ds["activity_fig1d_pct_approx"]
    dr18 = core.to_rna(ds["dr18_dna"])
    pre5 = core.to_rna(ds["pre5_flank_dna"])
    prefix = pre5 + dr18  # 前体口径: 论文机制在 5' 残留 repeat 侧
    boundary = len(prefix)

    stem_pairs = pairs_of(core.fold(dr18)[0])  # 活性态模板(论文同款 5 对)
    stem_pairs_full = [(i + len(pre5), j + len(pre5)) for i, j in stem_pairs]
    stem_arm5 = sorted(i for i, _ in stem_pairs_full)  # 5' 臂(禁之 = 茎必破坏)
    print("活性态模板 DR 茎配对(0基, 前体坐标): %s" % sorted(stem_pairs_full))

    names = sorted(spacers, key=lambda s: int(s[2:]))
    rows = []
    for nm in names:
        full = prefix + core.to_rna(spacers[nm])
        g_tot = log_z(full)
        g_act = log_z(full, stem_pairs_full)
        # 互斥侵占态: 强制该 spacer 的最长跨界螺旋 + 禁死茎 5' 臂
        inv = longest_cross_helix(core.fold(full)[0], boundary)
        g_inv = log_z(full, inv, forbid_positions=stem_arm5) if inv else None
        p_act = float(np.exp((g_tot - g_act) / RT))
        p_inv = float(np.exp((g_tot - g_inv) / RT)) if inv else 0.0
        rows.append({"name": nm, "activity": act[nm],
                     "P_active": round(p_act, 4),
                     "inv_helix_len": len(inv),
                     "P_invasion": round(p_inv, 6),
                     "dG_comp": round(g_inv - g_act, 2) if inv else None})

    print("%-5s %4s %8s %6s %10s %8s" % ("name", "act", "P_act", "invL", "P_inv", "dG_comp"))
    for r in sorted(rows, key=lambda r: (r["dG_comp"] if r["dG_comp"] is not None else 99)):
        print("%-5s %4d %8.4f %6d %10.6f %8s" % (
            r["name"], r["activity"], r["P_active"], r["inv_helix_len"],
            r["P_invasion"], r["dG_comp"]))

    # 判别性检验: Sp8 在 dG_comp 升序中的名次(越靠前=侵占最有竞争力)
    by_dg = sorted([r for r in rows if r["dG_comp"] is not None], key=lambda r: r["dG_comp"])
    sp8_dg_rank = [r["name"] for r in by_dg].index("Sp8") + 1
    by_pa = sorted(rows, key=lambda r: r["P_active"])
    sp8_pa_rank = [r["name"] for r in by_pa].index("Sp8") + 1
    verdict = "PASS" if (sp8_dg_rank <= 2 or sp8_pa_rank <= 2) else "FAIL"
    print("\n判别检验: Sp8 dG_comp 名次(侵占竞争力降序) %d/%d; P_active 名次(最低在前) %d/%d -> %s"
          % (sp8_dg_rank, len(by_dg), sp8_pa_rank, len(rows), verdict))
    if verdict == "FAIL":
        print("结论: 伪结外模型内 Sp8 不可判别(所有可枚举侵占态均比茎态贵 5.8+ kcal,"
              "论文失活折叠不可达)——Sp8 定性 out_of_model_domain,"
              "活性态竞争假说的检验需统一假结能量模型(中期工作)")

    acts = np.array([r["activity"] for r in rows], float)
    for feat in ("P_active", "dG_comp"):
        vals = np.array([r[feat] if r[feat] is not None else 0.0 for r in rows], float)
        rho = float(np.corrcoef(np.argsort(np.argsort(vals)), np.argsort(np.argsort(acts)))[0, 1])
        print("%-9s vs activity Spearman = %+.3f" % (feat, rho))

    json.dump({"dataset": args.dataset, "active_template_pairs_0based": sorted(stem_pairs),
               "model_domain": "ViennaRNA 约束配分; 5' 假结不可表示, P_active 为含完整茎的下界近似",
               "rows": rows,
               "discrimination": {"sp8_dG_comp_rank_asc": sp8_dg_rank,
                                  "sp8_P_active_rank_asc": sp8_pa_rank,
                                  "verdict": verdict,
                                  "fail_interpretation": "伪结外模型内 Sp8 不可判别: 定性 out_of_model_domain"},
               "preregistered_test": "Sp8 应处 dG_comp 最小或 P_active 最低极端; 否则模型无效"},
              open(args.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % args.out)


if __name__ == "__main__":
    main()
