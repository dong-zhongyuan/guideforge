"""MYCg1(型1)补偿突变实验设计(联审方案 #5, 因果检验; v2 去混杂, 2026-09)。

逻辑: t1_MYCg1 的 WT 游离态存在 10 对 DR-spacer 侵占配对(最长连续段 6 对,
Liao 2018 型; 非 10 对连续)。设计三臂:
  A break      : 侵占螺旋 3 个 spacer 位点换成不配对碱基(尽量保 GC) -> 预测活性恢复
  B break+comp : A 基础上再突变对应 DR 位点恢复配对 -> 预测活性再降
  C WT         : 原样对照
关键控制: spacer 改动会改变靶向——每条 spacer 突变构建必须配对其同源靶序列
(target = revcomp(spacer), 附 MYC 原始 PFS 口径注释), 以分离折叠效应与靶向效应。

v2 去混杂(round-3 评审 common concern 6 / bravo #4): v1 的 B 臂补偿位点
(DR 11/14/17, 1基)中 14/17 落在天然茎 3' 臂(茎对 (9,14)/(6,17)), 导致 B 臂
DR-only 折叠全毁(MFE 0.0, p_fold 0.0, bp_dist 5)——A vs B 同时改变"侵占配对"
与"天然茎完整性"两个变量, 无法归因单一折叠竞争。v2 改为枚举设计:
  补偿位点只允许落在 WT DR-only 折叠中不配对的侵占位点(环区/尾端);
  枚举所有 (3 位点组合 x spacer 替换碱基) 并保留满足全部约束的可行解:
    1) B 臂 DR-only 茎完整概率 >= STEM_PFOLD_MIN_RATIO * WT 参照;
    2) B 臂 DR-only MFE 结构与 WT 的 bp_distance <= MAX_STEM_BP_DIST;
    3) B 臂全长 inv_max_run 恢复至 WT 水平, A 臂低于 WT。
  评分取茎完整概率最高者(同分取位点沿螺旋跨度最大者)。
若枚举不到可行解则回退 v1 位点选择, 面板仍产出, 但 design_mode 记为
"joint_stem_x_invasion", prediction 如实重述为"茎x侵占联合效应", 不再作
单一折叠竞争归因。
脚本内自动验证(断言): A/B/C 各臂 inv_max_run 应分别为 下降/恢复/原值;
deconfounded 模式下另断言 B 臂茎完整性(p_fold 比率与 bp_dist), 不达则报错。

输出: data/ivt_compensatory_panel.json + .fasta (构建 + 同源靶, 可直接合成)
运行: python scripts/crrna_compensatory_design.py
"""
import argparse
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
PAIRABLE = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"),
            ("G", "T"), ("T", "G")}  # 含 G-T 摆动(U 的 DNA 字母)
SAME_GROUP = {"A": "T", "T": "A", "G": "C", "C": "G"}

STEM_PFOLD_MIN_RATIO = 0.8   # B 臂 DR-only 茎完整概率不得低于 WT 参照的 80%
MAX_STEM_BP_DIST = 1         # B 臂 DR-only MFE 结构与 WT 允许的碱基对距离上限


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def cross_helix_pairs(dr, sp):
    full = core.to_rna(dr + sp)
    ss, _ = core.fold(full)
    stack, cross = [], []
    for i, ch in enumerate(ss):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            a, b = sorted((i, j))
            if (a < len(dr)) != (b < len(dr)):
                dr_pos = a if a < len(dr) else b
                sp_pos = b if a < len(dr) else a
                cross.append((dr_pos, sp_pos))
    cross.sort()
    return cross, full


def spacer_candidates(drb, spb):
    """spacer 突变候选: 换为不与 drb 配对的碱基(优先保 GC 的同组互换)。"""
    return [c for c in (SAME_GROUP[spb], COMP[spb], COMP[SAME_GROUP[spb]])
            if c != spb and (drb, c) not in PAIRABLE]


def evaluate_combo(dr, sp, combo, pick, wt_dr_ss, stem_pairs, wt_run):
    """对一组 (位点组合, spacer 替换) 构建 A/B 两臂并计算去混杂指标。"""
    sp_a = list(sp)
    dr_b = list(dr)
    for di, sj0, cand in pick:
        sp_a[sj0] = cand
        dr_b[di] = COMP[cand]
    sp_a, dr_b = "".join(sp_a), "".join(dr_b)
    dr_b_ss, dr_b_mfe = core.fold(core.to_rna(dr_b))
    p_fold_b = core.stem_intact_prob(core.to_rna(dr_b), stem_pairs)
    bp_dist_b = RNA.bp_distance(wt_dr_ss, dr_b_ss)
    run_a = core.inv_max_run(core.to_rna(dr + sp_a), len(dr))
    run_b = core.inv_max_run(core.to_rna(dr_b + sp_a), len(dr_b))
    return {"sp_a": sp_a, "dr_b": dr_b, "dr_b_mfe": dr_b_mfe,
            "p_fold_b_dr_only": p_fold_b, "bp_dist_b_dr_only": bp_dist_b,
            "run_a": run_a, "run_b": run_b}


def enumerate_b_arm(dr, sp, cross, wt_dr_ss, stem_pairs, wt_dr_pfold, wt_run):
    """枚举茎安全补偿设计; 返回 (best, n_feasible)。无可行解返回 (None, 0)。"""
    stem_pos = {p for pair in stem_pairs for p in pair}
    safe = [c for c in cross if c[0] not in stem_pos]
    best, n_feasible = None, 0
    seen = set()
    for combo in itertools.combinations(safe, 3):
        cand_lists = []
        for di, sj in combo:
            opts = spacer_candidates(dr[di], sp[sj - len(dr)])
            if not opts:
                break
            cand_lists.append([(di, sj - len(dr), c) for c in opts])
        else:
            for pick in itertools.product(*cand_lists):
                key = tuple(sorted((di, c) for di, _, c in pick))
                if key in seen:
                    continue
                seen.add(key)
                r = evaluate_combo(dr, sp, combo, pick, wt_dr_ss, stem_pairs, wt_run)
                ok = (r["p_fold_b_dr_only"] >= STEM_PFOLD_MIN_RATIO * wt_dr_pfold
                      and r["bp_dist_b_dr_only"] <= MAX_STEM_BP_DIST
                      and r["run_b"] >= wt_run and r["run_a"] < wt_run)
                if not ok:
                    continue
                n_feasible += 1
                # 评分: 茎完整概率优先, 同分取位点沿侵占螺旋跨度最大者
                span = max(di for di, _ in combo) - min(di for di, _ in combo)
                score = (r["p_fold_b_dr_only"], span, -r["bp_dist_b_dr_only"])
                if best is None or score > best[0]:
                    best = (score, combo, pick, r)
    return (best[1:] if best else None), n_feasible


def legacy_pick(dr, sp, cross):
    """v1 位点选择(螺旋中部间隔 3 对)——仅作枚举失败时的回退, 不去混杂。"""
    idx = sorted(set([len(cross) // 4, len(cross) // 2, 3 * len(cross) // 4]))[:3]
    combo = [cross[k] for k in idx]
    pick = []
    for di, sj in combo:
        sj0 = sj - len(dr)
        cand = spacer_candidates(dr[di], sp[sj0])[0]
        pick.append((di, sj0, cand))
    return combo, pick


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top-json", default=os.path.join(ROOT, "data", "cross_spacer.t1_MYCg1.top.json"))
    ap.add_argument("--out-prefix", default=os.path.join(ROOT, "data", "ivt_compensatory_panel"))
    args = ap.parse_args()

    d = json.load(open(args.top_json, encoding="utf-8"))
    wt_entry = d["top"][0]
    dr = wt_entry["dr_seq"]
    sp = d["spacer_fixed_dna"]
    cross, full = cross_helix_pairs(dr, sp)
    assert len(cross) >= 6, "WT 侵占螺旋不足 6 对, 与诊断不符"

    wt_dr_ss, wt_dr_mfe = core.fold(core.to_rna(dr))
    stem_pairs = core.stem_pairs_of(wt_dr_ss)
    wt_dr_pfold = core.stem_intact_prob(core.to_rna(dr), stem_pairs)
    wt_run = core.inv_max_run(full, len(dr))

    result, n_feasible = enumerate_b_arm(dr, sp, cross, wt_dr_ss, stem_pairs,
                                         wt_dr_pfold, wt_run)
    if result is not None:
        combo, pick, ev = result
        design_mode = "deconfounded"
    else:
        combo, pick = legacy_pick(dr, sp, cross)
        ev = evaluate_combo(dr, sp, combo, pick, wt_dr_ss, stem_pairs, wt_run)
        design_mode = "joint_stem_x_invasion"
        print("[警告] 茎安全补偿位点枚举无可行解, 回退 v1 位点; "
              "本面板只能作 茎x侵占 联合效应解释")

    notes = []
    for di, sj0, cand in pick:
        notes.append({"dr_pos_1based": di + 1, "dr_wt": dr[di],
                      "spacer_pos_1based": sj0 + 1, "spacer_wt": sp[sj0],
                      "spacer_mut": cand, "dr_compensatory": COMP[cand],
                      "dr_pos_in_wt_stem": any(di in p for p in stem_pairs)})
    sp_a, dr_b = ev["sp_a"], ev["dr_b"]

    arms = {
        "C_WT": (dr, sp),
        "A_break": (dr, sp_a),
        "B_break_compensate": (dr_b, sp_a),
    }
    print("%-22s %5s %5s %6s %6s %8s %7s" % ("arm", "GC_dr", "GC_sp", "x_nt",
                                             "inv_run", "DR_pfold", "DR_bpdist"))
    out = []
    for name, (dri, spi) in arms.items():
        cx, _ = cross_helix_pairs(dri, spi)
        run = core.inv_max_run(core.to_rna(dri + spi), len(dri))
        dri_ss, dri_mfe = core.fold(core.to_rna(dri))
        pfold = core.stem_intact_prob(core.to_rna(dri), stem_pairs)
        bpd = RNA.bp_distance(wt_dr_ss, dri_ss)
        print("%-22s %5d %5d %6d %6d %8.3f %7d" % (name, gc(dri), gc(spi),
                                                  len(cx) * 2, run, pfold, bpd))
        out.append({"arm": name, "dr_dna": dri, "spacer_dna": spi,
                    "construct_dna": dri + spi,
                    "target_rna_dna": revcomp(spi),  # 同源靶(spacer 互补)
                    "cross_nt": len(cx) * 2, "inv_max_run": run,
                    "gc_dr": gc(dri), "gc_sp": gc(spi),
                    "dr_only_mfe_kcal": round(dri_mfe, 2),
                    "dr_only_stem_pfold": round(pfold, 4),
                    "dr_only_bp_dist_vs_wt": bpd})
    runs = {o["arm"]: o["inv_max_run"] for o in out}
    assert runs["A_break"] < wt_run, "break 臂未压低侵占螺旋"
    assert runs["B_break_compensate"] >= runs["A_break"], "补偿臂未恢复配对"
    if design_mode == "deconfounded":
        b = next(o for o in out if o["arm"] == "B_break_compensate")
        assert runs["B_break_compensate"] >= wt_run, "补偿臂未恢复至 WT 侵占水平"
        assert b["dr_only_stem_pfold"] >= STEM_PFOLD_MIN_RATIO * wt_dr_pfold, \
            "B 臂天然茎完整概率低于 WT 参照的 %.0f%%, 茎完整性不达标" % (STEM_PFOLD_MIN_RATIO * 100)
        assert b["dr_only_bp_dist_vs_wt"] <= MAX_STEM_BP_DIST, \
            "B 臂 DR-only 结构偏离 WT(bp_dist > %d), 茎完整性不达标" % MAX_STEM_BP_DIST

    if design_mode == "deconfounded":
        prediction = ("A(破坏侵占, 靶向配对保持, 天然茎完整)活性应高于 WT; "
                      "B(恢复侵占配对, 且补偿位点全部避开天然茎, DR-only 茎完整概率 "
                      "%.2f vs WT %.2f, bp_dist %d)活性应回落; 两臂天然茎均完整, "
                      "A vs B 差异可归于折叠竞争单一变量; 每臂须用各自同源靶。 "
                      "局限: 茎完整性为 MFE/配分函数计算口径, IVT 验证前为计算预测"
                      % (ev["p_fold_b_dr_only"], wt_dr_pfold, ev["bp_dist_b_dr_only"]))
    else:
        prediction = ("A(破坏侵占)活性应高于 WT; B(恢复侵占配对)活性应回落。 "
                      "注意: 本面板 B 臂补偿同时破坏天然茎(DR-only p_fold %.2f vs WT %.2f, "
                      "bp_dist %d), A vs B 只能作 茎x侵占 联合效应解释, "
                      "不能归因于单一折叠竞争; 每臂须用各自同源靶"
                      % (ev["p_fold_b_dr_only"], wt_dr_pfold, ev["bp_dist_b_dr_only"]))

    payload = {"source": os.path.basename(args.top_json),
               "design_mode": design_mode,
               "design_notes": notes,
               "wt_dr_only": {"mfe_struct": wt_dr_ss,
                              "mfe_kcal": round(wt_dr_mfe, 2),
                              "stem_pfold": round(wt_dr_pfold, 4)},
               "stem_integrity_gate": {"p_fold_min_ratio_of_wt": STEM_PFOLD_MIN_RATIO,
                                       "max_bp_dist_vs_wt": MAX_STEM_BP_DIST,
                                       "n_feasible_designs": n_feasible},
               "prediction": prediction,
               "claim_scope": ("deconfounded 模式: B 臂补偿位点全部落在 WT DR-only 折叠的"
                               "非茎区(环区/尾端), 两臂天然茎均完整, A vs B 支持折叠竞争的"
                               "单一变量归因(计算口径, 待 IVT)" if design_mode == "deconfounded"
                               else "joint 模式: A vs B 为 茎x侵占 联合效应, 不可作单一归因"),
               "pfs_note": "MYC 原靶 PFS 口径沿用原设计; 突变臂靶序列为 spacer 的互补链",
               "arms": out}
    json.dump(payload, open(args.out_prefix + ".json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(args.out_prefix + ".fasta", "w") as f:
        for o in out:
            f.write(">%s|crRNA_DR+spacer|inv_run=%d\n%s\n" % (o["arm"], o["inv_max_run"], o["construct_dna"]))
            f.write(">%s|cognate_target|revcomp(spacer)\n%s\n" % (o["arm"], o["target_rna_dna"]))
    print("输出 -> %s.json / .fasta (构建+同源靶各 3 条)" % args.out_prefix)


def gc(s):
    return s.count("G") + s.count("C")


if __name__ == "__main__":
    main()
