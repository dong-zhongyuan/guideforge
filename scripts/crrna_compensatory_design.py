"""MYCg1(型1)补偿突变实验设计(联审方案 #5, 因果检验, 2026-09-01)。

逻辑: t1_MYCg1 的 WT 游离态存在 10 对连续 DR-spacer 侵占螺旋(Liao 2018 型)。
设计三臂:
  A break      : 侵占螺旋中部 3 个 spacer 位点换成不配对碱基(尽量保 GC) -> 预测活性恢复
  B break+comp : A 基础上再突变对应 DR 位点恢复配对 -> 预测活性再降
  C WT         : 原样对照
关键控制: spacer 改动会改变靶向——每条 spacer 突变构建必须配对其同源靶序列
(target = revcomp(spacer), 附 MYC 原始 PFS 口径注释), 以分离折叠效应与靶向效应。
脚本内自动验证: 各臂 inv_max_run 应分别为 下降/恢复/原值, 不达则报错。

输出: data/ivt_compensatory_panel.json + .fasta (构建 + 同源靶, 可直接合成)
运行: python scripts/crrna_compensatory_design.py
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402

COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
PAIRABLE = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"),
            ("G", "T"), ("T", "G")}  # 含 G-T 摆动(U 的 DNA 字母)
SAME_GROUP = {"A": "T", "T": "A", "G": "C", "C": "G"}


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

    # 取螺旋中部间隔 3 对(索引 1,4,7 或居中三点)
    idx = sorted(set([len(cross) // 4, len(cross) // 2, 3 * len(cross) // 4]))[:3]
    chosen = [cross[k] for k in idx]

    def gc(s):
        return s.count("G") + s.count("C")

    sp_a = list(sp)
    dr_b = list(dr)
    notes = []
    for di, sj in chosen:  # di=DR 0基, sj=spacer 0基(全长坐标)
        sj0 = sj - len(dr)
        drb, spb = dr[di], sp_a[sj0]
        assert (drb, spb) in PAIRABLE
        # A: spacer 换成不与 drb 配对的碱基(优先保 GC 的同组互换)
        for cand in (SAME_GROUP[spb], COMP[spb], COMP[SAME_GROUP[spb]]):
            if cand != spb and (drb, cand) not in PAIRABLE:
                break
        sp_a[sj0] = cand
        # B: DR 突变恢复与 cand 的配对
        dr_b[di] = COMP[cand]
        notes.append({"dr_pos_1based": di + 1, "dr_wt": drb,
                      "spacer_pos_1based": sj0 + 1, "spacer_wt": spb,
                      "spacer_mut": cand, "dr_compensatory": COMP[cand]})
    sp_a = "".join(sp_a)
    dr_b = "".join(dr_b)

    arms = {
        "C_WT": (dr, sp),
        "A_break": (dr, sp_a),
        "B_break_compensate": (dr_b, sp_a),
    }
    print("%-22s %5s %5s %6s %6s" % ("arm", "GC_dr", "GC_sp", "x_nt", "inv_run"))
    out = []
    for name, (dri, spi) in arms.items():
        cx, _ = cross_helix_pairs(dri, spi)
        run = core.inv_max_run(core.to_rna(dri + spi), len(dri))
        print("%-22s %5d %5d %6d %6d" % (name, gc(dri), gc(spi), len(cx) * 2, run))
        out.append({"arm": name, "dr_dna": dri, "spacer_dna": spi,
                    "construct_dna": dri + spi,
                    "target_rna_dna": revcomp(spi),  # 同源靶(spacer 互补)
                    "cross_nt": len(cx) * 2, "inv_max_run": run,
                    "gc_dr": gc(dri), "gc_sp": gc(spi)})
    runs = {o["arm"]: o["inv_max_run"] for o in out}
    wt_run = runs["C_WT"]
    assert runs["A_break"] < wt_run, "break 臂未压低侵占螺旋"
    assert runs["B_break_compensate"] >= runs["A_break"], "补偿臂未恢复配对"

    payload = {"source": os.path.basename(args.top_json), "design_notes": notes,
               "prediction": "A(破坏侵占, 靶向配对保持)活性应高于 WT; B(恢复侵占配对)活性应回落;"
                             "A vs B 差异 = 折叠竞争的因果证据; 每臂须用各自同源靶",
               "pfs_note": "MYC 原靶 PFS 口径沿用原设计; 突变臂靶序列为 spacer 的互补链",
               "arms": out}
    json.dump(payload, open(args.out_prefix + ".json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(args.out_prefix + ".fasta", "w") as f:
        for o in out:
            f.write(">%s|crRNA_DR+spacer|inv_run=%d\n%s\n" % (o["arm"], o["inv_max_run"], o["construct_dna"]))
            f.write(">%s|cognate_target|revcomp(spacer)\n%s\n" % (o["arm"], o["target_rna_dna"]))
    print("输出 -> %s.json / .fasta (构建+同源靶各 3 条)" % args.out_prefix)


if __name__ == "__main__":
    main()
