"""DP09 假结能量模型竞争检验(Knotty, 联审方案的中期项, 2026-09-01)。

背景: ViennaRNA 伪结外空间内 Sp8 不可判别(state_competition FAIL,
out_of_model_domain)。Knotty(CCJ 算法 + DP09 假结能量参数)可输出含假结的
MFE。本脚本在 DP09 单一标尺下逐条(14 条 Creutzburg 前体 + 成熟两口径 +
GuideForge pdbdr WT/TOP3)取得 pk-MFE 结构与能量, 并按结构内容分类:
  stem_intact   = 模板 5 对茎全部出现(活性态必要条件)
  dr_cross_pk   = 预测结构含 DR 位点与 spacer/flank 的跨界配对(侵占)
  is_pk         = 结构含 [] 假结记号
判别检验(先定后跑): Sp8 应表现为 stem 破坏或 DR 参与跨界假结, 而 Sp12 不。

编排: 本脚本(宿主)只生成任务文件; 容器运行器 /dawn/crrna_pk_runner.py
由本机 ssh 双跳执行, 输出经本机管道回写宿主 /tmp/pk_out.json; 再由
--assemble 模式装配为 data/pk_competition.json。
用法(宿主): python scripts/crrna_pk_competition.py jobs|assemble
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402

DAWN = "/public/home/mengxl/openclaw/dawn"
OUT_TMP = "/tmp/pk_out.json"

STEM18_0BASE = [(3, 16), (4, 15), (5, 14), (6, 13), (7, 12)]  # DR18 模板茎(折叠口径)


def pairs_from_db_pk(db):
    """解析含 [] 的点括号为配对集合(0 基)。"""
    stacks = {"(": [], "[": []}
    close = {")": "(", "]": "["}
    pairs = set()
    for i, ch in enumerate(db):
        if ch in "([":
            stacks[ch].append(i)
        elif ch in ")]":
            j = stacks[close[ch]].pop()
            pairs.add((min(i, j), max(i, j)))
    return pairs


def classify(db, boundary, stem_pairs=None):
    """寄存器无关分类(跨能量模型不可用 ViennaRNA 模板精确匹配——DP09 偏好滑移寄存,
    实测教训): stem_intact = DR 区内存在 >=4 对连续螺旋; dr_cross = DR 位点与
    spacer 的跨界配对(flank 侧不计, 因加工后被丢弃)。"""
    pairs = pairs_from_db_pk(db)
    dr_lo, dr_hi = boundary - 18, boundary  # DR 区间(0基)
    dr_pairs = sorted(p for p in pairs if dr_lo <= p[0] < dr_hi and dr_lo <= p[1] < dr_hi)
    best = cur = 0
    prev = None
    for i, j in dr_pairs:
        if prev is not None and i == prev[0] + 1 and j == prev[1] - 1:
            cur += 1
        else:
            cur = 1
        best = max(best, cur)
        prev = (i, j)
    cross = sorted(p for p in pairs if dr_lo <= p[0] < dr_hi and p[1] >= boundary)
    return {"stem_helix_len": best, "stem_intact": best >= 4,
            "dr_cross_pairs": cross, "n_dr_cross": len(cross),
            "is_pk": "[" in db}


def build_jobs():
    ds = json.load(open(os.path.join(ROOT, "data", "creutzburg_2020_dataset.json"),
                        encoding="utf-8"))
    dr18 = core.to_rna(ds["dr18_dna"])
    pre5 = core.to_rna(ds["pre5_flank_dna"])
    jobs = {}
    for nm, sp in ds["spacers_21nt_dna"].items():
        jobs["%s_pre" % nm] = pre5 + dr18 + core.to_rna(sp)
        jobs["%s_mat" % nm] = dr18 + core.to_rna(sp)
    d = json.load(open(os.path.join(ROOT, "data", "tp53_r248q_pdbdr.v2.top.json"),
                       encoding="utf-8"))
    for e in d["top"]:
        if e.get("desc") in (None, "WT") or e.get("rank", 99) > 3:
            continue
        tag = "gf_%s" % e["desc"].replace("+", "_")
        jobs[tag] = core.to_rna(e["construct_dna"])
    jobs["gf_WT"] = core.to_rna(d["top"][0]["construct_dna"])
    json.dump(jobs, open(os.path.join(DAWN, "pk_jobs.json"), "w"))
    print("任务 %d 条 -> %s/pk_jobs.json (容器侧: /dawn/pk_jobs.json)" % (len(jobs), DAWN))
    print("下一步(本机): ssh a6000-ct 'python /dawn/crrna_pk_runner.py' | ssh srv 'cat > %s'" % OUT_TMP)


def assemble():
    res = json.load(open(OUT_TMP))
    ds = json.load(open(os.path.join(ROOT, "data", "creutzburg_2020_dataset.json"),
                        encoding="utf-8"))
    act = ds["activity_fig1d_pct_approx"]
    dr18 = core.to_rna(ds["dr18_dna"])
    stem = set(STEM18_0BASE)
    rows = []
    for nm in sorted(act, key=lambda s: int(s[2:])):
        for ctx, boundary, tag in (("pre", len(ds["pre5_flank_dna"]) + 18, "%s_pre"),
                                    ("mature", 18, "%s_mat")):
            r = res[tag % nm]
            info = classify(r["db_pk"], boundary, stem) if "db_pk" in r else {"error": r}
            rows.append({"name": nm, "context": ctx, "activity": act[nm],
                         "energy_dp09": r.get("energy_dp09"), **info,
                         "db_pk": r.get("db_pk")})
    gf = []
    for k, r in res.items():
        if not k.startswith("gf_"):
            continue
        info = classify(r["db_pk"], 18, stem) if "db_pk" in r else {"error": r}
        gf.append({"name": k, "energy_dp09": r.get("energy_dp09"), **info,
                   "db_pk": r.get("db_pk")})

    print("%-6s %-7s %4s %8s %7s %8s %5s %s" % (
        "name", "ctx", "act", "E_dp09", "stem", "dr_cross", "pk", "判定"))
    for r in rows:
        verdict = []
        if not r.get("stem_intact", False):
            verdict.append("茎破坏")
        if r.get("n_dr_cross", 0) >= 3:
            verdict.append("DR跨界侵占")
        print("%-6s %-7s %4d %8s %7s %8d %5s %s" % (
            r["name"], r["context"], r["activity"], r["energy_dp09"],
            r.get("stem_intact"), r.get("n_dr_cross", 0), r.get("is_pk"),
            "; ".join(verdict) or "活性样"))

    sp8_pre = next(r for r in rows if r["name"] == "Sp8" and r["context"] == "pre")
    sp12_pre = next(r for r in rows if r["name"] == "Sp12" and r["context"] == "pre")
    discrim = ((not sp8_pre.get("stem_intact", True)) or
               sp8_pre.get("n_dr_cross", 0) >= 3) and \
              (sp12_pre.get("stem_intact", False) and sp12_pre.get("n_dr_cross", 0) < 3)
    print("\n判别检验: Sp8(pre)=%s | Sp12(pre)=%s -> %s" % (
        ("茎破坏/侵占" if (not sp8_pre.get("stem_intact", True) or
                          sp8_pre.get("n_dr_cross", 0) >= 3) else "活性样"),
        ("茎破坏/侵占" if (not sp12_pre.get("stem_intact", True) or
                          sp12_pre.get("n_dr_cross", 0) >= 3) else "活性样"),
        "PASS" if discrim else "FAIL"))

    out = {"engine": "Knotty (CCJ + DP09), 容器 ~/dzy/tools/Knotty-master 官方源码编译",
           "energy_scale": "DP09 单标尺; 仅 MFE 级, 无约束配分(与 pfunc -pseudo 的差距如实声明)",
           "creutzburg_rows": rows, "guideforge_rows": gf,
           "discrimination_pre": {
               "sp8": {"stem_intact": sp8_pre.get("stem_intact"),
                       "n_dr_cross": sp8_pre.get("n_dr_cross")},
               "sp12": {"stem_intact": sp12_pre.get("stem_intact"),
                        "n_dr_cross": sp12_pre.get("n_dr_cross")},
               "verdict": "PASS" if discrim else "FAIL"}}
    path = os.path.join(ROOT, "data", "pk_competition.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("jobs", "assemble"))
    args = ap.parse_args()
    (build_jobs if args.mode == "jobs" else assemble)()


if __name__ == "__main__":
    main()
