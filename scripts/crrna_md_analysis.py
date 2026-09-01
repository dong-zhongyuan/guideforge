"""8D4A 结合态 MD 预注册分析(证据链补强 #7 收尾; 计划先于轨迹固定, 2026-09-01)。

预注册指标(主判据两项, 阴性结果同样报告):
  M1 茎区 RMSD: DR 茎 5 对碱基(0基 (3,16)(4,15)(5,14)(6,13)(7,12), 即 1 基 4-8 与 13-17)
     共 10 个 C1' 原子, 以蛋白链 A Cα 拟合后, 对最小化后构型(em.gro)的 RMSD;
     统计 2-10 ns 窗口中位数。判据: 候选 <= WT + 0.5 A 为"茎几何不劣于 WT"。
  M2 DR-蛋白接触保持: 8D4A 口径 5 A 接触表(data/8D4A_dr_contacts_*.json)所列
     DR 位点对应残基的全部重原子, 与蛋白链 A 全部重原子的最小距离逐帧(每 100 ps)
     判存活(<4.5 A)。判据: 候选存活率 >= WT - 10 个百分点为"接触不劣于 WT"。
  M3(次要) DR-蛋白氢键数: DR 区极性原子(N/O)与蛋白极性原子距离 <3.5 A 计一条,
     逐帧计数(近似距离口径, 不含角度)。参考性输出, 不作判据。
  解释性(不判据): 蛋白 Cα RMSD、体系回旋半径、温度/压强(md.log 摘录)。

输入: 容器 ~/md_runs/<tag>/{sys.gro, md.xtc, em.gro}; 宿主接触表经 /dawn 传入。
用法(容器, ~/dzy/envs/gmx/bin/python): python crrna_md_analysis.py --run-dir ~/md_runs/wt \
  --contact-json /dawn/md/8D4A_dr_contacts_zeng.json --ref-dir ...
输出: <run-dir>/analysis.json + 控制台 WT vs 候选对照表(由汇总入口聚合)。
"""
import argparse
import json
import os
import sys

import MDAnalysis as mda
import numpy as np

STEM_PAIRS_0BASE = [(3, 16), (4, 15), (5, 14), (6, 13), (7, 12)]  # DR 18nt, 5bp 茎
CONTACT_DIST_A = 4.5
HBOND_DIST_A = 3.5
RMSD_TOL_A = 0.5
CONTACT_TOL_PP = 10.0
START_NS = 2.0


def heavy(ag):
    return ag.select_atoms("not name H*")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--contact-json", required=True,
                    help="8D4A DR 接触表 json(取其 DR 位点列表)")
    ap.add_argument("--dr-first-resid", type=int, required=True,
                    help="crRNA 链(41nt)在体系拓扑中的第一个残基号(sys.gro 重编号后)")
    ap.add_argument("--protein-name", default="protein")
    args = ap.parse_args()
    d = args.run_dir

    u = mda.Universe(os.path.join(d, "sys.gro"), os.path.join(d, "md.xtc"))
    ref = mda.Universe(os.path.join(d, "em.gro"))

    # crRNA 链 = 41 连续 RNA 残基(名称 A/U/G/C)
    rna = u.select_atoms("resname A U G C")
    resid_offset = args.dr_first_resid
    dr = u.select_atoms("resid %d-%d and resname A U G C"
                        % (resid_offset, resid_offset + 17))
    prot = u.select_atoms("protein")

    stem_sel = " or ".join(
        "resid %d name C1' or resid %d name C1'" % (resid_offset + i, resid_offset + j)
        for i, j in STEM_PAIRS_0BASE)
    stem = u.select_atoms(stem_sel)
    stem_ref = ref.select_atoms(stem_sel.replace("resid", "resid"))  # 同序
    prot_ca = u.select_atoms("protein and name CA")
    prot_ca_ref = ref.select_atoms("protein and name CA")

    # M2 接触位点: 接触表中 counts>0 的 DR 位点(1 基)
    cj = json.load(open(args.contact_json, encoding="utf-8"))
    sites = sorted(int(k) for k, v in cj["contacts"].items() if v.get("n_contact_residues", 0) > 0)
    contact_sel = " or ".join("resid %d and not name H*"
                              % (resid_offset + s - 1) for s in sites)
    contact_atoms = u.select_atoms(contact_sel)
    prot_heavy = heavy(prot)
    dr_polar = dr.select_atoms("name N* O* and not name H*")

    dt = u.trajectory.dt
    rmsds, surv, hbs, frames = [], [], [], 0
    for ts in u.trajectory:
        t_ns = ts.time / 1000.0
        frames += 1
        if t_ns < START_NS:
            continue
        # M1: Kabsch 拟合蛋白 Cα 后量茎 C1' RMSD
        P, Q = prot_ca.positions - prot_ca.centroid(), prot_ca_ref.positions - prot_ca_ref.centroid()
        cov = P.T @ Q
        v, s, wt_ = np.linalg.svd(cov)
        if np.linalg.det(v @ wt_) < 0:
            wt_[-1] *= -1
        rot = v @ wt_
        sp = (stem.positions - prot_ca.centroid()) @ rot + prot_ca_ref.centroid()
        rmsds[-1] = float(np.sqrt(np.mean(np.sum((sp - stem_ref.positions) ** 2, axis=1))))
        # M2: 接触存活
        dmat = mda.lib.distances.distance_array(contact_atoms.positions,
                                                prot_heavy.positions)
        surv.append(float((dmat.min(axis=1) < CONTACT_DIST_A).mean()))
        # M3: 氢键(距离口径)
        dmat2 = mda.lib.distances.distance_array(dr_polar.positions,
                                                 prot_heavy.select_atoms("name N* O*").positions)
        hbs.append(int((dmat2 < HBOND_DIST_A).sum()))

    out = {"run_dir": d, "n_frames_total": frames,
           "n_frames_scored": len(rmsds),
           "window_ns": [START_NS, float(u.trajectory.time / 1000)],
           "M1_stem_c1p_rmsd_A_median": round(float(np.median(rmsds)), 3),
           "M1_stem_c1p_rmsd_A_p90": round(float(np.percentile(rmsds, 90)), 3),
           "M2_contact_survival_frac": round(float(np.mean(surv)), 3),
           "M3_dr_protein_polar_contacts": round(float(np.mean(hbs)), 1),
           "criteria": {"M1_pass_delta_A": RMSD_TOL_A, "M2_pass_delta_pp": CONTACT_TOL_PP,
                        "contact_sites_1based": sites},
           "preregistered": "2026-09-01 先于轨迹固定; 阴性结果同样报告"}
    path = os.path.join(d, "analysis.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False)[:400])
    print("输出 -> %s" % path)


if __name__ == "__main__":
    main()
