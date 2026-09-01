"""8D4A 结合态 MD 预注册分析(指标先于轨迹固定, 见 data/md_analysis_plan.json)。

实现约束(用户 2026-09-01 铁律: 有官方工具不用自编译实现):
  预处理   gmx trjconv -pbc whole(GROMACS 官方, 分子完整性);
  M1 茎RMSD MDAnalysis.analysis.rms.RMSD(官方, 蛋白 CA 拟合);
  M2 接触  MDAnalysis.lib.distances.distance_array(官方距离原语);
  M3 氢键  MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis(官方, 距离+角度);
  茎配对/阈值/判据: 读 analysis_input.json(crrna_md_manifest.py 生成), 不硬编码;
  crRNA 首残基自动探测(--dr-first-resid 可覆盖)。

用法(容器, ~/dzy/envs/gmx/bin/python):
  python crrna_md_analysis.py --run-dir ~/md_runs/wt \
      --analysis-input /dawn/md/analysis_input.wt.json \
      --contact-json /dawn/md/contacts.wt.json
输出: <run-dir>/analysis.json
"""
import argparse
import json
import os
import shutil
import subprocess

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rms import RMSD

CONTACT_DIST_A = 4.5
START_NS = 2.0
RNA_DONOR_H = "H61 H62 H41 H42 H21 H22 HO2'"


def trjconv_whole(run_dir, gmx="gmx"):
    """官方 GROMACS: 分子完整性化(幂等, 生成 md_whole.xtc)。"""
    out = os.path.join(run_dir, "md_whole.xtc")
    if not os.path.isfile(out):
        subprocess.run([gmx, "trjconv", "-f", os.path.join(run_dir, "md.xtc"),
                        "-s", os.path.join(run_dir, "md.tpr"),
                        "-o", out, "-pbc", "whole"],
                       input="0\n", text=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out


def find_crna_first_resid(u):
    """残基序 = 原拓扑顺序(蛋白→crRNA 41nt→靶RNA→DNA); crRNA = 首个 RNA 残基。"""
    for r in u.residues:
        if r.resname in ("A", "U", "G", "C"):
            return r.resid
    raise SystemExit("未找到 RNA 残基")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--analysis-input", required=True)
    ap.add_argument("--contact-json", required=True)
    ap.add_argument("--dr-first-resid", type=int, default=None)
    ap.add_argument("--gmx", default="gmx")
    args = ap.parse_args()
    d = args.run_dir

    inp = json.load(open(args.analysis_input, encoding="utf-8"))
    stem_pairs = [tuple(p) for p in inp["dr_stem_pairs_0based"]]
    criteria = inp["criteria"]

    xtc = trjconv_whole(d, args.gmx)
    u = mda.Universe(os.path.join(d, "sys.gro"), xtc)
    ref = mda.Universe(os.path.join(d, "em.gro"))

    first_resid = args.dr_first_resid or find_crna_first_resid(u)
    print("crRNA 首残基 resid=%d, DR %dnt" % (first_resid, inp["dr_len"]))

    stem_sel = " or ".join("(resid %d name C1')" % (first_resid + i)
                           for pair in stem_pairs for i in pair)
    cj = json.load(open(args.contact_json, encoding="utf-8"))
    sites = sorted(int(k) for k, v in cj["contacts"].items()
                   if v.get("n_contact_residues", 0) > 0)
    contact_sel = " or ".join("(resid %d and not name H*)"
                              % (first_resid + s - 1) for s in sites)
    contact_atoms = u.select_atoms(contact_sel)
    prot_heavy = u.select_atoms("protein and not name H*")

    start_frame = int(START_NS * 1000 / u.trajectory.dt)

    # M1: 官方 RMSD, 蛋白 CA 拟合, 茎 C1' 为量测组
    rms = RMSD(u, reference=ref, select="protein and name CA",
               groupselections=[stem_sel])
    rms.run(start=start_frame)
    rmsd_vals = rms.results["rmsd"][:, 2]

    # M2: 逐 100 ps 接触存活
    step = max(1, int(round(100.0 / u.trajectory.dt)))
    surv = []
    for ts in u.trajectory[start_frame::step]:
        dmat = mda.lib.distances.distance_array(contact_atoms.positions,
                                                prot_heavy.positions)
        surv.append(float((dmat.min(axis=1) < CONTACT_DIST_A).mean()))

    # M3: 官方氢键分析(DR 供体氢 vs 蛋白受体 N/O)
    hb = HydrogenBondAnalysis(
        u,
        hydrogens_sel="resid %d-%d and name %s" % (
            first_resid, first_resid + inp["dr_len"] - 1, RNA_DONOR_H),
        acceptors_sel="protein and name O* N*",
        d_h_cutoff=1.2, d_a_cutoff=3.5, angle_cutoff=130.0)
    hb.run(start=start_frame)
    counts = np.zeros(len(u.trajectory))
    for row in hb.results["hbonds"]:
        counts[int(row[0])] += 1
    hb_mean = float(counts[start_frame:].mean())

    out = {"run_dir": d, "analysis_input": args.analysis_input,
           "n_frames_scored": int(len(rmsd_vals)),
           "window_ns": [START_NS, round(float(u.trajectory.time / 1000.0), 2)],
           "M1_stem_c1p_rmsd_A_median": round(float(np.median(rmsd_vals)), 3),
           "M1_stem_c1p_rmsd_A_p90": round(float(np.percentile(rmsd_vals, 90)), 3),
           "M2_contact_survival_frac": round(float(np.mean(surv)), 3),
           "M3_dr_donor_hbond_mean": round(hb_mean, 1),
           "contact_sites_1based": sites,
           "criteria_pass_thresholds": {
               "M1_candidate_leq_WT_plus_A": criteria.get("M1_pass_delta_A"),
               "M2_candidate_geq_WT_minus_pp": criteria.get("M2_pass_delta_pp")},
           "preregistered": criteria.get("preregistered", ""),
           "note": "M1/M2 主判据, M3 次要; 阴性结果同样报告; 全部指标官方库实现"}
    path = os.path.join(d, "analysis.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False)[:500])
    print("输出 -> %s" % path)


if __name__ == "__main__":
    main()
