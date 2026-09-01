"""OpenMM 体系 → GROMACS 输入转换(证据链补强 #7 配套)。

背景: 宿主机驱动 550 (CUDA 12.4) 无法 JIT openmm 任何现构建的 PTX(8.1-8.6 全试);
GROMACS 2026.3 conda 构建带 sm_86 预编译 SASS, 在容器 ~/dzy/envs/gmx 跑 GPU。
体系参数化(amber14 全套, TIP3P, 0.15M NaCl)已在 OpenMM 侧完成(prep 阶段),
本脚本用 ParmEd 将 system.pdb/system.xml 原样导出为 GROMACS sys.top/sys.gro,
并生成四阶段 mdp(最陡下降 EM → NVT 100ps → NPT 200ps → 产出 10ns)与容器侧
启动脚本 run_gmx.sh。

用法(宿主, gpumd env):
  python scripts/crrna_md_togmx.py --sys-dir data/md/wt --gpu 0
  python scripts/crrna_md_togmx.py --sys-dir data/md/A7C_U14G --gpu 0
输出: sys.top / sys.gro / em.mdp / nvt.mdp / npt.mdp / md.mdp / run_gmx.sh
"""
import argparse
import os

import parmed
from openmm import XmlSerializer, unit
from openmm.app import PDBFile, ForceField, PME

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

MDP = {
"em.mdp": """integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = none
""",
"nvt.mdp": """integrator  = md
dt          = 0.002
nsteps      = 50000
nstxout     = 5000
nstlog      = 5000
nstenergy   = 5000
nstxout-compressed = 5000
continuation    = no
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 1.0
ref_t           = 300
pcoupl          = no
constraints     = h-bonds
DispCorr        = EnerPres
gen_vel         = yes
gen_temp        = 300
gen_seed        = 20260901
""",
"npt.mdp": """integrator  = md
dt          = 0.002
nsteps      = 100000
nstxout     = 5000
nstlog      = 5000
nstenergy   = 5000
nstxout-compressed = 5000
continuation    = yes
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 1.0
ref_t           = 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 5.0
ref_p           = 1.0
compressibility = 4.5e-5
constraints     = h-bonds
DispCorr        = EnerPres
gen_vel         = no
""",
"md.mdp": """integrator  = md
dt          = 0.002
nsteps      = 5000000
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstlog      = 50000
nstenergy   = 50000
nstxout-compressed = 10000
nstcheckpoint    = 25000
continuation    = yes
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 1.0
ref_t           = 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 5.0
ref_p           = 1.0
compressibility = 4.5e-5
constraints     = h-bonds
DispCorr        = EnerPres
gen_vel         = no
""",
}

RUN_SH = """#!/bin/bash
# 容器内执行: bash run_gmx.sh  (GPU {gpu})
set -e
GMX={gmx}
$GMX grompp -f em.mdp -c sys.gro -r sys.gro -p sys.top -o em.tpr -maxwarn 2 > grompp_em.log 2>&1
$GMX mdrun -deffnm em -nb gpu -pme gpu -bonded gpu -gpu_id {gpu} -ntomp 8 > mdrun_em.log 2>&1
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p sys.top -o nvt.tpr -maxwarn 2 > grompp_nvt.log 2>&1
$GMX mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -gpu_id {gpu} -ntomp 8 > mdrun_nvt.log 2>&1
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p sys.top -o npt.tpr -maxwarn 2 > grompp_npt.log 2>&1
$GMX mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -gpu_id {gpu} -ntomp 8 > mdrun_npt.log 2>&1
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p sys.top -o md.tpr -maxwarn 2 > grompp_md.log 2>&1
$GMX mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -gpu_id {gpu} -ntomp 8 -maxh 590 > mdrun_md.log 2>&1
echo DONE > done.flag
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sys-dir", required=True)
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--gmx", default="/home/node/dzy/envs/gmx/bin/gmx")
    args = ap.parse_args()
    d = args.sys_dir

    pdb = PDBFile(os.path.join(d, "system.pdb"))
    # parmed 无法导出 openmm 约束(SETTLE 几何判定失败), 且约束键在 system.xml 中无键参数;
    # 改为: 用同拓扑重建无约束体系(全部键显式参数化), 约束交给 grompp constraints=h-bonds。
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=None, rigidWater=False)
    struct = parmed.openmm.load_topology(pdb.topology, system, xyz=pdb.positions)
    # parmed 对水残基强制走 SETTLE 分支且几何判定失败; 改名后按普通分子写键/约束
    for res in struct.residues:
        if res.name in ("HOH", "WAT", "SOL"):
            res.name = "TP3"
    top = os.path.join(d, "sys.top")
    gro = os.path.join(d, "sys.gro")
    struct.save(top)
    struct.save(gro)
    print("导出 %s (%d 原子) / %s" % (top, len(struct.atoms), gro))

    for name, content in MDP.items():
        with open(os.path.join(d, name), "w") as f:
            f.write(content)
    with open(os.path.join(d, "run_gmx.sh"), "w") as f:
        f.write(RUN_SH.format(gpu=args.gpu, gmx=args.gmx))
    print("mdp + run_gmx.sh 就绪 (GPU %d)" % args.gpu)


if __name__ == "__main__":
    main()
