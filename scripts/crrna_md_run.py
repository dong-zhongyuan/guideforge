"""8D4A 结合态 MD 运行(证据链补强 #7): 最小化 → NVT → NPT → 产出 MD。

用法:
  CUDA_VISIBLE_DEVICES=0 python scripts/crrna_md_run.py --sys-dir data/md/wt --ns 10
断点续跑: 同目录 production.chk 存在则从检查点继续。
输出: min.pdb(最小化后), equil.dcd, production.dcd, md.log, production.chk, done.json
"""
import argparse
import json
import os
import time

from openmm import (app, unit, XmlSerializer, LangevinMiddleIntegrator,
                    MonteCarloBarostat, Platform, Context, LocalEnergyMinimizer)
from openmm.app import PDBFile, DCDReporter, StateDataReporter, CheckpointReporter

DT_PS = 0.002  # 2 fs


def context_for(system, device):
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                          DT_PS * unit.picoseconds)
    platform = Platform.getPlatformByName("CUDA")
    return Context(system, integrator, platform, {"DeviceIndex": str(device)}), integrator


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sys-dir", required=True)
    ap.add_argument("--ns", type=float, default=10.0)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--nvt-ps", type=float, default=100.0)
    ap.add_argument("--npt-ps", type=float, default=200.0)
    args = ap.parse_args()
    d = args.sys_dir

    pdb = PDBFile(os.path.join(d, "system.pdb"))
    system = XmlSerializer.deserialize(open(os.path.join(d, "system.xml")).read())

    # NVT
    ctx, integ = context_for(system, args.device)
    ctx.setPositions(pdb.positions)
    LocalEnergyMinimizer.minimize(ctx, tolerance=10 * unit.kilojoule_per_mole / unit.nanometer)
    state = ctx.getState(getPositions=True)
    with open(os.path.join(d, "min.pdb"), "w") as f:
        PDBFile.writeFile(pdb.topology, state.getPositions(), f, keepIds=True)
    ctx.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    print("[nvt] %g ps" % args.nvt_ps, flush=True)
    integ.step(int(round(args.nvt_ps / DT_PS)))

    # NPT
    system.addForce(MonteCarloBarostat(1 * unit.atmosphere, 300 * unit.kelvin, 25))
    ctx2, integ2 = context_for(system, args.device)
    state = ctx.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    ctx2.setPositions(state.getPositions())
    ctx2.setVelocities(state.getVelocities())
    ctx2.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    print("[npt] %g ps" % args.npt_ps, flush=True)
    integ2.step(int(round(args.npt_ps / DT_PS)))

    # production
    prod_steps = int(round(args.ns * 1000 / DT_PS))
    ctx3, integ3 = context_for(system, args.device)
    st2 = ctx2.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    ctx3.setPositions(st2.getPositions())
    ctx3.setVelocities(st2.getVelocities())
    ctx3.setPeriodicBoxVectors(*st2.getPeriodicBoxVectors())

    dcd = DCDReporter(os.path.join(d, "production.dcd"), 10000, append=True)  # 20 ps/帧
    rep = StateDataReporter(open(os.path.join(d, "md.log"), "a"), 50000, step=True,
                            time=True, speed=True, progress=True,
                            totalSteps=prod_steps, remainingTime=True,
                            potentialEnergy=True, temperature=True, separator=" | ")
    chk = CheckpointReporter(os.path.join(d, "production.chk"), 250000)

    chk_file = os.path.join(d, "production.chk")
    start_step = 0
    if os.path.isfile(chk_file):
        with open(chk_file, "rb") as f:
            ctx3.loadCheckpoint(f.read())
        start_step = integ3.getStep()
        print("[resume] from step %d / %d" % (start_step, prod_steps), flush=True)

    t0 = time.time()
    chunk = 50000
    while integ3.getStep() < prod_steps:
        integ3.step(min(chunk, prod_steps - integ3.getStep()))
        dcd.report(ctx3, integ3)
        rep.report(ctx3, integ3)
        chk.report(ctx3, integ3)
    dt_h = (time.time() - t0) / 3600
    done = {"sys_dir": d, "ns": args.ns, "steps": prod_steps, "device": args.device,
            "wall_hours_fresh": round(dt_h, 2) if start_step == 0 else None,
            "resumed_from": start_step, "completed": True}
    json.dump(done, open(os.path.join(d, "done.json"), "w"), indent=1)
    print(json.dumps(done), flush=True)


if __name__ == "__main__":
    main()
