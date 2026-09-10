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

    # production(OpenMM 官方 Simulation API: 报告器由 Simulation 托管初始化
    # 与按间隔触发; 手写 report 循环绕过初始化会在 8.4 抛 _out AttributeError,
    # 2026-09-10 重写)。注意: 生产段只建 integrator, Context 由 Simulation
    # 自建——integrator 一旦绑定 Context 不可复用(8.4 抛 already bound)。
    prod_steps = int(round(args.ns * 1000 / DT_PS))
    st2 = ctx2.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)

    integ3 = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                      DT_PS * unit.picoseconds)
    sim = app.Simulation(pdb.topology, system, integ3,
                         Platform.getPlatformByName("CUDA"),
                         {"DeviceIndex": str(args.device)})
    sim.context.setPositions(st2.getPositions())
    # DCDReporter(append=True) 以 r+b 打开, 文件不存在直接 FileNotFoundError
    # (2026-09-10 根因; 此前的 _out AttributeError 只是该主错误的 __del__ 次生噪音)
    dcd_path = os.path.join(d, "production.dcd")
    sim.reporters.append(
        DCDReporter(dcd_path, 10000, append=os.path.isfile(dcd_path)))  # 20 ps/帧
    sim.reporters.append(
        StateDataReporter(open(os.path.join(d, "md.log"), "a"), 50000, step=True,
                          time=True, speed=True, progress=True,
                          totalSteps=prod_steps, remainingTime=True,
                          potentialEnergy=True, temperature=True,
                          separator=" | "))
    sim.reporters.append(
        CheckpointReporter(os.path.join(d, "production.chk"), 250000))

    chk_file = os.path.join(d, "production.chk")
    start_step = 0
    if os.path.isfile(chk_file):
        sim.loadCheckpoint(chk_file)
        start_step = sim.currentStep
        print("[resume] from step %d / %d" % (start_step, prod_steps), flush=True)
    else:
        sim.context.setVelocities(st2.getVelocities())
        sim.context.setPeriodicBoxVectors(*st2.getPeriodicBoxVectors())

    t0 = time.time()
    sim.step(prod_steps - start_step)
    dt_h = (time.time() - t0) / 3600
    done = {"sys_dir": d, "ns": args.ns, "steps": prod_steps, "device": args.device,
            "wall_hours_fresh": round(dt_h, 2) if start_step == 0 else None,
            "resumed_from": start_step, "completed": True}
    json.dump(done, open(os.path.join(d, "done.json"), "w"), indent=1)
    print(json.dumps(done), flush=True)


if __name__ == "__main__":
    main()
