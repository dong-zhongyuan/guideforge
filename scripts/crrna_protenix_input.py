"""Protenix 输入生成器(2026-09-05, AF3 权重不可得的开源替代通道)。

背景: 策划案 V3 表1 的 "AF3/Chai-1" 共折叠层需要独立于 Chai-1 的
template-free 对照。DeepMind AF3 权重 gated 申请未果; Protenix-v1
(ByteDance, Apache-2.0 代码+权重全开放, 基准在同等约束下达到/超过 AF3,
bioRxiv 2026.02.05.703733) 是合法的 AF3 级替代引擎。

本脚本把 data/af3_inputs/GF_*.json(alphafoldserver 方言: 顶层任务列表,
蛋白链 useStructureTemplate=false) 转换为 Protenix 方言(顶层任务列表,
链 {proteinChain|rnaSequence: {sequence, count}}; 不含 templates 字段
= template-free, 与 AF3 输入口径一致)。

运行: python scripts/crrna_protenix_input.py
输出: data/protenix_inputs/PX_*.json + manifest.json
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "data", "af3_inputs")
DST = os.path.join(ROOT, "data", "protenix_inputs")


def convert_job(job):
    seqs = []
    for entry in job["sequences"]:
        for kind, v in entry.items():
            out = {"sequence": v["sequence"], "count": 1}
            seqs.append({kind: out})
    return {"name": job["name"], "covalent_bonds": [], "sequences": seqs}


def main():
    os.makedirs(DST, exist_ok=True)
    manifest = []
    for f in sorted(glob.glob(os.path.join(SRC, "GF_*.json"))):
        if f.endswith("GF_all13_alphafoldserver.json"):
            continue
        jobs = json.load(open(f, encoding="utf-8"))
        assert isinstance(jobs, list) and len(jobs) == 1, f
        job = convert_job(jobs[0])
        n_res = sum(len(list(s.values())[0]["sequence"]) for s in job["sequences"])
        out = os.path.join(DST, "PX_" + os.path.basename(f)[3:])
        json.dump([job], open(out, "w", encoding="utf-8"), indent=1)
        manifest.append({"file": os.path.basename(out), "name": job["name"],
                         "n_chains": len(job["sequences"]), "n_residues": n_res})
        print(f"  {os.path.basename(out)}: {len(job['sequences'])} 链 {n_res} 残基")
    json.dump({"generated_by": "scripts/crrna_protenix_input.py (2026-09-05)",
               "purpose": "AF3 级开源引擎(Protenix-v1) template-free 共折叠输入; "
                          "替代 gated AF3, 作 Chai-1 自模板口径的独立对照层",
               "model": "protenix_base_default_v1.0.0 (368M, 训练截止 2021-09-30 对齐 AF3)",
               "license": "Apache-2.0 代码+权重",
               "template_free": "输入不含 templates 字段即 template-free",
               "jobs": manifest},
              open(os.path.join(DST, "manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"[protenix] {len(manifest)} 套输入 -> {DST}")


if __name__ == "__main__":
    main()
