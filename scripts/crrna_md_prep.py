"""8D4A 结合态 MD 系统搭建(证据链补强 #7, 用户 2026-09-01 指定)。

体系: SuCas12a2(链A) + crRNA(链B, 前18nt=DR) + 靶 RNA(链C) + dsDNA(D/E)
      + ZN/MG 离子; TIP3P 水, 1.0 nm padding, 0.15 M NaCl, amber14 全套
      (蛋白 ff14SB + RNA.OL3 + DNA.OL15 + tip3p)。
突变: crRNA DR 上的点突变(如 A7C+U14G)——残基改名 + 删旧碱基原子,
      PDBFixer.addMissingAtoms 按模板重建碱基重原子, 再加氢。
内部缺失环: 不自动补建(保持实验构型, 缺口记录进 manifest); 只补缺失原子/氢。

用法:
  python scripts/crrna_md_prep.py --mut WT --out-dir data/md/wt
  python scripts/crrna_md_prep.py --mut A7C+U14G --out-dir data/md/A7C_U14G
输出: system.pdb / system.xml / prep.json
"""
import argparse
import json
import os
import sys
import time

from openmm import app, unit, XmlSerializer, LangevinMiddleIntegrator, Platform
from openmm.app import Modeller, ForceField, PDBxFile, PDBFile, PME, HBonds
from pdbfixer import PDBFixer

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
CIF = os.path.join(ROOT, "data", "8D4A.cif")

BASE_ATOMS = {
    "A": {"N9", "C8", "N7", "C5", "C6", "N1", "N6", "C2", "N3"},
    "G": {"N9", "C8", "N7", "C5", "C6", "N1", "O6", "C2", "N2", "N3"},
    "C": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "U": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
}
SUGAR = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "H5'", "H5''", "HO5'"}


def parse_mut(mut):
    out = []
    if mut.upper() != "WT":
        for tok in mut.split("+"):
            wt_b, pos, new_b = tok[0], int(tok[1:-1]), tok[-1]
            assert wt_b in "AUGC" and new_b in "AUGC", tok
            out.append((pos, wt_b, new_b))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mut", default="WT", help="如 A7C+U14G (链B 1-based) 或 WT")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--padding-nm", type=float, default=1.0)
    ap.add_argument("--drop-zn", action="store_true", help="amber 模板缺 ZN 时手动剔除")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    t0 = time.time()
    src = PDBxFile(CIF)
    muts = parse_mut(args.mut)

    mod = Modeller(src.topology, src.positions)

    # 定位 crRNA 链(全 RNA 且最长)
    rna_chains = [c for c in mod.topology.chains()
                  if len(list(c.residues())) > 0
                  and all(r.name in ("A", "U", "G", "C") for r in c.residues())]
    cr_chain = max(rna_chains, key=lambda c: len(list(c.residues())))
    cr_res = list(cr_chain.residues())
    wt_dr = "".join(r.name for r in cr_res[:18])
    assert wt_dr == "AUUUCUACUAUUGUAGAU", "crRNA 链 5' 18nt 非 PDB 口径 DR: %s" % wt_dr

    # 删杂 + 突变残基碱基原子 + 各核酸链首残基 5' 磷酸(amber 终端模板无 5'P)
    del_atoms = [a for a in mod.topology.atoms()
                 if a.residue.name == "HOH"
                 or (args.drop_zn and a.residue.name == "ZN")]
    nuc = set("A U G C DA DT DG DC".split())
    for ch in mod.topology.chains():
        res = list(ch.residues())
        if res and res[0].name in nuc:
            del_atoms += [a for a in res[0].atoms() if a.name in ("P", "OP1", "OP2")]
    mut_log = []
    for pos, wt_b, new_b in muts:
        r = cr_res[pos - 1]
        assert r.name == wt_b, "位点 %d 现为 %s 非 %s" % (pos, r.name, wt_b)
        del_atoms += [a for a in r.atoms() if a.name in BASE_ATOMS[r.name]]
        mut_log.append({"pos": pos, "wt": wt_b, "mut": new_b})
        r.name = new_b
    if del_atoms:
        mod.delete(del_atoms)

    # 本版 PDBFixer 只收文件: 中间 PDB(保留链 ID)
    mut_pdb = os.path.join(args.out_dir, "mutated.pdb")
    with open(mut_pdb, "w") as f:
        PDBFile.writeFile(mod.topology, mod.positions, f, keepIds=True)
    fixer = PDBFixer(filename=mut_pdb)
    fixer.findMissingResidues()
    internal_gaps = []
    for (ci, ri), seq in list(fixer.missingResidues.items()):
        ch = list(fixer.topology.chains())[ci]
        n = ch.getNumResidues()
        if ri not in (0, n):  # 内部缺口
            internal_gaps.append({"chain": ch.id, "after_res_index": ri, "n_missing": len(seq)})
    fixer.missingResidues = {}  # 不补建环
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values()) + \
        sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    fixer.addMissingHydrogens(forcefield=ff)
    mod2 = Modeller(fixer.topology, fixer.positions)
    mod2.addSolvent(ff, model="tip3p", padding=args.padding_nm * unit.nanometer,
                    ionicStrength=0.15 * unit.molar,
                    positiveIon="Na+", negativeIon="Cl-", neutralize=True)

    system = ff.createSystem(mod2.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=HBonds, rigidWater=True)

    pdb_out = os.path.join(args.out_dir, "system.pdb")
    with open(pdb_out, "w") as f:
        PDBFile.writeFile(mod2.topology, mod2.positions, f, keepIds=True)
    xml = XmlSerializer.serialize(system)
    with open(os.path.join(args.out_dir, "system.xml"), "w") as f:
        f.write(xml)

    man = {"mut": args.mut, "mutations": mut_log, "wt_dr_rna": wt_dr,
           "internal_gaps_not_built": internal_gaps,
           "missing_atoms_rebuilt": n_missing,
           "n_atoms_total": mod2.topology.getNumAtoms(),
           "n_residues": mod2.topology.getNumResidues(),
           "padding_nm": args.padding_nm,
           "forcefields": ["amber14-all.xml", "amber14/tip3p.xml"],
           "prep_seconds": round(time.time() - t0, 1),
           "cif_source": "data/8D4A.cif (SuCas12a2-crRNA-靶RNA 复合物, 结合态口径)"}
    with open(os.path.join(args.out_dir, "prep.json"), "w") as f:
        json.dump(man, f, ensure_ascii=False, indent=1)
    print(json.dumps(man, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
