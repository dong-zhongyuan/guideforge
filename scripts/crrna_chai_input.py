"""Chai-1 共折叠输入 FASTA 生成器(8D4A -> 三元复合物输入, 2026-09-02)。

从 data/8D4A.cif 提取 SuCas12a2 蛋白链(1207aa)与两条 RNA 链(crRNA 41nt /
靶 RNA 28nt), 写成 Chai 头部约定(>protein|name=... / >rna|name=...)的 FASTA。
此前 data/chai_input_WT.fasta 由内联代码生成(纪律违规), 本脚本为其正式出处,
幂等可再生。
运行(gpumd env): python scripts/crrna_chai_input.py [--out data/chai_input_WT.fasta]
"""
import argparse
import os

from openmm.app import PDBxFile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cif", default=os.path.join(ROOT, "data", "8D4A.cif"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "chai_input_WT.fasta"))
    args = ap.parse_args()

    f = PDBxFile(args.cif)
    chains = list(f.topology.chains())
    prot = max((c for c in chains
                if all(r.name in THREE2ONE for r in c.residues())),
               key=lambda c: len(list(c.residues())))
    rna = [c for c in chains
           if len(list(c.residues())) > 0
           and all(r.name in ("A", "U", "G", "C") for r in c.residues())]
    crna = max(rna, key=lambda c: len(list(c.residues())))
    target = next(c for c in rna if c is not crna)

    seq_prot = "".join(THREE2ONE[r.name] for r in prot.residues())
    seq_crna = "".join(r.name for r in crna.residues())
    seq_tgt = "".join(r.name for r in target.residues())
    with open(args.out, "w") as out:
        out.write(">protein|name=SuCas12a2\n%s\n" % seq_prot)
        out.write(">rna|name=crRNA\n%s\n" % seq_crna)
        out.write(">rna|name=target\n%s\n" % seq_tgt)
    print("蛋白 %daa + crRNA %dnt + 靶 %dnt -> %s" % (
        len(seq_prot), len(seq_crna), len(seq_tgt), args.out))


if __name__ == "__main__":
    main()
