"""Chai-1 增量折叠输入生成器(2026-09-05, 矩阵补齐: APC 列 + 新补偿臂)。

背景: data/chai_matrix_4t.json 是旧四靶(含 KRAS_G12C)+旧补偿臂(round-3 前
退役分子)时代的产物。V3 对齐后矩阵应为 R248Q/G12D/R273H/APC × 8 骨架。
增量策略: 旧矩阵 3 靶 × 7 非补偿骨架 = 21 行协议一致可复用; 需新跑
  APC × 8 骨架(含 WT) = 8 折
  B_break_compensate(新臂 zengDR+T10A/T12G/T19G) × 旧 3 靶 = 3 折
共 11 折, 与旧 21 行拼成新 4×8=32 矩阵。协议与旧矩阵完全一致
(8D4A 自模板 m8 + ESM), 保证行间可比; template-free 独立对照由
Protenix/AF3 层承担, 不在本增量内。

输入构造(与 crrna_af3_input.py 同源口径):
  蛋白 = SuCas12a2 1207aa(data/chai_input_WT.fasta)
  crRNA = DR 变体 19nt + 靶标 spacer 24nt (RNA)
  靶 RNA = 突变转录本上 protospacer 窗口 24nt + 3' 下游 PFS 5nt(原位提取,
           不经手抄; R248Q 重构序列与 data/af3_inputs/GF_WT.json 的 29nt
           靶链逐字符一致, 作构造正确性的对拍证据)

运行: python scripts/crrna_chai_increment_input.py
输出: data/chai_increment_inputs/CHAI_<靶标>_<骨架>.fasta × 11 + manifest.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "chai_increment_inputs")

COMP = str.maketrans("ACGTU", "TGCAA")

TARGETS = {  # 靶标 -> (突变转录本 fasta, spacer DNA)
    "TP53_R248Q": ("tp53_r248q_mrna_NM000546.fa", "GTTCATGCCGCCCATGCAGGAACT"),
    "KRAS_G12D": ("agent/kras_g12d.fa", "CAGCTCCAACTACCACAAGTTTAT"),
    "TP53_R273H": ("agent/tp53_r273h.fa", "CACCTCAAAGCTGTTCCGTCCCAG"),
    "APC_Q1312x": ("agent/apc_q1338x.fa", "CTTCCTGTGTCGTCTGATTACATC"),
}

SCAFFOLDS = {  # 8 员族 DR (DNA)
    "WT": "AATTTCTACTGTTGTAGAT",
    "A8C_U15G": "AATTTCTCCTGTTGGAGAT",
    "A1G_U3A_A8G_U15C": "GAATTCTGCTGTTGCAGAT",
    "A1C": "CATTTCTACTGTTGTAGAT",
    "A1U": "TATTTCTACTGTTGTAGAT",
    "U5G_A18C": "AATTGCTACTGTTGTAGCT",
    "A1U_A2C": "TCTTTCTACTGTTGTAGAT",
    "B_break_compensate": "AATTTCTACAGGTGTAGAG",
}


def to_rna(s):
    return s.upper().replace("T", "U")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper().replace("U", "T")


def target29(mut_fasta, spacer):
    """突变转录本上 protospacer 窗口 + 3' PFS 5nt (sense, RNA)。"""
    tx = read_fasta(mut_fasta)
    win = revcomp(spacer)
    i = tx.index(win)
    return to_rna(tx[i:i + len(win) + 5])


def main():
    os.makedirs(OUT, exist_ok=True)
    prot = None
    with open(os.path.join(DATA, "chai_input_WT.fasta")) as fh:
        lines = fh.read().splitlines()
    prot = lines[1]  # >protein|name=SuCas12a2 下一行

    # 对拍: R248Q 重构靶链须等于 AF3 输入的 29nt 靶链
    gf = json.load(open(os.path.join(DATA, "af3_inputs", "GF_WT.json"),
                        encoding="utf-8"))[0]
    af3_tgt = [list(s.values())[0]["sequence"] for s in gf["sequences"]
               if "rnaSequence" in s][1]
    r248q_tgt = target29(os.path.join(DATA, TARGETS["TP53_R248Q"][0]),
                         TARGETS["TP53_R248Q"][1])
    assert r248q_tgt == af3_tgt, (
        f"R248Q 靶链重构与 AF3 输入不符: {r248q_tgt} != {af3_tgt}")
    print("[对拍] R248Q 靶链 29nt 与 GF_WT.json 逐字符一致 ✓")

    jobs = []
    for tgt, (fa, spacer) in TARGETS.items():
        t29 = target29(os.path.join(DATA, fa), spacer)
        scafs = SCAFFOLDS if tgt == "APC_Q1312x" else {
            "B_break_compensate": SCAFFOLDS["B_break_compensate"]}
        for sc, dr in scafs.items():
            crna = to_rna(dr + spacer)
            name = f"CHAI_{tgt}_{sc}"
            with open(os.path.join(OUT, name + ".fasta"), "w") as out:
                out.write(">protein|name=SuCas12a2\n%s\n" % prot)
                out.write(">rna|name=crRNA\n%s\n" % crna)
                out.write(">rna|name=target\n%s\n" % t29)
            jobs.append({"file": name + ".fasta", "target": tgt,
                         "scaffold": sc, "crRNA": crna, "target_rna": t29})
    json.dump({"generated_by": "scripts/crrna_chai_increment_input.py (2026-09-05)",
               "protocol": "与旧矩阵一致: 8D4A 自模板 m8(data/chai_template_hits.m8) + ESM; "
                           "template-free 对照由 Protenix/AF3 层承担",
               "increment": "APC×8 + 新补偿臂×3 = 11 折, 拼旧矩阵 21 行成新 4×8=32",
               "crosscheck": "R248Q 靶链重构与 GF_WT.json 逐字符一致",
               "jobs": jobs},
              open(os.path.join(OUT, "manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"[chai增量] {len(jobs)} 套输入 -> {OUT}")


if __name__ == "__main__":
    main()
