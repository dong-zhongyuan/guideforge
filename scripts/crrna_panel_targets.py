"""策划案 V3 靶标 panel 计算案例补齐(2026-09-01, 对标 V3 表2; 当日修正 ORF 定位 bug)。

panel: TP53-R248Q(已有, 已核验 ATG+742 密码子248 正确) / KRAS G12C(修复) /
KRAS G12D / TP53 R273H / APC 无义(MCR 区首个 CAG->TAG)。

定位规则(教训: 早期 KRAS 案例把 CDS 坐标当 mRNA 坐标用, 突变落在 5'UTR):
  一律以转录本首个 ATG 为 ORF 原点, 按"密码子号+密码子内编辑"定位,
  引入前断言 WT 密码子与预期氨基酸, 杜绝坐标错位。

运行: python scripts/crrna_panel_targets.py [--only kras_g12d,tp53_r273h,apc_q1338x,kras_g12c]
"""
import argparse
import json
import os
import subprocess
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")
AGENT = os.path.join(DATA, "agent")
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id=%s&rettype=fasta&retmode=text")

AA = {}
_b = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for x in _b:
    for y in _b:
        for z in _b:
            AA[x + y + z] = _aas[_i]
            _i += 1

# (key, accession, 密码子号, WT密码子, 突变密码子, 预期WT氨基酸, 标签)
TARGETS = [
    ("kras_g12c", "NM_033360", 12, "GGT", "TGT", "G", "KRAS_G12C_c.34G>T"),
    ("kras_g12d", "NM_033360", 12, "GGT", "GAT", "G", "KRAS_G12D_c.35G>A"),
    ("tp53_r273h", "NM_000546", 273, "CGT", "CAT", "R", "TP53_R273H_c.818G>A"),
]


def fetch(acc, cache):
    if os.path.isfile(cache):
        return open(cache).read()
    with urllib.request.urlopen(EFETCH % acc, timeout=120) as r:
        txt = r.read().decode()
    open(cache, "w").write(txt)
    return txt


def read_seq(txt):
    return "".join(l.strip() for l in txt.splitlines()
                   if not l.startswith(">")).upper()


def mutate(seq, codon_no, wt_codon, mut_codon, label, key):
    atg = seq.find("ATG")
    assert atg >= 0, "%s 无 ATG" % key
    start = atg + (codon_no - 1) * 3
    obs = seq[start:start + 3]
    assert obs == wt_codon, "%s 密码子 %d 为 %s(%s) 非预期 %s" % (
        key, codon_no, obs, AA.get(obs, "?"), wt_codon)
    edited = [i for i in range(3) if wt_codon[i] != mut_codon[i]]
    assert len(edited) == 1, "%s 编辑须单碱基" % key
    pos = start + edited[0] + 1  # mRNA 1 基
    out = seq[:start] + mut_codon + seq[start + 3:]
    return out, pos, wt_codon[edited[0]], mut_codon[edited[0]]


def apc_design(seq):
    """APC 无义: MCR(密码子~1300-1450)内首个 CAG->TAG; 记录实际密码子号。"""
    atg = seq.find("ATG")
    for codon_no in range(1300, 1451):
        c = seq[atg + (codon_no - 1) * 3: atg + codon_no * 3]
        if c == "CAG":
            out = seq[:atg + (codon_no - 1) * 3] + "TAG" + seq[atg + codon_no * 3:]
            return out, atg + (codon_no - 1) * 3 + 1, codon_no, "C", "T"
    raise SystemExit("APC MCR 区未找到 CAG")


def run_agent(key, mut_fa, wt_fa):
    out = os.path.join(AGENT, key)
    cmd = [sys.executable, os.path.join(HERE, "crrna_design_agent.py"),
           "--mut-fasta", mut_fa, "--wt-fasta", wt_fa,
           "--effector", "cas12a2_zeng2026",
           "--clusters", os.path.join(DATA, "context_typing", "typing.clusters.json"),
           "--spacers", os.path.join(DATA, "context_typing", "spacers.txt"),
           "--out-prefix", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    print("\n".join(r.stdout.strip().splitlines()[-4:]))
    if r.returncode != 0:
        print("[错误] %s: %s" % (key, r.stderr[-300:]))
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    os.makedirs(AGENT, exist_ok=True)
    only = set(args.only.split(",")) if args.only else None

    jobs = []
    for key, acc, codon, wtc, mutc, aa, label in TARGETS:
        if only and key not in only:
            continue
        wt_fa = os.path.join(AGENT, "%s_wt.fa" % key)
        seq = read_seq(fetch(acc, wt_fa))
        mseq, pos, wb, mb = mutate(seq, codon, wtc, mutc, label, key)
        mut_fa = os.path.join(AGENT, "%s.fa" % key)
        open(mut_fa, "w").write(">%s\n%s\n" % (label, mseq))
        if not os.path.isfile(wt_fa):
            open(wt_fa, "w").write(">%s_WT\n%s\n" % (key, seq))
        print("[%s] 密码子 %d %s->%s (mRNA 第 %d 位 %s>%s) ORF 定位核验通过" % (
            key, codon, wtc, mutc, pos, wb, mb))
        jobs.append(key)

    if (not only or "apc_q1338x" in only):
        key = "apc_q1338x"
        wt_fa = os.path.join(AGENT, "%s_wt.fa" % key)
        seq = read_seq(fetch("NM_000038", wt_fa))
        mseq, pos, codon, wb, mb = apc_design(seq)
        label = "APC_NM000038_Q%dx_c.%dC>T(MCR区设计无义,非特定COSMIC条目)" % (
            codon, (codon - 1) * 3 + 1)
        mut_fa = os.path.join(AGENT, "%s.fa" % key)
        open(mut_fa, "w").write(">%s\n%s\n" % (label, mseq))
        if not os.path.isfile(wt_fa):
            open(wt_fa, "w").write(">%s_WT\n%s\n" % (key, seq))
        print("[%s] MCR 区密码子 %d CAG->TAG (mRNA 第 %d 位 C>T)" % (key, codon, pos))
        jobs.append(key)

    for key in jobs:
        wt_fa = os.path.join(AGENT, "%s_wt.fa" % key)
        mut_fa = os.path.join(AGENT, "%s.fa" % key)
        run_agent(key, mut_fa, wt_fa)


if __name__ == "__main__":
    main()
