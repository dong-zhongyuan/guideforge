"""策划案 V3 靶标 panel 计算案例补齐(2026-09-01, 对标 V3 表2; 当日修正 ORF 定位 bug)。

panel: TP53-R248Q(已有, 已核验 ATG+742 密码子248 正确) / KRAS G12C(修复) /
KRAS G12D / TP53 R273H / APC 无义(MCR 区首个 CAG->TAG, 正式口径 Q1328*)。

定位规则(教训: 早期 KRAS 案例把 CDS 坐标当 mRNA 坐标用, 突变落在 5'UTR):
  一律以转录本首个 ATG 为 ORF 原点, 按"密码子号+密码子内编辑"定位,
  引入前断言 WT 密码子与预期氨基酸, 杜绝坐标错位。
  APC 例外与最终口径(2026-09-08 拍板): 任务⑩附带核查发现 APC NM_000038.6
  的 5'UTR 有上游 ATG(13-15 位), "首个 ATG"不是 CDS 起点(真实 CDS=最长
  ORF, 60..8591, 2843 aa)——首个 ATG 锚定的"密码子 1312"(转录本 3946 位
  C>T)在真实阅读框是 A1296V 错义, 并非无义; 且真实密码子 1312 为
  GGA(Gly), "Q1312 无义"在真实阅读框物理不存在。真实阅读框 MCR
  1300-1450 首个 CAG = 密码子 1328(c.3982C>T, 转录本 4041 位)。策划案
  旧口径"Q1312x"曾为湿实验矩阵/演示契约所锚定而短暂锁定, 经用户拍板
  正式整体切换为 Q1328*: apc_design 采用最长 ORF 锚点
  (crrna_mutation_typing.longest_orf), panel 键名同步更名 apc_q1328x,
  下游(湿实验 CSV/演示/Chai 增量/选型特征)全部重生成。TP53/KRAS 的
  首个 ATG 经核验即 CDS 起点, 不受影响。

运行: python scripts/crrna_panel_targets.py [--only kras_g12d,tp53_r273h,apc_q1328x,kras_g12c]
"""
import argparse
import json
import os
import subprocess
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
DATA = os.path.join(ROOT, "data")
AGENT = os.path.join(DATA, "agent")

import crrna_target_abundance as tabund  # noqa: E402
import crrna_mutation_typing as mt  # noqa: E402
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
    """APC 无义: MCR(密码子~1300-1450)内首个 CAG->TAG; 记录实际密码子号。

    锚点 = 最长 ORF(crrna_mutation_typing.longest_orf), 不是 seq.find("ATG"):
    NM_000038.6 5'UTR 有上游 ATG(13-15 位), 首个 ATG 锚定的"密码子 1312"
    (3946 位)在真实阅读框为 A1296V 错义(任务⑩核查发现); 真实锚点下 MCR
    首个 CAG = 密码子 1328(4041 位)——2026-09-08 起为正式口径(Q1328*)。"""
    orf = mt.longest_orf(seq)
    assert orf, "APC 无最长 ORF"
    atg = orf["start0"]
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
    # 丰度档维只对 V3 表2 panel 靶标配发(kras_g12c 为干实验附加, 无 panel 上下文)
    if key in tabund.PANEL_TARGETS:
        cmd += ["--target-key", key]
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

    if (not only or "apc_q1328x" in only):
        key = "apc_q1328x"
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
