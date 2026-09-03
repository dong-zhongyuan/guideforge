"""全转录组脱靶扫描（Cas12a2, RNA 靶向; 精确算法非种子近似）

与 crrna_specificity_scan.py 的区别: 那个脚本扫少数几条转录本(逐窗口纯 Python),
本脚本面向 GENCODE 全转录组(~20 万条), 用 numpy 向量化在拼接序列上对所有
24nt 窗口精确计算错配数(无种子近似、无假阴性), 并按 PFS 相似度分级。

判定口径(与 crrna_specificity_scan.py / PDB 8D4A 一致):
  - 扫描目标 = revcomp(spacer), 窗口位于转录本 sense 链
  - PFS 位于窗口 3' 下游; Cas12a2 PFS 松散(Zeng 2026: R248Q 靶点 PFS=CAGAG 为
    GAAAG 的 1 错配), 故按 PFS 错配数分级: 0mm=高危, 1mm=中危, >=2mm=低危
  - 种子区 = spacer 3' 端 7nt(Bravo 2023) ↔ 窗口 5' 端 7nt(revcomp 方向);
    种子区错配权重高于远端(单独报告 n_seed_mm)

用法:
  python crrna_transcriptome_scan.py --spacer GTTCATGCCGCCCATGCAGGAACT \
      --fasta ../data/transcriptome/gencode.v49.pc_transcripts.fa \
      --fasta2 ../data/transcriptome/gencode.v49.lncRNA_transcripts.fa \
      --out ../data/tp53_r248q_transcriptome_scan
输出: <out>.sites.csv / <out>.summary.json / <out>.png
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)

from crrna_specificity_scan import revcomp, normalize  # noqa: E402
from scaffold_registry import get_pam  # noqa: E402

B2I = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
SEED_LEN = 7  # 窗口 5' 端 7nt ↔ spacer 3' 种子区(Bravo 2023)


def read_fasta_annot(path):
    """读 FASTA, 返回 [(tid, gene, seq)]; GENCODE 头为竖线分隔, 第 6 段为基因名。"""
    records, name, gene, chunks = [], None, '', []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    records.append((name, gene, ''.join(chunks)))
                parts = line[1:].strip().split()[0].split('|')
                name = parts[0]
                gene = parts[5] if len(parts) > 5 else ''
                chunks = []
            elif line.strip():
                chunks.append(line.strip().upper().replace('U', 'T'))
    if name is not None:
        records.append((name, gene, ''.join(chunks)))
    return records


def encode_concatenated(records):
    """拼接全部转录本为 int8 数组(N=4 作分隔), 记录每条转录本偏移。"""
    parts, offsets, pos = [], [], 0
    for tid, gene, seq in records:
        offsets.append((tid, gene, pos, len(seq)))
        parts.append(seq)
        pos += len(seq) + 1
    blob = 'N'.join(parts)
    arr = np.frombuffer(blob.encode('ascii'), dtype=np.uint8)
    lut = np.full(256, 4, dtype=np.int8)
    for b, i in B2I.items():
        lut[ord(b)] = i
    return lut[arr], offsets


def scan(records, target, pfs, max_mm):
    """向量化精确扫描。返回命中列表(dict)。"""
    S, offsets = encode_concatenated(records)
    n = len(S)
    L, lp = len(target), len(pfs)
    tcode = np.array([B2I[c] for c in target], dtype=np.int8)
    pcode = np.array([B2I[c] for c in pfs], dtype=np.int8)

    valid = n - L - lp + 1
    mm = np.zeros(valid, dtype=np.int8)
    for j in range(L):
        mm += (S[j:j + valid] != tcode[j])
    pmm = np.zeros(valid, dtype=np.int8)
    for j in range(lp):
        pmm += (S[L + j:L + j + valid] != pcode[j])
    # 含分隔符 N 的窗口/PFS 无效: N 编码为 4, 用累计和标记
    is_n = (S == 4).astype(np.int32)
    csum = np.concatenate(([0], np.cumsum(is_n)))
    n_in_win = csum[L:valid + L] - csum[:valid]
    n_in_pfs = csum[L + lp:valid + L + lp] - csum[L:valid + L - 0] if lp else 0
    ok = (mm <= max_mm) & (n_in_win == 0) & (np.asarray(n_in_pfs) == 0)
    hits = np.where(ok)[0]

    # 定位到转录本
    starts = np.array([o[2] for o in offsets])
    rec_idx = np.searchsorted(starts, hits, side='right') - 1
    sites = []
    tcode_set_seed = set(range(SEED_LEN))
    for h, ri in zip(hits, rec_idx):
        tid, gene, off, ln = offsets[ri]
        w0 = S[h:h + L]
        mmpos = [j for j in range(L) if w0[j] != tcode[j]]
        sites.append({
            'transcript': tid, 'gene': gene, 'pos': int(h - off),
            'mismatches': int(mm[h]),
            'seed_mm': sum(1 for j in mmpos if j in tcode_set_seed),
            'distal_mm': sum(1 for j in mmpos if j not in tcode_set_seed),
            'window': ''.join('ACGT'[x] for x in w0),
            'pfs': ''.join('ACGT'[x] for x in S[h + L:h + L + lp]),
            'pfs_mm': int(pmm[h]),
            'mm_positions': ';'.join(str(j + 1) for j in mmpos),
        })
    return sites


def risk_tier(s):
    """风险分级: 0-1 错配+种子区完整+PFS 0-1mm = 高危。"""
    if s['mismatches'] <= 1 and s['seed_mm'] == 0 and s['pfs_mm'] <= 1:
        return 'high'
    if s['mismatches'] <= 2 and s['seed_mm'] <= 1 and s['pfs_mm'] <= 1:
        return 'medium'
    return 'low'


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacer', required=True)
    ap.add_argument('--fasta', required=True, help='主 FASTA(如 pc_transcripts)')
    ap.add_argument('--fasta2', default=None, help='追加 FASTA(如 lncRNA)')
    ap.add_argument('--pfs', default=None,
                    help='PFS 规则(默认取注册表 cas12a2 条目 pam.rule; 换体系用 --pfs 覆盖)')
    ap.add_argument('--max-mismatch', type=int, default=4)
    ap.add_argument('--out', default='transcriptome_scan')
    args = ap.parse_args()

    spacer = normalize(args.spacer)
    if not 17 <= len(spacer) <= 25 or set(spacer) - set('ACGT'):
        ap.error('--spacer 必须为 17-25nt ACGT/U')
    target = revcomp(spacer)
    pfs = normalize(args.pfs or get_pam('cas12a2')['rule'])

    records = read_fasta_annot(args.fasta)
    n1 = len(records)
    if args.fasta2:
        records += read_fasta_annot(args.fasta2)
    print(f"转录本: {len(records)} 条(主编码 {n1} + lncRNA {len(records) - n1})")
    print(f"扫描目标(revcomp): {target}  PFS: {pfs}(按错配数分级)")

    sites = scan(records, target, pfs, args.max_mismatch)
    for s in sites:
        s['risk'] = risk_tier(s)
    sites.sort(key=lambda s: (s['mismatches'], s['pfs_mm'], s['seed_mm'], s['gene']))

    by_mm, by_risk = {}, {}
    for s in sites:
        by_mm[s['mismatches']] = by_mm.get(s['mismatches'], 0) + 1
        by_risk[s['risk']] = by_risk.get(s['risk'], 0) + 1
    print(f"命中位点(<= {args.max_mismatch} 错配): {len(sites)}")
    print(f"  按错配数: {dict(sorted(by_mm.items()))}  按风险: {by_risk}")

    out_csv = args.out + '.sites.csv'
    fields = ['transcript', 'gene', 'pos', 'mismatches', 'seed_mm', 'distal_mm',
              'window', 'pfs', 'pfs_mm', 'mm_positions', 'risk']
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sites)
    payload = {'spacer_dna': spacer, 'target_scanned': target, 'pfs_rule': pfs,
               'pfs_tiers': 'pfs_mm 0=高危 1=中危 >=2=低危(Cas12a2 PFS 松散, Zeng 2026)',
               'seed_len_window5p': SEED_LEN,
               'max_mismatch': args.max_mismatch,
               'n_transcripts': len(records), 'n_sites': len(sites),
               'by_mismatch': by_mm, 'by_risk': by_risk,
               'high_risk_sites': [s for s in sites if s['risk'] == 'high'],
               'medium_risk_sites': [s for s in sites if s['risk'] == 'medium'],
               'disclaimer': '纯序列近似; 激活与否取决于表达量/结构可及性/实验条件, '
                             '最终特异性以湿实验为准'}
    out_json = args.out + '.summary.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f"\n输出: {out_csv} / {out_json}")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        ax = axes[0]
        mms = sorted(by_mm)
        ax.bar([str(m) for m in mms], [by_mm[m] for m in mms], color='#4878CF')
        ax.set_xlabel('mismatches')
        ax.set_ylabel('sites')
        ax.set_title(f'Off-target sites by mismatch (n={len(sites)})')
        ax = axes[1]
        tiers = ['high', 'medium', 'low']
        ax.bar(tiers, [by_risk.get(t, 0) for t in tiers],
               color=['#C05630', '#E0B070', '#8FBF9F'])
        ax.set_ylabel('sites')
        ax.set_title('Risk tiers (mismatch + seed + PFS)')
        fig.suptitle(f'Transcriptome off-target scan ({len(records)} transcripts)')
        fig.tight_layout()
        out_png = args.out + '.png'
        fig.savefig(out_png, dpi=150)
        print(f'图: {out_png}')
    except ImportError:
        print('(matplotlib 不可用, 跳过图)')

    hi = [s for s in sites if s['risk'] == 'high']
    print('\n高危位点(0-1 错配 + 种子区完整 + PFS≤1mm):')
    for s in hi[:20]:
        print(f"  {s['gene'] or s['transcript']:<16} {s['mismatches']}mm "
              f"(seed {s['seed_mm']}/distal {s['distal_mm']}) pfs {s['pfs']}({s['pfs_mm']}mm)")


if __name__ == '__main__':
    main()
