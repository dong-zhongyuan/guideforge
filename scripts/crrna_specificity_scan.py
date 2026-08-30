"""crRNA 特异性扫描（Cas12a2, RNA 靶向；体系切换后 crRNA 链）

定位：Cas12a2 由靶 RNA 激活，"脱靶"= 转录组中与 spacer 部分错配、可能错误激活的
其他 RNA。本模块在转录本序列（FASTA，sense 链）中扫描 spacer 互补区的近似位点，
并检查位点 3' 侧的 PFS（注册表 cas12a2 条目: GAAAG，位于靶 RNA 上）。

判定口径（与 PDB 8D4A 复合物一致）:
  - crRNA spacer 与靶 RNA 反向互补配对 → 扫描目标 = revcomp(spacer)
  - PFS 5'-GAAAG-3' 位于靶 RNA 互补区的 3' 下游（8D4A 靶 RNA 实测排布）
  - 错配数 <= --max-mismatch 的窗口计入位点表；0 错配+PFS = 预期靶点，
    1 错配位点尤其关键（对应 TP53 野生型 vs R248Q 突变型的单碱基鉴别场景）

用法:
  python crrna_specificity_scan.py --spacer <23nt> --fasta transcripts.fa
  python crrna_specificity_scan.py --spacer <23nt> --gene TP53   # 经 NCBI 拉 mRNA
输出: <out>.sites.csv + <out>.summary.json + 屏幕汇总
"""
import argparse
import csv
import json
import os
from urllib.parse import quote
from urllib.request import urlopen

def _build_comp():
    table = str.maketrans({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'})
    return table


COMP = _build_comp()


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def normalize(seq):
    return str(seq).upper().replace('U', 'T')


def fetch_gene_seq(gene, *, organism='Homo sapiens', db='nuccore', timeout=30.0,
                   base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils'):
    """从 NCBI 拉取基因 mRNA 序列（与旧链 design.py 同口径, 纯 urllib）。"""
    term = quote(f'{gene}[Gene] AND {organism}[Organism] AND mRNA[Filter]')
    url = f'{base_url}/esearch.fcgi?db={db}&term={term}&retmax=1&retmode=json'
    with urlopen(url, timeout=timeout) as resp:
        ids = json.loads(resp.read().decode('utf-8')).get('esearchresult', {}).get('idlist', [])
    if not ids:
        raise RuntimeError(f'NCBI 未找到基因 {gene} 的 mRNA 记录')
    url2 = f'{base_url}/efetch.fcgi?db={db}&id={ids[0]}&rettype=fasta&retmode=text'
    with urlopen(url2, timeout=timeout) as resp:
        lines = resp.read().decode('utf-8').splitlines()
    return ''.join(l.strip() for l in lines if not l.startswith('>'))


def read_fasta(path):
    """读多记录 FASTA, 返回 [(name, seq)]（DNA 字母表, U→T）。"""
    records, name, chunks = [], None, []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    records.append((name, normalize(''.join(chunks))))
                name, chunks = line[1:].split()[0] or f'record{len(records)}', []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, normalize(''.join(chunks))))
    if not records:
        raise ValueError(f'{path} 中没有 FASTA 记录')
    return records


def scan_record(name, seq, target, pfs_rule, max_mm):
    """在单条转录本上找 target 的 <=max_mm 错配窗口, 并读 3' 下游 PFS。"""
    L, lp = len(target), len(pfs_rule)
    sites = []
    for i in range(0, len(seq) - L - lp + 1):
        window = seq[i:i + L]
        mm = sum(a != b for a, b in zip(window, target))
        if mm <= max_mm:
            pfs = seq[i + L:i + L + lp]
            sites.append({'transcript': name, 'pos': i, 'mismatches': mm,
                          'window': window, 'pfs': pfs,
                          'pfs_match': pfs == pfs_rule})
    return sites


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacer', required=True, help='固定 spacer(17-25nt ACGT/U)')
    ap.add_argument('--fasta', default=None, help='转录本 FASTA(可多条)')
    ap.add_argument('--gene', default=None, help='基因名(经 NCBI 拉 mRNA, 与 --fasta 二选一)')
    ap.add_argument('--pfs', default='GAAAG', help='PFS 规则(默认注册表 cas12a2 口径 GAAAG)')
    ap.add_argument('--max-mismatch', type=int, default=4)
    ap.add_argument('--out', default='crrna_specificity', help='输出前缀')
    args = ap.parse_args()

    spacer = normalize(args.spacer)
    if not 17 <= len(spacer) <= 25 or set(spacer) - set('ACGT'):
        ap.error('--spacer 必须为 17-25nt ACGT/U')
    pfs_rule = normalize(args.pfs)
    target = revcomp(spacer)  # spacer 与靶 RNA 反向互补

    if args.fasta:
        records = read_fasta(args.fasta)
    elif args.gene:
        seq = fetch_gene_seq(args.gene)
        records = [(args.gene, normalize(seq))]
    else:
        ap.error('需要 --fasta 或 --gene')

    sites = []
    for name, seq in records:
        sites.extend(scan_record(name, seq, target, pfs_rule, args.max_mismatch))

    sites.sort(key=lambda s: (s['mismatches'], not s['pfs_match'], s['transcript'], s['pos']))
    by_mm = {}
    for s in sites:
        key = f"{s['mismatches']}mm" + ('+PFS' if s['pfs_match'] else '')
        by_mm[key] = by_mm.get(key, 0) + 1

    print(f'spacer({len(spacer)}nt): {spacer}')
    print(f'扫描目标(revcomp): {target}  PFS: {pfs_rule}(位点 3\' 下游)')
    print(f'转录本: {len(records)} 条  命中位点(<= {args.max_mismatch} 错配): {len(sites)}')
    for key in sorted(by_mm):
        print(f'  {key}: {by_mm[key]}')
    print('\n位点明细(错配升序):')
    print(f"{'转录本':<16}{'位置':>8}{'错配':>5}{'PFS':>8}{'窗口'}")
    for s in sites[:50]:
        print(f"{s['transcript']:<16}{s['pos']:>8}{s['mismatches']:>5}"
              f"{s['pfs'] + ('✓' if s['pfs_match'] else '✗'):>8}{s['window']}")

    out_csv = args.out + '.sites.csv'
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(sites[0].keys()) if sites else
                                ['transcript', 'pos', 'mismatches', 'window', 'pfs', 'pfs_match'])
        writer.writeheader()
        writer.writerows(sites)
    out_json = args.out + '.summary.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump({'spacer_dna': spacer, 'target_scanned': target, 'pfs_rule': pfs_rule,
                   'max_mismatch': args.max_mismatch, 'n_transcripts': len(records),
                   'n_sites': len(sites), 'by_mismatch': by_mm, 'sites': sites},
                  fh, ensure_ascii=False, indent=1)
    print(f'\n输出: {out_csv} / {out_json}')
    print('说明: 0 错配+PFS = 预期靶点; 1 错配位点对应单碱基鉴别压力(如 TP53 WT vs 突变);')
    print('      本扫描为纯序列近似评估, 激活与否还取决于实验条件, 最终特异性以湿实验为准。')


if __name__ == '__main__':
    main()
