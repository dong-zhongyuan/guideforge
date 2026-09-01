"""合成优先级表生成（v2 候选 → 湿实验初筛清单）

从 crrna_directed_design.py 的产出中按"假设驱动"挑选首批合成对象:
  1. WT（对照, 必含）
  2. 茎稳定化假设: ddG_dr 最小且 hbond_loss=0 的共变变体（RNP 组装增强方向）
  3. 种子区释放假设: d_seed 最大的变体（验证该轴在 DR-only 空间内是否可忽略）
  4. 总分 TOP-N（综合排序参考）
每条标注假设类别与预期优势论述, 输出 CSV + FASTA（DNA 字母表, 可直接送合成）。

用法:
  python make_synthesis_priority.py --run-prefix ../data/tp53_r248q_zengdr_v2 \
      --n-stab 3 --n-seed 2 --n-top 4
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--run-prefix', required=True, help='v2 产出的 out-prefix')
    ap.add_argument('--n-stab', type=int, default=3, help='茎稳定化假设候选数')
    ap.add_argument('--n-seed', type=int, default=2, help='种子区释放假设候选数')
    ap.add_argument('--n-top', type=int, default=4, help='总分 TOP 参考数')
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.run_prefix + '.variants.csv', encoding='utf-8')))
    wt, variants = rows[0], rows[1:]
    passed = [r for r in variants if r['passed'] == 'True']

    picks, seen = [], set()

    def pick(row, hypothesis):
        if row['dr_seq'] in seen:
            return
        seen.add(row['dr_seq'])
        picks.append({'hypothesis': hypothesis, **row})

    pick(wt, 'WT 对照')
    stab = sorted([r for r in passed if float(r['ddG_dr']) < -0.1
                   and float(r['hbond_loss']) == 0.0],
                  key=lambda r: float(r['ddG_dr']))
    for r in stab[:args.n_stab]:
        pick(r, '茎稳定化(RNP组装增强)')
    seed = sorted(passed, key=lambda r: -float(r['d_seed']))
    for r in seed[:args.n_seed]:
        pick(r, "种子区释放(3'端可及性)")
    for r in sorted(passed, key=lambda r: -float(r['score']))[:args.n_top]:
        pick(r, '综合评分 TOP')

    out_csv = args.run_prefix + '.synthesis_priority.csv'
    fields = ['hypothesis', 'desc', 'dr_seq', 'construct_dna', 'ddG_dr', 'd_seed',
              'bp_dist', 'contact_sum', 'hbond_loss', 'score', 'advantage', 'source']
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for p in picks:
            writer.writerow({k: p[k] for k in fields})
    out_fasta = args.run_prefix + '.synthesis_priority.fasta'
    with open(out_fasta, 'w', encoding='utf-8') as fh:
        for p in picks:
            fh.write(f">{p['desc']}|{p['hypothesis']}\n{p['construct_dna']}\n")

    print(f"首批合成清单 {len(picks)} 条 → {out_csv} / {out_fasta}")
    for p in picks:
        print(f"  [{p['hypothesis']}] {p['desc']:<24} ddG_dr {p['ddG_dr']:>6} "
              f"d_seed {p['d_seed']:>6} score {p['score']:>8}  {p['advantage']}")


if __name__ == '__main__':
    main()
