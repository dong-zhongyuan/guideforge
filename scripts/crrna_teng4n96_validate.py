"""文献提活变体回头验证——"4n96 检验"(Teng et al. 2019, Genome Biology, 10.1186/s13059-019-1620-8)

依据(序列来源已核对): 该文 Additional file 5(pUC19-U6-crRNA backbone partial
sequence, pUC19-U6-crRNA1 / pUC19-U6-crRNA4n96 两个构建, 同一位点提取):
  crRNA1 (LbCas12a WT 型, 19nt): AATTTCTACTATTGTAGAT
  crRNA4n96 (提活变体, 19nt):    AATTTCTACTATGGTAGAT  = crRNA1 + 位置13 T->G(loop 区)
正文口径: 4n96(loop UAUG) 编辑效率高于原始 scaffold。

检验逻辑(同论文内部对照, 避开注册表口径错位):
  以 crRNA1 为 WT 基准、4n96 为变体, 逐项跑本管线的结构口径指标——
  若管线对"文献实证提活"的单点变体给出非惩罚性评分, 则打分体系不误伤已知
  正向变体(与 DR 互换检验互补: 那个验功能保持, 这个验提活方向)。

判据(相对分位数, 不用绝对阈值): 4n96 的相对分(<=0)在参照分布
(默认 zengdr.v2 passed 变体 score 分布)中的百分位 >= 25% 即 PASS。

诚实边界: 4n96 为 Cas12a/哺乳编辑体系, Cas12a2 为 RNA 靶向远缘体系;
本检验只证明"打分不与文献正向证据矛盾", 不证明 Cas12a2 上提活。
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
sys.path.insert(0, os.path.join(ROOT, 'src'))

import RNA  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

WT_DR = 'AATTTCTACTATTGTAGAT'
V4N96_DR = 'AATTTCTACTATGGTAGAT'
PROVENANCE = {
    'paper': 'Teng et al. 2019, Genome Biology, doi:10.1186/s13059-019-1620-8',
    'source': 'Additional file 5 backbone: pUC19-U6-crRNA1 与 pUC19-U6-crRNA4n96 同一位点提取',
    'wt_dr_dna': WT_DR, 'v4n96_dr_dna': V4N96_DR,
    'difference': 'position 13 T->G (loop region; 正文口径 loop UAUG)',
    'note': 'Teng crRNA1 与注册表 cas12a2 PDB 口径 DR(18nt) 仅差 5\' 端 1nt, 故不用注册表 cas12a 占位条目作基准',
}


def percentile_of(x, ref):
    return 100.0 * sum(1 for v in ref if v <= x) / len(ref)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacer', default='GTTCATGCCGCCCATGCAGGAACT')
    ap.add_argument('--reference-csv', default=os.path.join(
        ROOT, 'data', 'tp53_r248q_zengdr.v2.variants.csv'),
        help='参照分数分布(passed 变体的 score 列)')
    ap.add_argument('--contact-dr', default=None,
                    help='接触表锚定 DR(默认直接按 WT_DR 计算, 对齐转移)')
    ap.add_argument('--out', default=os.path.join(ROOT, 'data', 'teng4n96_validation.json'))
    args = ap.parse_args()

    spacer = core.to_rna(args.spacer)
    wt_dr, v_dr = core.to_rna(WT_DR), core.to_rna(V4N96_DR)
    n = len(wt_dr)
    assert n == len(v_dr) and sum(a != b for a, b in zip(wt_dr, v_dr)) == 1

    # WT(crRNA1) 基准
    full_wt = wt_dr + spacer
    ss_wt, mfe_wt = core.fold(full_wt)
    _, up_wt = core.pf_stats(full_wt, n)
    dr_ss_wt, _ = core.fold(wt_dr)

    # 变体(4n96)
    full_v = v_dr + spacer
    ss_v, mfe_v = core.fold(full_v)
    _, up_v = core.pf_stats(full_v, n)
    bp_dist = RNA.bp_distance(ss_wt, ss_v)
    ddg = mfe_v - mfe_wt

    # 接触/氢键(按 WT_DR 锚定的 8D4A 接触表, 对齐转移)
    contact = core.load_contacts(wt_dr, enabled=True)
    diff_pos = [i + 1 for i in range(n) if wt_dr[i] != v_dr[i]]
    csum = core.contact_sum(diff_pos, contact)
    hloss, hfrac = core.hbond_loss(wt_dr, v_dr, diff_pos, contact)
    contact_note = ('对齐转移' if contact and contact.get('transfer') else '精确锚定')

    # 相对打分(同主管线权重, 3'保守窗: 位置13 不在末5nt窗内)
    w = dict(w_bp=0.3, w_ddg=0.1, w_contact=1.0, w_ens=0.05, w_hbond=0.5, w_cons3=0.3)
    cons3_n = sum(1 for p in diff_pos if p > n - 5)
    rel_score = -(w['w_bp'] * bp_dist + w['w_ddg'] * max(ddg, 0) + w['w_contact'] * csum
                  + w['w_hbond'] * hloss + w['w_cons3'] * cons3_n / 5)

    # 硬过滤口径(以 crRNA1 为 WT)
    ok = (bp_dist <= 4 and up_v >= up_wt - 0.10
          and 'TTTT' not in core.to_dna(v_dr) and 'GGGG' not in core.to_dna(v_dr))

    # 相对分位数判据
    ref = []
    if os.path.isfile(args.reference_csv):
        with open(args.reference_csv, encoding='utf-8') as fh:
            for row in csv.DictReader(fh):
                if row.get('passed') == 'True':
                    ref.append(float(row['score']))
    pct = percentile_of(rel_score, ref) if ref else None
    verdict = ('PASS(不误伤): 分数百分位 >= 25%' if pct is not None and pct >= 25
               else ('FAIL: 落入参照分布下四分位以下' if pct is not None else '无参照分布, 仅报告数值'))

    out = {
        'date': '2026-08-31', 'provenance': PROVENANCE,
        'wt_full_mfe_kcal': round(mfe_wt, 2), 'variant_full_mfe_kcal': round(mfe_v, 2),
        'metrics': {'bp_dist': bp_dist, 'ddG': round(ddg, 2),
                    'spacer_unpaired_wt': round(up_wt, 3),
                    'spacer_unpaired_variant': round(up_v, 3),
                    'contact_sum': round(csum, 3), 'hbond_loss': round(hloss, 3),
                    'cons3_hits': cons3_n, 'contact_mode': contact_note},
        'relative_score_at_pipeline_weights': round(rel_score, 4),
        'hard_filters_pass': bool(ok),
        'reference_distribution': {'file': os.path.basename(args.reference_csv),
                                   'n_passed': len(ref),
                                   'percentile_of_variant': round(pct, 1) if pct is not None else None},
        'verdict': verdict,
        'disclaimer': 'Cas12a/哺乳体系证据; 检验对象是打分体系与文献正向证据的一致性, 非 Cas12a2 提活证明',
    }
    with open(args.out, 'w', encoding='utf-8') as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    print(f'[4n96] WT(crRNA1)={core.to_dna(wt_dr)}  4n96={core.to_dna(v_dr)}  差异位置 {diff_pos}')
    print(f'[4n96] bp_dist={bp_dist}  ddG={ddg:+.2f}  spacer_up {up_wt:.3f}->{up_v:.3f}  '
          f'contact={csum:.3f}({contact_note})  hbond_loss={hloss:.3f}  cons3命中={cons3_n}')
    print(f'[4n96] 硬过滤(以 crRNA1 为基准): {"通过" if ok else "未通过"}')
    print(f'[4n96] 相对分(主管线权重) = {rel_score:.4f}  参照分布百分位 = '
          f'{"%.1f" % pct if pct is not None else "N/A"}')
    print(f'[4n96] 判定: {verdict}')
    print(f'[4n96] 输出 -> {args.out}')


if __name__ == '__main__':
    main()
