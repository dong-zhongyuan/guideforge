# -*- coding: utf-8 -*-
"""DeWeirdt 2020 大库赢家 -> SuCas12a2 注册表映射与语境评估(2026-09-16).

背景(用户质询"为什么WT一直是最好/读文献了吗"的核查结论):
- "WT 是活性上限"此前的证据基础 = Han 2025 工具箱 7 骨架(0 条 trans 超 WT)
  + 选型器在我们点突变空间的预测(0/1300)——范围有限, 不得作为普遍结论;
- DeWeirdt 2020(Nat Biotechnol) 35,883 条 AsCas12a 替代 DR 实测扫描:
  64 条(0.2%)敲低深于 WT(lfc127 < WT-0.1), 最强 -4.88 vs WT -3.41;
- Dmytrenko 2023(Nature 613:588, ED Fig.2c): Cas12a<->Cas12a2 DR 互换功能保持
  ——跨酶移植有文献通路(边界: 互换实验是功能保持, 非"变强"证明)。

本脚本:
1) 解析 64 条赢家, 统计获胜突变模式(位置/方向/频次);
2) 映射到我们 19nt 注册表(As-WT 与 Su-WT 仅第 11 位差 C/G, 映射=保 Su 第11位);
3) 在 JD12 双 spacer 语境评 cross_nt(硬过滤口径 <=WT)/p_fold/ddG_dr/接触代价;
4) 输出 data/deweirdt_winner_map.json。

运行: PYTHONUTF8=1 python scripts/crrna_deweirdt_winner_map.py
"""
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_scaffold_design as core  # noqa: E402

DATA = os.path.join(ROOT, 'data')
OURS = 'AAUUUCUACUGUUGUAGAU'
CONTACTS = {i + 1: v for i, v in enumerate(
    [0, 8, 6, 11, 7, 5, 6, 6, 7, 7, 7, 10, 7, 7, 3, 5, 9, 7, 6])}
CTX = {'JD12-sp1': 'ACAGGCACAAACATGCACCTCAA',
       'JD12-sp2': 'AGGACAGGCACAAACATGCACCT'}


def main():
    d = json.load(open(os.path.join(DATA, 'raw', 'deweirdt2020_dr_scan.json'),
                       encoding='utf-8'))
    rows = d['rows']
    wt = next(r for r in rows if r['cls'] == 'wildtype')
    winners = sorted((r for r in rows
                      if r['lfc127'] is not None
                      and r['lfc127'] < wt['lfc127'] - 0.1),
                     key=lambda r: r['lfc127'])
    cnt = Counter()
    for r in winners:
        for i, (a, b) in enumerate(zip(wt['dr_rna'], r['dr_rna'])):
            if a != b:
                cnt[(i + 1, a, b)] += 1
    pattern = [{'pos': k[0], 'from': k[1], 'to': k[2], 'n': v}
               for k, v in cnt.most_common(12)]

    cands = []
    for r in winners[:12]:
        dr19 = r['dr_rna'][1:]
        mapped = dr19[:10] + OURS[10] + dr19[11:]
        muts = [i + 1 for i, (a, b) in enumerate(zip(OURS, mapped)) if a != b]
        rec = {'tag': r['tag'], 'as_dr_rna': dr19, 'su_mapped_rna': mapped,
               'as_lfc127': r['lfc127'], 'n_mut_vs_su': len(muts),
               'mut_positions': muts,
               'mut_contact_sum': sum(CONTACTS[p] for p in muts),
               'contexts': {}}
        for ctx, sp in CTX.items():
            spacer = core.to_rna(sp)
            wt_full_ss, _ = core.fold(OURS + spacer)
            _, wt_mfe = core.fold(OURS)
            wt_stem = core.stem_pairs_of(core.fold(OURS)[0])
            wt_cross, _ = core.cross_pairs(OURS + spacer, 19)
            wt_pfold = core.stem_intact_prob(OURS + spacer, wt_stem)
            full = mapped + spacer
            ss, _ = core.fold(full)
            _, dr_mfe = core.fold(mapped)
            cross, _ = core.cross_pairs(full, 19)
            p_fold = core.stem_intact_prob(full, wt_stem)
            rec['contexts'][ctx] = {
                'wt_cross': wt_cross, 'cross_nt': cross,
                'cross_pass': cross <= wt_cross,
                'wt_p_fold': round(wt_pfold, 3), 'p_fold': round(p_fold, 3),
                'ddG_dr': round(dr_mfe - wt_mfe, 2)}
        cands.append(rec)

    out = {'generated_by': 'scripts/crrna_deweirdt_winner_map.py',
           'claim_correction': ('"WT 是活性上限"仅在 Han 2025 七骨架与我们'
                               '点突变预测空间内成立; DeWeirdt 大库 64 条实测'
                               '胜 WT, 文献口径下 DR 可以变强(胜率 0.2%, 需筛选)'),
           'wt_as_lfc127': wt['lfc127'], 'n_winners': len(winners),
           'best_as_lfc127': winners[0]['lfc127'],
           'winner_mutation_pattern': pattern,
           'transfer_basis': 'Dmytrenko 2023 Nature ED Fig.2c '
                             '(Cas12a<->Cas12a2 DR 互换功能保持; 边界: 保持非增强)',
           'top_candidates': cands}
    dst = os.path.join(DATA, 'deweirdt_winner_map.json')
    json.dump(out, open(dst, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('winners=%d best=%.2f (WT %.2f) -> %s'
          % (len(winners), winners[0]['lfc127'], wt['lfc127'], dst))
    for c in cands[:6]:
        p = {k: (v['cross_pass'], v['p_fold']) for k, v in c['contexts'].items()}
        print('  %-22s lfc=%.2f muts=%s cross/pfold=%s'
              % (c['tag'][-19:], c['as_lfc127'], c['mut_positions'], p))


if __name__ == '__main__':
    main()
