# -*- coding: utf-8 -*-
"""赢家方法论与 spacer 联系的完整分析(2026-09-16, 用户两次方法学纠错的响应)。

纠错史(如实):
1. 首版 §J 规律为边际统计(合并语境), 用户指出"要和 spacer 联系";
2. 语境拆分后: 赢家个体迁移率仅 42%; 协变承载者(9,104 条)赢率 ~0.3%,
   中位 Δlfc 比 WT 差 +2.9 —— "赢家富集该模式" != "该模式致强"(逆概率陷阱),
   §J 先验的处方性质被本次分析实质性削弱(保留为排序先验, 不得表述为增强配方);
3. 真配方在上位效应层: 赢家载体 63% 携带 13 位环伴随 / 32% 携带 11 位伴随,
   输家载体携带 2/12/3/10 位杂伴随 —— 配方 = 一条远端对协变 + 正确的环伴随。

本脚本固化全部分析并输出判读; 载体: data/raw/deweirdt2020_dr_scan.json。
输出: data/winner_recipe_spacerlinked.json
运行: PYTHONUTF8=1 python scripts/crrna_winner_recipe.py
"""
import json
import os
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data')
TH = 0.1
STEM = {(5, 18), (6, 17), (7, 16), (8, 15), (9, 14)}


def main():
    d = json.load(open(os.path.join(DATA, 'raw', 'deweirdt2020_dr_scan.json'),
                       encoding='utf-8'))
    rows = d['rows']
    wt = next(r for r in rows if r['cls'] == 'wildtype')
    wt20 = wt['dr_rna']

    def muts(r):
        return {i: (a, b) for i, (a, b) in enumerate(zip(wt20, r['dr_rna']))
                if a != b and i > 0}

    def has_pair(m, pa, pb):
        return pa in m and pb in m

    w127 = {r['tag'] for r in rows
            if r['lfc127'] is not None and r['lfc127'] < wt['lfc127'] - TH}
    w128 = {r['tag'] for r in rows
            if r['lfc128'] is not None and r['lfc128'] < wt['lfc128'] - TH}
    both = w127 & w128

    carriers = [r for r in rows if r['lfc127'] is not None
                and (has_pair(muts(r), 8, 15) or has_pair(muts(r), 5, 18))]
    win_car = [r for r in carriers if r['tag'] in w127]
    lose_car = [r for r in carriers if r['tag'] not in w127]

    def companion(grp):
        cm = Counter()
        for r in grp:
            for p in muts(r):
                if p not in (5, 8, 15, 18):
                    cm[p] += 1
        n = max(len(grp), 1)
        return {str(p): round(100 * c / n, 1) for p, c in cm.most_common(6)}

    # 三连模式的双语境表现(迁移王候选)
    pat = {}
    for name, key in (('(8,13,15)', (8, 13, 15)), ('(5,13,18)', (5, 13, 18)),
                      ('(8,15)', (8, 15)), ('(5,18)', (5, 18))):
        sub = [r for r in rows if r['lfc127'] is not None
               and r['lfc128'] is not None
               and all(k in muts(r) for k in key)
               and not any(p in muts(r) for p in (5, 8, 13, 15, 18)
                           if p not in key)]
        if not sub:
            continue
        d127 = np.array([r['lfc127'] - wt['lfc127'] for r in sub])
        d128 = np.array([r['lfc128'] - wt['lfc128'] for r in sub])
        pat[name] = {'n': len(sub),
                     'win127_pct': round(100 * (d127 < -TH).mean(), 1),
                     'win128_pct': round(100 * (d128 < -TH).mean(), 1),
                     'median_d127': round(float(np.median(d127)), 2),
                     'median_d128': round(float(np.median(d128)), 2)}

    out = {
        'generated_by': 'scripts/crrna_winner_recipe.py',
        'correction_history': [
            'v1 §J 边际规律(已削弱: 承载者赢率~0.3%, 逆概率陷阱, 用户纠错)',
            'v2 语境拆分: 个体迁移率 42%, 规则处方性不成立',
            'v3 上位效应层配方(本版): 远端对协变 + 13/11 位环伴随'],
        'context_split': {'winners_127': len(w127), 'winners_128': len(w128),
                          'both': len(both),
                          'transfer_rate_127_to_128': round(
                              len(both) / max(len(w127), 1), 3)},
        'carriers': {'n': len(carriers), 'win_rate_127': round(
            len(win_car) / max(len(carriers), 1), 4)},
        'companion_positions': {'winners': companion(win_car),
                                'losers': companion(lose_car[:4000])},
        'pattern_dual_context': pat,
        'recipe': {
            'iron_rules': ['永不单边破对(0/64 赢家含错配)',
                           '只动一条远端对 (8,15) 或 (5,18)',
                           '不碰 (7,16)(0/64, 反偏好)',
                           '不做双对(全库 0 例)'],
            'companion': '13 位环伴随最优(63% 赢家载体), 11 位次之(32%)',
            'mapped_candidates': {
                '(8,13,15)': 'A8C+U15G+U13C = AAUUUCUCCUGUCGGAGAU',
                '(5,13,18)': 'U5C+A18G+U13C = AAUUCCUACUGUCGUAGGU'},
            'boundary': 'As 体系实测(敲低读数, 激活上游共用论证迁移到旁切); '
                        '迁移通路 Dmytrenko DR 互换(保持非增强); '
                        '个体优势语境绑定 -> 候选需逐语境装配验证'},
    }
    dst = os.path.join(DATA, 'winner_recipe_spacerlinked.json')
    json.dump(out, open(dst, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('[recipe] 迁移率 %.0f%% | 载体赢率 %.1f%% | 伴随: 赢 %s vs 输 %s'
          % (100 * out['context_split']['transfer_rate_127_to_128'],
             100 * out['carriers']['win_rate_127'],
             out['companion_positions']['winners'],
             out['companion_positions']['losers']))
    for k, v in pat.items():
        print('[pattern] %-9s n=%4d win127=%.0f%% win128=%.0f%% med=%+.2f/%+.2f'
              % (k, v['n'], v['win127_pct'], v['win128_pct'],
                 v['median_d127'], v['median_d128']))
    print('->', dst)


if __name__ == '__main__':
    main()
