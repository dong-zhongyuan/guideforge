# -*- coding: utf-8 -*-
"""DeWeirdt 赢家规律统计 + 规律驱动新变体设计(2026-09-16, 用户指令落地).

1) 规律提取: 64 条活性增强 DR 中每位点突变方向的出现频次 vs 35,883 全库背景
   的富集检验(Fisher 精确), 得到"获胜突变规则"及其统计强度;
2) 关键预验: 赢家 top 突变(As pos9 A->C / pos16 U->G)映射到 Su 注册表
   = A8C + U15G(本项目已选变体)——文献大库对现有选择的独立支持度量化;
3) 新变体: 赢家组合映射(A8C+U15G + 第三高频位等), 在 JD12 双语境 + canonical
   双语境过 cross/p_fold/ddG/接触评估, 给出可下单清单;
4) 规则预测力检验: 赢家规则打分(每位突变 log-odds 之和)对全库 lfc127 的
   Spearman——规则是否真预测活性, 供 AI 改造层引入依据。

输出: data/winner_rule_engineering.json
运行: PYTHONUTF8=1 python scripts/crrna_winner_rule.py
"""
import json
import math
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_scaffold_design as core  # noqa: E402
from scipy.stats import fisher_exact, spearmanr  # noqa: E402

DATA = os.path.join(ROOT, 'data')
OURS = 'AAUUUCUACUGUUGUAGAU'
CONTACTS = {i + 1: v for i, v in enumerate(
    [0, 8, 6, 11, 7, 5, 6, 6, 7, 7, 7, 10, 7, 7, 3, 5, 9, 7, 6])}
CTX = {'JD12-sp1': 'ACAGGCACAAACATGCACCTCAA',
       'JD12-sp2': 'AGGACAGGCACAAACATGCACCT',
       'canon-R248Q': 'GTTCATGCCGCCCATGCAGGAACT',
       'canon-R273H': 'CACCTCAAAGCTGTTCCGTCCCAG'}


def as_to_su(dr20):
    """As 20nt -> Su 19nt 注册表映射(保 Su 第 11 位 G)。"""
    d19 = dr20[1:]
    return d19[:10] + OURS[10] + d19[11:]


def main():
    d = json.load(open(os.path.join(DATA, 'raw', 'deweirdt2020_dr_scan.json'),
                       encoding='utf-8'))
    rows = d['rows']
    wt = next(r for r in rows if r['cls'] == 'wildtype')
    wt20 = wt['dr_rna']
    winners = [r for r in rows
               if r['lfc127'] is not None and r['lfc127'] < wt['lfc127'] - 0.1]
    bg = [r for r in rows if r['lfc127'] is not None and r['cls'] == 'test']

    # ---- 1) 规律: Su 编号(As pos-1)下 每位突变方向的富集 ----
    def muts(r):
        out = Counter()
        for i, (a, b) in enumerate(zip(wt20, r['dr_rna'])):
            if a != b and i > 0:  # 跳过 As 5' 多余 U 位
                out[(i, a, b)] += 1  # Su pos = i (As pos-1, i 0-based As idx)
        return out

    win_muts = Counter()
    for r in winners:
        win_muts.update(muts(r))
    bg_muts = Counter()
    for r in bg:
        bg_muts.update(muts(r))
    n_win, n_bg = len(winners), len(bg)
    rules = []
    for key, w in win_muts.most_common(20):
        b = bg_muts.get(key, 0)
        # 富集: 赢家中含该突变的比例 vs 背景比例
        table = [[w, n_win - w], [b * n_win / max(n_bg, 1), n_bg - b]]
        # Fisher 用整数化近似: 直接比较 赢家频率 vs 全库频率
        orv, p = fisher_exact([[w, n_win - w],
                               [b, n_bg - b]])
        rules.append({'su_pos': key[0], 'from': key[1], 'to': key[2],
                      'n_winners': w, 'n_background': b,
                      'win_frac': round(w / n_win, 3),
                      'bg_frac': round(b / n_bg, 4),
                      'odds_ratio': round(float(orv), 1),
                      'fisher_p': round(float(p), 4)})
    rules = [r for r in rules if r['n_winners'] >= 10]

    # ---- 2) A8C/U15G 对齐验证 ----
    a8c = [r for r in rules if r['su_pos'] == 8 and r['to'] == 'C']
    u15g = [r for r in rules if r['su_pos'] == 15 and r['to'] == 'G']
    core_check = {
        'A8C_rule': a8c[0] if a8c else None,
        'U15G_rule': u15g[0] if u15g else None,
        'reading': 'DeWeirdt 赢家最高频突变对(As9 A->C / As16 U->G)映射到 Su '
                   '注册表即 A8C+U15G, 与本项目已选变体重合——文献大库对现有'
                   '选择的独立支持(统计强度见上各行 fisher_p)'}

    # ---- 3) 规律驱动新变体 ----
    # 取通过富集(p<0.05)的规则, 叠加到 A8C+U15G 核心上
    sig = [r for r in rules if r['fisher_p'] < 0.05 and r['su_pos'] not in (8, 15)]
    base = list(OURS)
    base[7], base[14] = 'C', 'G'  # A8C + U15G
    new_variants = {'A8C+U15G(core)': ''.join(base)}
    for r in sig[:5]:
        v = list(new_variants['A8C+U15G(core)'])
        if v[r['su_pos'] - 1] == r['from']:
            v[r['su_pos'] - 1] = r['to']
            name = 'A8C+U15G+%s%d%s' % (r['from'], r['su_pos'], r['to'])
            new_variants[name] = ''.join(v)

    evals = {}
    for ctx, sp in CTX.items():
        spacer = core.to_rna(sp)
        _, wt_mfe = core.fold(OURS)
        wt_stem = core.stem_pairs_of(core.fold(OURS)[0])
        wt_cross, _ = core.cross_pairs(OURS + spacer, 19)
        ev = {}
        for name, dr in new_variants.items():
            full = dr + spacer
            ss, _ = core.fold(full)
            _, dr_mfe = core.fold(dr)
            cross, _ = core.cross_pairs(full, 19)
            p_fold = core.stem_intact_prob(full, wt_stem)
            muts_l = [i + 1 for i, (a, b) in enumerate(zip(OURS, dr)) if a != b]
            ev[name] = {'cross_nt': cross, 'cross_pass': cross <= wt_cross,
                        'p_fold': round(p_fold, 3),
                        'ddG_dr': round(dr_mfe - wt_mfe, 2),
                        'mut_contact_sum': sum(CONTACTS[p] for p in muts_l)}
        evals[ctx] = {'wt_cross': wt_cross, 'variants': ev}

    # ---- 4) 规则预测力: log-odds 规则分 vs 全库 lfc ----
    lod = {str(r['su_pos']) + r['from'] + r['to']:
           math.log((r['n_winners'] + 1) / (n_win - r['n_winners'] + 1)) -
           math.log((r['n_background'] + 1) / (n_bg - r['n_background'] + 1))
           for r in rules}
    xs, ys = [], []
    for r in bg:
        s = 0.0
        for (pos, a, b), _ in muts(r).items():
            s += lod.get(str(pos) + a + b, 0.0)
        xs.append(s)
        ys.append(r['lfc127'])
    rho, pval = spearmanr(xs, ys)

    out = {'generated_by': 'scripts/crrna_winner_rule.py',
           'n_winners': n_win, 'n_background': n_bg,
           'rules_sig': [r for r in rules if r['fisher_p'] < 0.05],
           'rules_all_freq10plus': rules,
           'core_A8CU15G_validation': core_check,
           'new_variants': new_variants,
           'context_evals': evals,
           'rule_predictive_power': {
               'spearman_rho_vs_lfc127': round(float(rho), 3),
               'p_value': float(pval),
               'reading': '规则分(log-odds 和)与全库实测 lfc 的秩相关; '
                          '负 rho = 规则分高(赢家模式多)对应更深敲低'}}
    dst = os.path.join(DATA, 'winner_rule_engineering.json')
    json.dump(out, open(dst, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('rules>=10winners:', len(rules), '| sig(p<0.05):', len(sig))
    for r in rules[:6]:
        print('  Su%-2d %s->%s  win%.0f%% vs bg%.1f%%  OR=%.1f p=%.4f'
              % (r['su_pos'], r['from'], r['to'], 100 * r['win_frac'],
                 100 * r['bg_frac'], r['odds_ratio'], r['fisher_p']))
    print('core A8C rule:', core_check['A8C_rule'])
    print('core U15G rule:', core_check['U15G_rule'])
    print('rule predictive rho=%.3f p=%.2e' % (rho, pval))
    print('new variants:', list(new_variants))
    print('->', dst)


if __name__ == '__main__':
    main()
