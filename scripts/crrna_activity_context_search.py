# -*- coding: utf-8 -*-
"""跨 spacer 语境的 DR 变体活性搜索(2026-09-16, 用户需求: 提高活性).

问题: 选型器(fig1g 归一化 GFP, 低=抑制强=活性强)特征随 spacer 语境变化
(分型学说: 每型各有最优 DR)。R248Q 语境下 0/204 变体胜 WT
(data/selector_ranked_candidates.csv)。本脚本在其他语境重算:
对每个语境跑全变体枚举管线 -> 按该语境 WT 为参照算 5 维特征 ->
复刻选型器(DecisionTree max_leaf=4, random_state=42, 16 样本)预测 ->
找 pred < WT_pred(即预测活性更强)的变体。

口径边界(如实): 选型器为同源小样本决策树(16 对, 叶=4), 预测为同源外推;
"胜 WT"指预测抑制读数更低, 不等于细胞杀伤更强; 任何命中仍需过
装配/耐受门槛并经湿实验裁决。

语境: canonical R248Q/R273H + JD12_sp2 + §A3 三型代表 9 条 = 12 个。
输出: data/activity_context_search.json
运行: PYTHONUTF8=1 python scripts/crrna_activity_context_search.py
"""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import RNA                                          # noqa: E402
from sklearn.tree import DecisionTreeRegressor      # noqa: E402

import crrna_scaffold_design as core                # noqa: E402
import crrna_train_selector as sel                  # noqa: E402

DATA = os.path.join(ROOT, 'data')
TMP = os.path.join(ROOT, 'tmp_actsearch')
OUT = os.path.join(DATA, 'activity_context_search.json')

CANON = {
    'canonical_R248Q': 'GTTCATGCCGCCCATGCAGGAACT',
    'canonical_R273H': 'CACCTCAAAGCTGTTCCGTCCCAG',
    'JD12_sp2': 'AGGACAGGCACAAACATGCACCT',
}


def load_rep_spacers():
    """§A3 三型代表 spacer(名 -> 序列), 从分型队列文本取。"""
    names = {}
    ref = json.load(open(os.path.join(
        DATA, 'context_typing', 'typing_v47_expansion.json'), encoding='utf-8'))
    pool = {}
    for fn in ('spacers.txt', 'spacers_v47_expansion.txt'):
        p = os.path.join(DATA, 'context_typing', fn)
        if os.path.exists(p):
            for line in open(p, encoding='utf-8'):
                if line.strip():
                    n, s = line.rstrip('\n').split('\t')
                    pool[n] = s
    for t, reps in ref['a3_enriched']['representatives'].items():
        for rn in reps:
            if rn in pool:
                names['t%s_%s' % (t, rn.replace(' ', '_'))] = pool[rn]
    return names


def train_selector():
    X, y, _, _, _ = sel.load_han()
    dt = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                               random_state=42)
    dt.fit(X, y)
    return dt


def features_for_context(spacer_dna, tmp_prefix):
    """跑管线取全变体, 按本语境 WT 参照算特征, 返回 (rows, wt_pred行)。"""
    subprocess.run([sys.executable, os.path.join(HERE, 'crrna_scaffold_design.py'),
                    '--effector', 'cas12a2_zeng2026', '--spacer', spacer_dna,
                    '--topk', '8', '--out-prefix', tmp_prefix],
                   capture_output=True, text=True, check=True)
    top = json.load(open(tmp_prefix + '.top.json', encoding='utf-8'))
    wt_dr = core.to_rna(top['top'][0]['dr_seq'])
    spacer = core.to_rna(spacer_dna)
    wt_full_ss, _ = core.fold(wt_dr + spacer)
    _, wt_dr_mfe = core.fold(wt_dr)
    wt_stem = core.stem_pairs_of(core.fold(wt_dr)[0])
    out = []
    with open(tmp_prefix + '.variants.csv', newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['passed'] not in ('True', 'true', '1'):
                continue
            dr = core.to_rna(r['dr_seq'])
            full = dr + spacer
            ss, mfe = core.fold(full)
            _, dr_mfe = core.fold(dr)
            cross_nt, _ = core.cross_pairs(full, len(dr))
            p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
            _, sp_up, _ = core.pf_stats(full, len(dr), 7)
            out.append({
                'desc': r['desc'], 'score': float(r['score']),
                'ddG_dr': round(dr_mfe - wt_dr_mfe, 2),
                'bp_dist': RNA.bp_distance(wt_full_ss, ss),
                'cross_nt': cross_nt,
                'p_fold': round(p_fold, 5),
                'spacer_up': round(sp_up, 3)})
    return out


def main():
    os.makedirs(TMP, exist_ok=True)
    dt = train_selector()
    contexts = dict(CANON)
    contexts.update(load_rep_spacers())
    print('[search] %d 个语境' % len(contexts))
    report = {'generated_by': 'scripts/crrna_activity_context_search.py',
              'direction': 'fig1g 低=抑制强; hit = pred < 该语境 WT pred',
              'model': 'DecisionTree(max_leaf=4, rs=42) 复刻自 '
                       'crrna_train_selector(16 对同源)',
              'contexts': {}}
    any_hit = False
    for cname, sp in sorted(contexts.items()):
        try:
            rows = features_for_context(sp, os.path.join(TMP, 'c_' + cname))
        except Exception as e:  # noqa: BLE001
            report['contexts'][cname] = {'error': str(e)[:120]}
            continue
        wt = next((r for r in rows if r['desc'] == 'WT'), None)
        wt_pred = float(dt.predict([[wt['ddG_dr'], wt['bp_dist'],
                                     wt['cross_nt'], wt['p_fold'],
                                     wt['spacer_up']]])[0]) if wt else None
        for r in rows:
            r['pred'] = round(float(dt.predict(
                [[r['ddG_dr'], r['bp_dist'], r['cross_nt'],
                  r['p_fold'], r['spacer_up']]])[0]), 4)
        hits = sorted((r for r in rows
                       if wt_pred is not None and r['desc'] != 'WT'
                       and r['pred'] < wt_pred - 1e-9),
                      key=lambda r: r['pred'])
        report['contexts'][cname] = {
            'spacer_dna': sp, 'n_variants': len(rows),
            'wt_pred': round(wt_pred, 4) if wt_pred is not None else None,
            'n_hits': len(hits),
            'hits': [{'desc': h['desc'], 'pred': h['pred'],
                      'gain_vs_wt': round(h['pred'] - wt_pred, 4),
                      'ddG_dr': h['ddG_dr'], 'pipeline_score': h['score']}
                     for h in hits[:12]]}
        if hits:
            any_hit = True
        print('[search] %-28s n=%3d WT_pred=%.4f hits=%d%s'
              % (cname, len(rows), wt_pred or -1, len(hits),
                 (' 最优: %s %.4f (%+.4f)' %
                  (hits[0]['desc'], hits[0]['pred'],
                   hits[0]['pred'] - wt_pred)) if hits else ''))
    report['any_context_has_hit'] = any_hit
    json.dump(report, open(OUT, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print('->', OUT, '| any_hit =', any_hit)


if __name__ == '__main__':
    main()
