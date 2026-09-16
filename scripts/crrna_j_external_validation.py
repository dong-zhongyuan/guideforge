# -*- coding: utf-8 -*-
"""§J 赢家规则的 Cas12a2 体系外部一致性检验 + §L 二轮提议协议预登记载体(2026-09-16).

外部检验(§J-ext): Han 2025 工具箱 7 骨架(Cas12a2 体系, 实测 fig1g 低=抑制强)
对 Su 注册表做最优偏移对齐后, 计算各骨架的 §J 规则匹配数/得分,
与 fig1g 的 Spearman —— 规则在独立体系的方向一致性。
n=7 功效极低, 预登记口径为"方向一致性检查", 任何结果不构成替换/否决依据。

§L 预登记(登记先于湿数据): BayesOpt 二轮提议协议参数固化, 见
docs/preregistration.md §L 与本脚本 L_PROTOCOL 字典(与 crrna_bayesopt.py
实现一致性由 tests 保证)。

输出: data/j_rule_external_validation.json
运行: PYTHONUTF8=1 python scripts/crrna_j_external_validation.py
"""
import json
import os

from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data')

SU_WT = 'AAUUUCUACUGUUGUAGAU'

# §L 二轮提议协议(预登记, 先于任何湿实验数据回填)
L_PROTOCOL = {
    'trigger': 'v6 panel 6 条构建的 n>=5 细胞存活数据经 --ingest-cell 回填后',
    'gp': 'RBF 核, alpha=1e-2, normalize_input=True(与 crrna_bayesopt.py 现实现一致)',
    'y': '每构建 z(delta_kill) = z(WT存活 - 变体存活), 每靶内标准化',
    'acquisition': 'EI(xi=0.01), 在全变体枚举空间(硬过滤通过者)上取 top3',
    'proposal_rules': 'top3 中排除已测构建; 若 EI top3 全部低于已测最优的 EI, '
                      '停止二轮(提前停止条款)',
    'reporting': '二轮提议与理由落盘 data/bayesopt_round2_proposals.json, '
                 '提议先于二轮实验登记',
}


def best_align(dr):
    """Han DR(20-21nt) 对 Su-WT19 的最优偏移对齐, 返回 (offset, 对齐19mer, 恒等分)。"""
    best = None
    for off in range(-2, 3):
        if off >= 0:
            seg = dr[off:off + 19]
        else:
            seg = dr[:19] if len(dr) >= 19 else dr
        seg = (seg + 'N' * 19)[:19]
        ident = sum(a == b for a, b in zip(seg, SU_WT))
        if best is None or ident > best[2]:
            best = (off, seg, ident)
    return best


def main():
    han = json.load(open(os.path.join(DATA, 'han2025_dataset.json'),
                         encoding='utf-8'))
    wr = json.load(open(os.path.join(DATA, 'winner_rule_engineering.json'),
                        encoding='utf-8'))
    rules = wr['rules_sig']
    fig = han['toolbox_activity_fig1g']

    rows = []
    for tag, full in sorted(han['toolbox_sequences'].items()):
        dr = full.replace('T', 'U')
        off, seg, ident = best_align(dr)
        matches = []
        for r in rules:
            pos, frm, to = r['su_pos'], r['from'], r['to']
            # 注册表差异容忍: 只要求对齐后该位为目标碱基(方向), from 不强求
            if seg[pos - 1] == to and SU_WT[pos - 1] != to:
                matches.append('%d%s' % (pos, to))
        rows.append({'tag': tag, 'aligned': seg, 'offset': off,
                     'identity_to_su_wt': ident,
                     'rule_matches': matches, 'n_matches': len(matches),
                     'fig1g': fig.get(tag)})
    valid = [r for r in rows if r['fig1g'] is not None]
    rho, p = spearmanr([r['n_matches'] for r in valid],
                       [r['fig1g'] for r in valid])

    out = {
        'generated_by': 'scripts/crrna_j_external_validation.py',
        'protocol': '§J-ext(方向一致性检查, n=7 功效极低, 不构成替换/否决依据)',
        'alignment_note': 'Han DR(20-21nt)对 Su-WT19 最优偏移对齐; '
                          '规则匹配容忍注册表 from 差异(只要求目标碱基落位)',
        'rows': rows,
        'spearman_nmatches_vs_fig1g': {'rho': round(float(rho), 3),
                                       'p': float(p)},
        'reading': ('期望方向: 规则匹配多 -> fig1g 低(抑制强) -> 负 rho; '
                    'n=7 仅作方向参考, 已在 §J-ext 口径中声明低功效'),
        'L_protocol_prereg': L_PROTOCOL,
    }
    dst = os.path.join(DATA, 'j_rule_external_validation.json')
    json.dump(out, open(dst, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('[J-ext] 对齐与规则匹配:')
    for r in rows:
        print('  %-4s ident=%2d/19 matches=%s fig1g=%s'
              % (r['tag'], r['identity_to_su_wt'],
                 r['rule_matches'] or '-', r['fig1g']))
    print('[J-ext] Spearman(n_matches vs fig1g) = %.3f (p=%.3f)' % (rho, p))
    print('->', dst)


if __name__ == '__main__':
    main()
