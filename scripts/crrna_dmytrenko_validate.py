"""crRNA 骨架管线的 Dmytrenko 2023 文献验证（零湿实验成本）

依据（已逐条核对原文, 2026-08-31）:
  1. ED Fig.2c + 正文: Cas12a 与 Cas12a2 可互换向导而不损失免疫功能
     ("Cas12a and Cas12a2 can interchange guides without impairing immunity")
     → 检验A: 注册表 cas12a DR 作为互换变体过本管线结构口径, 应被判为合格
       （过滤器不误杀功能保持的大改 DR; 两 DR 长度不同, 结构距离不可定义,
        改用茎完整性 + spacer 游离 + 序列纪律三项判定, 如实注明口径）
  2. Fig.1c: 跨家族 DR 比对 3' 端高度保守、loop 区可变（图注原文
     "CRISPR repeats with a conserved 3' end"; "The loop of the hairpin is variable"）
     → 检验B: 单点扫描的位点耐受度 vs 保守性分区(3'窗/茎-悬垂/loop)的
       Spearman 相关 —— 打分逻辑是否与天然进化数据一致
  3. ED Fig.3: 加工切点位于 Cas12a 切点下游 1 nt, 成熟边界 = 结构口径 18nt DR
     → 已实现为主管线硬过滤(--proc-window), 本脚本不重复

用法:
  python crrna_dmytrenko_validate.py --variants-csv ../data/tp53_r248q_zengdr.variants.csv \
      --effector cas12a2_zeng2026 --spacer GTTCATGCCGCCCATGCAGGAACT \
      --out-prefix ../data/tp53_r248q_zengdr
输出: <prefix>.dmytrenko_validation.json + <prefix>.conservation_validation.png
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
from scaffold_registry import get_scaffold  # noqa: E402


def check_dr_swap(spacer_dna, effector_wt, effector_swap, margin=0.10):
    """检验A: 互换 DR 过结构口径（不误杀判定）。"""
    dr_wt = core.to_rna(get_scaffold(effector_wt))
    dr_sw = core.to_rna(get_scaffold(effector_swap))
    spacer = core.to_rna(spacer_dna)
    wt = core.wt_reference(dr_wt, spacer)

    ss_sw, mfe_sw = core.fold(dr_sw + spacer)
    _, spacer_up_sw = core.pf_stats(dr_sw + spacer, len(dr_sw))
    ss_dr_sw, dr_mfe_sw = core.fold(dr_sw)

    pairs_wt = sum(1 for c in wt['dr_only_struct'] if c in '()')
    pairs_sw = sum(1 for c in ss_dr_sw if c in '()')
    checks = {
        'stem_pairs_kept': {'wt_pairs': pairs_wt, 'swap_pairs': pairs_sw,
                            'ok': bool(pairs_sw >= pairs_wt - 2)},
        'spacer_unpaired': {'wt': wt['spacer_mean_unpaired'],
                            'swap': round(spacer_up_sw, 3),
                            'ok': bool(spacer_up_sw >= wt['spacer_mean_unpaired'] - margin)},
        'sequence_discipline': {'ok': 'TTTT' not in core.to_dna(dr_sw)
                                and 'GGGG' not in core.to_dna(dr_sw)},
    }
    verdict = all(v['ok'] for v in checks.values())
    note = ('两 DR 长度不同(%d vs %d nt), 与 WT 的 bp_distance 不可定义; '
            '改用茎完整性(配对数 >= WT-2)+spacer 游离+序列纪律三项, 口径如实注明'
            % (len(dr_sw), len(dr_wt)))
    return {'swap_effector': effector_swap, 'swap_dr_rna': dr_sw,
            'swap_dr_only_struct': ss_dr_sw, 'swap_dr_only_mfe': dr_mfe_sw,
            'full_mfe_kcal': round(mfe_sw, 2), 'checks': checks,
            'verdict': 'PASS(未误杀)' if verdict else 'FAIL(被过滤)',
            'note': note}


def conservation_partition(dr, cons3_window=5):
    """位置保守性分区(Fig.1c 定性口径): 3'保守窗=2, 茎/悬垂=1, loop=0。"""
    ss = core.fold(dr)[0]
    n = len(dr)
    part = {}
    loop_pos = set()
    # 点括号直接给配对; loop = 茎环之间的未配对段(取连续未配对且两侧为配对的最长段)
    unpaired_runs, run = [], []
    for i, ch in enumerate(ss):
        if ch == '.':
            run.append(i)
        else:
            if run:
                unpaired_runs.append(run)
                run = []
    if run:
        unpaired_runs.append(run)
    if unpaired_runs:
        loop_pos = set(max(unpaired_runs, key=len))
    for i in range(n):
        pos = i + 1
        if pos > n - cons3_window:
            part[pos] = 2          # 3' 保守窗(跨家族保守)
        elif i in loop_pos:
            part[pos] = 0          # loop(Fig.1c 灰色标注, 可变)
        else:
            part[pos] = 1          # 5' 悬垂/茎区
    return part, sorted(loop_pos)


def check_conservation_correlation(variants_csv, effector, cons3_window=5, cons3_weight=0.3):
    """检验B: 位点耐受度 × 保守性分区 Spearman + 置换检验。"""
    dr = core.to_rna(get_scaffold(effector))
    part, loop_idx = conservation_partition(dr, cons3_window)

    by_pos = {}
    with open(variants_csv, encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if int(row['n_mut']) != 1:
                continue
            p = int(row['mut_positions'][1:-1])
            # score_adj: 还原 3'保守窗惩罚后的原始指标分(排除 w_cons3 自证的循环性)
            adj = float(row['score']) + float(row.get('cons3_frac', 0) or 0) * cons3_weight
            by_pos.setdefault(p, []).append(
                (row['passed'] == 'True', float(row['score']), adj))

    rows = []
    for p, items in sorted(by_pos.items()):
        tol_pass = sum(1 for ok, _, _ in items if ok) / len(items)
        mean_score = float(np.mean([s for _, s, _ in items]))
        mean_adj = float(np.mean([a for _, _, a in items]))
        rows.append({'pos': p, 'cons_class': part.get(p, 1),
                     'tolerance_pass': round(tol_pass, 3),
                     'mean_score': round(mean_score, 4),
                     'mean_score_adj': round(mean_adj, 4), 'n_alt': len(items)})

    cons = np.array([r['cons_class'] for r in rows])
    tol = np.array([r['tolerance_pass'] for r in rows])
    msc = np.array([r['mean_score'] for r in rows])
    madj = np.array([r['mean_score_adj'] for r in rows])
    rho_tol = float(_spearman(cons, tol))
    rho_msc = float(_spearman(cons, msc))
    rho_adj = float(_spearman(cons, madj))
    p_tol = _perm_p(cons, tol, rho_tol)
    p_msc = _perm_p(cons, msc, rho_msc)
    p_adj = _perm_p(cons, madj, rho_adj)
    return {'rows': rows, 'loop_positions_1based': [i + 1 for i in loop_idx],
            'cons_partition': {str(k): v for k, v in part.items()},
            'spearman_cons_vs_tolerance': {'rho': round(rho_tol, 3), 'perm_p': round(p_tol, 4)},
            'spearman_cons_vs_meanscore': {'rho': round(rho_msc, 3), 'perm_p': round(p_msc, 4)},
            'spearman_cons_vs_scoreadj': {'rho': round(rho_adj, 3), 'perm_p': round(p_adj, 4),
                                           'note': 'score_adj 剔除 w_cons3 惩罚项, '
                                                   '为非循环的相关口径(仅 ddG/bp_dist/contact/hbond/ens 分量)'},
            'note': '预期方向为负(保守度高→耐受低); rho<0 且 p<0.05 即与天然进化数据一致; '
                    'n=19 位点, 置换检验功效有限, 未达显著时如实报告方向'}


def _spearman(x, y):
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return np.corrcoef(rx, ry)[0, 1]


def _perm_p(x, y, observed, n_perm=20000, seed=0):
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        if abs(_spearman(rng.permutation(x), y)) >= abs(observed):
            cnt += 1
    return (cnt + 1) / (n_perm + 1)


def make_figure(res_corr, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['WenQuanYi Zen Hei', 'Noto Sans CJK SC',
                                       'Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    rows = res_corr['rows']
    pos = [r['pos'] for r in rows]
    tol = [r['tolerance_pass'] for r in rows]
    msc = [r['mean_score'] for r in rows]
    cls = [r['cons_class'] for r in rows]
    cmap = {0: '#8FBF9F', 1: '#E0B070', 2: '#C05630'}
    labels = {0: 'loop(可变)', 1: '茎/悬垂', 2: "3'保守窗"}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.bar(pos, tol, color=[cmap[c] for c in cls], edgecolor='white')
    ax.set_xlabel('DR 位置 (nt)')
    ax.set_ylabel('单点耐受度 (通过硬过滤比例)')
    ax.set_title('位点敏感性 × 保守性分区')
    ax.set_ylim(0, 1.05)
    handles = [plt.Rectangle((0, 0), 1, 1, color=cmap[k]) for k in (2, 1, 0)]
    ax.legend(handles, [labels[k] for k in (2, 1, 0)], fontsize=8, loc='lower left')

    ax = axes[1]
    for c in (0, 1, 2):
        ys = [m for m, k in zip(msc, cls) if k == c]
        if ys:
            xs = np.random.default_rng(c).normal(c, 0.06, len(ys))
            ax.scatter(xs, ys, color=cmap[c], s=45, alpha=0.85, edgecolor='white',
                       label=labels[c], zorder=3)
            ax.scatter([c], [np.mean(ys)], color=cmap[c], marker='D', s=90,
                       edgecolor='black', linewidth=0.8, zorder=4)
    r1 = res_corr['spearman_cons_vs_meanscore']
    r2 = res_corr['spearman_cons_vs_tolerance']
    r3 = res_corr['spearman_cons_vs_scoreadj']
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([labels[k] for k in (0, 1, 2)], fontsize=9)
    ax.set_ylabel('单点平均 score')
    ax.set_title("全分 ρ=%s(p=%s)  去惩罚分 ρ=%s(p=%s)\n耐受度口径 ρ=%s (p=%s)"
                 % (r1['rho'], r1['perm_p'], r3['rho'], r3['perm_p'],
                    r2['rho'], r2['perm_p']), fontsize=9.5)
    ax.legend(fontsize=8, loc='lower left')

    fig.suptitle('打分逻辑的天然进化数据验证（Dmytrenko 2023 Fig.1c 口径）', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_png, dpi=180)
    print(f'[validate] 验证图 -> {out_png}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--variants-csv', required=True, help='单点扫描 variants.csv(zeng2026 口径)')
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--swap-effector', default='cas12a', help='互换检验用的 DR 来源条目')
    ap.add_argument('--spacer', default='GTTCATGCCGCCCATGCAGGAACT')
    ap.add_argument('--cons3-window', type=int, default=5)
    ap.add_argument('--w-cons3', type=float, default=0.3,
                    help='与主管线一致的 3\' 窗惩罚权重(用于 score_adj 还原)')
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    swap = check_dr_swap(args.spacer, args.effector, args.swap_effector)
    print('[检验A] DR 互换:', swap['verdict'])
    for k, v in swap['checks'].items():
        print('   ', k, v)

    corr = check_conservation_correlation(args.variants_csv, args.effector,
                                          args.cons3_window, args.w_cons3)
    print('[检验B] 保守性 vs 耐受度:', corr['spearman_cons_vs_tolerance'])
    print('[检验B] 保守性 vs 平均分:', corr['spearman_cons_vs_meanscore'])
    print('[检验B] 保守性 vs 去惩罚分:', corr['spearman_cons_vs_scoreadj'])

    payload = {'date': '2026-08-31', 'source': 'Dmytrenko et al. 2023 Nature (10.1038/s41586-022-05559-3)',
               'A_dr_swap': swap, 'B_conservation_corr': corr,
               'disclaimer': '检验A为过滤器不误杀检验(非活性预测); 检验B为计算内部一致性验证; 均不替代湿实验'}
    out_json = args.out_prefix + '.dmytrenko_validation.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f'[validate] 结果 -> {out_json}')
    make_figure(corr, args.out_prefix + '.conservation_validation.png')


if __name__ == '__main__':
    main()
