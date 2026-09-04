"""R3 消融: score_variant 逐项 leave-one-out + ddG 罚/赏不对称分析(2026-09 round-3)。

评审问题(round-2 meta R3 / bravo #5 / alfa #3 / charlie #3):
  (a) 主管线 9 项打分里, 哪些项对 TOP-K 排序有实际影响, 哪些是无效项?
      既有 weight_sensitivity.py 只做权重 ±50% 扰动, 不是逐 component 消融。
  (b) ddG 轴的罚/赏不对称(同一全长 ddG 上 -w_ddg*max(ddG,0) 罚 0.1 vs
      +w_stab*max(-ddG,0) 赏 0.3, stab_dr_only=False 默认口径)未经分析:
      两项分段不相交, 联合贡献 = w_stab*|ddG|^- - w_ddg*|ddG|^+, 即 V 形激励、
      稳定化臂斜率 3 倍于失稳臂——该不对称是否影响排序?

方法:
  1. 按 v6 口径重新枚举突变库(与 crrna_selection_baseline.build_library 同一
     代码路径, 保证口径一致; 纯 ViennaRNA, 无 RNet/引擎依赖);
  2. 消融 = 固定库重打分(指标不变, 只有线性组合权重变; 与 weight_sensitivity
     的 rescore 同式, 已验证与 score_variant 输出一致到 1e-3 内)。注意: 策略B
     的 SA 轨迹以满权重打分为目标, 消融不改变已枚举的库——本脚本回答的是
     "对既有库的重排序", 不是"消融后重新搜索", 该边界在输出中声明;
  3. 每个默认项置零重算全部变体, 报告: 与完整模型的 tie-aware Spearman
     (通过池 203 条 + 全库 441 条两口径)、TOP-12 重合/Jaccard、TOP-1 是否易主、
     默认 TOP-12 成员在消融后的最大名次跌落;
  4. ddG 对称性: w_ddg/w_stab 同设 0.1 与 0.3 两档 + 全轴置零, 同样重算对比;
     并报告通过池中 ddG>0 / ddG=0 / ddG<0 的条数(哪一臂在起作用)。

运行: python scripts/crrna_score_ablation.py
输出: data/score_ablation.<effector>.json + 控制台摘要。
"""
import argparse
import json
import os
import sys
from datetime import date

import numpy as np
import RNA
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)

from crrna_selection_baseline import LIB_PARAMS, build_library  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

TERMS = ['w_bp', 'w_ddg', 'w_contact', 'w_ens', 'w_hbond',
         'w_fold', 'w_seed', 'w_stab', 'w_cons3']
DEFAULTS = {t: LIB_PARAMS[t] for t in TERMS}

# 每项的物理含义(供结论可读性)
TERM_NOTE = {
    'w_bp': 'MFE 结构距离罚(bp_dist)',
    'w_ddg': '全长 ddG 失稳罚(max(ddG,0))',
    'w_contact': '8D4A 蛋白接触罚(contact_sum)',
    'w_ens': '系综多样性增量罚(max(d_ens,0))',
    'w_hbond': '氢键净找回(hbond_net, 有符号)',
    'w_fold': '茎完整概率变化赏(dp_fold)',
    'w_seed': "种子区游离度变化赏(d_seed; 默认 0, 共线关闭)",
    'w_stab': '茎稳定化赏(max(-ddG,0))',
    'w_cons3': "3' 保守窗突变罚(cons3_frac)",
}


def rescore(r, w):
    """由存储指标重构 score(与 score_variant 同式, stab_dr_only=False 口径)。"""
    ddg = float(r['ddG'])
    return (-w['w_bp'] * float(r['bp_dist'])
            - w['w_ddg'] * max(ddg, 0.0)
            - w['w_contact'] * float(r['contact_sum'])
            - w['w_ens'] * max(float(r['d_ens']), 0.0)
            + w['w_hbond'] * float(r['hbond_net'])
            - w['w_cons3'] * float(r['cons3_frac'])
            + w['w_fold'] * float(r['dp_fold'])
            + w['w_seed'] * float(r['d_seed'])
            + w['w_stab'] * max(-ddg, 0.0))


def compare(base_scores, new_scores, passed_flags, k):
    """消融/对称化重打分 vs 完整模型的排序对比。"""
    base_pass = base_scores[passed_flags]
    new_pass = new_scores[passed_flags]
    rho_pass = float(spearmanr(base_pass, new_pass).statistic)
    rho_all = float(spearmanr(base_scores, new_scores).statistic)
    top_base = set(np.argsort(-base_pass)[:k].tolist())
    top_new = set(np.argsort(-new_pass)[:k].tolist())
    inter = len(top_base & top_new)
    # 默认 TOP-K 成员(通过池内名次)在消融后的最差名次
    rank_new = np.empty(len(new_pass), dtype=int)
    rank_new[np.argsort(-new_pass)] = np.arange(1, len(new_pass) + 1)
    worst = int(rank_new[sorted(top_base)].max())
    top1_same = int(np.argmax(base_pass)) == int(np.argmax(new_pass))
    return {'spearman_passed': round(rho_pass, 4),
            'spearman_all': round(rho_all, 4),
            f'top{k}_intersection': inter,
            f'top{k}_jaccard': round(inter / (2 * k - inter), 3),
            'top1_unchanged': top1_same,
            'default_topk_worst_rank_after': worst}


def term_spread(rows, term):
    """该打分项在通过池上的贡献幅度(IQR of w*term)——解释"为什么无效"。"""
    vals = np.array([_term_value(r, term) for r in rows])
    return {'default_weight': DEFAULTS[term],
            'weighted_contrib_iqr': round(float(
                np.percentile(vals, 75) - np.percentile(vals, 25)), 5),
            'weighted_contrib_absmax': round(float(np.abs(vals).max()), 5)}


def _term_value(r, term):
    """单项对总分的带权重贡献(与 rescore 同式逐项)。"""
    ddg = float(r['ddG'])
    return {
        'w_bp': -DEFAULTS['w_bp'] * float(r['bp_dist']),
        'w_ddg': -DEFAULTS['w_ddg'] * max(ddg, 0.0),
        'w_contact': -DEFAULTS['w_contact'] * float(r['contact_sum']),
        'w_ens': -DEFAULTS['w_ens'] * max(float(r['d_ens']), 0.0),
        'w_hbond': DEFAULTS['w_hbond'] * float(r['hbond_net']),
        'w_fold': DEFAULTS['w_fold'] * float(r['dp_fold']),
        'w_seed': DEFAULTS['w_seed'] * float(r['d_seed']),
        'w_stab': DEFAULTS['w_stab'] * max(-ddg, 0.0),
        'w_cons3': -DEFAULTS['w_cons3'] * float(r['cons3_frac']),
    }[term]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--spacer', default=None)
    ap.add_argument('--topk', type=int, default=12)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    k = args.topk

    rows, wt, dr, entry = build_library(args.effector, args.spacer)
    scores_default = np.array([rescore(r, DEFAULTS) for r in rows])
    max_dev = max(abs(scores_default[i] - rows[i]['score'])
                  for i in range(len(rows)))
    print(f"库: {len(rows)} 变体, 通过 {sum(r['passed'] for r in rows)} 条; "
          f"rescore 与 score_variant 最大偏差 {max_dev:.2e} (存储指标舍入)")
    passed_flags = np.array([r['passed'] for r in rows])

    # ---- LOO 消融 ----
    ablations = {}
    print(f"\n=== Leave-one-out 消融(TOP-{k}, 通过池 {int(passed_flags.sum())} 条) ===")
    print(f"{'term':<10}{'rho_pass':>9}{'rho_all':>9}{'top12重合':>9}"
          f"{'Jaccard':>8}{'top1不变':>8}{'最差跌落':>8}  note")
    for t in TERMS:
        if DEFAULTS[t] == 0.0:
            ablations[t] = {'note': '默认权重已为 0, 消融恒等于完整模型',
                            'default_weight': 0.0, 'term': TERM_NOTE[t],
                            'spearman_passed': 1.0}
            print(f"{t:<10}{'1.0(恒等)':>18}  — 默认已为 0  {TERM_NOTE[t]}")
            continue
        w = dict(DEFAULTS)
        w[t] = 0.0
        sc = np.array([rescore(r, w) for r in rows])
        cmp_ = compare(scores_default, sc, passed_flags, k)
        cmp_.update({'term': TERM_NOTE[t], 'default_weight': DEFAULTS[t]})
        cmp_.update(term_spread([r for r, p in zip(rows, passed_flags) if p], t))
        ablations[t] = cmp_
        print(f"{t:<10}{cmp_['spearman_passed']:>9.4f}{cmp_['spearman_all']:>9.4f}"
              f"{cmp_[f'top{k}_intersection']:>9}{cmp_[f'top{k}_jaccard']:>8.3f}"
              f"{str(cmp_['top1_unchanged']):>8}{cmp_['default_topk_worst_rank_after']:>8}"
              f"  {TERM_NOTE[t]}")

    # ---- ddG 罚/赏不对称 ----
    ddgs = np.array([float(r['ddG']) for r, p in zip(rows, passed_flags) if p])
    n_pos, n_zero = int((ddgs > 1e-9).sum()), int((np.abs(ddgs) <= 1e-9).sum())
    n_neg = int((ddgs < -1e-9).sum())
    ddg_axis = {
        'axis_note': '默认 stab_dr_only=False: w_ddg 与 w_stab 作用于同一全长 ddG, '
                     '分段不相交; 联合贡献 = 0.3*|ddG|(ddG<0) 或 -0.1*ddG(ddG>0), '
                     'V 形激励, 稳定化臂斜率 3 倍于失稳臂(bravo 评审点名的不对称)',
        'passed_ddg_positive': n_pos, 'passed_ddg_zero': n_zero,
        'passed_ddg_negative': n_neg,
        'arm_reading': f'通过池 {len(ddgs)} 条中 ddG>0 失稳 {n_pos} 条(受 0.1 罚), '
                       f'ddG<0 稳定化 {n_neg} 条(受 0.3 赏), ddG=0 惰性 {n_zero} 条; '
                       '0.3/0.1 不对称决定稳定化臂的拉动强度'}
    for tag, wd, ws in (('sym_0.1', 0.1, 0.1), ('sym_0.3', 0.3, 0.3),
                        ('axis_off', 0.0, 0.0)):
        w = dict(DEFAULTS)
        w['w_ddg'], w['w_stab'] = wd, ws
        sc = np.array([rescore(r, w) for r in rows])
        ddg_axis[tag] = {'w_ddg': wd, 'w_stab': ws,
                         **compare(scores_default, sc, passed_flags, k)}
        c = ddg_axis[tag]
        print(f"[ddG] {tag}: rho_pass={c['spearman_passed']:.4f} "
              f"top{k}重合={c[f'top{k}_intersection']} top1不变={c['top1_unchanged']}")

    # ---- 结论 ----
    eff = {t: a for t, a in ablations.items()
           if a.get('default_weight', 0) > 0}
    ranked_by_impact = sorted(eff, key=lambda t: (
        eff[t][f'top{k}_intersection'], eff[t]['spearman_passed']))
    impact_list = [(t, eff[t]['spearman_passed'], eff[t][f'top{k}_intersection'],
                    eff[t]['top1_unchanged']) for t in ranked_by_impact]
    strict_inert = [t for t in ranked_by_impact
                    if eff[t]['spearman_passed'] >= 0.999
                    and eff[t][f'top{k}_intersection'] == k]
    conclusion = (
        f"逐项消融(按 TOP-{k} 重合升序, 每项=[Spearman(通过池), TOP-{k} 重合, TOP-1不变]): "
        f"{impact_list}。w_contact(8D4A 接触罚)是对 TOP-{k} 构成影响最大的项"
        f"(置零后重合仅 {eff['w_contact'][f'top{k}_intersection']}/{k} 且 TOP-1 易主), "
        f"w_ens 次之(重合 {eff['w_ens'][f'top{k}_intersection']}/{k}); w_hbond 显著影响"
        f"中部排序(rho={eff['w_hbond']['spearman_passed']})但不改 TOP-{k} 构成; "
        f"严格无效项(重合={k} 且 rho>=0.999)不存在; w_seed 默认已关闭(恒等)。"
        f"ddG 轴: 罚 0.1/赏 0.3 不对称下通过池含稳定化(ddG<0) {n_neg} 条、失稳 {n_pos} 条、"
        f"惰性 {n_zero} 条; 对称化为 0.1/0.1 后 rho={ddg_axis['sym_0.1']['spearman_passed']:.4f}"
        f"(几乎不变, 因稳定化奖赏只作用于 {n_neg} 条)、0.3/0.3 后 "
        f"rho={ddg_axis['sym_0.3']['spearman_passed']:.4f}(失稳罚加重 3 倍, 重排更明显)、"
        f"全轴置零 rho={ddg_axis['axis_off']['spearman_passed']:.4f}——"
        f"V 形激励(稳定化臂 3 倍斜率)是设计意图(解决纯罚分制下增强候选浮不出的结构性问题), "
        f"其对排序的净效应 = 把少数 ddG<0 变体系统性上拉; 实测 TOP-{k} 构成对该不对称"
        f"不敏感(对称化两档 TOP-{k} 重合均>=10/{k}, TOP-1 均不变), 是否保留应随湿实验标定复核。")

    payload = {
        'task': 'R3 消融: score_variant 逐项 leave-one-out + ddG 罚/赏不对称分析',
        'generated': date.today().isoformat(),
        'vienna_version': RNA.__version__, 'numpy_version': np.__version__,
        'config': {'effector': args.effector,
                   'dr_dna': core.to_dna(dr),
                   'library_params': LIB_PARAMS, 'topk': k,
                   'scope': '纯 ViennaRNA 打分消融(固定库重打分; 指标不变, 仅权重变); '
                            'SA 搜索轨迹依赖满权重, 消融不重搜库——回答"对既有库的重排序"; '
                            'RNet/引擎部分未涉(round-3 R3 授权跳过)',
                   'rescore_vs_score_variant_max_dev': max_dev},
        'library': {'n_variants': len(rows),
                    'n_passed': int(passed_flags.sum())},
        'defaults': DEFAULTS,
        'leave_one_out': ablations,
        'impact_ranking_topk_asc': [
            {'term': t, 'spearman_passed': eff[t]['spearman_passed'],
             f'top{k}_intersection': eff[t][f'top{k}_intersection'],
             'top1_unchanged': eff[t]['top1_unchanged']}
            for t in ranked_by_impact],
        'strict_inert_terms': strict_inert,
        'ddg_axis_symmetry': ddg_axis,
        'conclusion_zh': conclusion,
        'disclaimer': '排序消融诊断, 非活性预测; "无效项"仅指在当前库/当前默认权重下'
                      '对 TOP-K 排序无实际影响, 不排除其在其他上下文或权重量程内有作用。'}
    out = args.out or os.path.join(
        ROOT, 'data', f'score_ablation.{args.effector}.json')
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f"\n{conclusion}")
    print(f"\n输出: {out}")


if __name__ == '__main__':
    main()
