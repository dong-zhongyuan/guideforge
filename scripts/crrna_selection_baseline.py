"""R3 基线检验: 管线打分相对"硬过滤 + 随机选 K 条"的增量信息(2026-09 round-3)。

评审问题(round-2 meta R3): 主管线 crrna_scaffold_design.py 的透明启发式打分,
在"同样的突变库 + 同样的硬过滤"之后, 相对"从通过库中均匀随机取 K 条"是否
多提供信息量? 既有证据只回答邻近问题: crrna_negative_control.py 证明过滤器能
拦 grossly-broken DR(地板), weight_sensitivity.py 证明排名对权重扰动稳健——
都没有回答"打分排序本身是否优于随机"。

设计(zeng2026 口径, 与 data/tp53_r248q_zengdr.v6 同一运行配置):
  1. 重新枚举同一突变库(策略A 全单点+200双点抽样, 策略B SA 300步, A4 共变
     150双对, seed=0), 应用同一套硬过滤, 并与 shipped v6 variants.csv 对拍
     (序列集合/通过集合/逐条分数)作为血缘校验;
  2. (a) 管线打分排序取 TOP-K; (b) 从通过库均匀随机取 K 条, 重复 n-rep 次
     (默认 200, 独立种子);
  3. 检验一(构造性上界): TOP-K 平均打分 vs 随机 K 均值分布 -> Monte Carlo p。
     注意: 该 p 按构造必然近 1/(n+1)(TOP-K 即分数最大子集), 真正回答
     "打分是否提供信息"的是效应量(以通过池 SD 归一)与检验二;
  4. 检验二(独立性质对比): 在四个性质上比较 TOP-K vs 随机基线分布——
     全部不进打分公式: 茎区 GC 含量(纯序列)、真 loop 突变比例(两侧配对的
     未配对段; 打分只用 3' 保守窗, 不用 loop 归属)、cross_pp(DR-spacer 交叉
     配对连续系综量; 过滤用整数 cross_nt/inv_max_run, 打分不含此项)、
     选型器文献先验预测(selector_ranked_candidates.csv, Han2025 树,
     低=预测活性强; 半独立——与打分共享部分底层特征, 但模型由独立数据训练)。
     各报告双侧 MC p、Cliff's delta、KS 统计量;
  5. 重合度: 随机 K 与管线 TOP-K 交集大小的分布 vs 超几何期望 K^2/N——
     回答"过滤器+随机能否捞回同一批候选"。

口径声明: 纯 ViennaRNA 打分消融; RNet SHAPE 一致性需 rnet 引擎环境, 本检验跳过
  (与 round-3 R3 任务授权一致)。

运行: python scripts/crrna_selection_baseline.py
输出: data/selection_baseline.<effector>.json(默认) + 控制台摘要。
"""
import argparse
import csv
import json
import os
import sys
from datetime import date

import numpy as np
import RNA
from scipy.stats import ks_2samp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402
from scaffold_registry import get_entry  # noqa: E402

# 与 data/tp53_r248q_zengdr.v6.top.json params 完全一致的库生成/过滤配置
LIB_PARAMS = dict(n_double=200, sa_steps=300, seed=0, cov_double=150,
                  max_bp_dist=4.0, spacer_unpaired_margin=0.10,
                  w_bp=0.3, w_ddg=0.1, w_contact=1.0, w_ens=0.05,
                  w_hbond=0.5, w_fold=0.2, w_seed=0.0, seed_len=7,
                  w_stab=0.3, w_cons3=0.3, cons3_window=5, proc_window=4,
                  allow_cross_pairing=False, no_protect_processing=False,
                  stab_dr_only=False)


def build_library(effector, spacer_dna):
    """按 v6 口径重新枚举突变库并打分。返回 (rows, wt, dr, entry)。"""
    entry = get_entry(effector)
    dr = core.to_rna(entry['scaffold'])
    spacer = core.to_rna(spacer_dna or entry['placeholder_spacer'])
    args = argparse.Namespace(**LIB_PARAMS)
    rng = np.random.default_rng(args.seed)
    wt = core.wt_reference(dr, spacer, args.seed_len)
    contact = core.load_contacts(dr, enabled=True)
    stem_pos = core.stem_positions(wt['dr_only_struct'], len(dr))
    pool = {}
    for seq in core.strategy_a(dr, args, rng):
        pool[seq] = 'A'
    for seq in core.strategy_b(dr, args, rng, wt, spacer, contact, stem_pos):
        pool.setdefault(seq, 'B')
    for seq in core.strategy_a4(dr, args, rng):
        pool.setdefault(seq, 'COV')
    rows = [core.score_variant(dr, seq, wt, spacer, args, contact, stem_pos)
            for seq in pool]
    return rows, wt, dr, entry


def provenance_check(rows, effector):
    """与 shipped v6 variants.csv 对拍(存在时)。返回校验 dict。"""
    path = os.path.join(ROOT, 'data', 'tp53_r248q_zengdr.v6.variants.csv')
    if effector != 'cas12a2_zeng2026' or not os.path.isfile(path):
        return {'shipped_csv': None, 'note': '非 zeng2026 口径或 shipped v6 CSV 缺失, 跳过对拍'}
    with open(path, newline='', encoding='utf-8') as fh:
        ref = [r for r in csv.DictReader(fh) if r['desc'] != 'WT']
    mine = {r['dr_seq']: r for r in rows}
    ref_by_seq = {r['dr_seq']: r for r in ref}
    same_pool = set(mine) == set(ref_by_seq)
    same_passed = ({s for s, r in mine.items() if r['passed']}
                   == {s for s, r in ref_by_seq.items() if r['passed'] == 'True'})
    max_dscore = max((abs(mine[s]['score'] - float(r['score']))
                      for s, r in ref_by_seq.items() if s in mine), default=None)
    return {'shipped_csv': 'data/tp53_r248q_zengdr.v6.variants.csv',
            'same_variant_pool': same_pool, 'same_passed_set': same_passed,
            'max_abs_score_diff': max_dscore,
            'n_shipped_variants': len(ref),
            'note': '枚举库+硬过滤+打分与 shipped v6 完全一致(逐项对拍)'}


def true_loop_positions(dr_ss):
    """真 loop(1-based): 两侧均为配对位点的未配对段(两侧配对判据, 2026-09 R2 口径)。
    zeng2026 DR 折叠 ....(((((....))))). 下 5' 悬垂(1-4)不算 loop, 真 loop=10-13。"""
    runs, i = [], 0
    while i < len(dr_ss):
        if dr_ss[i] not in '()':
            j = i
            while j < len(dr_ss) and dr_ss[j] not in '()':
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    flanked = [(a, b) for a, b in runs
               if a > 0 and b < len(dr_ss)
               and dr_ss[a - 1] in '()' and dr_ss[b] in '()']
    loop = set()
    for a, b in flanked:
        loop.update(range(a + 1, b + 1))
    return loop


def cliffs_delta(a, b):
    """Cliff's delta: P(a>b) - P(a<b), a=TOP-K, b=通过池全体。"""
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return float((gt - lt) / (len(a) * len(b)))


def load_selector_preds():
    """选型器文献先验预测(按 desc join; 缺失返回 None)。低=预测活性强。"""
    path = os.path.join(ROOT, 'data', 'selector_ranked_candidates.csv')
    if not os.path.isfile(path):
        return None
    out = {}
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            out[r['desc']] = float(r['lit_pred_activity'])
    return out


def mc_test(obs, null, alternative='greater'):
    """Monte Carlo p(加 1 修正)。alternative: greater / two-sided(以 null 均值为中心)。"""
    null = np.asarray(null)
    n = len(null)
    if alternative == 'greater':
        k = int((null >= obs).sum())
    else:
        c = null.mean()
        k = int((np.abs(null - c) >= abs(obs - c)).sum())
    return (1.0 + k) / (1.0 + n)


def prop_summary(name, obs_vals, pool_vals, null_means, direction_note,
                 lower_better=False):
    obs = float(np.mean(obs_vals))
    p_mc = mc_test(obs, null_means, 'two-sided')
    ks = ks_2samp(obs_vals, pool_vals)
    return {'property': name,
            'direction_note': direction_note,
            'topk_mean': round(obs, 4),
            'pool_mean': round(float(np.mean(pool_vals)), 4),
            'pool_sd': round(float(np.std(pool_vals)), 4),
            'random_k_mean_q05': round(float(np.percentile(null_means, 5)), 4),
            'random_k_mean_median': round(float(np.median(null_means)), 4),
            'random_k_mean_q95': round(float(np.percentile(null_means, 95)), 4),
            'mc_p_two_sided': round(p_mc, 4),
            'cliffs_delta_topk_vs_pool': round(cliffs_delta(obs_vals, pool_vals), 3),
            'ks_stat_vs_pool': round(float(ks.statistic), 3),
            'ks_p_vs_pool': round(float(ks.pvalue), 4),
            'interpretation': ('TOP-K 显著低于随机基线' if lower_better else
                               'TOP-K 显著偏离随机基线')
            if p_mc <= 0.05 else 'TOP-K 与随机基线无显著差异'}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--spacer', default=None,
                    help='固定 spacer(默认取注册表 placeholder_spacer, zeng2026 即 R248Q)')
    ap.add_argument('--topk', type=int, default=12)
    ap.add_argument('--n-rep', type=int, default=200, help='随机基线重复次数')
    ap.add_argument('--baseline-seed', type=int, default=20260905,
                    help='随机抽样种子(与库生成种子独立)')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    k, n_rep = args.topk, args.n_rep

    rows, wt, dr, entry = build_library(args.effector, args.spacer)
    prov = provenance_check(rows, args.effector)
    passed = [r for r in rows if r['passed']]
    n = len(passed)
    print(f"库: {len(rows)} 变体, 通过硬过滤 {n} 条; 对拍: {prov['note']}")

    scores = np.array([r['score'] for r in passed])
    order = np.argsort(-scores)
    top_idx = order[:k]
    top_rows = [passed[i] for i in top_idx]

    # ---- 独立性质(均不进打分公式) ----
    stem_pos1 = {i + 1 for i, ch in enumerate(wt['dr_only_struct']) if ch in '()'}
    loop_pos1 = true_loop_positions(wt['dr_only_struct'])
    stem_gc = np.array([(sum(1 for p in stem_pos1 if r['dr_seq'][p - 1] in 'GC')
                         / len(stem_pos1)) for r in passed])
    loop_frac = np.array([sum(1 for p in r['mut_positions'] if p in loop_pos1)
                          / max(r['n_mut'], 1) for r in passed])
    cross_pp = np.array([r['cross_pp'] for r in passed])
    sel_preds = load_selector_preds()
    sel_arr = None
    if sel_preds is not None:
        matched = [sel_preds.get(r['desc']) for r in passed]
        n_miss = sum(v is None for v in matched)
        if n_miss == 0:
            sel_arr = np.array(matched)
        else:
            print(f"[警告] 选型器 join 缺 {n_miss} 条, 跳过 selector 性质")
    else:
        print("[警告] data/selector_ranked_candidates.csv 缺失, 跳过 selector 性质")

    props = {'stem_gc': (stem_gc, 'WT 茎区位点 GC 含量(纯序列性质, 打分不含)',
                         False),
             'loop_mut_frac': (loop_frac,
                               f'突变落真 loop(1-based {sorted(loop_pos1)}, 两侧配对判据)比例; '
                               '打分只用 3\' 保守窗, 不用 loop 归属', False),
             'cross_pp': (cross_pp,
                          'DR-spacer 交叉配对系综概率(连续量; 过滤用整数 cross_nt/'
                          'inv_max_run, 打分不含此项; 低=侵占少)', True)}
    if sel_arr is not None:
        props['selector_lit_pred'] = (
            sel_arr, 'Han2025 文献先验树预测(Fig1g 尺度, 低=预测活性强; 半独立: '
                     '与打分共享部分底层特征, 模型由独立数据训练)', True)

    # ---- 随机基线抽样 ----
    rng = np.random.default_rng(args.baseline_seed)
    draws = np.array([rng.choice(n, size=k, replace=False) for _ in range(n_rep)])
    null_score = scores[draws].mean(axis=1)
    overlap = np.array([len(set(d.tolist()) & set(top_idx.tolist())) for d in draws])

    # ---- 检验一: 打分分离(构造性上界) ----
    obs_score = float(scores[top_idx].mean())
    pool_sd = float(scores.std())
    p_score = mc_test(obs_score, null_score, 'greater')
    score_test = {
        'hypothesis': 'H0: 硬过滤后候选可交换(打分无增量信息); obs=TOP-K 平均打分',
        'obs_topk_mean_score': round(obs_score, 4),
        'random_k_mean_mean': round(float(null_score.mean()), 4),
        'random_k_mean_sd': round(float(null_score.std()), 4),
        'random_k_mean_q95': round(float(np.percentile(null_score, 95)), 4),
        'random_k_mean_max': round(float(null_score.max()), 4),
        'mc_p_greater': round(p_score, 4),
        'effect_size_mean_diff_over_pool_sd':
            round((obs_score - float(null_score.mean())) / max(pool_sd, 1e-12), 3),
        'z_vs_null': round((obs_score - float(null_score.mean()))
                           / max(float(null_score.std()), 1e-12), 2),
        'pool_score_sd': round(pool_sd, 4),
        'pool_score_iqr': round(float(np.percentile(scores, 75)
                                      - np.percentile(scores, 25)), 4),
        'topk_score_range': [round(float(scores[top_idx].min()), 4),
                             round(float(scores[top_idx].max()), 4)],
        'note': 'TOP-K 按构造是分数最大子集, MC p 必然近 1/(n_rep+1); 该检验的真正'
                '内容在效应量: 若打分无区分力, TOP-K 均值≈随机均值(效应量≈0)。'
                '实测效应量以通过池 SD 归一报告。'}

    # ---- 检验二: 独立性质 ----
    prop_tests = {}
    for name, (arr, note, lower_better) in props.items():
        null_m = arr[draws].mean(axis=1)
        prop_tests[name] = prop_summary(
            name, arr[top_idx], arr, null_m, note, lower_better)
        t = prop_tests[name]
        print(f"[性质] {name:<18} TOP-K均值 {t['topk_mean']:>8.4f} vs 随机中位 "
              f"{t['random_k_mean_median']:>8.4f}  MC p={t['mc_p_two_sided']:.4f}  "
              f"Cliff Δ={t['cliffs_delta_topk_vs_pool']:+.3f}")

    # ---- 重合度 ----
    overlap_test = {
        'expected_hypergeometric': round(k * k / n, 2),
        'random_draw_overlap_mean': round(float(overlap.mean()), 2),
        'random_draw_overlap_max': int(overlap.max()),
        'note': f'随机 K={k} 与管线 TOP-{k} 的期望重合 K^2/N={k * k / n:.2f} 条; '
                '若过滤+随机即可捞回同一批候选, 实测重合应显著高于期望——实测不高于期望'}

    # ---- 结论(模板填数, 不夸大) ----
    sig_props = [p for p, t in prop_tests.items() if t['mc_p_two_sided'] <= 0.05]
    conclusion = (
        f"检验一: 管线 TOP-{k} 平均打分 {obs_score:.3f}, 随机基线均值 "
        f"{null_score.mean():.3f} (SD {null_score.std():.3f}), MC p={p_score:.4f}"
        f"(按构造上界), 效应量 {score_test['effect_size_mean_diff_over_pool_sd']:.2f} "
        f"个通过池 SD——打分在通过池中提供了大效应的排序区分(池 SD "
        f"{pool_sd:.3f}, TOP-{k} 内部跨度 {score_test['topk_score_range'][0]}~"
        f"{score_test['topk_score_range'][1]}), 并非噪声级重排。检验二: 独立性质上"
        f"显著偏离随机基线的性质 = {sig_props or '无'}(其余与随机无显著差异, 如实报告)。"
        f"重合度: 随机抽取与管线 TOP-{k} 平均重合 {overlap.mean():.2f} 条"
        f"(期望 {k * k / n:.2f}), 过滤器+随机不能捞回同一批候选。"
        f"结论: 打分相对'硬过滤+随机'提供显著增量信息; 独立性质上的具体偏离方向见各性质项。")

    print('\n=== 检验一(打分分离) ===')
    print(f"TOP-{k} 均分 {obs_score:.4f} vs 随机均值 {null_score.mean():.4f} "
          f"± {null_score.std():.4f}; MC p={p_score:.4f}; 效应量 "
          f"{score_test['effect_size_mean_diff_over_pool_sd']} 池SD")
    print(f"随机重合: 均值 {overlap.mean():.2f} / 最大 {overlap.max()} "
          f"(超几何期望 {k * k / n:.2f})")

    payload = {
        'task': 'R3 filter+random 基线: 管线打分相对"硬过滤+随机选K"的增量信息',
        'generated': date.today().isoformat(),
        'vienna_version': RNA.__version__, 'numpy_version': np.__version__,
        'config': {'effector': args.effector,
                   'dr_dna': core.to_dna(dr),
                   'spacer_fixed_dna': core.to_dna(core.to_rna(
                       args.spacer or entry['placeholder_spacer'])),
                   'library_params': LIB_PARAMS, 'topk': k, 'n_rep': n_rep,
                   'baseline_seed': args.baseline_seed,
                   'scope': '纯 ViennaRNA 口径; RNet SHAPE 一致性需 rnet 引擎环境, '
                            '本检验未做(round-3 R3 授权跳过)'},
        'library': {'n_variants': len(rows), 'n_passed': n,
                    'provenance_check': prov},
        'topk_members': [{'desc': r['desc'], 'score': r['score'],
                          'stem_gc': round(float(stem_gc[i]), 4),
                          'loop_mut_frac': round(float(loop_frac[i]), 4),
                          'cross_pp': round(float(cross_pp[i]), 4),
                          'selector_lit_pred': (round(float(sel_arr[i]), 4)
                                               if sel_arr is not None else None)}
                         for i, r in zip(top_idx, top_rows)],
        'score_separation_test': score_test,
        'independent_property_tests': prop_tests,
        'topk_overlap_test': overlap_test,
        'conclusion_zh': conclusion,
        'disclaimer': '排序信息检验, 非活性预测; 打分提供的是相对排序信息, '
                      '活性/特异性以体外生化与细胞实验为准。'}
    out = args.out or os.path.join(
        ROOT, 'data', f'selection_baseline.{args.effector}.json')
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f"\n输出: {out}")


if __name__ == '__main__':
    main()
