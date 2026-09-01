"""v2 奖赏制权重灵敏度分析（Monte Carlo / 拉丁超立方）

回答"奖赏权重是不是拍脑袋": 在默认权重 ±50% 范围内拉丁超立方采样 N 组,
对 v2 候选库重打分(指标已在 variants.csv 中, 无需重新折叠), 统计:
  - 默认 TOP-K 候选在扰动下留在 TOP-K 的频率(稳健性)
  - 全排名与默认排名的 Spearman 相关分布
  - 稳健核心集(>= --robust-frac 频率留在 TOP-K 的变体)

输出: <run-prefix>.sensitivity.json / .png(柱状图) / 控制台摘要。

用法:
  python weight_sensitivity.py --run-prefix ../data/tp53_r248q_zengdr_v2
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

WEIGHTS = ['w_seed', 'w_fold', 'w_stab', 'w_bp', 'w_ddg',
           'w_contact', 'w_ens', 'w_hbond', 'w_cons3']
DEFAULTS = dict(w_seed=0.0, w_fold=0.2, w_stab=0.3, w_bp=0.3,
                w_ddg=0.1, w_contact=1.0, w_ens=0.05, w_hbond=0.5, w_cons3=0.3)
# w_seed 默认关闭(与 spacer_up 共线 0.835 超预注册线, 见主管线帮助); 灵敏度分析中
# 按 [0, 0.4] 均匀采样检验其潜在影响, 其余权重在默认值 ±range 内拉丁超立方采样


def rescore(r, w):
    """由 variants.csv 的存储指标重构主线 v1.4 score（与 score_variant 同式）。"""
    return (w['w_seed'] * float(r['d_seed'])
            + w['w_fold'] * float(r['dp_fold'])
            + w['w_stab'] * max(-float(r['ddG']), 0.0)
            - w['w_bp'] * float(r['bp_dist'])
            - w['w_ddg'] * max(float(r['ddG']), 0.0)
            - w['w_contact'] * float(r['contact_sum'])
            - w['w_ens'] * max(float(r['d_ens']), 0.0)
            + w['w_hbond'] * float(r['hbond_net'])
            - w['w_cons3'] * float(r['cons3_frac']))


def lhs(n, dim, rng, lo=0.5, hi=1.5):
    """拉丁超立方采样, 每维均匀分层于 [lo, hi]（相对默认权重）。"""
    u = np.zeros((n, dim))
    for d in range(dim):
        perm = rng.permutation(n)
        u[:, d] = (perm + rng.random(n)) / n
    return lo + (hi - lo) * u


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--run-prefix', required=True)
    ap.add_argument('--n-samples', type=int, default=256)
    ap.add_argument('--topk', type=int, default=12)
    ap.add_argument('--range', type=float, default=0.5, help='权重扰动范围(±比例)')
    ap.add_argument('--robust-frac', type=float, default=0.8,
                    help='稳健核心集的 TOP-K 保留频率阈值')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.run_prefix + '.variants.csv', encoding='utf-8')))
    variants = rows[1:]
    passed = [r for r in variants if r['passed'] == 'True']
    print(f"通过硬过滤变体 {len(passed)} 条, 采样 {args.n_samples} 组权重 "
          f"(非零默认 ±{int(args.range * 100)}% 拉丁超立方; w_seed 默认0按[0,0.4]采样)")

    rng = np.random.default_rng(args.seed)
    base_w = np.array([DEFAULTS[k] for k in WEIGHTS])
    samples = lhs(args.n_samples, len(WEIGHTS), rng,
                  lo=1 - args.range, hi=1 + args.range) * base_w
    seed_idx = WEIGHTS.index('w_seed')
    samples[:, seed_idx] = rng.random(args.n_samples) * 0.4  # 默认0的权重: 加性采样

    base_scores = np.array([rescore(r, DEFAULTS) for r in passed])
    base_top = set(np.argsort(-base_scores)[:args.topk].tolist())
    in_count = np.zeros(len(passed))
    rhos = []
    for s in samples:
        w = dict(zip(WEIGHTS, s))
        sc = np.array([rescore(r, w) for r in passed])
        top = set(np.argsort(-sc)[:args.topk].tolist())
        for k in top:
            in_count[k] += 1
        rhos.append(spearman(base_scores, sc))
    freq = in_count / args.n_samples

    print(f"Spearman 排名相关: 中位 {np.median(rhos):.3f} "
          f"[{np.min(rhos):.3f}, {np.max(rhos):.3f}]")
    robust = [i for i in range(len(passed)) if freq[i] >= args.robust_frac]
    print(f"稳健核心集(>= {args.robust_frac:.0%} 频率留在 TOP-{args.topk}): {len(robust)} 条")
    print(f"\n{'desc':<26}{'freq':>6}{'base_score':>11}  默认TOP{args.topk}?  advantage")
    for i in sorted(range(len(passed)), key=lambda i: -freq[i])[:args.topk + 4]:
        r = passed[i]
        print(f"{r['desc']:<26}{freq[i]:>6.2f}{base_scores[i]:>11.3f}  "
              f"{'Y' if i in base_top else '-':>5}    {r['advantage']}")

    payload = {
        'method': f'拉丁超立方 {args.n_samples} 组, 默认权重 ±{args.range:.0%}',
        'weights': DEFAULTS, 'topk': args.topk,
        'spearman_median': round(float(np.median(rhos)), 3),
        'spearman_min': round(float(np.min(rhos)), 3),
        'robust_frac_threshold': args.robust_frac,
        'robust_core': [{'desc': passed[i]['desc'], 'dr_seq': passed[i]['dr_seq'],
                         'freq': round(float(freq[i]), 3),
                         'base_score': round(float(base_scores[i]), 4)}
                        for i in sorted(robust, key=lambda i: -freq[i])],
        'topk_stability': [{'desc': passed[i]['desc'], 'freq': round(float(freq[i]), 3)}
                           for i in sorted(base_top, key=lambda i: -base_scores[i])],
    }
    out_json = args.run_prefix + '.sensitivity.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        idx = sorted(range(len(passed)), key=lambda i: -freq[i])[:20]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh([passed[i]['desc'] for i in idx][::-1],
                [freq[i] for i in idx][::-1], color='#4878CF')
        ax.axvline(args.robust_frac, color='red', ls='--',
                   label=f'robust threshold {args.robust_frac}')
        ax.set_xlabel(f'Frequency in TOP-{args.topk} under weight perturbation '
                      f'(n={args.n_samples}, ±{int(args.range * 100)}%)')
        ax.set_title('Weight sensitivity of v2 directed scoring')
        ax.legend()
        fig.tight_layout()
        out_png = args.run_prefix + '.sensitivity.png'
        fig.savefig(out_png, dpi=150)
        print(f"\n输出: {out_json} / {out_png}")
    except ImportError:
        print(f"\n输出: {out_json} (matplotlib 不可用, 跳过图)")


if __name__ == '__main__':
    main()
