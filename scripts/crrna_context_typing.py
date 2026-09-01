"""crRNA 骨架的 spacer 上下文分型（go/no-go 试验）

假说: spacer(A 部分)序列千差万别, 与 DR(B 部分)同处一条 RNA, 通过全长折叠
(交界配对/自互补竞争)影响 DR 的最优选择 —— 因此骨架不应只有一种, 而应按
spacer 上下文分 3-4 型, 每型各有最优 DR。

本模块两步(go/no-go):
  1) 对一组 spacer(Zeng 2026 MOESM3 全部 86 条)计算上下文特征, KMeans 聚类
     (k=2..6 按轮廓系数选型, 手写实现零额外依赖);
  2) 每型取最靠型心的代表, 子进程重跑 crrna_scaffold_design(参数一致),
     对比各型 TOP DR 的重合度 —— 型间 TOP 显著不同 = 分型有依据(阳性);
     高度重合 = spacer 上下文效应弱于 DR 内在偏好(阴性, 同样如实报告)。

特征(全部来自全长折叠, 复用主管线函数):
  GC%、spacer 单独折叠 MFE(自互补强度)、全长 spacer 平均未配对概率、
  DR 3'端-spacer 5'端交界配对对数、spacer 区内部配对对数。

用法:
  python crrna_context_typing.py --spacers ../data/context_typing/spacers.txt \
      --effector cas12a2_zeng2026 --out-prefix ../data/context_typing/typing
输出: typing.clusters.json + typing.report.txt + 各代表 <prefix>.<name>.* 管线产出
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_scaffold_design as core  # noqa: E402
from scaffold_registry import get_scaffold  # noqa: E402


def spacer_features(dr, name, seq):
    """单条 spacer 的上下文特征向量(5 维)。"""
    full = dr + seq
    ss_full, _ = core.fold(full)
    _, spacer_up = core.pf_stats(full, len(dr))
    _, self_mfe = core.fold(seq)
    n = len(full)
    stack, pairs = [], []
    for i, ch in enumerate(ss_full):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            pairs.append((stack.pop(), i))
    dl = len(dr)
    junction = sum(1 for i, j in pairs if i < dl <= j)
    internal = sum(1 for i, j in pairs if i >= dl)
    gc = (seq.count('G') + seq.count('C')) / len(seq)
    return [gc, self_mfe, spacer_up, junction, internal]


def kmeans(X, k, seeds):
    """手写 KMeans(多起点), 返回 (labels, centers, inertia)。"""
    best = None
    for sd in seeds:
        rng = np.random.default_rng(sd)
        C = X[rng.choice(len(X), k, replace=False)].copy()
        for _ in range(100):
            d = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
            lab = d.argmin(1)
            newC = np.array([X[lab == c].mean(0) if (lab == c).any() else C[c]
                             for c in range(k)])
            if np.allclose(newC, C):
                break
            C = newC
        inertia = ((X - C[lab]) ** 2).sum()
        if best is None or inertia < best[2]:
            best = (lab, C, inertia)
    return best


def silhouette(X, lab, k):
    """手写轮廓系数(平均)。"""
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    s = []
    for i in range(len(X)):
        own = lab == lab[i]
        if own.sum() <= 1:
            s.append(0.0)
            continue
        a = D[i, own].sum() / (own.sum() - 1)
        b = min(D[i, lab == c].mean() for c in range(k) if c != lab[i] and (lab == c).any())
        s.append((b - a) / max(a, b) if b > 0 else 0.0)
    return float(np.mean(s))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacers', required=True, help='TSV: name\\tRNA spacer 序列')
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--max-k', type=int, default=6)
    ap.add_argument('--reps-per-type', type=int, default=3)
    ap.add_argument('--topk', type=int, default=8)
    ap.add_argument('--python', default=sys.executable, help='跑主管线用的 python')
    ap.add_argument('--skip-pipeline', action='store_true', help='只聚类不重跑管线')
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    dr = core.to_rna(get_scaffold(args.effector))
    names, seqs = [], []
    with open(args.spacers, encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                n, s = line.rstrip('\n').split('\t')
                names.append(n)
                seqs.append(core.to_rna(s))
    print(f'[typing] {len(seqs)} 条 spacer, DR={core.to_dna(dr)}({len(dr)}nt)')

    F = np.array([spacer_features(dr, n, s) for n, s in zip(names, seqs)])
    Z = (F - F.mean(0)) / F.std(0)
    feat_names = ['GC', 'self_mfe', 'spacer_up', 'junction_pairs', 'internal_pairs']

    cand = []
    for k in range(2, min(args.max_k, len(Z)) + 1):
        lab, C, inertia = kmeans(Z, k, seeds=range(10))
        sil = silhouette(Z, lab, k)
        cand.append({'k': k, 'silhouette': round(sil, 3), 'inertia': round(inertia, 2),
                     'sizes': [int((lab == c).sum()) for c in range(k)], 'labels': lab.tolist()})
        print(f'  k={k}  silhouette={sil:.3f}  sizes={[int((lab==c).sum()) for c in range(k)]}')
    best = max(cand, key=lambda c: c['silhouette'])
    k, lab = best['k'], np.array(best['labels'])
    print(f"[typing] 选定 k={k} (轮廓系数 {best['silhouette']})")

    types = {}
    for c in range(k):
        idx = np.where(lab == c)[0]
        d = np.sqrt(((Z[idx] - Z[idx].mean(0)) ** 2).sum(1))
        order = idx[np.argsort(d)]
        types[str(c)] = {
            'n': int(len(idx)),
            'feature_means': {fn: round(float(F[idx, i].mean()), 3)
                              for i, fn in enumerate(feat_names)},
            'members': [names[i] for i in idx],
            'representatives': [{'name': names[i], 'spacer': core.to_dna(seqs[i]),
                                 'dist_to_centroid': round(float(d[np.where(order == i)][0]) if i in order else 0, 3)}
                                for i in order[:args.reps_per_type]],
        }
        print(f"  型{c}: n={len(idx)} 特征均值={types[str(c)]['feature_means']}")

    report = {'effector': args.effector, 'dr_dna': core.to_dna(dr), 'n_spacers': len(seqs),
              'features': feat_names, 'k_scan': [{kk: v for kk, v in c.items() if kk != 'labels'} for c in cand],
              'chosen_k': k, 'types': types,
              'note': 'go/no-go: 各型代表重跑管线后比对 TOP DR 重合度; 高重合=上下文效应弱(阴性)'}

    tops = {}
    if not args.skip_pipeline:
        design = os.path.join(HERE, 'crrna_scaffold_design.py')
        for c in range(k):
            for rep in types[str(c)]['representatives']:
                prefix = f"{args.out_prefix}.t{c}_{rep['name'].replace(' ', '_')}"
                cmd = [args.python, design, '--effector', args.effector,
                       '--spacer', rep['spacer'], '--topk', str(args.topk),
                       '--use-covariation', '--out-prefix', prefix]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    print(f"[pipeline] {rep['name']} 失败: {r.stderr[-200:]}")
                    continue
                with open(prefix + '.top.json', encoding='utf-8') as fh:
                    top = json.load(fh)
                key = f't{c}:{rep["name"]}'
                tops[key] = [row['desc'] for row in top['top'][1:args.topk + 1]]
                print(f"[pipeline] 型{c} {rep['name']}: TOP{args.topk} = {tops[key]}")

        # 型间 TOP DR 重合度(两两 Jaccard + 全型公共子集)
        keys = list(tops)
        overlap = {}
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = set(tops[keys[i]]), set(tops[keys[j]])
                inter, union = len(a & b), len(a | b) if a | b else 1
                overlap[f'{keys[i]} vs {keys[j]}'] = f'{inter}/{union} J={inter / union:.2f}'
        per_type_sets = {}
        for key in keys:
            per_type_sets.setdefault(key.split(':')[0], set()).update(tops[key])
        tkeys = list(per_type_sets)
        common = set.intersection(*per_type_sets.values()) if len(per_type_sets) > 1 else set()
        report['top_overlaps'] = overlap
        report['top_per_type_union'] = {t: sorted(per_type_sets[t]) for t in tkeys}
        report['top_common_to_all_types'] = sorted(common)
        report['go_no_go'] = ('阳性(型间 TOP 差异显著): 各型并集不同且公共集小'
                              if len(common) <= args.topk // 2 and len(per_type_sets) > 1
                              else '阴性(上下文效应弱): 各型 TOP 高度重合, 通用骨架即可')
        print('[typing] 型间 TOP 公共集:', sorted(common))
        print('[typing] go/no-go:', report['go_no_go'])

    out_json = args.out_prefix + '.clusters.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(f'[typing] 结果 -> {out_json}')


if __name__ == '__main__':
    main()
