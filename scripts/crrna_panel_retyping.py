# -*- coding: utf-8 -*-
"""面板 spacer 上下文重分型(2026-09-16)。

背景: 用户指出——骨架分型学说(§A/§A2/§A3)建立在 DR×spacer 互作特征上,
"spacer 上下文分 3-4 型, 每型各有最优 DR"; 换 spacer(如 JD12-crRNA-2)
后, 已选骨架变体(A1C / A8C+U15G)必须在新上下文里重新分型核验,
不能沿用旧 spacer 的结论。

流程(全部确定性, 复用 §A3 口径):
1) 重建 v47 扩库队列(文献 85 + 确定性抽样, scripts/crrna_typing_v47_expansion
   同一套规则/常数), 计算 10 维互作特征(§A3: crrna_context_typing.
   spacer_features_v2), k=3(range(50) 种子) 重聚类, 对齐 §A3 记录的
   cluster_sizes=[137,56,119] 自检。
2) 三条面板 spacer 就近分型(同标准化):
     canonical_R248Q(24nt, GTTCATGCCGCCCATGCAGGAACT)
     canonical_R273H(24nt, CACCTCAAAGCTGTTCCGTCCCAG)
     JD12_sp2(23nt,    AGGACAGGCACAAACATGCACCT)
   注: 队列 spacer 为 23nt, 24nt 面板 spacer 属轻度外推, 如实记录。
3) 对分得型别的 §A3 代表 spacer 重跑主管线(topk=16, --use-covariation,
   与 §A2/A3 同口径), 取各代表 TOP 骨架集, 核验 A1C / A8C+U15G 的
   在/不在与排名。

输出: data/panel_retyping.json
运行: PYTHONUTF8=1 python scripts/crrna_panel_retyping.py
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_context_typing as typing_mod  # noqa: E402
import crrna_scaffold_design as core       # noqa: E402
import crrna_typing_v47_expansion as v47   # noqa: E402
import crrna_typing_robustness as rb       # noqa: E402
from scaffold_registry import get_scaffold  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'context_typing')
A3_REF = os.path.join(ROOT, 'data', 'context_typing', 'typing_v47_expansion.json')
OUT = os.path.join(ROOT, 'data', 'panel_retyping.json')

PANEL_SPACERS = {
    'canonical_R248Q': 'GTTCATGCCGCCCATGCAGGAACT',
    'canonical_R273H': 'CACCTCAAAGCTGTTCCGTCCCAG',
    'JD12_sp2': 'AGGACAGGCACAAACATGCACCT',
}
CHECK_SCAFFOLDS = ['A1C', 'A8C+U15G']
TOPK = 16


def main():
    dr = core.to_rna(get_scaffold(rb.EFFECTOR))

    lit = [l.strip().split('\t') for l in open(v47.LIT_SPACERS,
                                               encoding='utf-8') if l.strip()]
    new = v47.sample_spacers()
    names = [n for n, _ in lit] + [n for n, _, _ in new]
    seqs = [core.to_rna(s) for _, s in lit] + \
           [core.to_rna(d) for _, d, _ in new]
    print('[retype] 队列 %d 条(文献 %d + v47 %d)'
          % (len(names), len(lit), len(new)))

    feat_names = ['GC', 'self_mfe', 'spacer_up', 'junction_pairs',
                  'internal_pairs', 'seed5_up', 'cross_stem', 'cross_loop',
                  'junction_maxrun', 'coupling_dG']
    F = np.array([typing_mod.spacer_features_v2(dr, n, s)
                  for n, s in zip(names, seqs)], dtype=float)
    Z = (F - F.mean(0)) / F.std(0)
    lab, centers, _ = typing_mod.kmeans(Z, v47.K_FIXED, v47.SEEDS)
    sizes = [int((lab == c).sum()) for c in range(v47.K_FIXED)]

    ref = json.load(open(A3_REF, encoding='utf-8'))
    ref_sizes = ref['a3_enriched']['clustering']['cluster_sizes']
    print('[retype] k=3 复算簇大小 %s | §A3 记录 %s' % (sizes, ref_sizes))

    # 面板 spacer 分型(同标准化)
    mu, sd = F.mean(0), F.std(0)
    assigns = {}
    for pname, pdna in PANEL_SPACERS.items():
        f = np.array(typing_mod.spacer_features_v2(
            dr, pname, core.to_rna(pdna)), dtype=float)
        z = (f - mu) / sd
        d = ((centers - z) ** 2).sum(1)
        t = int(d.argmin())
        assigns[pname] = {
            'spacer_dna': pdna, 'spacer_len': len(pdna),
            'features': dict(zip(feat_names, [round(float(x), 3) for x in f])),
            'type': t, 'dist_to_centers': [round(float(x), 3) for x in d],
            'dist_margin_2nd': round(float(sorted(d)[1] - sorted(d)[0]), 3)}
        print('[retype] %-16s -> 型%d (次近距离差 %.3f)'
              % (pname, t, assigns[pname]['dist_margin_2nd']))

    # 各型 §A3 代表重跑管线取 TOP
    reps_ref = ref['a3_enriched']['representatives']
    seq_of = dict(zip(names, seqs))
    tops = {}
    for t in sorted(set(a['type'] for a in assigns.values())):
        tops[str(t)] = {}
        for rn in reps_ref[str(t)]:
            if rn not in seq_of:
                tops[str(t)][rn] = {'error': '代表不在复算队列'}
                continue
            prefix = os.path.join(DATA, 'retype_t%d_%s'
                                  % (t, rn.replace(' ', '_').replace('/', '_')))
            cmd = [sys.executable, os.path.join(HERE, 'crrna_scaffold_design.py'),
                   '--effector', rb.EFFECTOR, '--spacer',
                   core.to_dna(seq_of[rn]), '--topk', str(TOPK),
                   '--use-covariation', '--out-prefix', prefix]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                tops[str(t)][rn] = {'error': r.stderr[-160:]}
                continue
            top = json.load(open(prefix + '.top.json', encoding='utf-8'))
            ranks = {x['desc']: i + 1 for i, x in enumerate(
                top.get('top', top.get('designs', [])))}
            tops[str(t)][rn] = {
                'top_names': list(ranks)[:TOPK],
                'check_ranks': {s: ranks.get(s) for s in CHECK_SCAFFOLDS}}
            print('[retype] 型%d 代表 %-28s TOP16 含: %s'
                  % (t, rn,
                     {s: ranks.get(s) for s in CHECK_SCAFFOLDS}))

    # 判读(程序化)
    verdict = {}
    for pname, a in assigns.items():
        t = str(a['type'])
        rep_res = [v for v in tops.get(t, {}).values() if 'check_ranks' in v]
        if not rep_res:
            verdict[pname] = '型%d 无可用代表管线结果' % a['type']
            continue
        hits = {s: [r['check_ranks'][s] for r in rep_res
                    if r['check_ranks'][s] is not None] for s in CHECK_SCAFFOLDS}
        verdict[pname] = {
            'type': a['type'],
            'A1C': '代表 %d/%d 进 TOP16(最好第%s名)' % (
                len(hits['A1C']), len(rep_res),
                min(hits['A1C']) if hits['A1C'] else '-') if hits['A1C'] else
                '代表 0/%d 进 TOP16' % len(rep_res),
            'A8C+U15G': '代表 %d/%d 进 TOP16(最好第%s名)' % (
                len(hits['A8C+U15G']), len(rep_res),
                min(hits['A8C+U15G']) if hits['A8C+U15G'] else '-') if
                hits['A8C+U15G'] else '代表 0/%d 进 TOP16' % len(rep_res)}

    payload = {
        'generated_by': 'scripts/crrna_panel_retyping.py',
        'doctrine': '§A3 骨架分型(spacer 上下文 10 维互作特征, k=3, '
                    '各型各有最优 DR); 换 spacer 必须重分型核验',
        'cohort_n': len(names), 'cohort_sizes_recomputed': sizes,
        'a3_sizes_recorded': ref_sizes,
        'cohort_match': sizes == ref_sizes,
        'caveat': '队列 spacer 为 23nt; canonical 面板 spacer 为 24nt, '
                  '属轻度长度外推; JD12_sp2 为 23nt 同长',
        'panel_assignments': assigns, 'type_reps_tops': tops,
        'verdict': verdict}
    json.dump(payload, open(OUT, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print('->', OUT)
    for k, v in verdict.items():
        print('[verdict]', k, '|', v)


if __name__ == '__main__':
    main()
