# -*- coding: utf-8 -*-
"""骨架分型 go/no-go 试验的四项事后稳健性补救(post-hoc, 非预登记)。

对 crrna_context_typing.py 的擦线阳性(公共集 3 条, 阈 ≤4)做敏感性分析:
  补救1 置换检验: 固定 9 个代表 TOP 列表(k=8 切片, 来自 topk=16 重跑), 对
         9->3 组(3,3,3)的全部无序划分(280 种, 统计量组间对称故等价于 1680
         种标号分配)精确枚举公共交集大小, 得零分布与精确 p 值;
  补救2 敏感性双联: (a) topk=6/8/10/16 四档判据翻动; (b) 剔除全部含
         position-1 突变的变体后按 rank 重排 TOP-8 重算判据;
  补救3 聚类稳健性谱: (a) 手写 KMeans / Agglomerative(ward) / GaussianMixture
         三方法 k=2..6 轮廓系数谱; (b) 预登记方法(手写 KMeans)k=2..6 判据重算
         (新代表跑主管线 topk=16, 按 spacer 名缓存复用);
  补救4 n=9 小簇(型1)来源审查: 来源基因汇总 / 两两序列同一性 / 5 维特征对比 /
         可扩充候选。

口径声明: 本脚本为 post-hoc 敏感性分析, 非预登记; 所有结论与
docs/preregistration.md §A 判据并列阅读, 不替代原判据。
全部判读文字由数值按显式阈值规则生成(同 crrna_chai_matrix_summary.py 惯例),
程序内不手写结论方向。

用法:
  PYTHONUTF8=1 python scripts/crrna_typing_robustness.py
输出: data/context_typing/typing.robustness.json + 控制台报告
特征缓存: data/context_typing/typing.robustness.features.json
管线缓存: data/context_typing/robustness_runs/rb_<name>.{top.json,variants.csv,...}
"""
import itertools
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_context_typing as typing_mod  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402
from scaffold_registry import get_scaffold  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'context_typing')
RUNS = os.path.join(DATA, 'robustness_runs')
FEAT_CACHE = os.path.join(DATA, 'typing.robustness.features.json')
CLUSTERS = os.path.join(DATA, 'typing.clusters.json')
SPACERS = os.path.join(DATA, 'spacers.txt')
DESIGN = os.path.join(HERE, 'crrna_scaffold_design.py')

FEAT_NAMES = ['GC', 'self_mfe', 'spacer_up', 'junction_pairs', 'internal_pairs']
CRIT_TOPK = 8            # 预登记判据的 topk
PERM_TOPK_RUN = 16       # 统一重跑深度, k=8/6/10/16 均从此切片
ALPHA = 0.05             # 置换检验显著性阈
SIL_MARGIN = 0.01        # 轮廓系数"明显最优"所需领先幅度
IDENTITY_HIGH = 0.80     # 簇内平均同一性高于此值判"高度同源"
JUNCTION_EXPAND_MIN = 2  # 可扩充候选的 junction_pairs 下限
MUT_POS_RE = re.compile(r'[ACGU](\d+)[ACGU]')

SCOPE = ('post-hoc 敏感性分析, 非预登记; 与 docs/preregistration.md §A 判据'
         '并列阅读, 不替代原判据; 判读文字全部由数值按显式阈值规则生成')


def load_spacers():
    names, seqs = [], []
    with open(SPACERS, encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                n, s = line.strip().split('\t')
                names.append(n)
                seqs.append(core.to_rna(s))
    return names, seqs


def feature_matrix(dr, names, seqs):
    """86 条特征矩阵, 带 JSON 缓存(spacer 序列核对, 不一致则重算)。"""
    cache = {}
    if os.path.exists(FEAT_CACHE):
        with open(FEAT_CACHE, encoding='utf-8') as fh:
            cache = json.load(fh).get('features', {})
    F, dirty = [], False
    for n, s in zip(names, seqs):
        dna = core.to_dna(s)
        rec = cache.get(n)
        if rec is None or rec.get('spacer') != dna:
            rec = {'spacer': dna,
                   'features': typing_mod.spacer_features(dr, n, s)}
            cache[n] = rec
            dirty = True
        F.append(rec['features'])
    if dirty:
        with open(FEAT_CACHE, 'w', encoding='utf-8') as fh:
            json.dump({'dr_dna': core.to_dna(dr), 'feature_names': FEAT_NAMES,
                       'features': cache}, fh, ensure_ascii=False, indent=1)
    return np.array(F, dtype=float)


def san(name):
    return name.replace(' ', '_')


def run_pipeline(name, spacer_dna, topk=PERM_TOPK_RUN):
    """按 spacer 名缓存的主管线重跑(topk=16), 返回 out-prefix。"""
    os.makedirs(RUNS, exist_ok=True)
    prefix = os.path.join(RUNS, 'rb_' + san(name))
    top_path = prefix + '.top.json'
    if os.path.exists(top_path):
        with open(top_path, encoding='utf-8') as fh:
            if json.load(fh).get('spacer_fixed_dna') == spacer_dna:
                return prefix, False
    t0 = time.time()
    cmd = [sys.executable, DESIGN, '--effector', EFFECTOR,
           '--spacer', spacer_dna, '--topk', str(topk),
           '--use-covariation', '--out-prefix', prefix]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
    if r.returncode != 0:
        raise RuntimeError('管线失败 %s: %s' % (name, r.stderr[-300:]))
    return prefix, True


def top_descs(prefix, k):
    """top.json 的 TOP-k desc 列表(跳过第 0 行 WT)。"""
    with open(prefix + '.top.json', encoding='utf-8') as fh:
        top = json.load(fh)['top']
    return [row['desc'] for row in top[1:k + 1]]


def top_descs_no_pos1(prefix, k):
    """剔除全部含 position-1 突变的变体后, 按 rank 顺序取 TOP-k。"""
    import csv
    with open(prefix + '.variants.csv', newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    keep = []
    for r in sorted(rows, key=lambda r: int(r['rank'])):
        if r['desc'] == 'WT' or r['passed'] != 'True':
            continue
        if 1 in [int(p) for p in MUT_POS_RE.findall(r['desc'])]:
            continue
        keep.append(r['desc'])
        if len(keep) >= k:
            break
    return keep


def criterion(lists_by_type, topk):
    """预登记判据: 各型并集的全型公共子集大小 <= topk//2 判阳性。"""
    unions = {t: set().union(*lists) for t, lists in lists_by_type.items()}
    common = set.intersection(*unions.values()) if len(unions) > 1 else set()
    return {'common': sorted(common), 'n_common': len(common),
            'threshold': topk // 2, 'pass': len(common) <= topk // 2}


def reps_from_clusters(lab, Z, names, seqs, k, reps_per_type=3):
    """每型取最靠型心的 reps_per_type 个代表(不足取全部)。"""
    out = {}
    for c in range(k):
        idx = np.where(lab == c)[0]
        d = np.sqrt(((Z[idx] - Z[idx].mean(0)) ** 2).sum(1))
        order = idx[np.argsort(d)]
        out[str(c)] = [(names[i], core.to_dna(seqs[i])) for i in order[:reps_per_type]]
    return out


def levenshtein(a, b):
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[-1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def seq_identity(a, b):
    return 1.0 - levenshtein(a, b) / max(len(a), len(b))


# ---------- 判读规则(全部由数值生成, 不手写结论方向) ----------

def read_permutation(obs, null, p_obs, p_crit):
    med = float(np.median(null))
    q = {str(qq): round(float(np.percentile(null, qq)), 2)
         for qq in (5, 25, 75, 95)}
    if p_obs >= ALPHA:
        verdict = ('随机分组下公共集<=%d 是常见事件(p=%.4f>=%.2f), '
                   '判据通过率在随机分组下也常见(%.1f%%), 分型判据无判别力'
                   % (obs, p_obs, ALPHA, 100 * p_crit))
    else:
        verdict = ('观测公共集=%d 在零分布中罕见(p=%.4f<%.2f), '
                   '随机分组下判据通过率仅 %.1f%%, 判据具判别力'
                   % (obs, p_obs, ALPHA, 100 * p_crit))
    return {'observed': obs, 'null_median': med, 'null_quantiles': q,
            'p_observed': round(p_obs, 4), 'criterion_pass_rate_null': round(p_crit, 4),
            'reading': verdict}


def read_topk_sensitivity(per_k):
    passes = [v['pass'] for v in per_k.values()]
    commons = {k: v['n_common'] for k, v in per_k.items()}
    if len(set(passes)) == 1:
        verdict = ('判据在 topk=6/8/10/16 四档下结论一致(%s), 公共集大小 %s, '
                   '不随 topk 翻动' % ('全部通过' if passes[0] else '全部不通过',
                                       commons))
    else:
        flips = [k for k, v in per_k.items() if v['pass']]
        verdict = ('判据随 topk 翻动: 通过档 %s, 公共集大小 %s —— '
                   '擦线结论对 topk 选择敏感' % (flips, commons))
    return {'per_topk': per_k, 'reading': verdict}


def read_de_pos1(res, removed_counts):
    if not res['pass']:
        verdict = ('剔除 position-1 突变后公共集 %d 条 > 阈 %d, 判据不再通过 —— '
                   '分型判据的通过依赖 position-1 惰性位点'
                   % (res['n_common'], res['threshold']))
    else:
        verdict = ('剔除 position-1 突变后公共集 %d 条 <= 阈 %d, 判据仍通过, '
                   '不依赖 position-1 惰性位点' % (res['n_common'], res['threshold']))
    return {'criterion': res, 'removed_per_rep': removed_counts, 'reading': verdict}


def read_silhouette_spectrum(table):
    n_best3, margins = 0, {}
    for m, row in table.items():
        ks = sorted(row)
        best = max(ks, key=lambda k: row[k])
        second = max((row[k] for k in ks if k != best), default=None)
        margins[m] = {'best_k': best, 'silhouette': row[best],
                      'margin_over_runner_up': round(row[best] - second, 3)}
        if best == 3:
            n_best3 += 1
    n_m = len(table)
    if n_best3 == n_m and all(m['best_k'] == 3 and m['margin_over_runner_up'] >= SIL_MARGIN
                              for m in margins.values()):
        verdict = 'k=3 在全部 %d 种方法下均为轮廓系数最优(领先幅度均>=%.2f), 选型稳健' % (n_m, SIL_MARGIN)
    elif n_best3 == 0:
        verdict = ('k=3 在全部 %d 种方法下均非轮廓系数最优(%s), 原选型(k=3)不稳健'
                   % (n_m, {m: v['best_k'] for m, v in margins.items()}))
    else:
        verdict = ('k=3 仅在 %d/%d 种方法下最优(%s), 选型依赖聚类方法, 稳健性有限'
                   % (n_best3, n_m, {m: v['best_k'] for m, v in margins.items()}))
    return {'silhouette_by_method_k': table, 'best_per_method': margins,
            'reading': verdict}


def read_criterion_vs_k(per_k):
    passes = {k: v['pass'] for k, v in per_k.items()}
    if all(passes.values()):
        verdict = 'k=2..%d 各档判据全部通过, 分型阳性对 k 选择稳健' % max(map(int, per_k))
    elif not any(passes.values()):
        verdict = 'k=2..%d 各档判据全部不通过, 原 k=3 阳性不可复现' % max(map(int, per_k))
    else:
        verdict = ('判据随 k 翻动: 通过档 %s, 未通过档 %s —— k=3 阳性是特定选型下的结果'
                   % ([k for k, v in passes.items() if v],
                      [k for k, v in passes.items() if not v]))
    return {'per_k': per_k, 'reading': verdict}


def read_small_cluster(genes, ident_in, ident_bg, feat_diff, dominant, expand):
    dispersed = len(genes) >= 4
    homologous = ident_in['mean'] >= IDENTITY_HIGH
    if homologous or not dispersed:
        verdict = ('型1 来源集中在 %d 个基因(%s), 簇内平均同一性 %.3f%s —— '
                   '符合"同一来源/同源序列的批次效应"而非独立证据'
                   % (len(genes), sorted(genes), ident_in['mean'],
                      '(>=%.2f, 高度同源)' % IDENTITY_HIGH if homologous
                      else '(<%0.2f, 但来源过少)' % IDENTITY_HIGH))
    else:
        verdict = ('型1 横跨 %d 个独立基因(%s), 簇内平均同一性 %.3f(<%.2f, 非同源), '
                   '分离由 %s 维主导(标准化差 %.2f) —— 符合"来源分散但 '
                   'junction-pairing 特征一致的真实亚型"'
                   % (len(genes), sorted(genes), ident_in['mean'], IDENTITY_HIGH,
                      dominant, feat_diff[dominant]))
    verdict += '; 可主动扩充候选 %d 条' % len(expand)
    return {'n_members_genes': len(genes), 'genes': sorted(genes),
            'identity_within': ident_in, 'identity_background': ident_bg,
            'feature_std_diff_vs_rest': feat_diff, 'dominant_feature': dominant,
            'expansion_candidates': expand, 'reading': verdict}


def read_verdict(perm, topk_sens, de_pos1, sil, crit_k, small):
    flags = []
    if perm['p_observed'] >= ALPHA:
        flags.append('置换检验 p>=%.2f(判据无判别力)' % ALPHA)
    if len({v['pass'] for v in topk_sens['per_topk'].values()}) > 1:
        flags.append('判据随 topk 翻动')
    if not de_pos1['criterion']['pass']:
        flags.append('去 position-1 后判据不通过')
    n_best3 = sum(1 for v in sil['best_per_method'].values() if v['best_k'] == 3)
    if n_best3 < len(sil['best_per_method']) / 2:
        flags.append('k=3 在多数方法下非最优')
    if not all(v['pass'] for v in crit_k['per_k'].values()):
        flags.append('判据随 k 翻动')
    if small['identity_within']['mean'] >= IDENTITY_HIGH or small['n_members_genes'] < 4:
        flags.append('n=9 小簇有批次效应嫌疑')
    if len(flags) >= 3:
        verdict = '综合判定: 证据变弱'
    elif len(flags) >= 1:
        verdict = '综合判定: 证据不变(原"弱证据阳性"评级维持, 脆弱项如实并列)'
    else:
        verdict = '综合判定: 证据变强'
    return {'weakening_flags': flags, 'n_flags': len(flags),
            'reading': '%s —— 触发项: %s' % (verdict, '; '.join(flags) if flags else '无')}


# ---------- 主流程 ----------

EFFECTOR = 'cas12a2_zeng2026'


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', default=os.path.join(DATA, 'typing.robustness.json'))
    args = ap.parse_args()

    t_start = time.time()
    dr = core.to_rna(get_scaffold(EFFECTOR))
    names, seqs = load_spacers()
    with open(CLUSTERS, encoding='utf-8') as fh:
        clusters = json.load(fh)
    name2seq = {n: core.to_dna(s) for n, s in zip(names, seqs)}

    print('[robustness] 计算/载入 86 条特征矩阵 ...')
    F = feature_matrix(dr, names, seqs)
    Z = (F - F.mean(0)) / F.std(0)

    # ---- 补救1+2 的 9 代表 topk=16 重跑(k=3 档案代表) ----
    reps9 = {}
    for t, info in clusters['types'].items():
        for rep in info['representatives']:
            reps9[rep['name']] = t
    prefixes, n_new = {}, 0
    for nm in reps9:
        prefixes[nm], fresh = run_pipeline(nm, name2seq[nm])
        n_new += fresh
    print('[robustness] k=3 九代表 topk=16: 新跑 %d, 缓存 %d' % (n_new, 9 - n_new))

    # 切片一致性交叉核对(新跑 k=8 切片 vs 档案 topk=8 直接输出)
    crosscheck = {}
    for nm, t in reps9.items():
        arch = os.path.join(DATA, 'typing.t%s_%s.top.json' % (t, san(nm)))
        if os.path.exists(arch):
            with open(arch, encoding='utf-8') as fh:
                a8 = [r['desc'] for r in json.load(fh)['top'][1:9]]
            crosscheck[nm] = (top_descs(prefixes[nm], 8) == a8)
    n_match = sum(crosscheck.values())
    print('[robustness] 切片一致性: %d/9 与档案 topk=8 输出一致' % n_match)

    lists9 = {t: [top_descs(prefixes[nm], CRIT_TOPK)
                  for nm, tt in reps9.items() if tt == t]
              for t in sorted(set(reps9.values()))}

    # 不一致的成因量化: 档案 top-8 中多少 desc 已被当前代码的硬过滤(v1.6/v1.7,
    # 晚于档案生成)从突变库整体剔除(当前 variants.csv 中不存在)
    import csv as _csv
    archived_missing = {}
    for nm, t in reps9.items():
        arch = os.path.join(DATA, 'typing.t%s_%s.top.json' % (t, san(nm)))
        if not os.path.exists(arch):
            continue
        with open(arch, encoding='utf-8') as fh:
            a8 = [r['desc'] for r in json.load(fh)['top'][1:9]]
        with open(prefixes[nm] + '.variants.csv', newline='', encoding='utf-8') as fh:
            cur = {r['desc'] for r in _csv.DictReader(fh)}
        gone = [d for d in a8 if d not in cur]
        if gone:
            archived_missing[nm] = gone
    if n_match == 9:
        slice_note = '当前代码下 k=8 切片与档案 topk=8 输出全部一致'
    else:
        slice_note = ('当前代码下 k=8 切片与档案 topk=8 输出仅 %d/9 一致: 档案 top-8 中 '
                      '%d 条变体已被晚于档案生成的硬过滤(v1.6 交叉配对/v1.7 inv_max_run)'
                      '从突变库整体剔除(当前 variants.csv 不存在), 档案判据(公共集 3 条)在'
                      '当前代码下不可原样复现(现测公共集 %d 条); 本文件全部分析基于当前代码'
                      % (n_match, sum(len(v) for v in archived_missing.values()),
                         criterion(lists9, CRIT_TOPK)['n_common']))

    # ---- 补救1 置换检验 ----
    obs = criterion(lists9, CRIT_TOPK)
    lists_flat = [top_descs(prefixes[nm], CRIT_TOPK) for nm in reps9]
    null = []
    idx = list(range(9))
    for rest_a in itertools.combinations(idx[1:], 2):  # 元素0固定在A组, 破组间对称
        ga = (0,) + rest_a
        rem = [i for i in idx if i not in ga]
        for gb in itertools.combinations(rem, 3):
            gc = tuple(i for i in rem if i not in gb)
            if gb[0] > gc[0]:  # B/C 无序去重 -> 280 种无序划分
                continue
            unions = [set().union(*(lists_flat[i] for i in g)) for g in (ga, gb, gc)]
            null.append(len(set.intersection(*unions)))
    null = np.array(null)
    p_obs = float((null <= obs['n_common']).mean())
    p_crit = float((null <= obs['threshold']).mean())
    perm = read_permutation(obs['n_common'], null, p_obs, p_crit)
    perm['criterion'] = obs
    perm['n_partitions'] = int(len(null))
    perm['note'] = '280 种无序划分(统计量组间对称, 等价 1680 种标号分配)'
    print('[补救1] 观测公共集=%d, P(common<=obs)=%.4f, 判据零通过率=%.4f'
          % (obs['n_common'], p_obs, p_crit))

    # ---- 补救2a topk 敏感性 ----
    per_topk = {}
    for k in (6, 8, 10, 16):
        lt = {t: [top_descs(prefixes[nm], k) for nm, tt in reps9.items() if tt == t]
              for t in sorted(set(reps9.values()))}
        res = criterion(lt, k)
        per_topk[str(k)] = res
        print('[补救2a] topk=%d: 公共集 %d 条 %s -> %s'
              % (k, res['n_common'], res['common'], 'PASS' if res['pass'] else 'FAIL'))
    topk_sens = read_topk_sensitivity(per_topk)

    # ---- 补救2b 去 position-1 ----
    lt_np, removed = {}, {}
    for t in sorted(set(reps9.values())):
        lt_np[t] = []
        for nm, tt in reps9.items():
            if tt != t:
                continue
            before = top_descs(prefixes[nm], PERM_TOPK_RUN)
            after = top_descs_no_pos1(prefixes[nm], CRIT_TOPK)
            lt_np[t].append(after)
            removed[nm] = sum(1 for d in before
                              if 1 in [int(p) for p in MUT_POS_RE.findall(d)])
    de_pos1 = read_de_pos1(criterion(lt_np, CRIT_TOPK), removed)
    print('[补救2b] 去 pos-1 后公共集 %d 条 -> %s'
          % (de_pos1['criterion']['n_common'],
             'PASS' if de_pos1['criterion']['pass'] else 'FAIL'))

    # ---- 补救3a 轮廓系数谱 ----
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.mixture import GaussianMixture
    from sklearn.metrics import silhouette_score
    sil_table = {'kmeans_repo': {}, 'agglomerative_ward': {}, 'gaussian_mixture': {}}
    for k in range(2, 7):
        lab_km, _, _ = typing_mod.kmeans(Z, k, seeds=range(10))
        sil_table['kmeans_repo'][k] = round(typing_mod.silhouette(Z, lab_km, k), 3)
        lab_ag = AgglomerativeClustering(n_clusters=k, linkage='ward').fit_predict(Z)
        sil_table['agglomerative_ward'][k] = round(float(silhouette_score(Z, lab_ag)), 3)
        lab_gm = GaussianMixture(n_components=k, random_state=0, n_init=10).fit_predict(Z)
        sil_table['gaussian_mixture'][k] = round(float(silhouette_score(Z, lab_gm)), 3)
    sil = read_silhouette_spectrum(sil_table)
    print('[补救3a]', {m: v['best_k'] for m, v in sil['best_per_method'].items()})

    # ---- 补救3b 判据面 k=2..6(手写 KMeans, 预登记方法) ----
    needed = {}
    per_k_crit, k_labels = {}, {}
    for k in range(2, 7):
        lab, _, _ = typing_mod.kmeans(Z, k, seeds=range(10))
        k_labels[k] = lab
        for c, reps in reps_from_clusters(lab, Z, names, seqs, k).items():
            for nm, dna in reps:
                needed[nm] = dna
    n_fresh = 0
    for nm, dna in needed.items():
        if nm not in prefixes:
            prefixes[nm], fresh = run_pipeline(nm, dna)
            n_fresh += fresh
    print('[robustness] k=2..6 判据面: 代表去重 %d 名, 新跑 %d, 缓存 %d'
          % (len(needed), n_fresh, len(needed) - n_fresh))
    for k in range(2, 7):
        lab = k_labels[k]
        reps = reps_from_clusters(lab, Z, names, seqs, k)
        lt = {c: [top_descs(prefixes[nm], CRIT_TOPK) for nm, _ in rl]
              for c, rl in reps.items()}
        res = criterion(lt, CRIT_TOPK)
        res['sizes'] = [int((lab == c).sum()) for c in range(k)]
        res['silhouette'] = round(typing_mod.silhouette(Z, lab, k), 3)
        res['representatives'] = {c: [nm for nm, _ in rl] for c, rl in reps.items()}
        per_k_crit[str(k)] = res
        print('[补救3b] k=%d: sizes=%s 公共集 %d 条 -> %s'
              % (k, res['sizes'], res['n_common'], 'PASS' if res['pass'] else 'FAIL'))
    crit_k = read_criterion_vs_k(per_k_crit)

    # ---- 补救4 n=9 小簇来源审查 ----
    t1 = clusters['types']['1']
    members = t1['members']
    genes = {m.split()[0] for m in members}
    m_seqs = [name2seq[m].replace('T', 'U') for m in members]
    rest = [n for n in names if n not in set(members)]
    r_seqs = [name2seq[n].replace('T', 'U') for n in rest]

    def pair_stats(list_a, list_b):
        vals = [seq_identity(a, b) for a, b in
                (itertools.combinations(list_a, 2) if list_b is None
                 else itertools.product(list_a, list_b))]
        return {'mean': round(float(np.mean(vals)), 3),
                'median': round(float(np.median(vals)), 3),
                'max': round(float(np.max(vals)), 3), 'n_pairs': len(vals)}

    ident_in = pair_stats(m_seqs, None)
    ident_bg = pair_stats(r_seqs, None)
    ident_cross = pair_stats(m_seqs, r_seqs)

    i1 = [names.index(m) for m in members]
    iR = [names.index(n) for n in rest]
    feat_diff, dominant = {}, None
    for j, fn in enumerate(FEAT_NAMES):
        sd = F[iR, j].std()
        d = (F[i1, j].mean() - F[iR, j].mean()) / (sd if sd > 0 else 1.0)
        feat_diff[fn] = round(float(d), 2)
    dominant = max(feat_diff, key=lambda f: abs(feat_diff[f]))
    expand = sorted(
        [{'name': n, 'junction_pairs': round(float(F[names.index(n), 3]), 1)}
         for n in rest if F[names.index(n), 3] >= JUNCTION_EXPAND_MIN],
        key=lambda r: -r['junction_pairs'])
    small = read_small_cluster(genes, ident_in, ident_bg, feat_diff, dominant, expand)
    small['identity_cross'] = ident_cross
    small['members'] = members
    small['feature_means_in'] = {fn: round(float(F[i1, j].mean()), 3)
                                 for j, fn in enumerate(FEAT_NAMES)}
    small['feature_means_rest'] = {fn: round(float(F[iR, j].mean()), 3)
                                   for j, fn in enumerate(FEAT_NAMES)}
    print('[补救4]', small['reading'][:120])

    verdict = read_verdict(perm, topk_sens, de_pos1, sil, crit_k, small)
    print('[综合]', verdict['reading'])

    report = {
        'generated_by': 'scripts/crrna_typing_robustness.py',
        'date': time.strftime('%Y-%m-%d'),
        'scope': SCOPE,
        'effector': EFFECTOR,
        'runtime_sec': round(time.time() - t_start, 1),
        'pipeline_reruns': {'topk': PERM_TOPK_RUN,
                            'distinct_representatives': len(needed) + 9 - len(set(needed) & set(reps9)),
                            'fresh_runs': n_new + n_fresh},
        'slice_consistency_check': {'n_match_of_9': n_match, 'per_rep': crosscheck,
                                    'archived_top8_now_filtered': archived_missing,
                                    'note': slice_note},
        'permutation_test': perm,
        'topk_sensitivity': topk_sens,
        'de_pos1_sensitivity': de_pos1,
        'clustering_spectrum': sil,
        'criterion_vs_k': crit_k,
        'small_cluster_audit': small,
        'verdict': verdict,
    }
    with open(args.out, 'w', encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print('[robustness] -> %s (%.1fs)' % (args.out, time.time() - t_start))


if __name__ == '__main__':
    main()
