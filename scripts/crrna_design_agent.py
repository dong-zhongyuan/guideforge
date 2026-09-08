"""crRNA 设计智能体（ spacer 设计 × 骨架分型选择 组合入口）

输入目标突变转录本(+WT 对照), 输出 spacer × 推荐骨架 的组合候选:
  1) spacer 设计器(A 部分): 围绕突变位点 tilling 枚举全部覆盖窗口
     (突变置于 PFS 或 protospacer 不同位置, Zeng 2026 口径), 按鉴别档位
     (PFS 鉴别 > protospacer 单错配 > 其它)与结构代理(spacer 游离度/自互补)排序;
  2) 骨架选择器(B 部分): 新 spacer 的 5 维上下文特征, 与既有分型
     (data/context_typing/typing.clusters.json, 3 型)的型心比欧氏距离,
     给最近型 + 置信度(最近/次近距离比); 边界样本如实报两型;
  3) 各型代表骨架从分型试验的管线产出中取(每型首个非 5' 悬垂惰性位变体 + WT 备选);
  4) 模块一特征向量的丰度维(V3 §4.4): --target-key 指定 panel 靶标时, 从
     data/target_abundance_tiers.json(CCLE 18q3 实测 + 显式分档规则, 见
     scripts/crrna_target_abundance.py)取靶RNA丰度档位写入每条候选的
     target_features, 并经 activation_gate(显式阈值判读)进入模块三推荐输出;
  4b) 模块一特征向量的突变类型维(V3 §4.4 / 表2 APC 入选理由, 2026-09-07):
     提供 --mut-fasta + --wt-fasta 时以最长 ORF 法判定突变类型
     (错义/无义/移码, 规则与边界处理显式声明于
     scripts/crrna_mutation_typing.py docstring), mutation_type 与紧凑 ORF
     证据(ORF 蛋白长度/提前终止密码子/长度变化)写入每条候选的
     target_features, 完整 ORF 证据(WT/MUT 最长 ORF 起止、MUT 同起点对应
     ORF、判读文字)写入输出 JSON 的 mutation.typing 块; 判读文字全部由
     数值按显式规则生成;
  4c) 结构层界面特征差量附加层(V3 §4.1 末段 / 表1「选型特征」, 2026-09-09
     任务⑫): --target-key 指定 panel 靶标时, target_context 增挂
     interface_deltas 块(数据源 data/chai_interface_deltas.json, 由
     scripts/crrna_chai_interface_features.py 从 Chai 矩阵提取: 该靶标列
     全部非 WT 骨架相对同靶 WT 组合的界面特征差量 + WT 参照 + 局限注记);
     文献决策树不消费该维(文献训练行无结构数据), 作结构描述特征随
     design.json 归档并与 webapp 同源展示; 矩阵无该靶标列(如 APC 待增量
     折叠回填)时块内 available=False 并显式标注缺失原因;
  5) 阴性结果退化模式(策划案 V3 §5.3 末段 / 表5 风险1 对策, 2026-09-06):
     每次运行先按预登记 go/no-go 判据(docs/preregistration.md §A)从
     --clusters 产物重算"骨架最优解与靶标(spacer 上下文)特征是否显著相关"
     (interaction_assessment): 判据未通过(阴性, 各型 TOP 骨架高度重合)
     时输出从分型推荐退化为通用型单骨架推荐(universal_scaffold,
     既有数据口径的全场最优), 并在 JSON/控制台如实标注退化模式、交互弱的
     量化证据与"科学结论依然完整可交付"的交付口径; --degrade on/off 可强制
     演示/覆盖, 强制状态在输出中如实标注(forced 字段)。

边界(诚实声明): 排序为透明启发式, 非活性预测; 分型功能差异未实验验证;
骨架选择器为最近邻规则(样本量决定), 不是训练模型; 丰度档位为模块间特征
与判读门控输入, 不参与最近邻距离计算(分型模型为固定 5 维, 重聚类才可选型
面扩维, 见 typing.clusters.json features)。退化模式与既有"置信度不足报
两型"(boundary 置信度给 secondary_type)是两个层级: 后者是单条 spacer 在
分型面上的归属不确定(逐条样本级), 前者是分型假说本身未通过判据(全局
模式级); 退化模式下分型结果仍逐条计算但仅作诊断(typing_nearest_type),
不参与推荐。

用法:
  python crrna_design_agent.py --mut-fasta ../data/tp53_r248q_mrna_NM_000546.fa \
      --wt-fasta ../data/tp53_mrna_NM_000546.fa --effector cas12a2_zeng2026 \
      --clusters ../data/context_typing/typing.clusters.json \
      --spacers ../data/context_typing/spacers.txt \
      --out-prefix ../data/agent/tp53_r248q
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

import crrna_scaffold_design as core  # noqa: E402
import crrna_context_typing as typing_mod  # noqa: E402
import crrna_target_abundance as tabund  # noqa: E402
import crrna_mutation_typing as mtyping  # noqa: E402
import crrna_chai_interface_features as chai_if  # noqa: E402
from scaffold_registry import get_scaffold, resolve_pfs, pfs_mismatches  # noqa: E402

COMP = str.maketrans({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'})


def read_fasta_single(path):
    seq = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if not line.startswith('>'):
                seq.append(line.strip())
    return core.to_rna(''.join(seq))


def find_mutation(wt, mut):
    """单碱基差异定位(1-based); 多差异时取第一个并如实返回全部。"""
    diffs = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
    if not diffs and len(wt) != len(mut):
        raise ValueError('长度不同且无替换差异, 暂只支持单碱基替换')
    if len(diffs) > 1:
        print(f'[agent] 警告: 检出 {len(diffs)} 处差异, 取第一处(其余如实列出): '
              f'{[d + 1 for d in diffs][:6]}')
    return diffs[0] + 1, diffs


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def design_spacers(wt, mut, mpos, spacer_len=24, pfs_rule=None, pfs_tol=None,
                   effector='cas12a2_zeng2026', topn=10):
    """围绕突变 tilling: 突变可落在 PFS(窗口 3' 下游 1..5 位)或 protospacer 内。
    返回候选列表(按鉴别档位+结构代理排序)。

    PFS 规则默认取注册表 pfs 字段(round-3 R5 统一口径: 共识 + 容忍错配数);
    pfs_rule/pfs_tol 可显式覆盖, 如 pfs_rule='CAGAG'(Zeng 2026 R248Q 位点实测
    PFS, 与 GAAAG 共识差 2 个错配——第 1、3 位)。
    鉴别判定: 突变型 PFS 落在共识容忍范围内 且 WT 型落在范围外
    (突变"产生"可识别 PFS), 替代旧的首字符硬匹配。"""
    spec = resolve_pfs(effector, consensus=pfs_rule, tolerant_mismatches=pfs_tol)
    consensus, tol = spec['consensus'], spec['tolerant_mismatches']
    cands = []
    # 突变在 PFS: 窗口终点 = mpos - k (k=1..5), 即窗口 [mpos-k-L, mpos-k)
    for k in range(1, 6):
        end = mpos - k              # 1-based 排他终点
        start = end - spacer_len
        if start < 0:
            continue
        spacer = revcomp(mut[start:end])       # spacer = 靶 RNA 窗口的反向互补
        wt_window = revcomp(wt[start:end])
        pfs_mut = mut[end:end + 5]
        pfs_wt = wt[end:end + 5]
        mm = sum(a != b for a, b in zip(spacer, wt_window))
        if mm != 0:
            continue                # protospacer 必须与 WT 相同(鉴别全靠 PFS)
        pfs_mm_mut = pfs_mismatches(pfs_mut, consensus)
        pfs_mm_wt = pfs_mismatches(pfs_wt, consensus)
        discriminates = (pfs_mm_mut is not None and pfs_mm_wt is not None
                         and pfs_mm_mut <= tol and pfs_mm_wt > tol)
        cands.append({'spacer': spacer, 'mut_in': f'PFS+{k}', 'proto_mm': mm,
                      'pfs_mut': pfs_mut, 'pfs_wt': pfs_wt,
                      'pfs_mm_mut': pfs_mm_mut, 'pfs_mm_wt': pfs_mm_wt,
                      'pfs_discriminates': discriminates,
                      'tier': 2 if pfs_mut != pfs_wt else 0,
                      'start': start + 1})
    # 突变在 protospacer 内: 窗口覆盖 mpos
    for start in range(max(mpos - spacer_len, 0), mpos):
        spacer = revcomp(mut[start:start + spacer_len])
        wt_window = revcomp(wt[start:start + spacer_len])
        mm = sum(a != b for a, b in zip(spacer, wt_window))
        if mm == 1:
            pfs_mut = mut[start + spacer_len:start + spacer_len + 5]
            cands.append({'spacer': spacer, 'mut_in': f'protospacer@{mpos - start}',
                          'proto_mm': mm, 'pfs_mut': pfs_mut,
                          'pfs_wt': wt[start + spacer_len:start + spacer_len + 5],
                          'pfs_discriminates': False, 'tier': 1, 'start': start + 1})
    # 结构代理(同档内弱自互补优先: self_mfe 越接近 0 越游离)
    for c in cands:
        _, self_mfe = core.fold(c['spacer'])
        c['self_mfe'] = round(self_mfe, 2)
    cands.sort(key=lambda c: (-c['tier'], -c['self_mfe']))
    return cands


def load_type_models(clusters_json, spacers_txt, dr):
    """从分型产物重建: 标准化参数 + 各型型心 + 各型代表骨架。"""
    with open(clusters_json, encoding='utf-8') as fh:
        rep = json.load(fh)
    names, seqs = [], []
    with open(spacers_txt, encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                n, s = line.rstrip('\n').split('\t')
                names.append(n)
                seqs.append(core.to_rna(s))
    F = np.array([typing_mod.spacer_features(dr, n, s) for n, s in zip(names, seqs)])
    mu, sd = F.mean(0), F.std(0)
    Z = (F - mu) / sd
    name2z = {n: z for n, z in zip(names, Z)}
    centers = {}
    for t, info in rep['types'].items():
        zs = np.array([name2z[m] for m in info['members'] if m in name2z])
        centers[int(t)] = zs.mean(0)
    # 各型代表骨架: 分型试验 top.json 里首个非 5' 悬垂惰性位(位置1)变体 + WT
    type_dr = {}
    prefix_dir = os.path.dirname(clusters_json)
    for t in rep['types']:
        best, wt_dr = None, None
        for pj in sorted(glob.glob(os.path.join(prefix_dir, f'typing.t{t}_*.top.json'))):
            with open(pj, encoding='utf-8') as fh:
                top = json.load(fh)
            wt_dr = core.to_dna(dr)
            for row in top['top'][1:]:
                if row['mut_positions'] and min(row['mut_positions']) > 1 and row['passed']:
                    if best is None or row['score'] > best[1]:
                        best = (row['desc'], row['score'], row['dr_seq'])
        type_dr[int(t)] = {'representative_dr': best[2] if best else core.to_dna(dr),
                           'representative_desc': best[0] if best else 'WT',
                           'representative_score': best[1] if best else 0.0,
                           'wt_dr': wt_dr or core.to_dna(dr)}
    return {'centers': centers, 'mu': mu, 'sd': sd, 'feature_names': rep['features'],
            'type_profiles': {t: i['feature_means'] for t, i in rep['types'].items()},
            'type_dr': type_dr}


def select_scaffold_type(spacer_rna, dr, model):
    f = np.array(typing_mod.spacer_features(dr, 'query', spacer_rna))
    z = (f - model['mu']) / model['sd']
    d = {t: float(np.sqrt(((z - c) ** 2).sum())) for t, c in model['centers'].items()}
    rank = sorted(d, key=d.get)
    t1, t2 = rank[0], rank[1]
    ratio = d[t1] / d[t2] if d[t2] else 0.0
    conf = 'high' if ratio < 0.6 else ('medium' if ratio < 0.85 else 'boundary')
    return {'features': {n: round(float(v), 3) for n, v in zip(model['feature_names'], f)},
            'distances': {str(t): round(v, 3) for t, v in d.items()},
            'nearest_type': t1, 'confidence': conf,
            'secondary_type': t2 if conf == 'boundary' else None}


# ---------- 阴性结果退化模式(策划案 V3 §5.3 末段 / 表5 风险1 对策, 2026-09-06) ----------
#
# 判据本体 = 预登记 go/no-go 判据(docs/preregistration.md §A, 项目设定常数):
#   各型代表重跑管线的 TOP-k DR desc, 每型取并集, 全型公共子集大小
#   n_common <= topk//2(默认 topk=8, 阈=4) 且型数>1 -> 阳性(骨架最优解随
#   spacer 上下文变化, 维持分型推荐); 否则 -> 阴性(骨架最优解与靶标特征
#   无显著相关, 骨架-靶标交互弱), 智能体退化为通用型单骨架推荐。
# 本模块是该判据的下游消费点(不重定义判据, 阈值 topk//2 与
# crrna_context_typing.py:188-190 同口径); 判读文字全部由数值按显式规则
# 生成, 不手写结论方向。
DEGRADE_TOPK = 8        # 预登记判据 TOP 深度(§A; 与分型试验 --topk 默认一致)
DEGRADE_ALPHA = 0.05    # 事后置换检验显著性阈(与 crrna_typing_robustness.ALPHA 一致)


def gonogo_verdict(per_type_unions, topk=DEGRADE_TOPK):
    """预登记 §A go/no-go 判据的纯规则实现(与 crrna_context_typing.py 同口径)。

    per_type_unions: {型号: 该型各代表 TOP-k DR desc 的并集(set/list)}。
    返回 dict(significant, n_common, threshold, common, n_types);
    significant=False 即阴性 -> 退化模式。"""
    n_types = len(per_type_unions)
    common = (set.intersection(*[set(u) for u in per_type_unions.values()])
              if n_types > 1 else set())
    threshold = topk // 2
    significant = n_types > 1 and len(common) <= threshold
    return {'significant': significant, 'n_common': len(common),
            'threshold': threshold, 'common': sorted(common), 'n_types': n_types}


def per_type_top_unions(clusters_json, topk=DEGRADE_TOPK):
    """从分型试验产物重建各型 TOP-k DR desc 并集(判据输入)。

    数据源 = 与 load_type_models 同源的 typing.t<t>_*.top.json(每型各代表
    重跑管线的归档输出), 每文件取 top[1:topk+1](跳过 rank0 WT 行)的 desc,
    型内并集——与 crrna_context_typing.py 判据的 tops 构造逐行对应。
    返回 (unions, archived_common): archived_common 为 clusters.json 归档的
    top_common_to_all_types(无则 None), 供交叉核对。"""
    with open(clusters_json, encoding='utf-8') as fh:
        rep = json.load(fh)
    prefix_dir = os.path.dirname(clusters_json)
    unions = {}
    for t in rep['types']:
        u = set()
        for pj in sorted(glob.glob(os.path.join(
                prefix_dir, 'typing.t%s_*.top.json' % t))):
            with open(pj, encoding='utf-8') as fh:
                top = json.load(fh)['top']
            u.update(row['desc'] for row in top[1:topk + 1])
        unions[str(t)] = u
    return unions, rep.get('top_common_to_all_types')


def _verdict_reading(v, topk):
    """判读文字: 由判据数值按显式规则生成(两分支模板, 方向由 significant 决定)。"""
    if v['significant']:
        return ('各型代表管线 TOP-%d 最优骨架的全型公共子集 %d 条 <= 阈 %d'
                '(topk//2, 预登记判据 §A): 骨架最优解随靶标(spacer 上下文)'
                '特征变化, 分型判据通过(阳性), 维持分型推荐'
                % (topk, v['n_common'], v['threshold']))
    return ('各型代表管线 TOP-%d 最优骨架的全型公共子集 %d 条 > 阈 %d'
            '(topk//2, 预登记判据 §A): 各型最优骨架高度重合, 骨架最优解与'
            '靶标特征无显著相关(阴性), 智能体退化为通用型单骨架推荐'
            % (topk, v['n_common'], v['threshold']))


def _robustness_evidence(clusters_json):
    """交互弱/强的补充量化证据: 同源事后稳健性分析的置换检验(若存在)。

    只读 data/context_typing/typing.robustness.json 的 permutation_test
    (post-hoc, 非预登记, 不作触发判据); 判读由 p 值按 DEGRADE_ALPHA 规则生成。"""
    path = os.path.join(os.path.dirname(clusters_json), 'typing.robustness.json')
    if not os.path.isfile(path):
        return None
    with open(path, encoding='utf-8') as fh:
        perm = (json.load(fh).get('permutation_test') or {})
    p = perm.get('p_observed')
    if p is None:
        return None
    if p >= DEGRADE_ALPHA:
        reading = ('事后置换检验(post-hoc): P(零分布公共集<=观测 %d 条)=%.4f '
                   '>= %.2f, 观测公共集在随机分组下常见, 与"骨架-靶标交互弱"'
                   '一致' % (perm.get('observed', -1), p, DEGRADE_ALPHA))
    else:
        reading = ('事后置换检验(post-hoc): P(零分布公共集<=观测 %d 条)=%.4f '
                   '< %.2f, 观测公共集在随机分组下罕见, 支持分型判据具判别力'
                   % (perm.get('observed', -1), p, DEGRADE_ALPHA))
    return {'source': os.path.basename(path) + '(post-hoc, 非预登记, 不作触发判据)',
            'permutation_observed_common': perm.get('observed'),
            'permutation_p_observed': p, 'reading': reading}


def interaction_assessment(clusters_json, topk=DEGRADE_TOPK, force='auto'):
    """阴性退化判据评估: 从分型产物重算预登记 §A 判据, 决定推荐模式。

    返回的 mode 只有两种:
      'typed'              阳性, 维持分型推荐(现状);
      'degraded_universal' 阴性, 退化为通用型单骨架推荐。
    force: 'auto'(默认, 按判据数值)/'on'(强制退化, 阴性场景演示或未来湿实验
    阴性回流)/'off'(强制分型); 强制时判据数值仍如实计算并标注 forced 字段。
    判据不可评估(某型缺 top.json 产物)时 mode='typed' + evaluable=False,
    如实标注而非静默退化。"""
    unions, archived_common = per_type_top_unions(clusters_json, topk)
    evaluable = all(len(u) > 0 for u in unions.values()) and len(unions) > 1
    v = gonogo_verdict(unions, topk)
    ts = sorted(unions)
    jacs = []
    for i in range(len(ts)):
        for j in range(i + 1, len(ts)):
            a, b = unions[ts[i]], unions[ts[j]]
            jacs.append(len(a & b) / len(a | b) if a | b else 0.0)
    out = {'criterion_id': 'preregistration §A go/no-go',
           'criterion_rule': ('各型代表管线 TOP-%d DR desc 型内并集的全型公共子集 '
                              '<= %d(topk//2) 且型数>1 -> 阳性(分型推荐); 否则阴性 '
                              '-> 退化通用型单骨架推荐' % (topk, topk // 2)),
           'topk': topk, 'threshold': v['threshold'],
           'n_types': v['n_types'], 'evaluable': evaluable,
           'per_type_union': {t: sorted(unions[t]) for t in ts},
           'per_type_union_sizes': {t: len(unions[t]) for t in ts},
           'common': v['common'], 'n_common': v['n_common'],
           'pairwise_jaccard_mean': (round(float(np.mean(jacs)), 3)
                                     if jacs else None),
           'significant': v['significant'] if evaluable else None,
           'forced': None if force == 'auto' else force,
           'robustness_evidence': _robustness_evidence(clusters_json)}
    if archived_common is not None:
        out['archive_crosscheck'] = {
            'archived_common': sorted(archived_common),
            'recomputed_common': v['common'],
            'matches': sorted(archived_common) == v['common']}
    if not evaluable:
        out['mode'] = 'typed'
        out['reading'] = ('分型判据不可评估(型数 %d 或某型缺 top.json 产物), '
                          '维持分型推荐并如实标注' % v['n_types'])
    else:
        out['reading'] = _verdict_reading(v, topk)
        out['mode'] = 'typed' if v['significant'] else 'degraded_universal'
    if force == 'on' and out['mode'] != 'degraded_universal':
        out['mode'] = 'degraded_universal'
        out['reading'] += '; 本次 --degrade on 强制退化(判据实测为阳性, 强制仅作阴性场景演示)'
    elif force == 'off' and out['mode'] != 'typed':
        out['mode'] = 'typed'
        out['reading'] += '; 本次 --degrade off 强制分型(判据实测为阴性, 覆盖须审慎)'
    if out['mode'] == 'degraded_universal':
        rel = '<=' if out['n_common'] <= out['threshold'] else '>'
        out['deliverable_note'] = (
            '阴性结果交付口径(策划案 V3 §5.3 末段 / 表5 风险1 对策): 结论为'
            '"骨架-靶标交互弱, 通用型最优"; 通用型单骨架推荐保留系统价值, '
            '交互弱的量化证据见本块(公共集 %d 条 %s 阈 %d, 两两 Jaccard 均值 %s'
            '%s), 科学结论依然完整可交付'
            % (out['n_common'], rel, out['threshold'],
               out['pairwise_jaccard_mean'],
               '; 置换检验 p=%.4f' % out['robustness_evidence']['permutation_p_observed']
               if out['robustness_evidence'] else ''))
    return out


def universal_scaffold(model):
    """通用型单骨架(退化模式推荐): 既有数据口径的全场最优。

    候选池 = WT(管线分 0 基线) + 各型代表骨架(该型 spacer 上下文内管线分
    最高的非 position-1 变体, 与 load_type_models 同口径); 取分最高者,
    并列时 WT 优先、再按 desc 字典序(显式确定性规则)。
    跨型比较口径: 各代表得分均为其自身上下文内相对 WT(=0) 的管线分
    (score_variant; selector_model.json ranking_authority=pipeline_score)。"""
    cands = []
    wt_dr = None
    for t, v in sorted(model['type_dr'].items()):
        wt_dr = wt_dr or v['wt_dr']
        if v['representative_desc'] != 'WT':
            cands.append({'desc': v['representative_desc'],
                          'dr_dna': v['representative_dr'],
                          'score': v['representative_score'],
                          'source': 'type_%s_representative' % t})
    cands.append({'desc': 'WT', 'dr_dna': wt_dr, 'score': 0.0,
                  'source': 'wt_baseline'})
    cands.sort(key=lambda c: (-c['score'], c['desc'] != 'WT', c['desc']))
    best = dict(cands[0])
    best['rule'] = ('候选池=WT(0 基线)+各型代表骨架; 取管线分最高, 并列 WT '
                    '优先再按 desc 字典序')
    best['pool'] = cands
    return best


def target_context(target_key, abundance_path=None, interface_path=None):
    """模块一靶标解析的丰度维 + 模块三判读门控(同源 data/target_abundance_tiers.json)
    + 结构层界面特征差量附加层(同源 data/chai_interface_deltas.json, 任务⑫)。

    丰度档位不进入最近邻距离(分型模型固定 5 维), 作为靶标特征向量的一维
    随每条候选输出, 并由 activation_gate 按显式阈值生成细胞层判读。
    界面差量同样不进最近邻距离与文献决策树(文献训练行无结构数据),
    作结构描述特征(选型特征维)随 target_context.interface_deltas 归档与
    展示; 矩阵无该靶标列时 available=False 并显式标注缺失原因。"""
    if not target_key:
        return None
    path = abundance_path or os.path.join(ROOT, 'data',
                                          'target_abundance_tiers.json')
    if not os.path.isfile(path):
        raise SystemExit('缺丰度档位表 %s, 先运行 scripts/crrna_target_abundance.py'
                         % path)
    feat = tabund.target_feature(tabund.load_abundance(path), target_key)
    if feat is None:
        raise SystemExit('未知 panel 靶标 --target-key %s(可选: %s)'
                         % (target_key, sorted(tabund.PANEL_TARGETS)))
    ipath = interface_path or os.path.join(ROOT, 'data',
                                           'chai_interface_deltas.json')
    if not os.path.isfile(ipath):
        raise SystemExit('缺界面差量特征表 %s, 先运行 '
                         'scripts/crrna_chai_interface_features.py' % ipath)
    return {'feature': feat, 'gate': tabund.activation_gate(feat),
            'interface_deltas': chai_if.target_block(
                chai_if.load_deltas(ipath), target_key),
            'source': os.path.basename(path),
            'interface_source': os.path.basename(ipath)}


def _scaffold_fields(model, sel, degraded, uni):
    """每条候选的骨架推荐字段(分型推荐 vs 退化通用型, 两模式唯一分叉点)。

    分型模式(typed): 行为与退化模式引入前一致——最近型代表骨架; 置信度
    boundary 时 secondary_type 如实报次近型(样本级"置信度不足给多候选"
    的既有行为, 保留不变)。
    退化模式(degraded_universal): 全部候选统一取通用型单骨架(uni,
    universal_scaffold 的全场最优), scaffold_type=None; 分型结果仍逐条
    计算但降级为诊断(typing_nearest_type + confidence 键), 不参与
    推荐(全局模式级退化, 与样本级 boundary 多候选是两个层级, 见模块
    docstring 第 5 点)。"""
    if degraded:
        return {'scaffold_type': None,
                'typing_nearest_type': sel['nearest_type'],
                'dr_dna': uni['dr_dna'], 'dr_desc': uni['desc'],
                # WT 骨架对照构建: 各型 wt_dr 同源(= effector 注册表 DR)
                'wt_dr': model['type_dr'][sel['nearest_type']]['wt_dr']}
    dr_pick = model['type_dr'][sel['nearest_type']]
    return {'scaffold_type': sel['nearest_type'],
            'dr_dna': dr_pick['representative_dr'],
            'dr_desc': dr_pick['representative_desc'],
            'wt_dr': dr_pick['wt_dr']}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mut-fasta', default=None, help='突变转录本 FASTA(与 --wt-fasta 一起用于 tilling 设计 A)')
    ap.add_argument('--wt-fasta', default=None, help='WT 转录本 FASTA')
    ap.add_argument('--spacer', default=None,
                    help='直接输入 spacer(DNA/RNA, 17-25nt); 提供时跳过 tilling 设计, 只做骨架分型选择(项目本体口径)')
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--pfs', default=None,
                    help='PFS 共识序列(默认取注册表 pfs 字段; 如 CAGAG=Zeng 2026 R248Q 位点实测 PFS)')
    ap.add_argument('--pfs-tol', type=int, default=None,
                    help='PFS 容忍错配数(默认取注册表 pfs.tolerant_mismatches)')
    ap.add_argument('--spacer-len', type=int, default=24)
    ap.add_argument('--topn', type=int, default=10)
    ap.add_argument('--target-key', default=None,
                    help='panel 靶标键(tp53_r248q/kras_g12d/tp53_r273h/apc_q1328x); '
                         '提供时模块一特征向量附靶RNA丰度档维(data/target_abundance_tiers.json)')
    ap.add_argument('--abundance', default=None,
                    help='丰度档位表路径(默认 data/target_abundance_tiers.json)')
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--spacers', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--degrade', choices=('auto', 'on', 'off'), default='auto',
                    help='阴性退化模式(策划案 V3 §5.3 末段/表5 风险1): auto=按预登记 §A '
                         'go/no-go 判据从 --clusters 产物自动判定(默认); on=强制退化为'
                         '通用型单骨架推荐(阴性场景演示/未来湿实验阴性回流); off=强制'
                         '分型推荐(覆盖, 输出如实标注 forced)')
    args = ap.parse_args()

    dr = core.to_rna(get_scaffold(args.effector))
    model = load_type_models(args.clusters, args.spacers, dr)
    print(f'[agent] 骨架分型: {len(model["centers"])} 型; 各型代表骨架 '
          f'{ {t: v["representative_desc"] for t, v in model["type_dr"].items()} }')

    # 阴性退化判据(每次运行先评估): 模式 typed/degraded_universal + 量化证据
    assess = interaction_assessment(args.clusters, force=args.degrade)
    degraded = assess['mode'] == 'degraded_universal'
    uni = universal_scaffold(model) if degraded else None
    if degraded:
        assess['universal_scaffold'] = uni
    print(f"[agent] 交互判据(§A, topk={assess['topk']}): {assess['reading']}")
    if degraded:
        print(f"[agent] 退化模式: 通用型单骨架 = {uni['desc']}"
              f"(管线分 {uni['score']:.4f}, {uni['source']}); "
              f"分型结果仅作诊断不参与推荐")

    tctx = target_context(args.target_key, args.abundance)
    tfeat = tctx['feature'] if tctx else None
    if tctx:
        print(f"[agent] 靶标解析丰度维: {tfeat['gene']}@{tfeat['cell_line']} "
              f"total={tfeat['total_rpkm']:.2f} RPKM, 杂合场景 "
              f"{tfeat['mut_rpkm_het50']:.2f} -> 档位 "
              f"{tfeat['target_abundance_tier']}({tfeat['target_abundance_tier_label']})")
        print(f"[agent] 模块三门控: {tctx['gate']['note']}")

    if args.spacer:
        # 简版(项目本体口径): spacer 由实验侧输入, 只做骨架分型选择
        spacer_rna = core.to_rna(args.spacer)
        if not 17 <= len(spacer_rna) <= 25 or set(spacer_rna) - set('ACGU'):
            raise SystemExit('--spacer 必须为 17-25nt ACGT/U')
        sel = select_scaffold_type(spacer_rna, dr, model)
        sf = _scaffold_fields(model, sel, degraded, uni)
        row = {'spacer_dna': core.to_dna(spacer_rna), 'mut_in': 'input',
               'proto_mm': None, 'pfs_mut': '', 'pfs_wt': '', 'tier': None,
               'self_mfe': round(core.fold(spacer_rna)[1], 2),
               'scaffold_type': sf['scaffold_type'], 'confidence': sel['confidence'],
               'type_features': sel['features'], 'distances': sel['distances'],
               'target_features': tfeat,
               'dr_dna': sf['dr_dna'], 'dr_desc': sf['dr_desc'],
               'construct_dna': core.to_dna(sf['dr_dna'] + spacer_rna),
               'wt_dr_construct_dna': core.to_dna(sf['wt_dr'] + spacer_rna)}
        if degraded:
            # 分型置信度仍记于 confidence 键(语义=分型归属置信度, 诊断用);
            # 退化模式下推荐与分型解耦, 以 typing_nearest_type 标注诊断来源
            row['typing_nearest_type'] = sf['typing_nearest_type']
        out_rows = [row]
        if degraded:
            print(f"[agent] 输入 spacer {core.to_dna(spacer_rna)} -> 退化模式: "
                  f"通用型单骨架 {uni['desc']}(分型诊断: 最近型"
                  f"{sel['nearest_type']}({sel['confidence']}), 不参与推荐)")
        else:
            print(f"[agent] 输入 spacer {core.to_dna(spacer_rna)} -> 型{sel['nearest_type']}"
                  f"({sel['confidence']}) 推荐骨架 {sf['dr_desc']}")
        os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
        with open(args.out_prefix + '.design.json', 'w', encoding='utf-8') as fh:
            json.dump({'mode': 'spacer-input(骨架分型选择)', 'effector': args.effector,
                       'target_context': tctx,
                       'degradation': assess,
                       'designs': out_rows,
                       'disclaimer': '选择器为最近邻规则, 分型功能差异未实验验证'},
                      fh, ensure_ascii=False, indent=1)
        with open(args.out_prefix + '.design.fasta', 'w', encoding='utf-8') as fh:
            r = out_rows[0]
            tag = 'universal' if degraded else f"type{r['scaffold_type']}"
            fh.write(f">agent_input|{tag}|{r['dr_desc']}|{args.effector}\n"
                     f"{r['construct_dna']}\n>agent_input_WTdr|WT\n{r['wt_dr_construct_dna']}\n")
        print(f"[agent] 输出 -> {args.out_prefix}.design.json / .design.fasta")
        return

    if not (args.mut_fasta and args.wt_fasta):
        raise SystemExit('需要 --spacer 或 (--mut-fasta + --wt-fasta)')
    wt = read_fasta_single(args.wt_fasta)
    mut = read_fasta_single(args.mut_fasta)
    mpos, all_diffs = find_mutation(wt, mut)
    print(f'[agent] 突变定位: 1-based {mpos} ({wt[mpos-1]}>{mut[mpos-1]}), '
          f'转录本长度 WT={len(wt)} MUT={len(mut)}')

    # 模块一靶标解析: 突变类型判定(最长 ORF 法, 规则显式声明于
    # crrna_mutation_typing.py; 判读文字由数值按规则生成)
    typing = mtyping.classify(wt, mut)
    print(f'[agent] 突变类型判定(最长 ORF 法): {typing["reading"]}')
    # 特征向量: 丰度维(--target-key, 任务⑧) + 突变类型维(任务⑩)合并写入每行
    # target_features; 两者皆无时保持 None(与引入前口径一致)
    tfeat_row = dict(tfeat) if tfeat else {}
    tfeat_row.update(mtyping.feature_entry(typing))
    tfeat_row = tfeat_row or None

    cands = design_spacers(wt, mut, mpos, args.spacer_len, pfs_rule=args.pfs,
                           pfs_tol=args.pfs_tol, effector=args.effector,
                           topn=args.topn)
    print(f'[agent] tilling 候选 {len(cands)} 条, 取前 {args.topn} 条配骨架')
    out_rows = []
    for c in cands[:args.topn]:
        sel = select_scaffold_type(c['spacer'], dr, model)
        sf = _scaffold_fields(model, sel, degraded, uni)
        construct = core.to_dna(sf['dr_dna'] + c['spacer'])
        row = {'spacer_dna': core.to_dna(c['spacer']), 'mut_in': c['mut_in'],
               'proto_mm': c['proto_mm'], 'pfs_mut': c['pfs_mut'],
               'pfs_wt': c.get('pfs_wt', ''), 'tier': c['tier'],
               'self_mfe': c['self_mfe'],
               'scaffold_type': sf['scaffold_type'], 'confidence': sel['confidence'],
               'type_features': sel['features'],
               'target_features': tfeat_row,
               'dr_dna': sf['dr_dna'],
               'dr_desc': sf['dr_desc'],
               'construct_dna': construct,
               'wt_dr_construct_dna': core.to_dna(sf['wt_dr'] + c['spacer'])}
        if degraded:
            row['typing_nearest_type'] = sf['typing_nearest_type']
        out_rows.append(row)
        if degraded:
            print(f"  [{c['mut_in']:>14}] tier={c['tier']} self_mfe={c['self_mfe']:>6} "
                  f"-> 退化模式 通用骨架={uni['desc']} "
                  f"(分型诊断: 最近型{sel['nearest_type']}({sel['confidence']}), 不参与推荐)")
        else:
            print(f"  [{c['mut_in']:>14}] tier={c['tier']} self_mfe={c['self_mfe']:>6} "
                  f"-> 型{sel['nearest_type']}({sel['confidence']}) 骨架={sf['dr_desc']}")

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    pfs_out = resolve_pfs(args.effector, consensus=args.pfs,
                          tolerant_mismatches=args.pfs_tol)
    with open(args.out_prefix + '.design.json', 'w', encoding='utf-8') as fh:
        json.dump({'mutation': {'pos': mpos, 'wt_base': wt[mpos - 1], 'mut_base': mut[mpos - 1],
                                'all_diffs_1based': [d + 1 for d in all_diffs],
                                'typing': typing},
                   'effector': args.effector,
                   'target_context': tctx,
                   'degradation': assess,
                   'pfs': {'consensus': pfs_out['consensus'],
                           'tolerant_mismatches': pfs_out['tolerant_mismatches'],
                           'source': pfs_out['source']},
                   'n_candidates': len(cands), 'designs': out_rows,
                   'disclaimer': '排序为透明启发式(PFS 鉴别>单错配>结构代理), 非活性预测; '
                                 '骨架分型功能差异未实验验证; 选择器为最近邻规则'},
                  fh, ensure_ascii=False, indent=1)
    with open(args.out_prefix + '.design.fasta', 'w', encoding='utf-8') as fh:
        for i, r in enumerate(out_rows, 1):
            tag = 'universal' if degraded else f"type{r['scaffold_type']}"
            fh.write(f">agent_{i}|{r['mut_in']}|{tag}|{r['dr_desc']}|"
                     f"{args.effector}\n{r['construct_dna']}\n")
    print(f"[agent] 输出 -> {args.out_prefix}.design.json / .design.fasta")


if __name__ == '__main__':
    main()
