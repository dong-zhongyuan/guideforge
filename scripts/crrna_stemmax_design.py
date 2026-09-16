# -*- coding: utf-8 -*-
"""茎最大化(stem-max)骨架强化设计(2026-09-16, 用户指令: 骨架想办法变强).

设计逻辑(全部本地数据, 无新外部依赖):
1. 我们的 DR(AAUUUCUACUGUUGUAGAU, Zeng/Su 注册表) MFE 结构给出茎配对;
2. 8D4A 接触图谱(data/8D4A_dr_contacts_4d017c23.json)给出每位置蛋白接触数
   —— 界面安全位 = 接触 <=3 个残基;
3. 候选生成: 仅动"茎配对中双端均界面安全"的配对, 弱对(A:U / G:U)换强对
   (G:C 或 C:G, 取方向使 GC 分布均衡); 环区与高接触位一律不动;
   另含已有实测茎稳定化参照(A8C+U15G, U13C);
4. 评估: ViennaRNA ddG / 茎完整概率 / bp_dist / cross_nt / spacer_up
   (双 canonical spacer 语境) + 选型器复刻预测(DT max_leaf=4, rs=42)
   + 界面风险代理(改动位接触数之和, Chai 复验列为待办);
5. 排序输出 top 候选 + 风险标注。

边界(如实): 选型器为 16 对同源小样本树, 预测分辨率粗; "强化"主张的
最终判据仍是湿实验(H2 扩展臂); 跨酶(Han=GeCas12a2)证据迁移未证。

输出: data/stemmax_design.json
运行: PYTHONUTF8=1 python scripts/crrna_stemmax_design.py
"""
import json
import os
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
WT_DR = 'AAUUUCUACUGUUGUAGAU'
SAFE_CONTACT_MAX = 6          # 界面安全位阈值(残基数, 分层: <=3 严格/<=6 中等)
CONTEXTS = {'R248Q': 'GTTCATGCCGCCCATGCAGGAACT',
            'R273H': 'CACCTCAAAGCTGTTCCGTCCCAG'}
REF_SCAFFOLDS = {'A8C+U15G': 'AAUUUCUCCUGUUGGAGAU',
                 'U13C': 'AAUUUCUACUGUCGUAGAU'}  # 单点 13U->C(修正早版手写串错误)


def contacts():
    d = json.load(open(os.path.join(DATA, '8D4A_dr_contacts_4d017c23.json'),
                       encoding='utf-8'))
    return {int(k): v['n_contact_residues'] for k, v in d['contacts'].items()}


def stem_pairs(dr):
    ss, _ = core.fold(dr)
    stack, pairs = [], []
    for i, ch in enumerate(ss):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            pairs.append((stack.pop() + 1, i + 1))  # 1-based
    return pairs


def strengthen(dr, pairs, safe):
    """弱对换强对: 返回 (新DR, 改动列表)。仅动双端安全且当前为 AU/UA/GU/UG 的对。"""
    chars = list(dr)
    edits = []
    for a, b in sorted(pairs):
        if safe.get(a, 99) > SAFE_CONTACT_MAX or safe.get(b, 99) > SAFE_CONTACT_MAX:
            continue
        cur = dr[a - 1] + dr[b - 1]
        if cur in ('AU', 'UA', 'GU', 'UG'):
            # 方向: 5'端给 G 配 3'端 C, 或反之; 取 GC 含量均衡(轮流)
            use = 'GC' if (len(edits) % 2 == 0) else 'CG'
            if chars[a - 1] + chars[b - 1] != use:
                edits.append((a, b, cur, use))
                chars[a - 1], chars[b - 1] = use[0], use[1]
    return ''.join(chars), edits


def evaluate(dr, spacer, wt_dr, wt_full_ss, wt_dr_mfe, wt_stem, dt):
    full = dr + core.to_rna(spacer)
    ss, _ = core.fold(full)
    _, dr_mfe = core.fold(dr)
    cross_nt, _ = core.cross_pairs(full, len(dr))
    p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
    _, sp_up, _ = core.pf_stats(full, len(dr), 7)
    feats = [round(dr_mfe - wt_dr_mfe, 2),
             RNA.bp_distance(wt_full_ss, ss), cross_nt,
             round(p_fold, 5), round(sp_up, 3)]
    return {'ddG_dr': feats[0], 'bp_dist': feats[1], 'cross_nt': feats[2],
            'p_fold': feats[3], 'spacer_up': feats[4],
            'pred_fig1g': round(float(dt.predict([feats])[0]), 4)}


def main():
    safe = contacts()
    pairs = stem_pairs(WT_DR)
    print('[stemmax] WT DR 茎配对:', pairs)
    print('[stemmax] 接触数:', {i: safe[i] for i in range(1, 20)})

    new_dr, edits = strengthen(WT_DR, pairs, safe)
    print('[stemmax] 安全弱对强化:', edits or '无')
    cands = {'WT': WT_DR, 'stemmax': new_dr}
    cands.update(REF_SCAFFOLDS)

    X, y, _, _, _ = sel.load_han()
    dt = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                               random_state=42)
    dt.fit(X, y)

    out = {'generated_by': 'scripts/crrna_stemmax_design.py',
           'wt_dr': WT_DR, 'stem_pairs': [list(p) for p in pairs],
           'safe_contact_max': SAFE_CONTACT_MAX,
           'contact_counts': safe, 'stemmax_edits': [
               {'pair': [a, b], 'from': cur, 'to': use}
               for a, b, cur, use in edits],
           'candidates': {}}
    for cname, cdr in cands.items():
        n_mut = sum(a != b for a, b in zip(WT_DR, cdr))
        mut_pos = [i + 1 for i, (a, b) in enumerate(zip(WT_DR, cdr)) if a != b]
        risk = sum(safe.get(p, 0) for p in mut_pos)
        rec = {'dr_rna': cdr, 'n_mut': n_mut, 'mut_positions': mut_pos,
               'mut_contact_sum': risk}
        for ctx, sp in CONTEXTS.items():
            wt_full_ss, _ = core.fold(WT_DR + core.to_rna(sp))
            _, wt_dr_mfe = core.fold(WT_DR)
            rec[ctx] = evaluate(cdr, sp, WT_DR, wt_full_ss, wt_dr_mfe,
                                core.stem_pairs_of(core.fold(WT_DR)[0]), dt)
        out['candidates'][cname] = rec
        r248, r273 = rec['R248Q'], rec['R273H']
        print('[stemmax] %-10s n_mut=%d risk=%2d ddG=%+.2f/%+.2f pred=%.4f/%.4f'
              % (cname, n_mut, risk, r248['ddG_dr'], r273['ddG_dr'],
                 r248['pred_fig1g'], r273['pred_fig1g']))

    dst = os.path.join(DATA, 'stemmax_design.json')
    json.dump(out, open(dst, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('->', dst)


if __name__ == '__main__':
    main()
