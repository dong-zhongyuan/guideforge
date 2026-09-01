"""crRNA 设计智能体（ spacer 设计 × 骨架分型选择 组合入口）

输入目标突变转录本(+WT 对照), 输出 spacer × 推荐骨架 的组合候选:
  1) spacer 设计器(A 部分): 围绕突变位点 tilling 枚举全部覆盖窗口
     (突变置于 PFS 或 protospacer 不同位置, Zeng 2026 口径), 按鉴别档位
     (PFS 鉴别 > protospacer 单错配 > 其它)与结构代理(spacer 游离度/自互补)排序;
  2) 骨架选择器(B 部分): 新 spacer 的 5 维上下文特征, 与既有分型
     (data/context_typing/typing.clusters.json, 3 型)的型心比欧氏距离,
     给最近型 + 置信度(最近/次近距离比); 边界样本如实报两型;
  3) 各型代表骨架从分型试验的管线产出中取(每型首个非 5' 悬垂惰性位变体 + WT 备选)。

边界(诚实声明): 排序为透明启发式, 非活性预测; 分型功能差异未实验验证;
骨架选择器为最近邻规则(样本量决定), 不是训练模型。

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
from scaffold_registry import get_scaffold  # noqa: E402

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


def design_spacers(wt, mut, mpos, spacer_len=24, pfs_rule='CAGAG', topn=10):
    """围绕突变 tilling: 突变可落在 PFS(窗口 3' 下游 1..5 位)或 protospacer 内。
    返回候选列表(按鉴别档位+结构代理排序)。"""
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
        cands.append({'spacer': spacer, 'mut_in': f'PFS+{k}', 'proto_mm': mm,
                      'pfs_mut': pfs_mut, 'pfs_wt': pfs_wt,
                      'pfs_discriminates': pfs_mut.startswith(pfs_rule[:1]) and pfs_mut != pfs_wt,
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mut-fasta', default=None, help='突变转录本 FASTA(与 --wt-fasta 一起用于 tilling 设计 A)')
    ap.add_argument('--wt-fasta', default=None, help='WT 转录本 FASTA')
    ap.add_argument('--spacer', default=None,
                    help='直接输入 spacer(DNA/RNA, 17-25nt); 提供时跳过 tilling 设计, 只做骨架分型选择(项目本体口径)')
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--spacer-len', type=int, default=24)
    ap.add_argument('--topn', type=int, default=10)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--spacers', required=True)
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    dr = core.to_rna(get_scaffold(args.effector))
    model = load_type_models(args.clusters, args.spacers, dr)
    print(f'[agent] 骨架分型: {len(model["centers"])} 型; 各型代表骨架 '
          f'{ {t: v["representative_desc"] for t, v in model["type_dr"].items()} }')

    if args.spacer:
        # 简版(项目本体口径): spacer 由实验侧输入, 只做骨架分型选择
        spacer_rna = core.to_rna(args.spacer)
        if not 17 <= len(spacer_rna) <= 25 or set(spacer_rna) - set('ACGU'):
            raise SystemExit('--spacer 必须为 17-25nt ACGT/U')
        sel = select_scaffold_type(spacer_rna, dr, model)
        dr_pick = model['type_dr'][sel['nearest_type']]
        out_rows = [{'spacer_dna': core.to_dna(spacer_rna), 'mut_in': 'input',
                     'proto_mm': None, 'pfs_mut': '', 'pfs_wt': '', 'tier': None,
                     'self_mfe': round(core.fold(spacer_rna)[1], 2),
                     'scaffold_type': sel['nearest_type'], 'confidence': sel['confidence'],
                     'type_features': sel['features'], 'distances': sel['distances'],
                     'dr_dna': dr_pick['representative_dr'],
                     'dr_desc': dr_pick['representative_desc'],
                     'construct_dna': core.to_dna(dr_pick['representative_dr'] + spacer_rna),
                     'wt_dr_construct_dna': core.to_dna(dr_pick['wt_dr'] + spacer_rna)}]
        print(f"[agent] 输入 spacer {core.to_dna(spacer_rna)} -> 型{sel['nearest_type']}"
              f"({sel['confidence']}) 推荐骨架 {dr_pick['representative_desc']}")
        os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
        with open(args.out_prefix + '.design.json', 'w', encoding='utf-8') as fh:
            json.dump({'mode': 'spacer-input(骨架分型选择)', 'effector': args.effector,
                       'designs': out_rows,
                       'disclaimer': '选择器为最近邻规则, 分型功能差异未实验验证'},
                      fh, ensure_ascii=False, indent=1)
        with open(args.out_prefix + '.design.fasta', 'w', encoding='utf-8') as fh:
            r = out_rows[0]
            fh.write(f">agent_input|type{r['scaffold_type']}|{r['dr_desc']}|{args.effector}\n"
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

    cands = design_spacers(wt, mut, mpos, args.spacer_len, topn=args.topn)
    print(f'[agent] tilling 候选 {len(cands)} 条, 取前 {args.topn} 条配骨架')
    out_rows = []
    for c in cands[:args.topn]:
        sel = select_scaffold_type(c['spacer'], dr, model)
        dr_pick = model['type_dr'][sel['nearest_type']]
        construct = core.to_dna(dr_pick['representative_dr'] + c['spacer'])
        out_rows.append({'spacer_dna': core.to_dna(c['spacer']), 'mut_in': c['mut_in'],
                         'proto_mm': c['proto_mm'], 'pfs_mut': c['pfs_mut'],
                         'pfs_wt': c.get('pfs_wt', ''), 'tier': c['tier'],
                         'self_mfe': c['self_mfe'],
                         'scaffold_type': sel['nearest_type'], 'confidence': sel['confidence'],
                         'type_features': sel['features'],
                         'dr_dna': dr_pick['representative_dr'],
                         'dr_desc': dr_pick['representative_desc'],
                         'construct_dna': construct,
                         'wt_dr_construct_dna': core.to_dna(dr_pick['wt_dr'] + c['spacer'])})
        print(f"  [{c['mut_in']:>14}] tier={c['tier']} self_mfe={c['self_mfe']:>6} "
              f"-> 型{sel['nearest_type']}({sel['confidence']}) 骨架={dr_pick['representative_desc']}")

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    with open(args.out_prefix + '.design.json', 'w', encoding='utf-8') as fh:
        json.dump({'mutation': {'pos': mpos, 'wt_base': wt[mpos - 1], 'mut_base': mut[mpos - 1],
                                'all_diffs_1based': [d + 1 for d in all_diffs]},
                   'effector': args.effector, 'n_candidates': len(cands), 'designs': out_rows,
                   'disclaimer': '排序为透明启发式(PFS 鉴别>单错配>结构代理), 非活性预测; '
                                 '骨架分型功能差异未实验验证; 选择器为最近邻规则'},
                  fh, ensure_ascii=False, indent=1)
    with open(args.out_prefix + '.design.fasta', 'w', encoding='utf-8') as fh:
        for i, r in enumerate(out_rows, 1):
            fh.write(f">agent_{i}|{r['mut_in']}|type{r['scaffold_type']}|{r['dr_desc']}|"
                     f"{args.effector}\n{r['construct_dna']}\n")
    print(f"[agent] 输出 -> {args.out_prefix}.design.json / .design.fasta")


if __name__ == '__main__':
    main()
