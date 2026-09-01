"""crRNA 骨架定向优化 v2（Cas12a2 DR 增强设计 · 奖赏制打分 · 独立脚本不改动 v1）

定位:
  v1（crrna_scaffold_design.py）是纯罚分制"结构保持"管线，产出与 WT 高度相似的保守变体。
  本脚本面向另一个问题: DR 改变能否**增强** Cas12a2 杀伤力。打分改为奖赏制——
  变体可因定向改进得正分，TOP-K 中每个候选自带"预期优势"论述。

定向优化轴（文献依据见 README）:
  1. 释放 spacer 3' 端种子区: Cas12a2 种子区位于 spacer 3' 端 7nt
     (Bravo 2023 / Dmytrenko 2023 结构), WT 全长折叠中该区被配对困住
     (spacer_mean_unpaired ~0.43)。奖赏 seed 区未配对概率提升。
  2. DR 茎环正确预折叠: 预折叠的 crRNA 茎环诱导 Cas12a 构象变化、驱动 RNP 组装
     (Sudhakar 2023, PMC10200996; 因果方向为 RNA→蛋白, 骨架正确预折叠是组装前提,
     而非蛋白稳定 RNA 折叠)。奖赏 DR 单独折叠 ΔG 更负、系综多样性更低
     (折叠正确且均一, 而非无约束地追求稳定)。
  3. 加工位点保护(硬过滤): Cas12a2 在茎环基部切割 pre-crRNA
     (Dmytrenko 2023 ED Fig.3), 变体 DR 单独折叠必须保持最外层茎配对
     （默认 2 对, --protect-pairs 可调）, 否则成熟 crRNA 长度可能改变。

策略:
  A. 单点全扫描 + 抽样双点（复用 v1 strategy_a）
  C. 茎部共变枚举: 对 WT DR 茎区每个配对枚举规范碱基对替换
     (AU/GC/GU 及反向), 含单配对替换与抽样多配对联合替换。

打分（奖赏 + 罚分混合; 权重 CLI 可调）:
  score = + w_seed*(seed_unpaired - wt_seed_unpaired)      # 种子区释放奖赏
          + w_free*(spacer_unpaired - wt_spacer_unpaired)  # 整条 spacer 游离奖赏
          + w_stab*max(-ddG_dr, 0)                          # DR 茎稳定化奖赏
          - w_bp*bp_dist - w_ddg*max(ddG,0)                 # 结构/能量风险罚项
          - w_contact*contact_sum - w_ens*max(d_ens,0) - w_hbond*hbond_loss

输出: <out>.variants.csv / <out>.top.json / <out>.top.fasta，另对每个通过变体
生成 advantage 字段（优势论述）。Pareto 最优标记见 pareto 列（目标: seed 游离度↑、
DR 稳定性↑、结构距离↓、接触+氢键损失↓）。

用法:
  python crrna_directed_design.py --effector cas12a2_zeng2026 \
      --spacer GTTCATGCCGCCCATGCAGGAACT --topk 12 --out-prefix ../data/run_v2
"""
import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np
import RNA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crrna_scaffold_design import (  # noqa: E402  复用 v1 基础函数, 不改动原文件
    BASES, ROOT, fold, load_contacts, hbond_loss, contact_sum, mut_desc,
    strategy_a, to_rna, to_dna,
)
from scaffold_registry import get_entry  # noqa: E402

# Cas12a2 种子区: spacer 3' 端 7nt (Bravo 2023 Fig.1d, 溶剂暴露的靶 RNA 结合种子)
SEED_NT = 7
# 规范/摆动碱基对(共变替换候选)
CANON_PAIRS = [('G', 'C'), ('C', 'G'), ('A', 'U'), ('U', 'A'), ('G', 'U'), ('U', 'G')]


def pf_unpaired(seq):
    """配分函数逐位点未配对概率 + 系综多样性。返回 (unpaired[nt], diversity)。"""
    fc = RNA.fold_compound(seq)
    fc.pf()
    n = len(seq)
    P = np.array([list(r) for r in fc.bpp()])[:n + 1, :n + 1]
    paired = P.sum(axis=0) + P.sum(axis=1)
    return 1.0 - paired[1:], float(fc.mean_bp_distance())


def db_pairs(db):
    """dot-bracket → 配对列表 [(i, j), ...]（0-based, i<j）。仅处理 () 层。"""
    stack, pairs = [], []
    for i, c in enumerate(db):
        if c == '(':
            stack.append(i)
        elif c == ')':
            pairs.append((stack.pop(), i))
    return pairs


def protected_pairs(dr_struct, n_pairs):
    """最外层 n_pairs 对茎配对（加工位点位于茎环基部, Dmytrenko 2023 ED Fig.3）。"""
    pairs = db_pairs(dr_struct)
    pairs.sort(key=lambda p: p[1] - p[0], reverse=True)
    return sorted(pairs[:n_pairs])


def strategy_c(dr, wt_dr_struct, args, rng):
    """策略C: 茎部共变枚举。单配对全枚举 + 多配对联合抽样。"""
    pairs = db_pairs(wt_dr_struct)
    variants = set()
    for i, j in pairs:
        wt_pair = (dr[i], dr[j])
        for a, b in CANON_PAIRS:
            if (a, b) == wt_pair:
                continue
            seq = dr[:i] + a + dr[i + 1:j] + b + dr[j + 1:]
            variants.add(seq)
    # 多配对联合替换(抽样): 体现"同结构、不同序列"的实质改造
    n_pairs = len(pairs)
    for _ in range(args.n_joint):
        k = int(rng.integers(2, n_pairs + 1))
        idx = rng.choice(n_pairs, size=k, replace=False)
        seq = list(dr)
        changed = False
        for t in idx:
            i, j = pairs[t]
            a, b = CANON_PAIRS[rng.integers(len(CANON_PAIRS))]
            if (a, b) != (dr[i], dr[j]):
                changed = True
            seq[i], seq[j] = a, b
        if changed:
            variants.add(''.join(seq))
    variants.discard(dr)
    return sorted(variants)


def score_v2(dr, seq, wt, spacer, args, contact):
    """奖赏制打分。返回指标 dict（含 pass/score/advantage/pareto 占位）。"""
    full = seq + spacer
    ss, mfe = fold(full)
    bp_dist = RNA.bp_distance(wt['mfe_struct'], ss)
    unpaired, diversity = pf_unpaired(full)
    spacer_up = float(unpaired[len(dr):].mean())
    seed_up = float(unpaired[-min(SEED_NT, len(spacer)):].mean())
    dr_ss, dr_mfe = fold(seq)
    ddg = mfe - wt['mfe_kcal']
    ddg_dr = dr_mfe - wt['dr_only_mfe_kcal']
    d_ens = diversity - wt['ens_diversity']
    d_seed = seed_up - wt['seed_unpaired']
    d_free = spacer_up - wt['spacer_mean_unpaired']
    desc, mut_pos = mut_desc(dr, seq)
    csum = contact_sum(mut_pos, contact)
    hloss, hfrac = hbond_loss(dr, seq, mut_pos, contact)
    # 加工位点保护: 变体 DR 单独折叠须保持最外层茎配对
    var_pairs = set(db_pairs(dr_ss))
    proc_ok = all(p in var_pairs for p in wt['protected'])
    ok = (bp_dist <= args.max_bp_dist
          and spacer_up >= wt['spacer_mean_unpaired'] - args.spacer_unpaired_margin
          and proc_ok
          and 'TTTT' not in to_dna(seq) and 'GGGG' not in to_dna(seq))
    score = (args.w_seed * d_seed + args.w_free * d_free
             + args.w_stab * max(-ddg_dr, 0.0)
             - args.w_bp * bp_dist - args.w_ddg * max(ddg, 0.0)
             - args.w_contact * csum - args.w_ens * max(d_ens, 0.0)
             - args.w_hbond * hloss)
    adv = []
    if d_seed > 0.02:
        adv.append(f"种子区(3'端{SEED_NT}nt)游离度 {wt['seed_unpaired']:.2f}→{seed_up:.2f}")
    if d_free > 0.02:
        adv.append(f"spacer 游离度 {wt['spacer_mean_unpaired']:.2f}→{spacer_up:.2f}")
    if ddg_dr < -0.1:
        adv.append(f"DR 茎稳定化 {ddg_dr:+.1f} kcal/mol")
    if d_ens < -0.05:
        adv.append(f"折叠均一性提升(d_ens {d_ens:+.2f})")
    if bp_dist == 0 and hloss == 0.0:
        adv.append("结构与蛋白极性接触完全保持")
    return {'desc': desc, 'mut_positions': mut_pos, 'n_mut': len(mut_pos),
            'dr_seq': to_dna(seq), 'construct_dna': to_dna(full),
            'mfe_struct': ss, 'mfe_kcal': round(mfe, 2), 'ddG': round(ddg, 2),
            'dr_only_mfe_kcal': round(dr_mfe, 2), 'ddG_dr': round(ddg_dr, 2),
            'bp_dist': bp_dist, 'spacer_unpaired': round(spacer_up, 3),
            'seed_unpaired': round(seed_up, 3), 'd_seed': round(d_seed, 3),
            'ens_diversity': round(diversity, 2), 'd_ens': round(d_ens, 2),
            'contact_sum': round(csum, 3),
            'hbond_loss': round(hloss, 3), 'hbond_preserved': round(hfrac, 3),
            'proc_site_ok': bool(proc_ok),
            'advantage': '; '.join(adv) if adv else '与 WT 无显著差异',
            'pareto': None,
            'passed': bool(ok), 'score': round(score, 4)}


def pareto_flag(rows):
    """Pareto 最优标记。目标: seed 游离↑, DR 稳定(-ddG_dr)↑, bp_dist↓, 接触+氢键损失↓。"""
    passed = [r for r in rows if r['passed']]
    if not passed:
        return
    obj = np.array([[r['seed_unpaired'], -r['ddG_dr'],
                     -r['bp_dist'], -(r['contact_sum'] + r['hbond_loss'])]
                    for r in passed])
    for k, r in enumerate(passed):
        dominated = np.any(np.all(obj >= obj[k], axis=1)
                           & np.any(obj > obj[k], axis=1))
        r['pareto'] = not dominated


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacer', required=True, help='固定 spacer 序列(17-25nt ACGT/U), 全程不变')
    ap.add_argument('--effector', default='cas12a2_zeng2026')
    ap.add_argument('--registry', default=None, help='注册表路径(默认 configs/scaffold_registry.json)')
    ap.add_argument('--n-double', type=int, default=200, help='策略A 双点突变抽样数')
    ap.add_argument('--n-joint', type=int, default=100, help='策略C 多配对联合替换抽样数')
    ap.add_argument('--topk', type=int, default=12)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--max-bp-dist', type=float, default=4.0)
    ap.add_argument('--spacer-unpaired-margin', type=float, default=0.10)
    ap.add_argument('--protect-pairs', type=int, default=2,
                    help='加工位点保护的最外层茎配对数(Dmytrenko 2023 切割位于茎环基部)')
    ap.add_argument('--w-seed', type=float, default=1.0, help="种子区(spacer 3'端7nt)释放奖赏权重")
    ap.add_argument('--w-free', type=float, default=0.5, help='整条 spacer 游离奖赏权重')
    ap.add_argument('--w-stab', type=float, default=0.1, help='DR 茎稳定化奖赏权重(kcal/mol)')
    ap.add_argument('--w-bp', type=float, default=0.3)
    ap.add_argument('--w-ddg', type=float, default=0.1)
    ap.add_argument('--w-contact', type=float, default=1.0)
    ap.add_argument('--w-ens', type=float, default=0.05)
    ap.add_argument('--w-hbond', type=float, default=0.5)
    ap.add_argument('--no-contacts', action='store_true', help='关闭蛋白接触项(不推荐)')
    ap.add_argument('--out-prefix', default='crrna_directed_run')
    args = ap.parse_args()

    entry = get_entry(args.effector, args.registry)
    dr = to_rna(entry['scaffold'])
    spacer = to_rna(args.spacer)
    if not 17 <= len(spacer) <= 25 or set(spacer) - set(BASES):
        ap.error('--spacer 必须为 17-25nt ACGT/U')
    if entry.get('scaffold_side') != '5prime':
        ap.error(f"效应子 {args.effector} 骨架方位 {entry.get('scaffold_side')} 非 5prime, 暂不支持")
    rng = np.random.default_rng(args.seed)

    # WT 参考(含种子区与 DR 单独折叠、受保护配对)
    full = dr + spacer
    ss, mfe = fold(full)
    unpaired, diversity = pf_unpaired(full)
    dr_ss, dr_mfe = fold(dr)
    wt = {'full_len': len(full), 'mfe_struct': ss, 'mfe_kcal': round(mfe, 2),
          'ens_diversity': round(diversity, 2),
          'spacer_mean_unpaired': round(float(unpaired[len(dr):].mean()), 3),
          'seed_unpaired': round(float(unpaired[-SEED_NT:].mean()), 3),
          'dr_only_struct': dr_ss, 'dr_only_mfe_kcal': round(dr_mfe, 2),
          'protected': protected_pairs(dr_ss, args.protect_pairs)}
    print('=== WT 参考 (v2 奖赏制) ===')
    print(f"效应子: {entry['display_name']}")
    print(f"DR 骨架({len(dr)}nt): {to_dna(dr)}  spacer({len(spacer)}nt, 固定): {to_dna(spacer)}")
    print(f"全长 MFE: {wt['mfe_struct']}  {wt['mfe_kcal']} kcal/mol")
    print(f"DR 单独折叠: {wt['dr_only_struct']}  {wt['dr_only_mfe_kcal']} kcal/mol")
    print(f"spacer 平均未配对: {wt['spacer_mean_unpaired']}  "
          f"种子区(3'端{SEED_NT}nt)未配对: {wt['seed_unpaired']}")
    print(f"加工位点保护配对(0-based): {wt['protected']}")

    contact = None if args.no_contacts else load_contacts(dr, enabled=True)
    if contact is not None and contact.get('transfer'):
        t = contact['transfer']
        print(f"蛋白接触约束: {contact['source']} [对齐转移: 错配位点 "
              f"{t['mismatch_positions_1based']} 的极性接触已剔除]")
    else:
        print(f"蛋白接触约束: {'关闭' if args.no_contacts else contact['source']}")

    pool, source = {}, {}

    def _add(seq, tag):
        pool.setdefault(seq, None)
        source.setdefault(seq, set()).add(tag)

    for seq in strategy_a(dr, args, rng):
        _add(seq, 'A')
    c_variants = strategy_c(dr, wt['dr_only_struct'], args, rng)
    for seq in c_variants:
        _add(seq, 'C')
    print(f"策略A 突变扫描 + 策略C 茎部共变({len(c_variants)} 条)")

    rows = []
    for seq in pool:
        row = score_v2(dr, seq, wt, spacer, args, contact)
        row['source'] = '+'.join(sorted(source[seq]))
        rows.append(row)
    pareto_flag(rows)
    rows.sort(key=lambda r: (not r['passed'], -r['score']))
    for rank, row in enumerate(rows, 1):
        row['rank'] = rank

    wt_row = {'desc': 'WT', 'mut_positions': [], 'n_mut': 0, 'dr_seq': to_dna(dr),
              'construct_dna': to_dna(full), 'mfe_struct': wt['mfe_struct'],
              'mfe_kcal': wt['mfe_kcal'], 'ddG': 0.0,
              'dr_only_mfe_kcal': wt['dr_only_mfe_kcal'], 'ddG_dr': 0.0,
              'bp_dist': 0, 'spacer_unpaired': wt['spacer_mean_unpaired'],
              'seed_unpaired': wt['seed_unpaired'], 'd_seed': 0.0,
              'ens_diversity': wt['ens_diversity'], 'd_ens': 0.0,
              'contact_sum': 0.0, 'hbond_loss': 0.0, 'hbond_preserved': 1.0,
              'proc_site_ok': True, 'advantage': '参考株', 'pareto': None,
              'passed': True, 'score': 0.0, 'source': 'WT', 'rank': 0}

    n_pass = sum(r['passed'] for r in rows)
    n_pareto = sum(1 for r in rows if r['pareto'])
    print(f"\n=== 候选库 v2 ===  生成 {len(rows)} 条唯一变体, 通过硬过滤 {n_pass} 条, "
          f"Pareto 最优 {n_pareto} 条")
    print(f"{'rank':<5}{'desc':<22}{'ddG':>7}{'ddGdr':>7}{'bp_d':>5}{'seed':>6}"
          f"{'dseed':>7}{'hb':>6}{'P':>2}{'score':>8}  pass/src  advantage")
    for row in rows[:args.topk]:
        print(f"{row['rank']:<5}{row['desc']:<22}{row['ddG']:>7.2f}{row['ddG_dr']:>7.2f}"
              f"{row['bp_dist']:>5}{row['seed_unpaired']:>6.3f}{row['d_seed']:>7.3f}"
              f"{row['hbond_loss']:>6.2f}{'Y' if row['pareto'] else '-':>2}"
              f"{row['score']:>8.3f}  {'Y' if row['passed'] else 'N'}/{row['source']}"
              f"  {row['advantage']}")

    out_csv = args.out_prefix + '.variants.csv'
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(wt_row.keys()))
        writer.writeheader()
        writer.writerow(wt_row)
        writer.writerows(rows)
    top = [wt_row] + [r for r in rows if r['passed']][:args.topk]
    payload = {
        'disclaimer': '候选库压缩结果, 非活性预测; 活性/特异性以体外生化与细胞实验为准; '
                      '奖赏项为溶液态折叠预测, 结合态以 RNet SHAPE 一致性与实验为准',
        'scoring': 'v2 奖赏制: +种子区释放 +spacer游离 +DR茎稳定 -结构/接触/氢键罚项',
        'effector': args.effector, 'registry_entry_version': entry.get('version'),
        'spacer_fixed_dna': to_dna(spacer),
        'wt': {k: v for k, v in wt.items() if k != 'protected'},
        'params': {k: getattr(args, k) for k in
                   ('n_double', 'n_joint', 'seed', 'max_bp_dist',
                    'spacer_unpaired_margin', 'protect_pairs', 'w_seed', 'w_free',
                    'w_stab', 'w_bp', 'w_ddg', 'w_contact', 'w_ens', 'w_hbond',
                    'no_contacts', 'topk')},
        'n_variants': len(rows), 'n_passed': n_pass, 'n_pareto': n_pareto,
        'library_sha256': hashlib.sha256(
            json.dumps([r['dr_seq'] for r in rows], sort_keys=True).encode()).hexdigest(),
        'top': top,
    }
    out_json = args.out_prefix + '.top.json'
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    out_fasta = args.out_prefix + '.top.fasta'
    with open(out_fasta, 'w', encoding='utf-8') as fh:
        for row in top:
            fh.write(f">crRNA_v2_{row['rank']}_{row['desc']}|{args.effector}|spacer_fixed\n"
                     f"{row['construct_dna']}\n")
    print(f"\n输出: {out_csv} / {out_json} / {out_fasta}")
    print("说明: top.json 含 WT(第 0 名) + TOP-K 通过硬过滤变体; score>0 表示预测相对 WT "
          "有定向改进; advantage 列为该候选的优势论述。")


if __name__ == '__main__':
    main()
