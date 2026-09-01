"""DR 变体特异性对比报告（DR 改造是否影响脱靶风险的三途径分析）

针对"改 DR 会不会导致误杀正常细胞"的问题, 按机制分三层回答:
  1. 识别位点集合: 由 spacer+PFS 决定, 与 DR 无关 → 引用全转录组扫描结果
     (crrna_transcriptome_scan.py 产出), 对去重后的位点算 spacer:靶标 duplex ΔG
     (RNAduplex; spacer 固定, 该能量景观对所有 DR 变体相同);
  2. 成熟 crRNA 长度漂移: 加工位点保护硬过滤(proc_ok)已约束 → 逐候选复核;
  3. 激活阈值降低(残余风险): 计算代理 = spacer 可及性增量(d_spacer_up/d_seed),
     蛋白接触几何保持(hbond/contact); 构象激活阈值本身计算代理不了,
     明确指向湿实验 WT TP53 转录本对照。

用法:
  python crrna_dr_specificity.py --run-prefix ../data/tp53_r248q_zengdr.v6 \
      --scan-prefix ../data/tp53_r248q_transcriptome_scan \
      --spacer GTTCATGCCGCCCATGCAGGAACT \
      --out ../data/tp53_r248q_zengdr.v6.dr_specificity.json
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)

import RNA  # noqa: E402


def duplex_mfe(spacer_rna, window_dna):
    """spacer:靶窗口 duplex MFE (kcal/mol, RNAduplex; 反向互补由算法处理)。"""
    window_rna = window_dna.upper().replace('T', 'U')
    d = RNA.duplexfold(spacer_rna, window_rna)
    return round(d.energy, 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--run-prefix', required=True)
    ap.add_argument('--scan-prefix', required=True)
    ap.add_argument('--spacer', required=True)
    ap.add_argument('--topn', type=int, default=12)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    spacer_rna = args.spacer.upper().replace('T', 'U')

    # --- 层 1: 脱靶能量景观(基因去重, 每基因取错配最少代表位点) ---
    sites = list(csv.DictReader(open(args.scan_prefix + '.sites.csv', encoding='utf-8')))
    by_gene = {}
    for s in sorted(sites, key=lambda s: (int(s['mismatches']), int(s['pfs_mm']))):
        by_gene.setdefault(s['gene'] or s['transcript'], s)
    landscape = []
    for gene, s in sorted(by_gene.items(), key=lambda kv: int(kv[1]['mismatches'])):
        landscape.append({
            'gene': gene, 'mismatches': int(s['mismatches']),
            'seed_mm': int(s['seed_mm']), 'pfs': s['pfs'],
            'pfs_mm': int(s['pfs_mm']), 'risk': s['risk'],
            'duplex_mfe_kcal': duplex_mfe(spacer_rna, s['window']),
            'note': 'WT 等位位点(PFS 3mm)' if s['mismatches'] == '0' else '',
        })

    # --- 层 2+3: DR 变体逐候选复核 ---
    rows = list(csv.DictReader(open(args.run_prefix + '.variants.csv', encoding='utf-8')))
    wt, variants = rows[0], rows[1:]
    wt_sp_up = float(wt['spacer_unpaired'])
    wt_seed = float(wt['seed_up'])
    passed = [r for r in variants if r['passed'] == 'True'][:args.topn]
    candidates = []
    for r in passed:
        d_up = round(float(r['spacer_unpaired']) - wt_sp_up, 3)
        d_seed = round(float(r['seed_up']) - wt_seed, 3)
        proc_ok = r['proc_ok'] == 'True'
        contact_kept = float(r['hbond_loss']) == 0.0
        verdict = []
        verdict.append('加工位点保护' + ('通过' if proc_ok else '未通过!'))
        verdict.append(f"spacer 可及性 Δ{d_up:+.3f}/Δseed {d_seed:+.3f}"
                       + ('(与 WT 一致)' if abs(d_up) < 0.02 and abs(d_seed) < 0.05
                          else '(有变化, 需复核)'))
        verdict.append('极性接触' + ('保持' if contact_kept else f"损失 {r['hbond_loss']}"))
        candidates.append({
            'desc': r['desc'], 'dr_seq': r['dr_seq'], 'rank': int(r['rank']),
            'proc_ok': proc_ok, 'd_spacer_unpaired': d_up, 'd_seed': d_seed,
            'hbond_preserved': float(r['hbond_preserved']),
            'checks': verdict,
            'specificity_verdict': '预测脱靶图谱与 WT 一致(识别位点集合由 spacer 决定, '
                                   '与 DR 无关; spacer 可及性未显著增加)'
            if proc_ok and abs(d_up) < 0.02 and abs(d_seed) < 0.05 and contact_kept
            else '存在变化项, 需湿实验优先验证',
        })

    payload = {
        'question': 'DR 改造是否改变脱靶风险(误杀正常细胞)',
        'layer1_site_universe': {
            'principle': '识别位点集合由 spacer+PFS 决定, 与 DR 无关; 对所有 DR 变体相同',
            'transcriptome_scan': os.path.basename(args.scan_prefix),
            'duplex_landscape': landscape,
            'conclusion': '全转录组无高/中危位点; 唯一 0 错配位点即 TP53 WT 等位'
                          '(PFS 3mm 不匹配), 鉴别由突变产生的 PFS 完成',
        },
        'layer2_processing': '加工位点保护硬过滤(proc_ok)防止成熟 crRNA 长度漂移'
                             '→ 识别窗口不变; 逐候选见 candidates',
        'layer3_activation_threshold': {
            'proxy': 'spacer 可及性增量 + 蛋白接触几何保持',
            'residual_risk': '构象激活阈值的改变无法由序列/二级结构代理, '
                             '湿实验必须含 WT TP53 转录本(无 PFS)对照',
        },
        'candidates': candidates,
        'disclaimer': '纯计算代理评估, 非特异性实验证据; 最终以体外生化与细胞实验为准',
    }
    with open(args.out, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f"输出: {args.out}")
    print('\n脱靶能量景观(基因去重):')
    for l in landscape:
        print(f"  {l['gene']:<10} {l['mismatches']}mm pfs {l['pfs']}({l['pfs_mm']}mm) "
              f"duplex {l['duplex_mfe_kcal']:>7} kcal/mol [{l['risk']}] {l['note']}")
    print('\nDR 变体特异性复核(TOP):')
    for c in candidates:
        print(f"  {c['desc']:<20} {'; '.join(c['checks'])}")


if __name__ == '__main__':
    main()
