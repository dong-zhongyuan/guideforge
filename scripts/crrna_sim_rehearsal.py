# -*- coding: utf-8 -*-
"""湿实验分析链全流程彩排(2026-09-16, SIMULATED 数据, 先于真实数据)。

目的: 在真实细胞数据回填前, 用明确标注 SIMULATED 的模拟数据把
  --anova / --twin-check --cell-csv / --power / --ingest-cell
四段分析链完整跑通, 暴露接口/格式/逻辑问题, 并预生成报告骨架。

模拟情景(登记即情景, 非预测): H1 成立(A1C ≈ WT) + H2 阳性
(A8C+U15G 每靶 +12pp 杀伤), 噪声 sd=8pp, n=5 重复, 种子=42。
基线: HCT116/R248Q WT 存活 70.4%, SW480/R273H 45.0%(虚拟细胞 Scholz 标定中位)。

输出: data/sim_rehearsal/(模板CSV/SI表/各段stdout/汇总md), 全部带 SIMULATED 水印。
运行: PYTHONUTF8=1 python scripts/crrna_sim_rehearsal.py
"""
import csv
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data')
OUTD = os.path.join(DATA, 'sim_rehearsal')
PY = sys.executable

SCENARIO = {  # (target, scaffold) -> 存活中位%
    ('TP53-R248Q', 'WT'): 70.4, ('TP53-R248Q', 'A1C'): 69.0,
    ('TP53-R248Q', 'A8C+U15G'): 58.0,
    ('TP53-R273H', 'WT'): 45.0, ('TP53-R273H', 'A1C'): 44.0,
    ('TP53-R273H', 'A8C+U15G'): 33.0,
}
CELL = {'TP53-R248Q': 'HCT116(内源 R248Q 杂合)',
        'TP53-R273H': 'SW480(内源 R273H)'}


def main():
    os.makedirs(OUTD, exist_ok=True)
    rng = random.Random(42)

    # 1) 模拟回填 CSV(v6 模板格式, n=5)
    tpl = os.path.join(OUTD, 'sim_results_cell6.csv')
    with open(tpl, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['target', 'scaffold_desc'] +
                   ['survival_%d' % i for i in range(1, 6)] + ['note'])
        for (tgt, sca), mid in SCENARIO.items():
            vals = [round(max(0, min(100, rng.gauss(mid, 8))), 1)
                    for _ in range(5)]
            w.writerow([tgt, sca] + vals +
                       ['SIMULATED-REHEARSAL | %s | 情景: H1成立+H2阳性(+12pp), sd=8, seed=42'
                        % CELL[tgt]])

    # 2) 模拟 SI 表(scaffold_desc,target,SI_measured)
    si = os.path.join(OUTD, 'sim_si.csv')
    with open(si, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['scaffold_desc', 'target', 'SI_measured'])
        for (tgt, sca), mid in SCENARIO.items():
            wt_mid = SCENARIO[(tgt, 'WT')]
            si_v = round((100 - mid) / max(100 - wt_mid, 1e-9), 3)
            w.writerow([sca, tgt, si_v])

    # 3) 依次跑四段分析链, 收集 stdout
    runs = {
        'anova': [PY, os.path.join(HERE, 'crrna_ivt_template.py'),
                  '--anova', tpl],
        'twin': [PY, os.path.join(HERE, 'crrna_ivt_template.py'),
                 '--twin-check', tpl, '--cell-csv', si],
        'power': [PY, os.path.join(HERE, 'crrna_ivt_template.py'), '--power'],
        'bayes': [PY, os.path.join(HERE, 'crrna_bayesopt.py'),
                  '--ingest-cell', tpl, '--prior-han'],
    }
    results = {}
    for name, cmd in runs.items():
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        results[name] = {'rc': r.returncode,
                         'stdout': r.stdout[-3000:], 'stderr': r.stderr[-800:]}
        print('[rehearsal] %-6s rc=%d' % (name, r.returncode))
        if r.returncode != 0:
            print('   stderr:', r.stderr[-300:])

    # 4) 汇总报告
    md = ['# 湿实验分析链彩排报告(数据全部 SIMULATED, 不含任何真实实验数值)\n',
          '- 生成: scripts/crrna_sim_rehearsal.py, 2026-09-16, seed=42',
          '- 情景: H1 成立(A1C≈WT) + H2 阳性(A8C+U15G 每靶 +12pp), 噪声 sd=8pp, n=5',
          '- 基线: HCT116/R248Q WT 70.4% / SW480/R273H 45.0%(Scholz 标定中位)\n']
    for name, r in results.items():
        md.append('## %s (rc=%d)\n' % (name, r['rc']))
        md.append('```\n%s\n```\n' % r['stdout'].strip()[:1800])
    with open(os.path.join(OUTD, 'REHEARSAL_REPORT.md'), 'w',
              encoding='utf-8') as f:
        f.write('\n'.join(md))

    summary = {'generated_by': 'scripts/crrna_sim_rehearsal.py',
               'SIMULATED': True, 'scenario': 'H1+H2-positive, sd=8, n=5, seed=42',
               'baselines': {'TP53-R248Q_WT': 70.4, 'TP53-R273H_WT': 45.0},
               'runs': {k: v['rc'] for k, v in results.items()}}
    json.dump(summary, open(os.path.join(OUTD, 'summary.json'), 'w',
                            encoding='utf-8'), ensure_ascii=False, indent=1)
    ok = all(v['rc'] == 0 for v in results.values())
    print('[rehearsal] 全链 rc:', {k: v['rc'] for k, v in results.items()},
          '| ALL PASS' if ok else '| 有失败段, 见报告')
    print('->', OUTD)


if __name__ == '__main__':
    main()
