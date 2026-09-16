# -*- coding: utf-8 -*-
"""§J 文献赢家先验契约测试(2026-09-16).

覆盖:
1) 规则文件存在且含核心规则(A8C/U15G, DeWeirdt 64 赢家共识对);
2) load_lit_rules 产出正 log-odds 键;
3) 默认 w_lit=0 时管线打分不含先验项(回归口径: 默认行为不变的构造性检查);
4) 预测力边界如实携带(rho 弱, 只作先验)。
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

RULES_JSON = os.path.join(ROOT, 'data', 'winner_rule_engineering.json')


class TestLitPrior(unittest.TestCase):

    def test_rules_file_and_core_pair(self):
        d = json.load(open(RULES_JSON, encoding='utf-8'))
        self.assertGreaterEqual(d['n_winners'], 60)
        keys = {(r['su_pos'], r['from'], r['to']) for r in d['rules_sig']}
        self.assertIn((8, 'A', 'C'), keys, 'A8C 核心赢家规则缺失')
        self.assertIn((15, 'U', 'G'), keys, 'U15G 核心赢家规则缺失')
        for r in d['rules_sig']:
            self.assertLess(r['fisher_p'], 0.05)

    def test_load_lit_rules_positive(self):
        import crrna_scaffold_design as core
        rules = core.load_lit_rules(RULES_JSON)
        self.assertIn('8AC', rules)
        self.assertIn('15UG', rules)
        for v in rules.values():
            self.assertGreater(v, 0.0)

    def test_w_lit_behavioral_delta(self):
        """行为级契约: --w-lit 1.0 时 A8C+U15G 的分数增量 = 1.0 x (规则 log-odds 和)。

        基线取仓库既有无先验产物 data/jd12_sp1_direct.top.json(确定性管线);
        权重线性 => 增量在 w=0 处恰为 0, 覆盖默认行为不变的回归口径。
        运行成本 ~25s(枚举主导), 已接受。
        """
        import json
        import math
        import subprocess
        import sys

        import crrna_scaffold_design as core
        rules = core.load_lit_rules(RULES_JSON)
        expect = rules['8AC'] + rules['15UG']
        out = os.path.join(ROOT, 'data', 'tmp_test_lit.top')
        r = subprocess.run(
            [sys.executable, os.path.join(ROOT, 'scripts',
                                          'crrna_scaffold_design.py'),
             '--effector', 'cas12a2_zeng2026',
             '--spacer', 'ACAGGCACAAACATGCACCTCAA', '--topk', '8',
             '--use-covariation', '--w-lit', '1.0',
             '--out-prefix', out[:-4]],
            capture_output=True, text=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-300:])
        lit = json.load(open(out + '.json', encoding='utf-8'))['top']
        base = json.load(open(os.path.join(
            ROOT, 'data', 'jd12_sp1_direct.top.json'),
            encoding='utf-8'))['top']
        lv = next(x for x in lit if x['desc'] == 'A8C+U15G')
        bv = next(x for x in base if x['desc'] == 'A8C+U15G')
        self.assertAlmostEqual(lv['score'] - bv['score'], expect, places=3,
                               msg='先验增量不等于规则 log-odds 和')
        for f in (out + '.json', out + '.variants.csv', out + '.fasta'):
            if os.path.exists(f):
                os.remove(f)

    def test_predictive_power_boundary_recorded(self):
        d = json.load(open(RULES_JSON, encoding='utf-8'))
        pp = d['rule_predictive_power']
        self.assertIn('spearman_rho_vs_lfc127', pp)
        self.assertLess(abs(pp['spearman_rho_vs_lfc127']), 0.3,
                        '规则单独预测力必须如实标弱; 若变强须重审口径')


if __name__ == '__main__':
    unittest.main()
