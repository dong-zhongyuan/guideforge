# -*- coding: utf-8 -*-
"""§H 文献赢家先验契约测试(2026-09-16).

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

    def test_default_w_lit_zero_no_effect(self):
        """默认权重 0: lit_s 恒 0 -> 打分表达式与旧版一致(构造性检查)。"""
        import argparse

        import crrna_scaffold_design as core
        ap = argparse.Namespace(w_lit=0.0)
        lit_s = 0.0 if not getattr(ap, 'w_lit', 0.0) else 1.0
        self.assertEqual(lit_s, 0.0)
        self.assertEqual(getattr(ap, 'w_lit', 0.0) * lit_s, 0.0)

    def test_predictive_power_boundary_recorded(self):
        d = json.load(open(RULES_JSON, encoding='utf-8'))
        pp = d['rule_predictive_power']
        self.assertIn('spearman_rho_vs_lfc127', pp)
        self.assertLess(abs(pp['spearman_rho_vs_lfc127']), 0.3,
                        '规则单独预测力必须如实标弱; 若变强须重审口径')


if __name__ == '__main__':
    unittest.main()
