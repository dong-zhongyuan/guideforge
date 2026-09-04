# -*- coding: utf-8 -*-
"""注册表基线漂移测试（2026-09 round-3 R1 配套）。

对每个注册条目，用 scripts/register_scaffold.vienna_baseline 现算
ViennaRNA 侧基线，与 configs/scaffold_registry.json 中存储值比较：
容差内通过，漂移即失败。同时断言条目记录了 vienna_version 且与
当前 ViennaRNA 版本一致（版本升级会触发重算复核，而非静默漂移）。
"""
import json
import os
import sys
import unittest

import RNA

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

import register_scaffold  # noqa: E402

REGISTRY = os.path.join(ROOT, 'configs', 'scaffold_registry.json')

# 存储值按 round(x, 2)/round(x, 3) 落盘，容差取略大于舍入误差
TOL_2DP = 0.02
TOL_3DP = 0.005


class TestRegistryDrift(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(REGISTRY, encoding='utf-8') as f:
            cls.entries = json.load(f)['entries']

    def test_every_entry_records_vienna_version(self):
        for name, e in self.entries.items():
            with self.subTest(effector=name):
                self.assertIn('vienna_version', e['reference'],
                              f'{name} 缺 vienna_version 记录')
                self.assertEqual(e['reference']['vienna_version'], RNA.__version__,
                                 f'{name} 注册 ViennaRNA 版本与当前环境不一致，需重注册复核')

    def test_baselines_recompute_without_drift(self):
        for name, e in self.entries.items():
            with self.subTest(effector=name):
                scaffold = e['scaffold']
                spacer = e['placeholder_spacer']
                full_dna = (spacer + scaffold if e['scaffold_side'] == '3prime'
                            else scaffold + spacer)
                full_rna = full_dna.replace('T', 'U')
                base = register_scaffold.vienna_baseline(
                    full_rna, len(spacer), e['scaffold_side'])
                stored = e['reference']['baseline']
                for key, tol in (('vienna_mfe_kcal', TOL_2DP),
                                 ('vienna_ens_diversity', TOL_2DP),
                                 ('vienna_centroid_dist', TOL_2DP),
                                 ('spacer_mean_unpaired', TOL_3DP)):
                    self.assertIn(key, stored, f'{name} 缺基线字段 {key}')
                    self.assertAlmostEqual(
                        stored[key], base[key], delta=tol,
                        msg=(f'{name}.{key} 漂移: 注册值 {stored[key]} '
                             f'vs 现算 {base[key]}'))


if __name__ == '__main__':
    unittest.main()
