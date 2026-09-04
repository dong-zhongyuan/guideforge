# -*- coding: utf-8 -*-
"""R2 评审修复的回归测试(2026-09, round-2 R2 项)。

锁定两处统计口径:
  1. conservation_partition 的 loop 判据 = flanked-by-paired
     (连续未配对 run 且两侧紧邻位均配对); 5' 端悬垂不算 loop,
     即使它比真 loop 更长。
  2. _spearman 为 tie-aware midranks, 与 scipy.stats.spearmanr 一致
     (x 为 3 值离散变量, 并列大量存在, 手搓 argsort-of-argsort 不予平均)。
"""
import os
import sys
import unittest
from unittest import mock

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src"))

import crrna_dmytrenko_validate as dv  # noqa: E402


class TestConservationPartition(unittest.TestCase):
    def test_zeng2026_dr_loop_is_hairpin_not_5p_overhang(self):
        """真实 zeng2026 DR(19nt, fold ....(((((....))))).):
        5' 悬垂 1-4 与 loop 10-13 等长(4nt), loop 必须是 10-13。"""
        part, loop_idx = dv.conservation_partition("AAUUUCUACUGUUGUAGAU")
        self.assertEqual([i + 1 for i in loop_idx], [10, 11, 12, 13])
        for pos in (1, 2, 3, 4):
            self.assertEqual(part[pos], 1)          # 5' 悬垂 -> 茎/悬垂类
        for pos in (10, 11, 12, 13):
            self.assertEqual(part[pos], 0)          # loop 类
        for pos in (15, 16, 17, 18, 19):
            self.assertEqual(part[pos], 2)          # 3' 保守窗

    def _part_with_ss(self, ss, n):
        with mock.patch.object(dv.core, "fold", return_value=(ss, 0.0)):
            return dv.conservation_partition("A" * n, cons3_window=5)

    def test_longer_5p_overhang_still_not_loop(self):
        """5' 悬垂(6nt)比 loop(4nt)更长时, loop 仍取被夹住的那段。"""
        ss = "......((((....))))"                     # 18nt
        part, loop_idx = self._part_with_ss(ss, 18)
        self.assertEqual([i + 1 for i in loop_idx], [11, 12, 13, 14])
        self.assertTrue(all(part[p] == 1 for p in range(1, 7)))

    def test_3p_trailing_run_not_loop(self):
        ss = "(((((....)))))...."                     # 18nt
        _, loop_idx = self._part_with_ss(ss, 18)
        self.assertEqual([i + 1 for i in loop_idx], [6, 7, 8, 9])

    def test_no_flanked_run_means_no_loop(self):
        ss = "....(((())))"                           # 12nt, 仅 5' 悬垂
        _, loop_idx = self._part_with_ss(ss, 12)
        self.assertEqual(loop_idx, [])


class TestSpearmanTieAware(unittest.TestCase):
    def test_matches_scipy_on_tie_heavy_data(self):
        x = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2], float)
        y = np.array([1.0, 1.0, 0.5, 2.0, 1.0, 3.0, 3.0, 0.5, 2.0, 4.0,
                      1.0, 5.0, 5.0, 4.0, 0.5, 2.0])
        self.assertAlmostEqual(dv._spearman(x, y),
                               float(spearmanr(x, y).statistic), places=12)

    def test_differs_from_tie_blind_when_ties_present(self):
        """守护: 若退换回 argsort-of-argsort, 本断言会变红。"""
        x = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2], float)
        y = np.array([0.9, 0.1, 0.5, 0.4, 0.8, 0.2, 0.6, 0.3, 0.7, 1.0, 0.0, 0.55])
        naive = float(np.corrcoef(np.argsort(np.argsort(x)),
                                  np.argsort(np.argsort(y)))[0, 1])
        self.assertNotAlmostEqual(dv._spearman(x, y), naive, places=6)

    def test_no_ties_matches_naive(self):
        rng = np.random.default_rng(0)
        x, y = rng.normal(size=30), rng.normal(size=30)
        naive = float(np.corrcoef(np.argsort(np.argsort(x)),
                                  np.argsort(np.argsort(y)))[0, 1])
        self.assertAlmostEqual(dv._spearman(x, y), naive, places=12)


if __name__ == "__main__":
    unittest.main()
