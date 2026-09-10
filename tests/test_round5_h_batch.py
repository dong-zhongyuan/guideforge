# -*- coding: utf-8 -*-
"""Round-4 §H 批次契约测试(2026-09-11, 判据 docs/preregistration.md §H)。

沿用 §G 批次惯例: 不锁死具体数值, 锁定
1. 各产物 JSON 的 §H 块结构齐全;
2. 判读文字与数值按预登记规则同向(机器校验, 不手写结论方向);
3. 关键口径字段存在(随机斜率/方向门控/逐对胜率/SS 误差界/全量库加性)。
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def load(*parts):
    return json.load(open(os.path.join(ROOT, *parts), encoding="utf-8"))


class TestH1RandomSlopes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "hierarchical_pool.json")["h1_random_slopes"]

    def test_structure(self):
        for k in ("loeo_random_slope_without_deweirdt",
                  "loeo_random_slope_with_deweirdt",
                  "direction_conflicts_cos_neg", "cis_loeo_random_slope",
                  "cis_reference_random_intercept", "reading", "criterion"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_conflicts_are_negative_cos(self):
        """方向门控: 被列名的冲突终点其余弦必须确实为负。"""
        for g, c in self.d["direction_conflicts_cos_neg"].items():
            self.assertLess(c, 0.0, "%s 的 cos=%s 不为负" % (g, c))

    def test_reading_matches_rule(self):
        cis = self.d["cis_loeo_random_slope"]
        ref = self.d["cis_reference_random_intercept"]
        self.assertAlmostEqual(ref, -0.429, places=3)
        if cis is not None and cis > 0:
            self.assertIn("转正", self.d["reading"])
        elif cis is not None and cis > ref:
            self.assertIn("改善但仍负", self.d["reading"])
        else:
            self.assertIn("维持原值", self.d["reading"])

    def test_random_slope_loeo_ranges(self):
        for g, v in self.d["loeo_random_slope_with_deweirdt"].items():
            self.assertGreaterEqual(v, -1.0)
            self.assertLessEqual(v, 1.0)


class TestH2PairwiseStability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "selector_model.json")["rank_stability_g2"]["h2_pairwise"]

    def test_structure(self):
        for k in ("n_pairs", "n_robust", "n_lean", "n_indistinguishable",
                  "shah_samworth_fp_bound_pi06", "pairs", "reading", "criterion"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_pair_counts_consistent(self):
        self.assertEqual(self.d["n_pairs"], 28)  # 8 员族 C(8,2)
        n = (self.d["n_robust"] + self.d["n_lean"]
             + self.d["n_indistinguishable"])
        self.assertEqual(n, self.d["n_pairs"])
        self.assertEqual(len(self.d["pairs"]), self.d["n_pairs"])

    def test_pair_probabilities_valid(self):
        for p in self.d["pairs"]:
            self.assertGreaterEqual(p["P"], 0.5)  # 胜者在前, P>=0.5
            self.assertLessEqual(p["P"], 1.0)
            band = ("robust" if p["P"] >= 0.9 else
                    "lean" if p["P"] >= 0.6 else "indistinguishable")
            self.assertEqual(p["band"], band)

    def test_shah_samworth_bound_formula(self):
        q = self.d["n_robust"] + self.d["n_lean"]
        expect = q * q / ((2 * 0.6 - 1) * self.d["n_pairs"])
        self.assertAlmostEqual(self.d["shah_samworth_fp_bound_pi06"],
                               round(expect, 2), places=2)

    def test_reading_matches_rule(self):
        if self.d["n_robust"] >= 21:  # ceil(0.75*28)
            self.assertIn("逐对多数稳健", self.d["reading"])
        else:
            self.assertIn("维持弱稳定", self.d["reading"])


class TestH3FullLibraryAdditive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "deweirdt_position_sensitivity.json")
        cls.h = cls.d["h3_full_library_additive"]

    def test_structure(self):
        for k in ("additive_neg_lfc127", "additive_neg_lfc128"):
            self.assertIn(k, self.h["results"])
            self.assertIn("vs_identity", self.h["results"][k])
            self.assertIn("vs_class3", self.h["results"][k])
            self.assertGreater(self.h["results"][k]["n_library_rows"], 10000)
            self.assertGreater(self.h["results"][k]["n_positions"], 10)

    def test_reading_matches_rule(self):
        res = self.h["results"]
        ok = all(res[k]["vs_identity"]["rho"] < 0 and
                 res[k]["vs_identity"]["perm_p"] < 0.05 for k in res)
        neg = all(res[k]["vs_identity"]["rho"] < 0 for k in res)
        if ok:
            self.assertIn("升级", self.h["reading"])
        elif neg:
            self.assertIn("维持边缘", self.h["reading"])
        else:
            self.assertIn("方向分裂", self.h["reading"])

    def test_single_mutant_caliber_kept(self):
        """G4 单点口径保留并列, 不被 H3 覆盖。"""
        for k in ("sens_neglfc127", "sens_neglfc128"):
            self.assertIn(k, self.d["results"])


if __name__ == "__main__":
    unittest.main()
