# -*- coding: utf-8 -*-
"""§I 批次契约测试(2026-09-11, 判据 docs/preregistration.md §I)。

沿用 §G/§H 惯例: 不锁死具体数值, 锁定
1. data/transfer_domain.json 的 I1/I2 块结构齐全;
2. 判读文字与数值按预登记规则同向(机器校验, 不手写结论方向);
3. 关键口径字段存在(包络阈/越界清单/两两 Spearman/内部对照)。
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def load(*parts):
    return json.load(open(os.path.join(ROOT, *parts), encoding="utf-8"))


class TestI1DomainOfApplicability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "transfer_domain.json")["i1_domain_of_applicability"]

    def test_structure(self):
        for k in ("n_train", "n_candidates", "n_candidates_pass_v6",
                  "n_family", "envelope_threshold", "frac_inside",
                  "per_feature_minmax", "outside_candidates", "reading"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_counts_consistent(self):
        self.assertEqual(self.d["n_train"], 46)
        self.assertEqual(self.d["n_candidates"],
                         self.d["n_candidates_pass_v6"] + self.d["n_family"])
        self.assertEqual(self.d["n_family"], 8)
        self.assertIn("family_inside_count", self.d)
        self.assertEqual(self.d["family_inside_count"]
                         + len(self.d["family_outside"]), 8)

    def test_outside_list_bounded(self):
        self.assertEqual(len(self.d["outside_candidates"]),
                         self.d["n_candidates"] - self.d["n_inside"])

    def test_reading_matches_rule(self):
        if self.d["frac_inside"] >= 0.95:
            self.assertIn("内插", self.d["reading"])
        else:
            self.assertIn("适用域边界", self.d["reading"])


class TestI2FoldCongruence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "transfer_domain.json")["i2_fold_congruence"]

    def test_structure(self):
        for k in ("cas12a2_entries", "cas12a_entries", "cross_pairs",
                  "cross_median_rho", "internal_median_rho",
                  "profiles_paired_prob_per_column", "reading"):
            self.assertIn(k, self.d, "缺字段 %s" % k)
        self.assertEqual(len(self.d["cas12a2_entries"]), 2)
        self.assertEqual(len(self.d["cas12a_entries"]), 9)
        self.assertEqual(len(self.d["cross_pairs"]),
                         2 * len(self.d["cas12a_entries"]))

    def test_rho_ranges(self):
        for c in self.d["cross_pairs"]:
            self.assertGreaterEqual(c["rho"], -1.0)
            self.assertLessEqual(c["rho"], 1.0)

    def test_profiles_length(self):
        n_cols = self.d["cross_pairs"] and 18
        for name, prof in self.d["profiles_paired_prob_per_column"].items():
            self.assertEqual(len(prof), 18, "%s 谱长≠18" % name)
            for v in prof:
                if v is not None:
                    self.assertGreaterEqual(v, 0.0)
                    self.assertLessEqual(v, 1.0)

    def test_reading_matches_rule(self):
        med = self.d["cross_median_rho"]
        ref = self.d["internal_median_rho"]
        if med >= 0.7:
            self.assertIn("获支持", self.d["reading"])
        elif med >= 0.4 or med >= ref - 0.1:
            self.assertIn("部分支持", self.d["reading"])
        else:
            self.assertIn("同构不成立", self.d["reading"])


if __name__ == "__main__":
    unittest.main()
