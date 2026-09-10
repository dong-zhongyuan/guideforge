# -*- coding: utf-8 -*-
"""v47 扩库分型重评估契约测试(2026-09-10)。

锁定四件事:
1. typing_v47_expansion.json 结构完整(队列/聚类/两层 §A2 判定/判读);
2. 队列构成: n_total = n_literature + n_v47_new, 新增 >200 条,
   provenance 逐条携带 kind/gene/transcript/pos;
3. §A2 规则一致: k 固定 3, 两层 positive 标记与 p_obs < 0.05 同向,
   判读方向不手写;
4. 扩库结论如实: 当前数据下 A2c 为阴性(P=0.0571 擦线)——测试锁定
   reading 在阴性时必须含"阴性"字样(防未来口径悄悄翻案而不更新评审)。
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
PATH = os.path.join(ROOT, "data", "context_typing", "typing_v47_expansion.json")


class TestV47Expansion(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.d = json.load(open(PATH, encoding="utf-8"))

    def test_structure(self):
        for k in ("cohort", "provenance", "clustering", "representatives",
                  "a2b_permutation", "a2c_no_pos1", "reading",
                  "comparison_literature_cohort", "registered"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_cohort_composition(self):
        c = self.d["cohort"]
        self.assertEqual(c["n_total"], c["n_literature"] + c["n_v47_new"])
        self.assertGreater(c["n_v47_new"], 200)
        self.assertEqual(len(self.d["provenance"]), c["n_v47_new"])
        for n, p in self.d["provenance"].items():
            for f in ("kind", "gene", "transcript", "pos"):
                self.assertIn(f, p)
        # 聚类簇大小合计 = 总队列
        self.assertEqual(sum(self.d["clustering"]["cluster_sizes"]),
                         c["n_total"])

    def test_a2_rule_consistency(self):
        self.assertEqual(self.d["clustering"]["k_fixed"], 3)
        for layer in ("a2b_permutation", "a2c_no_pos1"):
            blk = self.d[layer]
            self.assertEqual(blk["threshold_p"], 0.05)
            self.assertEqual(blk["positive"], blk["p_obs"] < 0.05)
            self.assertEqual(blk["n_partitions"], 280)

    def test_reading_honest_direction(self):
        pos_c = self.d["a2c_no_pos1"]["positive"]
        pos_b = self.d["a2b_permutation"]["positive"]
        if pos_b and pos_c:
            self.assertIn("阳性", self.d["reading"])
            self.assertNotIn("阴性", self.d["reading"])
        elif not pos_c:
            self.assertIn("阴性", self.d["reading"])

    def test_a3_enriched_rule_consistency(self):
        """§A3 特征增强块: 结构与 §A3c 判读规则锁定。"""
        a3 = self.d["a3_enriched"]
        self.assertEqual(len(a3["feature_names"]), 10)
        self.assertEqual(a3["clustering"]["k_fixed"], 3)
        self.assertEqual(sum(a3["clustering"]["cluster_sizes"]),
                         self.d["cohort"]["n_total"])
        for layer in ("a3b_with_pos1", "a3c_no_pos1"):
            blk = a3[layer]
            self.assertEqual(blk["positive"], blk["p_obs"] < 0.05)
        pos_a3 = a3["a3c_no_pos1"]["positive"]
        pos_a2 = self.d["a2c_no_pos1"]["positive"]
        if pos_a3 and not pos_a2:
            self.assertIn("归因成立", a3["reading"])
        elif not pos_a3 and not pos_a2:
            self.assertIn("计算层不支持", a3["reading"])
        elif pos_a3 and pos_a2:
            self.assertIn("稳健", a3["reading"])
        else:
            self.assertIn("如实并列", a3["reading"])


if __name__ == "__main__":
    unittest.main()
