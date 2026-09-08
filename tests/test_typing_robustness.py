# -*- coding: utf-8 -*-
"""typing.robustness.json 契约测试(2026-09 骨架分型事后敏感性分析)。

锁定两件事:
1. JSON 结构键齐备(四项补救分析 + 综合判定 + 口径声明);
2. 判读文字与数值按脚本显式阈值规则一致(判读由数值生成, 不允许手写方向
   与数据脱节)——此处按 scripts/crrna_typing_robustness.py 的阈值常量重算。
"""
import json
import os
import re
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
JSON_PATH = os.path.join(ROOT, "data", "context_typing", "typing.robustness.json")

ALPHA = 0.05            # 与脚本 ALPHA 一致
IDENTITY_HIGH = 0.80    # 与脚本 IDENTITY_HIGH 一致
MUT_POS_RE = re.compile(r"[ACGU](\d+)[ACGU]")


def _load():
    with open(JSON_PATH, encoding="utf-8") as f:
        return json.load(f)


class TestTypingRobustness(unittest.TestCase):
    def test_structure_keys(self):
        d = _load()
        for key in ("generated_by", "date", "scope", "pipeline_reruns",
                    "slice_consistency_check", "permutation_test",
                    "topk_sensitivity", "de_pos1_sensitivity",
                    "clustering_spectrum", "criterion_vs_k",
                    "small_cluster_audit", "verdict"):
            self.assertIn(key, d)
        self.assertIn("post-hoc", d["scope"])

    def test_permutation_rule_consistency(self):
        p = _load()["permutation_test"]
        self.assertEqual(p["n_partitions"], 280)  # 9!/(3!3!3!)/6 无序划分
        self.assertTrue(0.0 <= p["p_observed"] <= 1.0)
        has_no_power = "无判别力" in p["reading"]
        self.assertEqual(has_no_power, p["p_observed"] >= ALPHA)
        self.assertEqual(p["criterion"]["n_common"], p["observed"])

    def test_topk_sensitivity_rule_consistency(self):
        s = _load()["topk_sensitivity"]
        self.assertEqual(set(s["per_topk"]), {"6", "8", "10", "16"})
        for k, v in s["per_topk"].items():
            self.assertEqual(v["n_common"], len(v["common"]))
            self.assertEqual(v["pass"], v["n_common"] <= int(k) // 2)
        uniform = len({v["pass"] for v in s["per_topk"].values()}) == 1
        self.assertEqual("不随 topk 翻动" in s["reading"], uniform)

    def test_de_pos1_rule_consistency(self):
        s = _load()["de_pos1_sensitivity"]
        crit = s["criterion"]
        self.assertEqual(crit["pass"], crit["n_common"] <= crit["threshold"])
        for desc in crit["common"]:
            self.assertNotIn(1, [int(p) for p in MUT_POS_RE.findall(desc)])
        self.assertEqual("通过依赖 position-1 惰性位点" in s["reading"],
                         not crit["pass"])

    def test_clustering_spectrum_argmax(self):
        s = _load()["clustering_spectrum"]
        self.assertEqual(set(s["silhouette_by_method_k"]),
                         {"kmeans_repo", "agglomerative_ward", "gaussian_mixture"})
        for m, row in s["silhouette_by_method_k"].items():
            self.assertEqual(set(row), {"2", "3", "4", "5", "6"})
            best = max(row, key=lambda k: row[k])
            self.assertEqual(s["best_per_method"][m]["best_k"], int(best))

    def test_criterion_vs_k_rule_consistency(self):
        s = _load()["criterion_vs_k"]
        self.assertEqual(set(s["per_k"]), {"2", "3", "4", "5", "6"})
        for k, v in s["per_k"].items():
            self.assertEqual(len(v["sizes"]), int(k))
            self.assertEqual(v["pass"], v["n_common"] <= v["threshold"])
        allpass = all(v["pass"] for v in s["per_k"].values())
        self.assertEqual("全部通过" in s["reading"], allpass)

    def test_small_cluster_rule_consistency(self):
        s = _load()["small_cluster_audit"]
        self.assertEqual(len(s["members"]), 9)
        self.assertEqual(s["n_members_genes"], len(s["genes"]))
        batch = (s["identity_within"]["mean"] >= IDENTITY_HIGH
                 or s["n_members_genes"] < 4)
        self.assertEqual("批次效应" in s["reading"], batch)

    def test_verdict_rule_consistency(self):
        v = _load()["verdict"]
        self.assertEqual(v["n_flags"], len(v["weakening_flags"]))
        if v["n_flags"] >= 3:
            self.assertIn("证据变弱", v["reading"])
        elif v["n_flags"] >= 1:
            self.assertIn("证据不变", v["reading"])
        else:
            self.assertIn("证据变强", v["reading"])


if __name__ == "__main__":
    unittest.main()
