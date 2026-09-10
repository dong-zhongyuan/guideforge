# -*- coding: utf-8 -*-
"""Round-4 §G 批次契约测试(2026-09-11, 判据 docs/preregistration.md §G)。

不锁死具体数值(置换检验带种子但仍属统计量), 锁定三件事:
1. 各产物 JSON 的 §G 块结构齐全;
2. 判读文字与数值按预登记规则同向(不手写结论方向的机器校验);
3. 关键口径字段存在(双控偏相关/分层置换/合并 p/LOCO/层级 LOEO 等)。
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def load(*parts):
    return json.load(open(os.path.join(ROOT, *parts), encoding="utf-8"))


class TestG1CrossDecompose(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "cross_decomposition.json")

    def test_dual_control_and_stratified_fields(self):
        f = self.d["contexts"]["mature18"]["features"]["raw_coreDR"]
        self.assertIn("partial_given_GC_A", f)
        self.assertIn("stratified_perm_p_GC", f)
        self.assertIsInstance(f["partial_given_GC_A"], float)
        self.assertGreaterEqual(f["stratified_perm_p_GC"], 0.0)
        self.assertLessEqual(f["stratified_perm_p_GC"], 1.0)

    def test_constant_column_guard(self):
        """mature18 无 flank -> raw_flank 常数列, 双控/分层必须为 None 而非
        NaN 或伪 p=0。"""
        f = self.d["contexts"]["mature18"]["features"]["raw_flank"]
        self.assertIsNone(f["partial_given_GC_A"])
        self.assertIsNone(f["stratified_perm_p_GC"])

    def test_g1c_reading_matches_rule(self):
        v = self.d["verdict"]["G1c_preregistered"]
        dual, p = v["partial_given_GC_A"], v["stratified_perm_p_GC"]
        if abs(dual) < 0.3 and p >= 0.05:
            self.assertIn("成立并强化", v["reading"])
        elif abs(dual) >= 0.3 and p < 0.05:
            self.assertIn("归因不完整", v["reading"])
        else:
            self.assertIn("部分支持", v["reading"])


class TestG2RankStability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "selector_model.json")["rank_stability_g2"]

    def test_structure(self):
        for k in ("top1_selection_frequency", "rank_ci95_median_rank",
                  "top1_single_tree", "top1_bootstrap_mode", "top1_borda",
                  "borda_order", "reading", "criterion"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_frequencies_valid(self):
        for n, f in self.d["top1_selection_frequency"].items():
            self.assertGreaterEqual(f, 0.0)
            self.assertLessEqual(f, 1.0)
        self.assertAlmostEqual(
            sum(self.d["top1_selection_frequency"].values()), 1.0, places=2)

    def test_reading_matches_rule(self):
        freq = self.d["top1_selection_frequency"][self.d["top1_bootstrap_mode"]]
        if freq >= 0.5:
            self.assertIn("首位稳定", self.d["reading"])
        else:
            self.assertIn("弱稳定", self.d["reading"])
        if self.d["top1_borda"] != self.d["top1_single_tree"]:
            self.assertIn("以 Borda 为报告口径", self.d["reading"])


class TestG3CombinedFisher(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "tp53_r248q_zengdr.v6.dmytrenko_validation.json")
        cls.c = cls.d["B_conservation_corr"]["combined_fisher_g3"]

    def test_structure_and_range(self):
        self.assertEqual(self.c["n_perm"], 20000)
        self.assertGreater(self.c["x2_obs"], 0)
        self.assertGreaterEqual(self.c["combined_perm_p"], 0.0)
        self.assertLessEqual(self.c["combined_perm_p"], 1.0)

    def test_reading_matches_rule(self):
        p = self.c["combined_perm_p"]
        if p < 0.05 and self.c["direction_all_negative"]:
            self.assertIn("升级", self.c["reading"])
        else:
            self.assertIn("维持", self.c["reading"])

    def test_marginal_calibers_still_reported(self):
        b = self.d["B_conservation_corr"]
        for k in ("spearman_cons_vs_tolerance", "spearman_cons_vs_meanscore",
                  "spearman_cons_vs_scoreadj"):
            self.assertIn("rho", b[k])
            self.assertIn("perm_p", b[k])


class TestG4DeweirdtPosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "deweirdt_position_sensitivity.json")

    def test_structure(self):
        self.assertGreater(self.d["n_positions"], 10)
        self.assertGreater(self.d["n_single_mutants_used"], 0)
        for k in ("sens_neglfc127", "sens_neglfc128"):
            self.assertIn("vs_identity", self.d["results"][k])
            self.assertIn("vs_class3", self.d["results"][k])

    def test_reading_matches_rule(self):
        res = self.d["results"]
        ok = all(res[k]["vs_identity"]["rho"] < 0 and
                 res[k]["vs_identity"]["perm_p"] < 0.05 for k in res)
        both_neg = all(res[k]["vs_identity"]["rho"] < 0 for k in res)
        if ok:
            self.assertIn("第二条独立证据链", self.d["reading"])
        elif both_neg:
            self.assertIn("方向一致为负", self.d["reading"])
        else:
            self.assertIn("方向分裂", self.d["reading"])

    def test_coordinate_offset_asserted(self):
        self.assertIn("offset_rule", self.d["data"])
        self.assertEqual(self.d["data"]["coordinate_offset"], 1)


class TestG5PkCaliber(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "creutzburg_pk_caliber.json")

    def test_structure(self):
        self.assertIn("metric_revision", self.d)
        self.assertEqual(self.d["stats"]["p_stem_paper"]["perm_p"] >= 0, True)
        for k in ("component_stem_unconstrained", "component_w_P1_state",
                  "component_dG_pk_cost"):
            self.assertIn(k, self.d["stats"])

    def test_reading_matches_rule(self):
        s = self.d["stats"]["p_stem_paper"]
        if s["spearman"] > 0 and s["perm_p"] < 0.05:
            self.assertIn("符号复现", self.d["reading"])
        else:
            self.assertIn("未达预登记阳性阈", self.d["reading"])

    def test_p1_pair_annotation(self):
        self.assertEqual(self.d["p1_pair"]["dr18_0based"], [1, 9])
        self.assertEqual(self.d["p1_pair"]["bases"], "A:U")


class TestG6RrsPositionModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "rrs_position_model.json")

    def test_structure(self):
        for k in ("position_only", "structure_only", "position_plus_structure"):
            self.assertIn(k, self.d["results"])
            self.assertEqual(len(self.d["results"][k]["per_fold"]), 8)

    def test_reading_matches_rule(self):
        rho = self.d["results"]["position_plus_structure"]["loco_spearman"]
        if rho >= 0.3:
            self.assertIn("升级", self.d["reading"])
            self.assertIn("Tian 同族 RRS 区", self.d["scope_limit"])
        else:
            self.assertIn("维持", self.d["reading"])

    def test_position_beats_structure(self):
        """论文自身结论(位置为主导)的数值化守护: 位置特征消融必须优于
        纯结构特征(否则升级理由不成立)。"""
        res = self.d["results"]
        if "升级" in self.d["reading"]:
            self.assertGreater(res["position_only"]["loco_spearman"],
                               res["structure_only"]["loco_spearman"])


class TestG7HierarchicalPool(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = load("data", "hierarchical_pool.json")

    def test_structure(self):
        for k in ("loeo_hierarchical_without_deweirdt",
                  "loeo_hierarchical_with_deweirdt",
                  "intercepts_full_fit", "shared_slopes_standardized",
                  "n_main5_positive_with_dw", "reading", "cis_note"):
            self.assertIn(k, self.d, "缺字段 %s" % k)

    def test_n_positive_recompute(self):
        lo = self.d["loeo_hierarchical_with_deweirdt"]
        n = sum(1 for g in ("fig1g", "cis_end", "trans_end", "rbs0", "rbs33")
                if lo.get(g, -1) > 0)
        self.assertEqual(n, self.d["n_main5_positive_with_dw"])

    def test_reading_matches_rule(self):
        if self.d["n_main5_positive_with_dw"] >= 3:
            self.assertIn("修订为「层级部分池化入池」", self.d["reading"])
        else:
            self.assertIn("维持拒入", self.d["reading"])


if __name__ == "__main__":
    unittest.main()
