# -*- coding: utf-8 -*-
"""fig1f 留出验证线性第二口径契约测试(2026-09-08)。

背景: 决策树按构造不可外推(预测恒在训练 y 值域内), fig1f 留出 9/9
标签超训练上界 → MAE 超阈。修复方案 = 并列第二口径(Ridge+训练折
StandardScaler, 可外推)。实测两口径 MAE 一致超阈(0.217 vs 0.216),
证伪「模型不可外推」假设、坐实图版间水平偏移为数据因素。

锁定四件事:
1. selector_model.json 的 fig1f_holdout_check 携带 linear_second_caliber
   块且关键字段齐全;
2. 数值自洽: spearman/mae 数值型, mae_within_threshold 与
   (mae < mae_threshold_trainfold) 一致;
3. 判读方向一致: reading 文案与 mae_within_threshold 同向(不手写);
4. 归因成立: 决策树与线性两口径 MAE 均超阈(模型因素证伪),
   extrapolation_note 保持图版偏移归因。
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MODEL = json.load(open(os.path.join(ROOT, "data", "selector_model.json"),
                       encoding="utf-8"))


class TestLinearSecondCaliber(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ho = MODEL["fig1f_holdout_check"]
        cls.lin = cls.ho["linear_second_caliber"]

    def test_block_present_with_keys(self):
        for k in ("model", "rationale", "spearman", "spearman_boot_ci95",
                  "mae", "mae_boot_ci95", "mae_threshold_trainfold",
                  "mae_within_threshold", "pred_range",
                  "n_pred_outside_train_range", "coef_standardized",
                  "reading"):
            self.assertIn(k, self.lin, "缺字段 %s" % k)

    def test_numeric_consistency(self):
        self.assertIsInstance(self.lin["spearman"], float)
        self.assertIsInstance(self.lin["mae"], float)
        self.assertGreater(self.lin["mae"], 0)
        thr = self.lin["mae_threshold_trainfold"]
        self.assertEqual(self.lin["mae_within_threshold"],
                         self.lin["mae"] < thr)
        self.assertAlmostEqual(thr, self.ho["mae_threshold_trainfold"],
                               places=4)
        self.assertEqual(len(self.lin["coef_standardized"]),
                         len(MODEL["features"]))

    def test_reading_direction_matches_flag(self):
        if self.lin["mae_within_threshold"]:
            self.assertIn("可采用", self.lin["reading"])
        else:
            self.assertIn("仍不采用", self.lin["reading"])

    def test_attribution_model_factor_falsified(self):
        # 两口径 MAE 均超阈 → 「树不可外推」假设被证伪, 归因为数据因素
        thr = self.ho["mae_threshold_trainfold"]
        self.assertGreaterEqual(self.ho["mae"], thr)
        self.assertGreaterEqual(self.lin["mae"], thr)
        self.assertIn("图版", self.ho["extrapolation_note"])


if __name__ == "__main__":
    unittest.main()
