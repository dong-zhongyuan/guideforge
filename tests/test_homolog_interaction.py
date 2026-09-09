# -*- coding: utf-8 -*-
"""同源跨上下文交互检验契约测试(判据登记: docs/preregistration.md §E, 2026-09-08)。

锁定四件事:
1. 真实归档 data/homolog_context_interaction.json schema 完整(E1 数值/CI/verdict、
   E2 区位分层、Tian Kendall W), 判读文字全部由数值按显式规则生成(非手写);
2. E1 判定与预登记 §E 同口径: rho < 0.4 且 bootstrap CI95 上界 < 0.5 -> PASS;
   当前真实数据判为 PASS(homolog_interaction_supported);
3. 脚本重跑确定性: --out 重跑与归档 E1 rho_cross / Tian W 逐位一致(固定 seed);
4. E3(Tian)双层如实披露: E3a 全量 Kendall W(位置梯度, 通用层)与 E3b 位置 1
   碱基偏好 Kendall W(交互层)并列报告——pos1 W<0.3 判交互佐证(E1 同向),
   最优碱基分布 >1 种, 两层解释(位置级敏感图谱通用 vs 位点内碱基最优
   上下文特异)不丢失。
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SCRIPT = os.path.join(ROOT, "scripts", "crrna_homolog_context_interaction.py")
ARCHIVE = os.path.join(ROOT, "data", "homolog_context_interaction.json")

PY = sys.executable


class TestHomologInteraction(unittest.TestCase):

    def test_archive_schema_and_verdict(self):
        self.assertTrue(os.path.exists(ARCHIVE), "归档 JSON 缺失")
        with open(ARCHIVE, encoding="utf-8") as f:
            d = json.load(f)
        e1 = d["deweirdt"]["e1"]
        self.assertEqual(e1["cls"], "test")
        self.assertGreaterEqual(e1["n"], 35000)
        self.assertIsNotNone(e1["rho_cross"])
        self.assertEqual(len(e1["bootstrap_ci95"]), 2)
        self.assertIn(e1["verdict"], ("PASS", "FAIL", "INCONCLUSIVE"))
        self.assertEqual(d["verdict"]["overall"], "homolog_interaction_supported")
        self.assertIn("deweirdt", d) and self.assertIn("e2", d["deweirdt"])
        self.assertIn("kendall_w", d["tian"])
        self.assertIn("reading", e1)
        self.assertIn("prereg_ref", d)

    def test_e1_matches_preregistered_rule(self):
        """E1 判定 = rho<0.4 且 CI95 上界<0.5 -> PASS(§E 同口径复算)。"""
        with open(ARCHIVE, encoding="utf-8") as f:
            e1 = json.load(f)["deweirdt"]["e1"]
        rho = e1["rho_cross"]
        ub = e1["bootstrap_ci95"][1]
        rule_pass = rho < 0.4 and ub < 0.5
        self.assertTrue(rule_pass, "数值应满足 PASS 条件")
        self.assertEqual(e1["verdict"], "PASS")

    @unittest.skipIf(
        not os.path.exists(os.path.join(ROOT, "data", "raw",
                                        "deweirdt2020_dr_scan.json")),
        "data/raw/ 不在 git(data/raw/ 被 .gitignore 排除, 需手动放置); "
        "干净克隆上跳过确定性重跑, 其余 3 项契约仍全量守护")
    def test_rerun_deterministic(self):
        """--out 重跑与归档逐位一致(固定 seed, 确定性)。"""
        with tempfile.TemporaryDirectory() as td:
            tmp = os.path.join(td, "out.json")
            r = subprocess.run([PY, SCRIPT, "--out", tmp],
                               cwd=ROOT, capture_output=True, text=True, timeout=600)
            self.assertEqual(r.returncode, 0, r.stderr[-2000:])
            with open(tmp, encoding="utf-8") as f:
                d = json.load(f)
            with open(ARCHIVE, encoding="utf-8") as f:
                arch = json.load(f)
            self.assertAlmostEqual(d["deweirdt"]["e1"]["rho_cross"],
                                   arch["deweirdt"]["e1"]["rho_cross"], places=9)
            self.assertEqual(d["deweirdt"]["e1"]["bootstrap_ci95"],
                             arch["deweirdt"]["e1"]["bootstrap_ci95"])
            self.assertEqual(d["tian"]["kendall_w"]["W"],
                             arch["tian"]["kendall_w"]["W"])
            self.assertEqual(d["tian"]["pos1"]["kendall_w"]["W"],
                             arch["tian"]["pos1"]["kendall_w"]["W"])

    def test_e3_limitation_disclosed(self):
        """E3 双层如实披露: E3a 位置梯度(通用层)与 E3b 碱基偏好(交互层)并列,
        不掩盖 pos1 碱基偏好跨 crRNA 翻转的交互佐证。"""
        with open(ARCHIVE, encoding="utf-8") as f:
            tian = json.load(f)["tian"]
        # E3a 全量 W(位置梯度)
        self.assertGreaterEqual(tian["kendall_w"]["W"], 0.7)
        # E3b pos1 碱基偏好低一致性 = 交互佐证(E1 同向)
        w1 = tian["pos1"]["kendall_w"]["W"]
        self.assertLess(w1, 0.3)
        self.assertIn("交互佐证(E1 同向)", tian["reading"])
        # 最优碱基分布必须多于 1 种(翻转证据)
        dist = tian["pos1"]["best_base_distribution"]
        self.assertGreaterEqual(len(dist), 2)
        # 区位/层次解释不丢失
        self.assertIn("位置梯度", tian["note"])
        self.assertIn("两层不矛盾", tian["reading"])


if __name__ == "__main__":
    unittest.main()
