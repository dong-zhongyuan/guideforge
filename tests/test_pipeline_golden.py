# -*- coding: utf-8 -*-
"""主管线黄金集回归测试(2026-09-02, 外部审计 P1 项)。

目的: 守护 crrna_scaffold_design.py 的打分/过滤行为——任何重构若静默改变
候选库结果(折叠口径/过滤阈值/打分公式/排序), 本测试的红灯先于结果漂移。

黄金值锁定口径(锁定日期 2026-09-02, ViennaRNA 2.7.2, py3.10):
  - 单元级: WT 构建(zeng2026 DR + R248Q spacer)的 MFE/结构/配分统计/交叉
    配对/侵占螺旋/茎完整概率;
  - 端到端: 策略A 确定性小跑(78 变体)的变体总数/通过数/TOP 排序首条。

变更黄金值的唯一合法理由: 有意修改口径(此时须同步更新本文件与锁定日期,
并在提交信息中说明口径变化)。
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

SPACER_DNA = "GTTCATGCCGCCCATGCAGGAACT"
DR_DNA = "AATTTCTACTGTTGTAGAT"          # cas12a2_zeng2026 注册表条目

GOLDEN_UNIT = {
    "mfe": -10.60,
    "ss": "....(((((....))))).((((.(((.......))).)))).",
    "div": 0.9015,
    "sp_up": 0.4267,
    "seed_up": 0.3066,
    "cross_nt": 0,
    "cross_pp": 0.0163,
    "inv_max_run": 0,
    "stem_pairs": [[4, 17], [5, 16], [6, 15], [7, 14], [8, 13]],
    "p_fold": 0.86603,
}
GOLDEN_E2E = {
    "n_variants": 78,
    "n_passed": 43,
    "top1_desc": "WT",
    "top1_score": 0.0,
    "top2_desc": "A1C",
}


class GoldenUnitTest(unittest.TestCase):
    """核心打分函数的数值锁定(相对 WT 构建的单元口径)"""

    def test_wt_construct_features(self):
        import crrna_scaffold_design as core
        dr = core.to_rna(DR_DNA)
        sp = core.to_rna(SPACER_DNA)
        full = dr + sp
        ss, mfe = core.fold(full)
        div, sp_up, seed_up = core.pf_stats(full, len(dr), 7)
        cross_nt, cross_pp = core.cross_pairs(full, len(dr))
        inv = core.inv_max_run(full, len(dr))
        pairs = sorted([list(p) for p in core.stem_pairs_of(core.fold(dr)[0])])
        pf = core.stem_intact_prob(full, [tuple(p) for p in pairs])
        self.assertEqual(ss, GOLDEN_UNIT["ss"])
        self.assertAlmostEqual(mfe, GOLDEN_UNIT["mfe"], places=2)
        self.assertAlmostEqual(div, GOLDEN_UNIT["div"], places=3)
        self.assertAlmostEqual(sp_up, GOLDEN_UNIT["sp_up"], places=3)
        self.assertAlmostEqual(seed_up, GOLDEN_UNIT["seed_up"], places=3)
        self.assertEqual(cross_nt, GOLDEN_UNIT["cross_nt"])
        self.assertAlmostEqual(cross_pp, GOLDEN_UNIT["cross_pp"], places=4)
        self.assertEqual(inv, GOLDEN_UNIT["inv_max_run"])
        self.assertEqual(pairs, [list(p) for p in GOLDEN_UNIT["stem_pairs"]])
        self.assertAlmostEqual(pf, GOLDEN_UNIT["p_fold"], places=5)


class GoldenE2ETest(unittest.TestCase):
    """策略A 确定性小跑的端到端锁定(变体数/过滤通过数/排序首条)"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="gf_golden_")
        cls.prefix = os.path.join(cls.tmp, "golden")
        cmd = [sys.executable,
               os.path.join(ROOT, "scripts", "crrna_scaffold_design.py"),
               "--effector", "cas12a2_zeng2026",
               "--spacer", SPACER_DNA,
               "--strategy", "A", "--n-double", "20", "--topk", "6",
               "--out-prefix", cls.prefix]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
        if r.returncode != 0:
            raise RuntimeError("管线小跑失败: %s" % r.stderr[-500:])

    def test_counts_and_top(self):
        with open(self.prefix + ".variants.csv", newline="",
                  encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        passed = [r for r in rows if r["passed"] in ("True", "true", "1")]
        top = json.load(open(self.prefix + ".top.json",
                             encoding="utf-8"))["top"]
        self.assertEqual(len(rows), GOLDEN_E2E["n_variants"])
        self.assertEqual(len(passed), GOLDEN_E2E["n_passed"])
        self.assertEqual(top[0]["desc"], GOLDEN_E2E["top1_desc"])
        self.assertAlmostEqual(top[0]["score"], GOLDEN_E2E["top1_score"],
                               places=4)
        self.assertEqual(top[1]["desc"], GOLDEN_E2E["top2_desc"])


if __name__ == "__main__":
    unittest.main()
