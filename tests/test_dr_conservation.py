# -*- coding: utf-8 -*-
"""任务⑪契约测试: 同源 DR 扩库 + 保守性先验 + 下游消费口径(2026-09-08)。

锁定四条契约:
  1. data/dr_homologs.json 每条序列必须有来源(provenance + sources 非空)——
     禁止编造序列的结构性护栏;
  2. 每条序列与 WT 口径可比对(锚 TTTCTACT/GTAGA 各唯一出现), 且注册表
     cas12a2 / cas12a2_zeng2026 的 DR 在库内、去 gap 后还原注册表序列;
  3. data/dr_conservation.json 的保守位点划分与数值自洽:
     identity 可由 composition 重算复现, class 与阈值规则一致,
     cons3_window_prior = strictly_conserved 3' run + trailing3 长度;
  4. 下游消费端 crrna_scaffold_design.load_cons3_window_prior 按 DR 序列
     读取新口径(= mapped_effectors 的 cons3_window_prior)。
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src"))

import crrna_dr_conservation as dc  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

HOMOLOGS = os.path.join(ROOT, "data", "dr_homologs.json")
CONSERVATION = os.path.join(ROOT, "data", "dr_conservation.json")


def load_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class TestHomologLibraryProvenance(unittest.TestCase):
    def setUp(self):
        self.entries = dc.load_library()

    def test_every_entry_has_provenance_and_sources(self):
        for e in self.entries:
            self.assertTrue(e.get("provenance"), e["name"])
            srcs = e.get("sources")
            self.assertIsInstance(srcs, list, e["name"])
            self.assertTrue(all(s.strip() for s in srcs), e["name"])

    def test_library_size_expanded(self):
        """扩库契约: 条目数达到两位数(round-2 批评 n=4 小样本)。"""
        self.assertGreaterEqual(len(self.entries), 10)

    def test_sequences_are_dna_with_unique_anchors(self):
        """每条序列可被比对规则消费: ACGT 字母 + 两锚各唯一出现。"""
        for e in self.entries:
            seq = e["dr_dna"]
            self.assertFalse(set(seq) - set("ACGT"), e["name"])
            self.assertEqual(seq.count(dc.ANCHOR_L), 1, e["name"])
            self.assertEqual(seq.count(dc.ANCHOR_R), 1, e["name"])
            seg = dc.segmented(e["name"], seq, e["provenance"])
            self.assertTrue(seg["middle"], e["name"])

    def test_registry_drs_present_and_recoverable(self):
        """WT 比对一致性: 注册表两条 Cas12a2 DR 在库内, 且分段拼接还原全长。"""
        registry = load_json(os.path.join(ROOT, "configs", "scaffold_registry.json"))["entries"]
        by_seq = {e["dr_dna"]: e for e in self.entries}
        for key in ("cas12a2", "cas12a2_zeng2026"):
            dr = registry[key]["scaffold"]
            self.assertIn(dr, by_seq, key)
            seg = dc.segmented("wt", dr, "test")
            rebuilt = seg["flank5"] + seg["anchor_l"] + seg["middle"] \
                + seg["anchor_r"] + seg["trailing3"]
            self.assertEqual(rebuilt, dr)


class TestConservationJsonInternalConsistency(unittest.TestCase):
    def setUp(self):
        self.d = load_json(CONSERVATION)

    def test_column_identity_recomputable_from_composition(self):
        """identity = 最高频非 gap 碱基数 / 全列条目数(gap 计入分母)。"""
        n = self.d["n_sequences"]
        for c in self.d["conservation"]:
            comp = c["composition"]
            self.assertEqual(sum(comp.values()), n)
            nongap = {b: v for b, v in comp.items() if b != "-"}
            expect = round(max(nongap.values()) / n, 3)
            self.assertEqual(c["identity"], expect, c)
            self.assertEqual(c["gap_n"], comp.get("-", 0))

    def test_class_matches_threshold_rule(self):
        for c in self.d["conservation"]:
            self.assertEqual(c["class"], dc.classify(c["identity"]), c)
            self.assertEqual(c["class_dedup"], dc.classify(c["identity_dedup"]), c)

    def test_cons3_window_prior_rule(self):
        """cons3_window_prior = strictly_conserved 3' run + trailing3 长度。"""
        strict, cons = dc.cons3_runs(self.d["conservation"])
        self.assertEqual(self.d["cons3_prior"]["strictly_conserved_3prime_run"], strict)
        self.assertEqual(self.d["cons3_prior"]["conserved_3prime_run_at_threshold"], cons)
        for key, m in self.d["mapped_effectors"].items():
            trailing = len(m["trailing3"])
            self.assertEqual(m["cons3_window_prior"], strict + trailing, key)

    def test_mapped_positions_consistent_with_columns(self):
        """逐位映射: core 位点的 identity/class 与对应列一致, base 还原 DR。"""
        for key, m in self.d["mapped_effectors"].items():
            dr = m["dr_dna"]
            for p in m["positions"]:
                self.assertEqual(dr[p["pos_1based"] - 1], p["base"])
                if p["region"] == "core":
                    c = self.d["conservation"][p["column"] - 1]
                    self.assertEqual(p["identity"], c["identity"])
                    self.assertEqual(p["class"], c["class"])
                else:
                    self.assertIsNone(p["column"])


class TestDownstreamConsumption(unittest.TestCase):
    def test_scaffold_design_resolves_window_from_prior(self):
        """下游消费端按 DR 序列读到新口径窗长。"""
        d = load_json(CONSERVATION)
        for key, m in d["mapped_effectors"].items():
            got = core.load_cons3_window_prior(m["dr_dna"])
            self.assertIsNotNone(got, key)
            self.assertEqual(got[0], m["cons3_window_prior"], key)
            self.assertIn(key, got[1])

    def test_unknown_dr_falls_back(self):
        self.assertIsNone(core.load_cons3_window_prior("ACGTACGTACGTACGT"))

    def test_argparse_default_is_prior_driven(self):
        """--cons3-window 默认 None(先验驱动), 不再硬编码 5。"""
        import argparse
        import inspect
        src = inspect.getsource(core.main)
        self.assertIn("'--cons3-window', type=int, default=None", src)


if __name__ == "__main__":
    unittest.main()
