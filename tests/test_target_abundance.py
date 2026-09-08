# -*- coding: utf-8 -*-
"""靶RNA丰度档位特征契约测试(策划案 V3 §4.4 对齐整改, 2026-09-05)。

锁定四件事:
1. data/target_abundance_tiers.json 的实测值与 CCLE gct 源数据逐一一致
   (独立重解析比对, 防 JSON 手写漂移);
2. 分档换算正确: 阈值锚定 Scholz 2026 标定 EC50 95%CI, 每条记录
   tier == tier_of(mut_rpkm_het50), 边界行为固定;
3. 智能体特征向量含丰度档维: agent.target_context 输出 + 四个 panel
   预计算设计(data/agent/*.design.json)每行 target_features 均携带该维;
4. webapp 与离线选型器特征列一致: crrna_agent_webapp 与
   crrna_train_selector 共用同一 FEATURES / TARGET_FEATURES 定义。
"""
import json
import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import crrna_target_abundance as tabund  # noqa: E402

ABN = json.load(open(os.path.join(DATA, "target_abundance_tiers.json"),
                     encoding="utf-8"))
PANEL_KEYS = ["tp53_r248q", "kras_g12d", "tp53_r273h", "apc_q1328x"]


def _gct_measurements():
    """独立重解析 gct(不复用被测模块的 parse_gct), 返回 {(gene, line): rpkm}"""
    path = os.path.join(DATA, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")
    out = {}
    with open(path, "rt", encoding="utf-8", errors="ignore") as f:
        header = None
        for line in f:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if header is None:
                if parts[0] == "Name":
                    header = [h.split(" (")[0].split("_")[0]
                              for h in parts]
                continue
            sym = parts[1] if len(parts) > 1 else ""
            if sym in tabund.GENES and (sym, "") not in out:
                for cl in tabund.CELL_LINES:
                    idx = header.index(cl)
                    out[(sym, cl)] = float(parts[idx])
    return out


class TestAbundanceJsonVsGct(unittest.TestCase):
    def test_values_match_gct_source(self):
        meas = _gct_measurements()
        self.assertEqual(len(meas), 6)
        for g in tabund.GENES:
            for cl in tabund.CELL_LINES:
                rec = ABN["genes"][g]["lines"][cl]
                self.assertEqual(rec["total_rpkm"], meas[(g, cl)],
                                 "%s@%s 与 gct 源数据不一致" % (g, cl))
                self.assertEqual(rec["mut_rpkm_het50"],
                                 round(meas[(g, cl)] * tabund.ALLELE_FACTOR_HET, 5))
        for key, p in ABN["panel_targets"].items():
            self.assertEqual(p["total_rpkm"],
                             meas[(p["gene"], p["cell_line"])])

    def test_thresholds_anchored_to_calibration(self):
        cal = json.load(open(
            os.path.join(DATA, "scholz2026_fig1h_calibration.json"),
            encoding="utf-8"))
        self.assertEqual(ABN["tier_rule"]["edges_mut_rpkm"], cal["ec50_ci95"])
        self.assertEqual((tabund.TIER_LO, tabund.TIER_HI),
                         tuple(cal["ec50_ci95"]))


class TestTierRule(unittest.TestCase):
    def test_boundaries(self):
        lo, hi = tabund.TIER_LO, tabund.TIER_HI
        self.assertEqual(tabund.tier_of(0.0), 0)
        self.assertEqual(tabund.tier_of(lo - 1e-9), 0)
        self.assertEqual(tabund.tier_of(lo), 1)
        self.assertEqual(tabund.tier_of(hi - 1e-9), 1)
        self.assertEqual(tabund.tier_of(hi), 2)

    def test_json_tiers_consistent_with_rule(self):
        for g in tabund.GENES:
            for cl in tabund.CELL_LINES:
                rec = ABN["genes"][g]["lines"][cl]
                self.assertEqual(rec["abundance_tier"],
                                 tabund.tier_of(rec["mut_rpkm_het50"]))
                self.assertEqual(rec["tier_label"],
                                 tabund.TIER_LABELS[rec["abundance_tier"]])
        for key, p in ABN["panel_targets"].items():
            rec = ABN["genes"][p["gene"]]["lines"][p["cell_line"]]
            self.assertEqual(p["abundance_tier"], rec["abundance_tier"])

    def test_gate_verdict_follows_tier(self):
        expect = {0: "below_ci_lo", 1: "within_ci", 2: "above_ci_hi"}
        for key, p in ABN["panel_targets"].items():
            feat = tabund.target_feature(ABN, key)
            gate = tabund.activation_gate(feat)
            self.assertEqual(gate["verdict"], expect[p["abundance_tier"]])
            self.assertIn("%.2f" % p["mut_rpkm_het50"], gate["note"])


class TestAgentFeatureVector(unittest.TestCase):
    def test_agent_context_has_abundance_dim(self):
        import crrna_design_agent as agent
        for key in PANEL_KEYS:
            tctx = agent.target_context(key)
            feat = tctx["feature"]
            self.assertEqual(feat["target_abundance_tier"],
                             ABN["panel_targets"][key]["abundance_tier"])
            self.assertIn("target_abundance_tier", feat)
            self.assertIn("mut_rpkm_het50", feat)
            self.assertEqual(tctx["gate"]["verdict"],
                             {0: "below_ci_lo", 1: "within_ci",
                              2: "above_ci_hi"}[feat["target_abundance_tier"]])

    def test_panel_design_artifacts_carry_tier(self):
        for key in PANEL_KEYS:
            d = json.load(open(os.path.join(DATA, "agent",
                                            key + ".design.json"),
                               encoding="utf-8"))
            self.assertIsNotNone(d.get("target_context"), key)
            self.assertTrue(d["designs"], key)
            for row in d["designs"]:
                tf = row.get("target_features")
                self.assertIsNotNone(tf, "%s 行缺 target_features" % key)
                self.assertEqual(
                    tf["target_abundance_tier"],
                    ABN["panel_targets"][key]["abundance_tier"])


class TestSelectorFeatureConsistency(unittest.TestCase):
    def test_webapp_shares_offline_feature_columns(self):
        import crrna_train_selector as train
        import crrna_agent_webapp as webapp
        # 同一对象(import 共用), 不是两份字面量
        self.assertIs(webapp.FEATURES, train.FEATURES)
        self.assertIs(webapp.TARGET_FEATURES, train.TARGET_FEATURES)
        self.assertIn("target_abundance_tier", train.TARGET_FEATURES)
        # webapp panel 四靶均有丰度上下文
        for key in webapp.PANEL.values():
            self.assertIn(key, tabund.PANEL_TARGETS)


if __name__ == "__main__":
    unittest.main()
