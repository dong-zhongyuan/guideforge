# -*- coding: utf-8 -*-
"""智能体阴性退化模式契约测试(策划案 V3 §5.3 末段 / 表5 风险1 对齐整改, 2026-09-06)。

锁定五件事:
1. 判据纯规则(gonogo_verdict)与预登记 §A(docs/preregistration.md)同口径:
   全型公共 TOP-k 子集 <= topk//2 且型数>1 -> 阳性(typed), 否则阴性
   (degraded_universal); 边界 4/5 与单型情形固定;
2. 真实归档产物(data/context_typing/)当前判为 typed(弱证据阳性, 公共集
   3 条), 重算公共集与归档 top_common_to_all_types 一致, 判读文字由数值
   按规则生成(含数值、方向与 significant 一致);
3. 端到端两场景(合成 fixture 子进程跑 CLI --spacer): 触发时全部候选退化
   为通用型单骨架(dr_desc=universal_scaffold 全场最优, scaffold_type=None,
   分型降级为 typing_nearest_type 诊断, 文案含退化标注+量化证据+"科学结论
   依然完整可交付"); 不触发时维持分型推荐(逐行 dr_desc=该型代表骨架);
4. --degrade on/off 强制路径如实标注 forced 字段;
5. webapp 载荷契约: /api/design 响应携带 degradation 块且模式与离线判据
   一致(当前真实数据 = typed)。
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import crrna_design_agent as agent  # noqa: E402

CLUSTERS = os.path.join(ROOT, "data", "context_typing", "typing.clusters.json")
SPACERS_REAL = os.path.join(ROOT, "data", "context_typing", "spacers.txt")
DR_DNA = "AATTTCTACTGTTGTAGAT"          # cas12a2_zeng2026 注册表条目
SPACER_DNA = "GTTCATGCCGCCCATGCAGGAACT"


def _variant_dr(i):
    """DR 的非 position-1 单点变体(第 3+i%3 位轮换), 仅供 fixture 使用。"""
    pos = 3 + i % 3
    return DR_DNA[:pos - 1] + "ACG"[i % 3] + DR_DNA[pos:]


def _write_fixture(dirpath, tops_per_type, top_scores):
    """合成分型产物 fixture: 2 型, 每型 3 spacer + 1 个代表 top.json。

    tops_per_type: {t: [desc x8]} —— 该型代表管线 TOP-8(跳过 rank0 WT),
    直接决定全型公共子集大小(判据输入); top_scores: {t: {desc: score}},
    每型代表骨架 = 该型分最高的非 pos-1 变体(与 load_type_models 同口径)。
    """
    with open(SPACERS_REAL, encoding="utf-8") as fh:
        pool = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    picks = [pool[0], pool[21], pool[42], pool[63], pool[7], pool[35]]
    names = ["T%s_%s" % (i % 2, p[0].replace(" ", "_")) for i, p in enumerate(picks)]
    with open(os.path.join(dirpath, "spacers.txt"), "w", encoding="utf-8") as fh:
        for n, (_, s) in zip(names, picks):
            fh.write("%s\t%s\n" % (n, s))
    types = {}
    for t in (0, 1):
        members = [n for n in names if n.startswith("T%d_" % t)]
        types[str(t)] = {"n": len(members), "members": members,
                         "feature_means": {}}
    with open(os.path.join(dirpath, "typing.clusters.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"effector": "cas12a2_zeng2026", "dr_dna": DR_DNA,
                   "features": ["GC", "self_mfe", "spacer_up",
                                "junction_pairs", "internal_pairs"],
                   "chosen_k": 2, "types": types}, fh, ensure_ascii=False)
    for t in (0, 1):
        rows = [{"rank": 0, "desc": "WT", "score": 0.0, "passed": True,
                 "mut_positions": [], "dr_seq": DR_DNA}]
        for i, desc in enumerate(tops_per_type[t], 1):
            rows.append({"rank": i, "desc": desc,
                         "score": top_scores[t].get(desc, -0.5),
                         "passed": True, "mut_positions": [3 + i % 3],
                         "dr_seq": _variant_dr(i)})
        with open(os.path.join(dirpath, "typing.t%d_REP.top.json" % t), "w",
                  encoding="utf-8") as fh:
            json.dump({"spacer_fixed_dna": picks[t][1].replace("U", "T"),
                       "top": rows}, fh, ensure_ascii=False)
    return os.path.join(dirpath, "typing.clusters.json")


def _run_agent(clusters_json, out_prefix, extra_args=()):
    env = dict(os.environ, PYTHONUTF8="1")
    cmd = [sys.executable, os.path.join(ROOT, "scripts",
                                        "crrna_design_agent.py"),
           "--spacer", SPACER_DNA, "--effector", "cas12a2_zeng2026",
           "--clusters", clusters_json,
           "--spacers", os.path.join(os.path.dirname(clusters_json),
                                     "spacers.txt"),
           "--out-prefix", out_prefix] + list(extra_args)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-500:]
    with open(out_prefix + ".design.json", encoding="utf-8") as fh:
        return json.load(fh), r.stdout


class TestGonogoRule(unittest.TestCase):
    """判据纯规则: 阈值/边界/型数(与预登记 §A 逐条对应)"""

    def test_boundary_common_4_vs_5(self):
        u4 = {"0": {"V1", "V2", "V3", "V4", "A"},
              "1": {"V1", "V2", "V3", "V4", "B"}}       # 公共 4 = 阈
        v = agent.gonogo_verdict(u4, topk=8)
        self.assertEqual(v["n_common"], 4)
        self.assertEqual(v["threshold"], 4)
        self.assertTrue(v["significant"])                # <= 阈 -> 阳性
        u5 = {"0": {"V1", "V2", "V3", "V4", "V5", "A"},
              "1": {"V1", "V2", "V3", "V4", "V5", "B"}}  # 公共 5 > 阈
        v = agent.gonogo_verdict(u5, topk=8)
        self.assertFalse(v["significant"])               # > 阈 -> 阴性(退化)

    def test_single_type_is_not_significant(self):
        v = agent.gonogo_verdict({"0": {"V1"}}, topk=8)   # §A 要求型数>1
        self.assertFalse(v["significant"])


class TestRealArchiveAssessment(unittest.TestCase):
    """真实归档产物: 当前数据判 typed(弱证据阳性), 与归档/判据口径一致"""

    def test_current_data_is_typed(self):
        a = agent.interaction_assessment(CLUSTERS)
        self.assertEqual(a["mode"], "typed")
        self.assertTrue(a["significant"])
        self.assertEqual(a["topk"], 8)
        self.assertEqual(a["threshold"], 4)
        self.assertEqual(a["n_types"], 3)
        self.assertEqual(a["n_common"], 3)               # {A1C, A1G, A1U}
        self.assertEqual(a["common"], ["A1C", "A1G", "A1U"])
        self.assertIsNone(a["forced"])
        self.assertNotIn("deliverable_note", a)

    def test_recompute_matches_archive(self):
        a = agent.interaction_assessment(CLUSTERS)
        self.assertTrue(a["archive_crosscheck"]["matches"])
        self.assertEqual(a["archive_crosscheck"]["archived_common"],
                         a["archive_crosscheck"]["recomputed_common"])

    def test_reading_follows_numbers(self):
        a = agent.interaction_assessment(CLUSTERS)
        self.assertIn("维持分型推荐", a["reading"])
        self.assertIn("3 条 <= 阈 4", a["reading"])
        # 事后置换检验证据(post-hoc): p=0.0357 < 0.05 -> 支持判别力
        ev = a["robustness_evidence"]
        self.assertIsNotNone(ev)
        self.assertAlmostEqual(ev["permutation_p_observed"], 0.0357)
        self.assertIn("罕见", ev["reading"])


class TestUniversalScaffold(unittest.TestCase):
    """通用型单骨架选取规则(既有数据口径: 管线分全场最优, WT=0 基线)"""

    def _model(self, scores):
        return {"type_dr": {t: {"representative_desc": d,
                                "representative_dr": _variant_dr(i),
                                "representative_score": s,
                                "wt_dr": DR_DNA}
                            for i, (t, (d, s)) in enumerate(scores.items())}}

    def test_best_positive_variant_wins(self):
        m = self._model({0: ("U3A", -0.5), 1: ("V9", 0.5), 2: ("U3C", -0.3)})
        best = agent.universal_scaffold(m)
        self.assertEqual(best["desc"], "V9")
        self.assertEqual(best["source"], "type_1_representative")

    def test_all_below_baseline_falls_back_to_wt(self):
        m = self._model({0: ("U3A", -0.5), 1: ("U15C", -0.3)})
        best = agent.universal_scaffold(m)
        self.assertEqual(best["desc"], "WT")
        self.assertEqual(best["score"], 0.0)

    def test_tie_prefers_wt(self):
        m = self._model({0: ("U3A", 0.0), 1: ("U15C", -0.3)})
        best = agent.universal_scaffold(m)
        self.assertEqual(best["desc"], "WT")             # 并列 WT 优先


class TestEndToEndFixtures(unittest.TestCase):
    """端到端两场景(合成 fixture 子进程): 触发退化 / 不触发维持分型"""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        base = cls._tmp.name
        # 触发场景: 两型 TOP-8 完全相同 -> 公共集 8 > 阈 4 -> 阴性
        cls.dir_deg = os.path.join(base, "deg")
        os.makedirs(cls.dir_deg)
        same = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]
        cls.clusters_deg = _write_fixture(
            cls.dir_deg, {0: same, 1: same},
            {0: {"V1": -0.2}, 1: {"V1": 0.5}})           # 型1 V1 正分 > WT 基线
        # 不触发场景: 两型 TOP-8 完全不交 -> 公共集 0 <= 阈 4 -> 阳性
        cls.dir_typ = os.path.join(base, "typ")
        os.makedirs(cls.dir_typ)
        cls.clusters_typ = _write_fixture(
            cls.dir_typ,
            {0: ["A%d" % i for i in range(1, 9)],
             1: ["B%d" % i for i in range(1, 9)]},
            {0: {"A1": -0.1}, 1: {"B1": -0.2}})

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_triggered_degradation(self):
        d, stdout = _run_agent(self.clusters_deg,
                               os.path.join(self.dir_deg, "out"))
        dg = d["degradation"]
        self.assertEqual(dg["mode"], "degraded_universal")
        self.assertFalse(dg["significant"])
        self.assertEqual(dg["n_common"], 8)
        self.assertIsNone(dg["forced"])
        # 判读/交付文案由数值生成: 含判据数值、退化方向与交付口径
        self.assertIn("8 条 > 阈 4", dg["reading"])
        self.assertIn("退化为通用型单骨架推荐", dg["reading"])
        self.assertIn("科学结论依然完整可交付", dg["deliverable_note"])
        self.assertIn("退化模式", stdout)
        # 通用型单骨架 = 全场最优(型1 代表 V1, 正分 0.5 > WT 0 基线)
        self.assertEqual(dg["universal_scaffold"]["desc"], "V1")
        self.assertEqual(dg["universal_scaffold"]["source"],
                         "type_1_representative")
        # 每条候选: 统一通用骨架, 分型仅作诊断不参与推荐
        for row in d["designs"]:
            self.assertIsNone(row["scaffold_type"])
            self.assertEqual(row["dr_desc"], "V1")
            self.assertIn("typing_nearest_type", row)
            self.assertIn("confidence", row)

    def test_not_triggered_typed(self):
        d, _ = _run_agent(self.clusters_typ,
                          os.path.join(self.dir_typ, "out"))
        dg = d["degradation"]
        self.assertEqual(dg["mode"], "typed")
        self.assertTrue(dg["significant"])
        self.assertEqual(dg["n_common"], 0)
        self.assertIn("维持分型推荐", dg["reading"])
        self.assertNotIn("deliverable_note", dg)
        rep = {0: "A1", 1: "B1"}                       # 各型分最高的非 pos-1 变体
        for row in d["designs"]:
            self.assertIn(row["scaffold_type"], (0, 1))
            self.assertEqual(row["dr_desc"], rep[row["scaffold_type"]])
            self.assertNotIn("typing_nearest_type", row)

    def test_forced_on_and_off_are_labeled(self):
        d, _ = _run_agent(self.clusters_typ,
                          os.path.join(self.dir_typ, "out_on"),
                          extra_args=("--degrade", "on"))
        dg = d["degradation"]
        self.assertEqual(dg["mode"], "degraded_universal")
        self.assertEqual(dg["forced"], "on")
        self.assertTrue(dg["significant"])             # 判据实测仍为阳性
        self.assertIn("强制退化", dg["reading"])
        for row in d["designs"]:
            self.assertIsNone(row["scaffold_type"])
            self.assertEqual(row["dr_desc"], "WT")     # 该池全部 <0 -> WT 基线
        d, _ = _run_agent(self.clusters_deg,
                          os.path.join(self.dir_deg, "out_off"),
                          extra_args=("--degrade", "off"))
        dg = d["degradation"]
        self.assertEqual(dg["mode"], "typed")
        self.assertEqual(dg["forced"], "off")
        self.assertFalse(dg["significant"])            # 判据实测仍为阴性
        self.assertIn("强制分型", dg["reading"])


class TestWebappContract(unittest.TestCase):
    def test_design_api_carries_degradation_block(self):
        import crrna_agent_webapp as webapp
        # webapp 与离线判据同源同结论(当前真实数据 = typed)
        self.assertEqual(webapp.DEGRADATION["mode"],
                         agent.interaction_assessment(CLUSTERS)["mode"])
        cli = webapp.app.test_client()
        r = cli.post("/api/design", json={"spacer": SPACER_DNA})
        self.assertEqual(r.status_code, 200)
        d = r.get_json()
        self.assertIn("degradation", d)
        self.assertEqual(d["degradation"]["criterion_id"],
                         "preregistration §A go/no-go")
        self.assertEqual(d["degradation"]["mode"], "typed")
        # typed 模式: 分型字段照常, 无退化诊断键
        self.assertIsNotNone(d["scaffold_type"])
        self.assertNotIn("typing_nearest_type", d)
        self.assertIn("recommended", d)

    def test_panel_api_carries_degradation_block(self):
        import crrna_agent_webapp as webapp
        cli = webapp.app.test_client()
        r = cli.get("/api/panel?target=TP53-R248Q")
        self.assertEqual(r.status_code, 200)
        d = r.get_json()
        self.assertEqual(d["degradation"]["mode"], "typed")


if __name__ == "__main__":
    unittest.main()
