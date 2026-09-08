# -*- coding: utf-8 -*-
"""策划案 V3 表2 面板与骨架族一致性契约测试(2026-09-04)。

锁定三件事(round-3 meta W6 + V3 对齐整改):
1. IVT 模板/订单表的 4 靶标 == V3 表2(TP53-R248Q / KRAS-G12D / TP53-R273H /
   APC-Q1328x); KRAS-G12C 只保留干实验证据, 不入湿实验矩阵;
2. 两文件的 8 骨架 == orientation_library 当前 best 去重集(含 round-3 去混杂
   补偿臂 DR AATTTCTACAGGTGTAGAG)——杜绝"订单表编码已退役分子"重演;
3. 订单表 RNA 序列 == 模板 DNA 的 T->U 转写, 32 行一一对应。
"""
import csv
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")

V3_TARGETS = {"TP53-R248Q", "KRAS-G12D", "TP53-R273H", "APC-Q1328x"}
NEW_COMP_ARM_DR = "AATTTCTACAGGTGTAGAG"  # zengDR+T10A/T12G/T19G (round-3 去混杂)
RETIRED_ARM_DR = "AATTTCTACTCTTCTACAT"   # round-3 前茎破坏混杂臂
APC_SPACER = "TGACACTGCTGGAACTTCGCTCAC"  # Q1328*(真实阅读框 MCR 首个 CAG, 4041 位)位点最优 spacer, 2026-09-08 口径切换


def _load_csv(name):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _library_scaffolds():
    lib = json.load(open(os.path.join(DATA, "orientation_library.json"),
                         encoding="utf-8"))
    seen, rows = {}, []
    for p in lib["orientations"]:
        b = p["best"]
        if b["desc"] in seen:
            continue
        seen[b["desc"]] = 1
        rows.append((b["desc"], b["construct_dna"][:19]))
    wt = json.load(open(os.path.join(DATA, "tp53_r248q_zengdr.v6.top.json"),
                        encoding="utf-8"))["top"][0]
    return [("WT", wt["dr_seq"])] + rows


class TestIvtPanelSync(unittest.TestCase):
    def test_targets_match_v3_table2(self):
        rows = _load_csv("ivt_round1_template.csv")
        self.assertEqual({r["target"] for r in rows}, V3_TARGETS)
        self.assertEqual(len(rows), 32)
        apc = [r for r in rows if r["target"] == "APC-Q1328x"]
        self.assertEqual(len(apc), 8)
        self.assertTrue(all(r["spacer_dna"] == APC_SPACER for r in apc))

    def test_scaffolds_match_library_incl_new_comp_arm(self):
        expect = dict(_library_scaffolds())
        self.assertEqual(len(expect), 8)
        for name in ("ivt_round1_template.csv", "ivt_round1_order_sheet.csv"):
            rows = _load_csv(name)
            self.assertEqual({r["scaffold_desc"] for r in rows}, set(expect))
        tmpl = _load_csv("ivt_round1_template.csv")
        for r in tmpl:
            self.assertTrue(
                r["construct_dna"].startswith(expect[r["scaffold_desc"]]),
                "%s 行的 DR 与当前库不一致" % r["combo_id"])
        arm = [r for r in tmpl if r["scaffold_desc"] == "B_break_compensate"]
        self.assertTrue(all(r["construct_dna"].startswith(NEW_COMP_ARM_DR)
                            for r in arm))
        self.assertFalse(any(r["construct_dna"].startswith(RETIRED_ARM_DR)
                             for r in tmpl))

    def test_order_sheet_rna_matches_template_dna(self):
        tmpl = {(r["target"], r["scaffold_desc"]): r
                for r in _load_csv("ivt_round1_template.csv")}
        sheet = _load_csv("ivt_round1_order_sheet.csv")
        self.assertEqual(len(sheet), 32)
        for r in sheet:
            key = (r["target"], r["scaffold_desc"])
            self.assertIn(key, tmpl)
            self.assertEqual(r["sequence_5to3_RNA"],
                             tmpl[key]["construct_dna"].replace("T", "U"))
            self.assertEqual(int(r["length_nt"]),
                             len(tmpl[key]["construct_dna"]))
            self.assertEqual(r["spacer_dna_ref"], tmpl[key]["spacer_dna"])

    def test_demo_panels_match_v3_table2(self):
        """演示层(webapp PANEL / 离线快照 PANEL)与湿实验矩阵同四靶;
        KRAS-G12C 退役后不得再出现在任何对外演示入口。"""
        import ast
        for rel, is_dict in (("scripts/crrna_agent_webapp.py", True),
                             ("scripts/crrna_demo_snapshot.py", False)):
            src = open(os.path.join(ROOT, rel), encoding="utf-8").read()
            tree = ast.parse(src)
            panel = None
            for node in ast.walk(tree):
                if (isinstance(node, ast.Assign)
                        and any(isinstance(t, ast.Name) and t.id == "PANEL"
                                for t in node.targets)):
                    panel = ast.literal_eval(node.value)
            self.assertIsNotNone(panel, rel + " 未找到 PANEL 定义")
            names = set(panel) if is_dict else {n for n, _ in panel}
            self.assertEqual(names, V3_TARGETS, rel + " panel 与 V3 表2 不符")
            self.assertNotIn("KRAS-G12C", names)


if __name__ == "__main__":
    unittest.main()
