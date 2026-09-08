# -*- coding: utf-8 -*-
"""突变类型判定(最长 ORF 法)契约测试(策划案 V3 §4.4 模块一 / 表2 APC 入选理由
对齐整改, 2026-09-07, 任务⑩)。

锁定五件事:
1. ORF 规则边界: 最长 ORF 选择的显式 tie-break(长度->起点->框序)、无框内
   终止子时的部分密码子处理(末端不足 3 nt 丢弃, terminated=False);
2. 四靶标 WT 序列 ORF 合理性: 最长 ORF 即各基因经典 CDS(TP53 393 aa /
   KRAS 189 aa / APC 2843 aa, 起 ATG 止终止子, 内部无终止子);
3. 四靶标类型判定: TP53-R248Q/KRAS-G12D/TP53-R273H=错义(ORF 长度不变、
   单氨基酸替换), APC-Q1328x=无义(提前终止子密码子 1328,
   蛋白 2843->1327 aa)——APC 序列经本模块核查后修订: 旧 Q1312x 口径
   (3946 位编辑)在真实阅读框为 A1296V 错义(系 panel 脚本以 5'UTR 上游
   ATG 锚定所致), 2026-09-08 拍板正式切换 Q1328*(见
   scripts/crrna_mutation_typing.py docstring);
4. 合成用例规则边界: 无义/移码(+1/-1)/同义/框内 indel/UTR/起始密码子
   丢失/复杂事件拒绝;
5. 接入契约: 四个 design.json 的 mutation.typing 与每行 target_features
   携带该维; webapp 与离线共用同一判定模块与 TARGET_FEATURES 定义,
   /api/panel 响应携带同源证据。
"""
import json
import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import crrna_mutation_typing as mt  # noqa: E402

# 靶标 -> (mut fasta, wt fasta); TP53-R248Q 的序列对在 data/ 根(历史布局),
# 其余在 data/agent/(panel 脚本产物)
PANEL = {
    "tp53_r248q": ("data/tp53_r248q_mrna_NM000546.fa", "data/tp53_mrna_NM000546.fa"),
    "kras_g12d": ("data/agent/kras_g12d.fa", "data/agent/kras_g12d_wt.fa"),
    "tp53_r273h": ("data/agent/tp53_r273h.fa", "data/agent/tp53_r273h_wt.fa"),
    "apc_q1328x": ("data/agent/apc_q1328x.fa", "data/agent/apc_q1328x_wt.fa"),
}
EXPECTED = {"tp53_r248q": ("missense", "R248Q", 248),
            "kras_g12d": ("missense", "G12D", 12),
            "tp53_r273h": ("missense", "R273H", 273),
            "apc_q1328x": ("nonsense", "Q1328*", 1328)}
# 各基因经典 CDS 蛋白长度(RefSeq 口径), WT ORF 合理性锚点
WT_PROTEIN_AA = {"tp53_r248q": 393, "kras_g12d": 189, "tp53_r273h": 393,
                 "apc_q1328x": 2843}


def _read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


class TestOrfRules(unittest.TestCase):
    """ORF 规则边界: tie-break / 部分密码子 / 无终止子"""

    def test_tie_break_length_then_start_then_frame(self):
        # 两条等长 ORF(6 nt): frame0 起点 0 与 frame1 起点 4 -> 起点靠前者胜
        seq = "ATGTAA" + "X".replace("X", "A") + "ATGTAA"  # 0: ATG TAA; 5: ATG TAA
        o = mt.longest_orf(seq)
        self.assertEqual(o["start0"], 0)
        # 更长 ORF 优先于更早起点
        seq2 = "ATGTAA" + "AA" + "ATGAAAAAATAA"  # 后者 12 nt > 前者 6 nt
        o2 = mt.longest_orf(seq2)
        self.assertEqual(o2["start0"], 8)
        self.assertEqual(o2["length_nt"], 12)     # ATG+AAA+AAA+TAA

    def test_unterminated_orf_discards_partial_codon(self):
        # 无框内终止子: ATG + 4 密码子 + 2 nt 残余 -> 末端部分密码子丢弃
        seq = "ATGAAAAAA" + "AAA" + "AA" + "AT"   # ATG + AAA AAA AAA AA + AT(残余)
        o = mt.longest_orf(seq)
        self.assertFalse(o["terminated"])
        self.assertIsNone(o["stop0"])
        self.assertEqual(o["length_nt"], 15)      # ATG + 4 完整密码子
        pub = mt._orf_public(o)
        self.assertEqual(pub["protein_aa"], 5)

    def test_no_atg_returns_none(self):
        self.assertIsNone(mt.longest_orf("AAAAAA"))
        with self.assertRaises(ValueError):
            mt.classify("AAAAAA", "AAATAA")


class TestPanelWTOrfSane(unittest.TestCase):
    """四靶标 WT 最长 ORF = 经典 CDS: M 起、终止子止、内部无终止子"""

    def test_wt_orf_is_canonical_cds(self):
        for key, (_, wtf) in PANEL.items():
            wt = _read(wtf)
            o = mt.longest_orf(wt)
            self.assertIsNotNone(o, key)
            prot = mt.translate(wt, o["start0"], o["end0"])
            self.assertTrue(prot.startswith("M"), key)
            self.assertTrue(prot.endswith("*"), key)
            self.assertNotIn("*", prot[:-1], key)   # 内部无提前终止子
            self.assertEqual(o["length_nt"] // 3 - 1, WT_PROTEIN_AA[key], key)


class TestPanelClassification(unittest.TestCase):
    """四靶标判定: 3 错义 + APC 无义(数值证据逐项锁定)"""

    def test_types_and_evidence(self):
        for key, (mutf, wtf) in PANEL.items():
            t = mt.classify(_read(wtf), _read(mutf))
            mtype, aa_change, codon = EXPECTED[key]
            self.assertEqual(t["mutation_type"], mtype, key)
            self.assertEqual(t["mutation_type_label"],
                             mt.TYPE_LABELS[mtype], key)
            self.assertEqual(t["aa_change"], aa_change, key)
            self.assertEqual(t["codon_no"], codon, key)
            self.assertEqual(t["wt_orf"]["protein_aa"], WT_PROTEIN_AA[key], key)
            self.assertIn(str(codon), t["reading"], key)
            if mtype == "missense":
                # ORF 起止与长度不变, 蛋白长度不变
                self.assertEqual(t["mut_orf_same_start"]["start_1based"],
                                 t["wt_orf"]["start_1based"], key)
                self.assertEqual(t["mut_orf_same_start"]["end_1based"],
                                 t["wt_orf"]["end_1based"], key)
                self.assertEqual(t["orf_length_change_aa"], 0, key)
                self.assertIsNone(t["premature_stop"], key)
            else:  # nonsense
                ps = t["premature_stop"]
                self.assertEqual(ps["codon_no"], 1328, key)
                self.assertEqual(ps["transcript_pos_1based"], 4041, key)
                self.assertEqual((ps["wt_codon"], ps["mut_codon"]),
                                 ("CAG", "TAG"), key)
                # ORF 变短: 2843 aa -> 1327 aa(提前终止), 截短 >50%
                self.assertEqual(t["mut_orf_same_start"]["protein_aa"], 1327, key)
                self.assertEqual(t["orf_length_change_aa"], 1327 - 2843, key)
                self.assertGreater(t["truncated_fraction"], 0.5, key)
                self.assertIn("无义", t["reading"], key)


class TestSyntheticBoundaries(unittest.TestCase):
    """合成序列验证规则边界(人为构造, 与真实 panel 无关)"""

    WT = "GCGGCC" + "ATG" + "AAA" * 20 + "CAG" + "AAA" * 10 + "TAA" + "GGG"

    def test_nonsense(self):
        t = mt.classify(self.WT, self.WT.replace("CAG", "TAG", 1))
        self.assertEqual(t["mutation_type"], "nonsense")
        self.assertEqual(t["premature_stop"]["codon_no"], 22)
        self.assertEqual(t["mut_orf_same_start"]["protein_aa"], 21)
        self.assertLess(t["orf_length_change_aa"], 0)

    def test_frameshift_insertion_and_deletion(self):
        for mut in (self.WT[:10] + "A" + self.WT[10:],       # +1 插入
                    self.WT[:12] + self.WT[13:]):            # -1 缺失
            t = mt.classify(self.WT, mut)
            self.assertEqual(t["mutation_type"], "frameshift")
            self.assertNotEqual((t["mut_orf_same_start"] or {}).get("stop_1based"),
                                t["wt_orf"]["stop_1based"])
            self.assertIn("移码", t["reading"])

    def test_synonymous(self):
        t = mt.classify(self.WT, self.WT.replace("CAG", "CAA", 1))  # Q>Q
        self.assertEqual(t["mutation_type"], "synonymous")
        self.assertEqual(t["orf_length_change_aa"], 0)

    def test_inframe_indel(self):
        t = mt.classify(self.WT, self.WT[:10] + "AAA" + self.WT[10:])  # +3
        self.assertEqual(t["mutation_type"], "inframe_indel")
        self.assertEqual(t["orf_length_change_aa"], 1)

    def test_utr_non_coding(self):
        t = mt.classify(self.WT, "T" + self.WT[1:])          # 5'UTR 替换
        self.assertEqual(t["mutation_type"], "non_coding")
        t2 = mt.classify(self.WT, self.WT[:-3] + "CCC")      # 3'UTR(终止子后)
        self.assertEqual(t2["mutation_type"], "non_coding")

    def test_start_lost(self):
        t = mt.classify(self.WT, self.WT[:6] + "TTG" + self.WT[9:])
        self.assertEqual(t["mutation_type"], "start_lost")
        self.assertIsNone(t["mut_orf_same_start"])

    def test_complex_event_rejected(self):
        with self.assertRaises(ValueError):
            mt.classify(self.WT, self.WT[:8] + "GC" + self.WT[12:])

    def test_reading_generated_from_numbers(self):
        # 判读文字含关键数值且方向与类型一致(不手写结论)
        t = mt.classify(self.WT, self.WT.replace("CAG", "TAG", 1))
        self.assertIn("22", t["reading"])
        self.assertIn("无义", t["reading"])
        self.assertIn("提前终止", t["reading"])


class TestDesignArtifacts(unittest.TestCase):
    """四个 panel design.json 携带 mutation.typing 与行级特征维"""

    def test_design_json_carry_typing(self):
        for key in PANEL:
            d = json.load(open(os.path.join(DATA, "agent",
                                            key + ".design.json"),
                               encoding="utf-8"))
            t = (d.get("mutation") or {}).get("typing")
            self.assertIsNotNone(t, key)
            self.assertEqual(t["mutation_type"], EXPECTED[key][0], key)
            self.assertTrue(d["designs"], key)
            for row in d["designs"]:
                tf = row.get("target_features")
                self.assertIsNotNone(tf, "%s 行缺 target_features" % key)
                self.assertEqual(tf["mutation_type"], EXPECTED[key][0], key)
                self.assertEqual(tf["orf_wt_protein_aa"], WT_PROTEIN_AA[key], key)
                if EXPECTED[key][0] == "nonsense":
                    self.assertEqual(tf["premature_stop_codon"], 1328, key)
                else:
                    self.assertIsNone(tf["premature_stop_codon"], key)


class TestSelectorWebappContract(unittest.TestCase):
    """选型契约与 webapp↔离线同源"""

    def test_target_features_contract(self):
        import crrna_train_selector as train
        import crrna_agent_webapp as webapp
        self.assertIn("mutation_type", train.TARGET_FEATURES)
        # 同一对象(import 共用), 不是两份字面量
        self.assertIs(webapp.TARGET_FEATURES, train.TARGET_FEATURES)
        # webapp 与离线共用同一判定模块(非各自实现)
        self.assertIs(webapp.mtyping, mt)

    def test_panel_api_carries_mutation_typing(self):
        import crrna_agent_webapp as webapp
        cli = webapp.app.test_client()
        for name, key in webapp.PANEL.items():
            r = cli.get("/api/panel?target=" + name)
            self.assertEqual(r.status_code, 200, name)
            d = r.get_json()
            mtp = d.get("mutation_typing")
            self.assertIsNotNone(mtp, name)
            self.assertEqual(mtp["mutation_type"], EXPECTED[key][0], name)
            # 与离线归档同源同结论: 响应证据 == design.json 的 mutation.typing
            arch = json.load(open(os.path.join(DATA, "agent",
                                               key + ".design.json"),
                                  encoding="utf-8"))["mutation"]["typing"]
            self.assertEqual(mtp["reading"], arch["reading"], name)


if __name__ == "__main__":
    unittest.main()
