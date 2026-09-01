"""同源 DR 位点保守性先验(策划案 V3 §4.1, 2026-09-01)。

序列集(全部为仓库已溯源条目, 不引入未核验硬编码):
  SuCas12a2  ATTTCTACTATTGTAGAT  注册表 cas12a2(8D4A 结构口径)
  Cas12a2    AATTTCTACTGTTGTAGAT 注册表 cas12a2_zeng2026(Zeng 2026)
  FnCas12a   AATTTCTACTGTTGTAGA  Creutzburg 2020 反向引物推得(数据集)
  LbCas12a   AATTTCTACTATTGTAGAT Teng 2019 crRNA1(teng4n96 脚本)
方法: mafft 不可用, DR 高度同源, 以保守锚 TTTCTACT 左对齐 + TAGAT 右对齐的
锚定对齐; 逐位一致性即保守性。诚实边界: n=4(3 直系同源), 小样本定量,
与 Dmytrenko 2023 Fig1c 跨家族定性结论(3' 保守/loop 可变)做一致性核对,
不替代其结论; --cons3-window 的窗宽仍为项目设定。
输出: data/dr_conservation.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

DRS = [(d["name"], d["dr_dna"], d["provenance"])
        for d in json.load(open(os.path.join(ROOT, "data", "dr_homologs.json"), encoding="utf-8"))]
ANCHOR_L = "TTTCTACT"
ANCHOR_R = "GTAGA"


def anchored(name, seq, prov):
    l = seq.index(ANCHOR_L)
    r = seq.rindex(ANCHOR_R)
    core = seq[l:r + len(ANCHOR_R)]
    return {"name": name, "full": seq, "core_start_in_full": l,
            "core": core, "provenance": prov}


def main():
    aln = [anchored(n, s, p) for n, s, p in DRS]
    width = min(len(a["core"]) for a in aln)
    cons = []
    for i in range(width):
        col = [a["core"][i] for a in aln]
        top = max(set(col), key=col.count)
        cons.append({"pos_1based": i + 1, "consensus": top,
                     "identity": round(col.count(top) / len(col), 2)})
    print("%-24s %s" % ("DR(锚定核心区)", "序列"))
    for a in aln:
        print("%-24s %s" % (a["name"], a["core"]))
    print("\n位点一致性:")
    print("  " + " ".join("%s%.0f" % (c["consensus"], c["identity"] * 10)
                          for c in cons))
    w3 = [c["identity"] for c in cons if c["pos_1based"] > width - 5]
    rest = [c["identity"] for c in cons if c["pos_1based"] <= width - 5]
    print("3' 末 5nt 平均一致性 %.2f vs 其余 %.2f" % (
        sum(w3) / len(w3), sum(rest) / len(rest)))

    out = os.path.join(ROOT, "data", "dr_conservation.json")
    json.dump({"alignment": aln, "conservation": cons,
               "n_sequences": len(DRS),
               "method": "TTTCTACT/GTAGAT 锚定对齐(mafft 不可用; DR 高度同源)",
               "limitation": "n=4 小样本定量, 仅为 cons3 窗的仓库内定量佐证, "
                             "不替代 Dmytrenko 2023 Fig1c 跨家族结论",
               "cons3_window_5nt_mean_identity": round(sum(w3) / len(w3), 3),
               "rest_mean_identity": round(sum(rest) / len(rest), 3)},
              open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
