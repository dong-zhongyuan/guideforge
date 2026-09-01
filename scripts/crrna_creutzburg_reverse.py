"""Creutzburg 2020 跨体系反向验证(FnCas12a, NAR gkz1240, PMC7102956)——证据汇聚"路径 A"。

数据来源(全部可溯源):
  1. 14 条 spacer(21nt): 补充材料 Table 1 引物 BG9822-BG9835, 取 GGTTT 后 21nt。
     位点核对(正文 0 起始编号): Sp8/Sp12 差异=位置 10,11,19,20(正文原句),
     Sp12/Sp14 差异=位置 5,7,9, 位置 19 为 K(G/U)——与提取序列逐一吻合。
  2. DR: 反向引物公共尾 revcomp = CTTTAAAT + AATTTCTACTGTTGTAGA(18nt 直接证据);
     19nt 版(3端补 T)按 Cas12a 家族成熟 DR 同源推断, 两口径都跑。
     该构建 DR 与注册表 cas12a2_zeng2026 条目(AATTTCTACTGTTGTAGAT)一致。
  3. 活性: Fig 1D 荧光损失 %(视觉读图, ±3%); Sp12 按正文 "75-80%" 锚定;
     正文确认 Sp1-Sp10 除 Sp8=0% 外均在 60-85%, Sp12 高, Sp13/Sp14 低。

检验问题: 游离态结构指标(spacer 游离度 / DR 茎完整概率 / 5端种子游离度)能否在
跨体系数据上区分好/坏 guide 并复现该文失活机制(Sp8 前 11nt 与 5端 repeat
互作形成失活折叠)。好/坏阈值 60%(正文区间下沿)。

诚实边界:
  - n=14, 活性为读图近似值, 非原始数字数据;
  - 只报单指标 Spearman 与 AUC, 不做复合打分权重拟合(14 点必然过拟合);
  - 体系为 FnCas12a 大肠杆菌质粒靶向, 与 Cas12a2 哺乳动物杀伤场景不同,
    相关性成立只说明指标方向跨体系一致, 不构成 Cas12a2 活性预测器。
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

SPACERS_DNA_21 = {
    "Sp1":  "ATCGATGGGAAACCTTACCCTC",
    "Sp2":  "AGACTTACCAATGAGCACGTG",
    "Sp3":  "AAATATTCCCCAGGTTCAGTG",
    "Sp4":  "AACACACTGCAATTCAGGTTG",
    "Sp5":  "ACCAGTCCATCATTGTAGTGG",
    "Sp6":  "AACGCCTACGAGATCATTGTG",
    "Sp7":  "ATGCTCCCGTAACACTAGCTG",
    "Sp8":  "ACCAGTCTAAATCGCTTACTG",
    "Sp9":  "AACTATATGCACCGGAGCTTG",
    "Sp10": "ACACAATTACCGGTCCCAATG",
    "Sp11": "ACCAGTCTAACCCGCTTACGT",
    "Sp12": "ACCAGTCTAACACGCTTACGT",
    "Sp13": "ACCAGCCAACCCCGCTTACGT",
    "Sp14": "ACCAGCCAACCACGCTTACGT",
}

ACTIVITY = {  # Fig 1D 荧光损失 % 近似值, Sp12 锚定正文 75-80
    "Sp1": 71, "Sp2": 73, "Sp3": 69, "Sp4": 72, "Sp5": 75, "Sp6": 78,
    "Sp7": 62, "Sp8": 0, "Sp9": 83, "Sp10": 69, "Sp11": 57, "Sp12": 79,
    "Sp13": 34, "Sp14": 29,
}

DR18_RNA = core.to_rna("AATTTCTACTGTTGTAGA")     # 引物直接证据(18nt)
DR19_RNA = core.to_rna("AATTTCTACTGTTGTAGAT")    # 家族同源补 3端T
PRE5_RNA = core.to_rna("TTTAAAT")                # 前体 crRNA 中成熟 DR 5侧残留 repeat
GOOD_TH = 60  # 好 guide 阈值: 正文 "60-85%" 区间下沿


def dr_stem_pairs(dr_rna):
    """DR 单独折叠的茎配对(0基 DR 内坐标)。"""
    db, _ = core.fold(dr_rna)
    stack, pairs = [], []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.append((stack.pop(), i))
    return pairs


def seed5_unpaired(full_rna, prefix_len):
    """spacer 5端 5nt(Cas12a 种子区, 论文口径 nucleotides 1-5)平均未配对概率。"""
    fc = RNA.fold_compound(full_rna)
    fc.pf()
    bpp = fc.bpp()
    probs = []
    for k in range(prefix_len, prefix_len + 5):
        p_pair = sum(bpp[k + 1][j + 1] for j in range(len(full_rna)) if j != k)
        probs.append(max(0.0, 1.0 - p_pair))
    return float(np.mean(probs))


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def perm_p(vals, acts, rho, n=20000, seed=0):
    rng = np.random.default_rng(seed)  # 只建一次, 每轮不同置换
    hits = 0
    for _ in range(n):
        if abs(spearman(rng.permutation(vals), acts)) >= abs(rho):
            hits += 1
    return (hits + 1) / (n + 1)


def auc(good_vals, bad_vals, higher_better=True):
    hits = 0
    for g in good_vals:
        for b in bad_vals:
            if (g > b) if higher_better else (g < b):
                hits += 1
    return hits / (len(good_vals) * len(bad_vals))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "creutzburg_reverse_validation.json"))
    args = ap.parse_args()

    names = sorted(SPACERS_DNA_21, key=lambda s: int(s[2:]))
    acts = np.array([ACTIVITY[n] for n in names], dtype=float)
    good = acts >= GOOD_TH

    contexts = {
        "mature18": (DR18_RNA, DR18_RNA),
        "mature19": (DR19_RNA, DR19_RNA),
        "pre24": (PRE5_RNA + DR18_RNA, DR18_RNA),  # 前体: 残留 repeat + 成熟 DR + spacer
    }

    report = {"source": "Creutzburg 2020 NAR gkz1240 PMC7102956", "n": len(names),
              "good_threshold_pct": GOOD_TH, "contexts": {}, "spacers": SPACERS_DNA_21,
              "activity_fig1d_approx": ACTIVITY,
              "evidence_checks": {
                  "Sp8_Sp12_diff_positions_0based": [10, 11, 19, 20],
                  "Sp12_Sp14_diff_positions_0based": [5, 7, 9],
                  "note": "与正文原句逐一吻合; Sp12 活性按正文 75-80% 锚定"}}

    print("== 数据核对 ==")
    d = lambda a, b: [k for k in range(21) if a[k] != b[k]]
    s8, s12, s14 = SPACERS_DNA_21["Sp8"], SPACERS_DNA_21["Sp12"], SPACERS_DNA_21["Sp14"]
    print("Sp8/Sp12 差异(0基):", d(s8, s12), "| 正文: 10,11,19,20")
    print("Sp12/Sp14 差异(0基):", d(s12, s14), "| 正文: 5,7,9")
    print("位置19: Sp8=%s Sp12=%s (正文: K=G/U)" % (s8[19], s12[19]))

    for ctx_name, (prefix, dr_for_pairs) in contexts.items():
        pairs = dr_stem_pairs(dr_for_pairs)
        rows = []
        for nm in names:
            full = prefix + core.to_rna(SPACERS_DNA_21[nm])
            _, sp_up, _ = core.pf_stats(full, len(prefix))
            stem_p = core.stem_intact_prob(full, [(len(prefix) - len(dr_for_pairs) + i,
                                                   len(prefix) - len(dr_for_pairs) + j)
                                                  for i, j in pairs])
            seed5 = seed5_unpaired(full, len(prefix))
            rows.append({"name": nm, "activity": ACTIVITY[nm],
                         "spacer_up": round(sp_up, 4),
                         "stem_prob": round(stem_p, 6),
                         "seed5_up": round(seed5, 4)})
        print("\n== 上下文 %s (DR茎配对 %s) ==" % (ctx_name, pairs))
        print("%-5s %4s  %7s %8s %7s" % ("name", "act", "sp_up", "stem_P", "seed5"))
        for r in rows:
            print("%-5s %4d  %7.4f %8.5f %7.4f" % (
                r["name"], r["activity"], r["spacer_up"], r["stem_prob"], r["seed5_up"]))

        stats = {}
        for feat in ("spacer_up", "stem_prob", "seed5_up"):
            vals = np.array([r[feat] for r in rows])
            rho = spearman(vals, acts)
            p = perm_p(vals, acts, rho)
            a = auc(vals[good].tolist(), vals[~good].tolist())
            stats[feat] = {"spearman": round(rho, 3), "perm_p": round(p, 4), "auc_good_bad": round(a, 3)}
            print("%-10s rho=%+.3f  perm_p=%.4f  AUC(good n=%d vs bad n=%d)=%.3f" % (
                feat, rho, p, good.sum(), (~good).sum(), a))
        report["contexts"][ctx_name] = {"rows": rows, "stats": stats}

    # 机制复现: Sp8 的 MFE 折叠中 spacer 与 DR 区的配对数(正文: 前 11nt 互作)
    print("\n== 机制复现 (mature19, Sp8 vs Sp12) ==")
    mech = {}
    for nm in ("Sp8", "Sp12"):
        full = DR19_RNA + core.to_rna(SPACERS_DNA_21[nm])
        db, mfe = core.fold(full)
        dr_len = len(DR19_RNA)
        cross = sum(1 for k in range(dr_len, len(full)) if db[k] != ".")
        dr_paired = sum(1 for k in range(dr_len) if db[k] != ".")
        print("%s: MFE=%+.1f  spacer配对nt=%d  DR区配对nt=%d" % (nm, mfe, cross, dr_paired))
        print("  %s" % db)
        mech[nm] = {"mfe": round(mfe, 2), "spacer_paired_nt": cross,
                    "dr_paired_nt": dr_paired, "db": db}
    report["mechanism_mfe"] = mech

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print("\n输出 -> %s" % args.out)


if __name__ == "__main__":
    main()
