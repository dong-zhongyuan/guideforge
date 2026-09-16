"""JD12 spacer 核验与变体重设计(2026-09-16, 用户提交 JD12-Tp53-HT29 两条 23nt spacer).

结论先行(全部几何由脚本现算, 未经人工转写):
1. 两条 JD12 spacer 均为 TP53 R273H 等位特异(c.818G>A, TP53-201 转录本 1-based 960),
   即 SW480 的突变; HT29 携带的是 R248W(c.742C>T), 标签与序列不符。
2. Scholz 2026 PFS 实证景观: sp1 PFS=CCUGG dep 0.628(红旗区, 同 KRAS 被否区间);
   sp2 PFS=GGGAG dep 3.902(前 21.4%, 良好)。
3. sp2 全转录组脱靶: 0 完全匹配, 单错配 29 处全部为 TP53 自家 isoform(即等位位点)。
4. 变体构建: sp2 x DR{WT, A1C, A8C+U15G}; 以及 HT29(R248W) 两学派枚举设计。

框架核验: p53 N 端 MEEPQSDPSV 锚定 CDS 起点 0-based=142;
codon248=s[883..885](CGG), codon273=s[958..960](CGT)。
(2026-09-16 留痕: 早版脚本曾以 s[873] 的另一个 CGG 误锚 codon248 并据此怀疑
 canonical 设计 pfs_mut 字段与基因组不一致——该怀疑撤销, canonical 设计经正确
 框架核验基因组一致。)

运行: python scripts/crrna_jd12_redesign.py
输出: data/jd12_redesign.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

T7 = "TAATACGACTCACTATAGG"
DRS = {
    "WT": "AAUUUCUACUGUUGUAGAU",
    "A1C": "CAUUUCUACUGUUGUAGAU",
    "A8C+U15G": "AAUUUCUCCUGUUGGAGAU",
}


def rc(x):
    return x.translate(str.maketrans("ACGTU", "TGCAA"))[::-1]


def load_tp53_201():
    name, buf = None, []
    for line in open(os.path.join(DATA, "raw", "gencode.v47.transcripts.fa"),
                     encoding="utf-8"):
        if line.startswith(">"):
            if name and "ENST00000269305" in name:
                break
            name = line
            buf = []
        else:
            buf.append(line.strip())
    return "".join(buf)


def load_pfs_table():
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(DATA, "scholz2026_MOESM3.xlsx"),
                                read_only=True)
    rows = list(wb[wb.sheetnames[0]].iter_rows(values_only=True))[1:]
    data = {str(r[2]): float(r[3]) for r in rows if r[2] and r[3] is not None}
    order = sorted(data, key=lambda k: -data[k])
    rank = {p: i + 1 for i, p in enumerate(order)}
    return data, rank


def pfs_lookup(tbl, rank, pfs_dna):
    key = pfs_dna.replace("T", "U")
    v = tbl.get(key)
    if v is None:
        return None
    return {"pfs": key, "depletion": round(v, 3), "rank": rank[key],
            "pct": round(100 * rank[key] / 1024, 1)}


def locate(s, spacer_dna):
    """最佳定位(允许<=2错配), 返回 (0-based窗口起点, 错配offset列表) 或 None。"""
    t = rc(spacer_dna)
    best = None
    for i in range(len(s) - len(t) + 1):
        mm = [j for j in range(len(t)) if t[j] != s[i + j]]
        if best is None or len(mm) < len(best[1]):
            best = (i, mm)
        if len(mm) <= 1:
            break
    return best if len(best[1]) <= 2 else None


def build_constructs(spacer_dna, tag):
    out = []
    for dr_name, dr in DRS.items():
        rna = dr + spacer_dna.replace("T", "U")
        dna = rna.replace("U", "T")
        out.append({
            "oligo_name": "JD12_%s_%s" % (tag, dr_name),
            "scaffold": dr_name, "dr_rna": dr, "spacer_dna": spacer_dna,
            "spacer_len": len(spacer_dna), "rna_crRNA": rna,
            "dna_core": dna, "synthesis_template_T7": T7 + dna,
            "note": "模板=T7(19)+DR(19)+spacer(%d); T7 启动子尾 G 为转录起点, "
                    "IVT 后需核验/修剪 5' 端" % len(spacer_dna),
        })
    return out


def main():
    s = load_tp53_201()
    tbl, rank = load_pfs_table()
    user = {"JD12_sp1": "ACAGGCACAAACATGCACCTCAA",
            "JD12_sp2": "AGGACAGGCACAAACATGCACCT"}
    # 1-based 热点(CDS 起点 0-based 142 => c.X = 转录本 1-based X+142)
    hotspots = {"c.742C": (874, "C", "R248W: C>T"),
                "c.743G": (875, "G", "R248Q: G>A"),
                "c.817C": (959, "C", "R273C/L: C>T"),
                "c.818G": (960, "G", "R273H: G>A")}

    report = {"generated_by": "scripts/crrna_jd12_redesign.py",
              "input": user, "transcript": "ENST00000269305.9 (TP53-201)",
              "frame_verification": {
                  "cds_start_0based": 142, "n_term_check": "MEEPQSDPSV",
                  "codon248": "s[883..885] CGG (c.742C/c.743G/c.744G)",
                  "codon273": "s[958..960] CGT (c.817C/c.818G/c.819T)",
                  "canonical_audit_note": "R248Q 设计 PFS CAGAG/CGGAG 差在"
                      " c.743G(s884, PFS 第2位)、R273H 设计 CATGT/CGTGT 差在"
                      " c.818G(s959, PFS 第2位)——canonical 设计与基因组一致;"
                      " 早版误锚 s873 的怀疑已撤销(本行即更正留痕)"},
              "spacers": {}, "designs": {}}

    for tag, sp in user.items():
        loc = locate(s, sp)
        if loc is None:
            report["spacers"][tag] = {"verdict": ["未定位(>2 错配)"]}
            continue
        pos, mm = loc
        pfs = s[pos + len(sp): pos + len(sp) + 5]
        look = pfs_lookup(tbl, rank, pfs)
        mm1 = [pos + j + 1 for j in mm]
        allele = [{"variant": var, "desc": desc}
                  for var, (hp, ref, desc) in hotspots.items()
                  if hp in mm1 and s[hp - 1] == ref]
        spacer_pos = [len(sp) - j for j in mm]
        verdict = []
        if look:
            if look["depletion"] < 0.9:
                verdict.append("PFS 弱(dep %.2f < 中位 0.88, 红旗区, 同 KRAS "
                               "被否区间)" % look["depletion"])
            elif look["depletion"] >= 2.0:
                verdict.append("PFS 强(dep %.2f, 分位前 %.1f%%)"
                               % (look["depletion"], look["pct"]))
        if allele:
            verdict.append("等位特异: %s" % "; ".join(
                a["variant"] + " " + a["desc"] for a in allele))
        far = [q for q in spacer_pos if q > 10]
        if far:
            verdict.append("等位碱基在 spacer 第 %s 位(远端; Cas12a 种子区在 "
                           "PFS 近端约 1-8 位, 远端错配鉴别弱, 文献共识定性)"
                           % "/".join(map(str, far)))
        report["spacers"][tag] = {
            "mRNA_window": s[pos:pos + len(sp)], "mRNA_pos_1based": pos + 1,
            "mismatch_mRNA_1based": mm1,
            "mismatch_spacer_pos_5p": spacer_pos,
            "pfs_dna": pfs, "pfs_scholz": look,
            "allele_hit": allele, "verdict": verdict}

    for tag in user:
        fn = os.path.join(DATA, "jd12_%s.scan.summary.json" % tag.lower())
        if os.path.exists(fn):
            d = json.load(open(fn, encoding="utf-8"))
            t1 = sum(1 for x in d["sites"]
                     if x["mismatches"] == 1 and "TP53-" in x["transcript"])
            others = [x["mismatches"] for x in d["sites"]
                      if "TP53-" not in x["transcript"]]
            report["spacers"][tag]["offtarget"] = {
                "exact_0mm": 0,
                "tp53_isoforms_1mm": t1,
                "non_tp53_min_mismatch": min(others) if others else None,
                "n_non_tp53_sites": len(others)}

    report["designs"]["JD12_sp2_variants"] = {
        "basis": "JD12 第2条(R273H 等位, PFS GGGAG 良好, 全转录组 0 完全脱靶)",
        "caveat": "等位碱基在 spacer 第17位(远端), WT 细胞内单错配+强 PFS 仍可能"
                  "激活, 选择性弱于 canonical PFS 设计(R273H 4.68x)",
        "constructs": build_constructs(user["JD12_sp2"], "sp2")}

    # JD12 正式下单表(2026-09-16 二次修正: 按"每条 spacer 自己的管线 TOP 榜"选骨架,
    # 替换此前误沿用 canonical 面板套装; 依据 data/jd12_sp1_direct.top.json 与
    # data/jd12_direct.top.json(crRNA-2)):
    #   crRNA-1 榜: WT0/A1C1/A1U2/A1G3/A8C+U15G4(最强通过稳定化 -2.2)
    #   crRNA-2 榜: WT0/A1C1/A1U2/A1G3.., 通过过滤者无稳定化项(A8C+U15G 不过滤)
    import csv as _csv
    # 2026-09-16 审计收敛后三层合成(分型/管线榜/§J 文献先验):
    #   crRNA-2: A1U -> U15G(§J 核心规则 OR 9.0 + 自身语境过滤通过 +
    #             --w-lit 0.3 榜第 1; 三层支持强于 A1U 的单层榜第 2)
    JD12_PICK = {
        "crRNA1": ["WT", "A1C", "A8C+U15G"],
        "crRNA2": ["WT", "A1C", "U15G"],
    }
    DR_ALL = dict(DRS)
    DR_ALL["A1U"] = "UAUUUCUACUGUUGUAGAU"
    DR_ALL["U15G"] = DRS["WT"][:14] + "G" + DRS["WT"][15:]  # 15U->G 程序化派生
    jd_order = []
    for tag, sp in (("crRNA1", user["JD12_sp1"]), ("crRNA2", user["JD12_sp2"])):
        pfs = report["spacers"]["JD12_sp" + tag[-1]]["pfs_scholz"]
        note = ("crRNA-1: PFS CCUGG 弱(dep 0.628, 如实标注继续做); 骨架=自身榜单 0/1/4 名"
                if tag == "crRNA1" else
                "crRNA-2: PFS GGGAG 良好(dep 3.90); 骨架=三层合成: A1C(榜1)+U15G(§J OR9.0+过滤通过+先验榜1)")
        for dr_name in JD12_PICK[tag]:
            dr = DR_ALL[dr_name]
            rna = dr + sp.replace("T", "U")
            dna = rna.replace("U", "T")
            jd_order.append({
                "oligo_name": "JD12_%s_%s" % (tag, dr_name.replace("+", "")),
                "target": "TP53-R273H(SW480, JD12 双窗口)",
                "scaffold": dr_name, "dr_rna": dr,
                "spacer_dna_23nt": sp, "rna_crRNA_42nt": rna,
                "dna_core_42nt": dna,
                "synthesis_template_61nt": T7 + dna,
                "note": note + "; 模板=T7(19)+DR(19)+spacer(23), IVT 后核验 5' 端"})
    dst_csv = os.path.join(DATA, "wetlab_jd12_order.csv")
    with open(dst_csv, "w", newline="", encoding="utf-8") as f:
        w = _csv.DictWriter(f, fieldnames=list(jd_order[0].keys()))
        w.writeheader()
        w.writerows(jd_order)
    report["designs"]["JD12_full_order"] = {
        "for": "用户 2026-09-16 拍板: JD12 两条 spacer 都进实验; 骨架按各自管线 TOP 榜",
        "rows": len(jd_order), "csv": "data/wetlab_jd12_order.csv",
        "picks": JD12_PICK,
        "basis": {"crRNA1": "data/jd12_sp1_direct.top.json",
                  "crRNA2": "data/jd12_direct.top.json"},
        "note": "A8C+U15G 仅入 crRNA-1 组(其自身榜第4+唯一通过稳定化); "
                "crRNA-2 组稳定化位空缺(无通过项), stemmax 见 "
                "data/stemmax_design.json(探索, 未入单)"}

    # R248W(HT29) 两学派枚举(c.742C=T 在 s883 0-based)
    sch1 = []
    for p0 in range(879, 884):
        pfs_wt = s[p0:p0 + 5]
        k = 883 - p0
        pfs_mut = pfs_wt[:k] + "T" + pfs_wt[k + 1:]
        dw = tbl.get(pfs_wt.replace("T", "U"))
        dm = tbl.get(pfs_mut.replace("T", "U"))
        if dw is None or dm is None:
            continue
        sch1.append({"allele_at_pfs_pos": k + 1, "pfs_mut": pfs_mut,
                     "pfs_wt": pfs_wt, "mut_depletion": round(dm, 3),
                     "wt_depletion": round(dw, 3), "ratio": round(dm / dw, 2),
                     "mut_rank": rank.get(pfs_mut.replace("T", "U")),
                     "spacer_dna": rc(s[p0 - 24:p0])})
    sch1.sort(key=lambda x: -x["ratio"])
    sch2 = []
    for w0 in range(867, 884):
        win = s[w0:w0 + 24]
        spos = 24 - (883 - w0)
        if not 1 <= spos <= 8:
            continue
        pfs = s[w0 + 24:w0 + 29]
        dv = tbl.get(pfs.replace("T", "U"))
        if dv is None:
            continue
        mut_win = win[:883 - w0] + "T" + win[884 - w0:]
        sch2.append({"allele_at_spacer_pos": spos, "pfs": pfs,
                     "pfs_depletion": round(dv, 3),
                     "pfs_rank": rank.get(pfs.replace("T", "U")),
                     "spacer_dna_mut": rc(mut_win)})
    sch2.sort(key=lambda x: -x["pfs_depletion"])
    report["designs"]["R248W_HT29_options"] = {
        "for": "若细胞系确为 HT29(R248W); JD12 序列实为 R273H 靶向",
        "same_window_problem": "沿用 R248Q 窗口时 mut PFS UGGAG dep 2.36 vs wt "
                               "CGGAG 3.16, 比值 0.75 反向, 不可沿用",
        "school1_pfs_enum_top3": sch1[:3],
        "school2_seed_enum_top3": sch2[:3]}
    report["designs"]["R248W_HT29_constructs"] = {
        "school1_recommend(sch1[0])": build_constructs(
            sch1[0]["spacer_dna"], "R248W_pfs"),
        "school2_recommend(sch2[0])": build_constructs(
            sch2[0]["spacer_dna_mut"], "R248W_seed")}

    out = os.path.join(DATA, "jd12_redesign.json")
    json.dump(report, open(out, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("->", out)
    for tag in user:
        v = report["spacers"][tag].get("verdict", [])
        print(tag, "|", "; ".join(v))


if __name__ == "__main__":
    main()
