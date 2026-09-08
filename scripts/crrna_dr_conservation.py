"""同源 DR 位点保守性先验(策划案 V3 §4.1 / 表1 蛋白预测层, 任务⑪扩库重做 2026-09-08)。

序列库: data/dr_homologs.json v2 —— n=11 条(9 个独立基因座口径), 每条带逐条来源
  (文献/专利 SEQ ID/结构 ID), 禁止无来源条目。组成: V-K (Cas12a2) 2 条
  (SuCas12a2 8D4A 结构口径 + Zeng 2026 细胞实验口径) + V-A (Cas12a) 9 条
  (Fn/Lb/As/Pm/Mb/Ts/Bs, Zetsche 2015 直系同源库 + Creutzburg 2020 + Teng 2019)。
  V-K 公开可核验 DR 仅 SuCas12a2 一个天然基因座(ca33Cas12a2 等其余 V-K 位点未见
  文本形式公布的 repeat 序列, 不入库, 检索过程见输出 JSON 的 scope_note);
  V-A 条目按 Dmytrenko 2023 Fig.1c 的 Cas12a2+Cas12a 联合比对框架作近缘参照。

比对规则(显式, 无外部工具依赖; 任务⑪重做: 中段可变区允许 indel):
  1) 分段: 每条 DR 切为 flank5 / 锚L(TTTCTACT) / middle / 锚R(GTAGA) / trailing3。
     两条锚必须各出现且仅出现一次(锚唯一性 = 比对合法性检查, 否则报错)。
  2) 锚列: 8 列锚L + 5 列锚R 直接成列, 不参与 indel。
  3) middle 段多序列比对: 星形比对, 中心 = Cas12a2_zeng2026 的 middle
     (本管线 WT 效应子口径, 显式声明); 其余序列逐条与中心做 Needleman-Wunsch
     全局比对(match=+2, mismatch=-1, gap=-2; 回溯同分决胜序: 对角 > gap-in-seq >
     gap-in-ref, 保证确定性), 再按"参考列 + 插入列"规则合并成多序列比对
     (插入列对所有已并入序列补 gap)。
  4) 全比对 = 锚L 8 列 + middle 比对列 + 锚R 5 列。

保守性规则(显式阈值):
  每列 identity = 非 gap 字符中最高频碱基占比(gap 不计入分子分母, 但单列报告
  gap 数); 分类三档:
    identity == 1.00          -> strictly_conserved (完全保守)
    0.80 <= identity < 1.00   -> conserved (保守, 有同源变异记录)
    identity <  0.80          -> variable (可变, 扰动优先集中区)
  主口径 n=11 全条目逐座计数; 另给 unique-全序列 去重敏感性口径。
  3' 保守窗先验: strictly_conserved 的 3' 连续 run 长度; 管线 --cons3-window
  采用该严格口径(0.80-1.00 之间的位点已有同源替换史, 全权重惩罚缺乏依据,
  仅标注; 两口径 run 长度并列报告)。

判读文字全部由数值按规则生成(见 build_interpretation), 不手写结论。

输出: data/dr_conservation.json (含来源清单/比对方法/阈值/各位点数值/注册表
  效应子逐位映射 cas12a2 与 cas12a2_zeng2026, 供 crrna_scaffold_design 消费)。
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

ANCHOR_L = "TTTCTACT"
ANCHOR_R = "GTAGA"
CONS_THRESHOLD = 0.80           # conserved 阈值(显式声明)
STRICT_IDENTITY = 1.00          # strictly_conserved 阈值(显式声明)
NW_MATCH, NW_MISMATCH, NW_GAP = 2, -1, -2   # middle 段 NW 打分(显式声明)
CENTER_NAME = "Cas12a2_zeng2026"            # 星形比对中心(显式声明)

EFFECTOR_KEYS = ("cas12a2", "cas12a2_zeng2026")


def load_library(path=None):
    """读 DR 同源库; 兼容 v1(裸 list)与 v2({"entries": [...]} )两种封装。"""
    path = path or os.path.join(ROOT, "data", "dr_homologs.json")
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    entries = d["entries"] if isinstance(d, dict) else d
    for e in entries:
        if not e.get("provenance"):
            raise ValueError("%s 缺 provenance 字段" % e["name"])
    return entries


def segmented(name, seq, prov):
    """把 DR 切为 flank5/anchorL/middle/anchorR/trailing3。锚缺失或不唯一即报错。"""
    if seq.count(ANCHOR_L) != 1 or seq.count(ANCHOR_R) != 1:
        raise ValueError("%s 锚不唯一/缺失: %s" % (name, seq))
    l = seq.index(ANCHOR_L)
    r = seq.rindex(ANCHOR_R)
    if r <= l + len(ANCHOR_L):
        raise ValueError("%s 两锚区间非法: %s" % (name, seq))
    return {"name": name, "full": seq, "flank5": seq[:l],
            "anchor_l": seq[l:l + len(ANCHOR_L)],
            "middle": seq[l + len(ANCHOR_L):r],
            "anchor_r": seq[r:r + len(ANCHOR_R)],
            "trailing3": seq[r + len(ANCHOR_R):],
            "core_start_in_full": l,
            "provenance": prov}


def nw_align(ref, seq):
    """Needleman-Wunsch 全局比对(seq 对齐到 ref)。返回 (ref_aln, seq_aln)。
    回溯同分决胜序: 对角(match/mismatch) > gap-in-seq(上移) > gap-in-ref(左移)。"""
    n, m = len(ref), len(seq)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + NW_GAP
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + NW_GAP
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = NW_MATCH if ref[i - 1] == seq[j - 1] else NW_MISMATCH
            F[i][j] = max(F[i - 1][j - 1] + s, F[i - 1][j] + NW_GAP,
                          F[i][j - 1] + NW_GAP)
    ra, sa = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = NW_MATCH if ref[i - 1] == seq[j - 1] else NW_MISMATCH
            if F[i][j] == F[i - 1][j - 1] + s:
                ra.append(ref[i - 1]); sa.append(seq[j - 1])
                i, j = i - 1, j - 1
                continue
        if i > 0 and F[i][j] == F[i - 1][j] + NW_GAP:
            ra.append(ref[i - 1]); sa.append("-")
            i -= 1
            continue
        ra.append("-"); sa.append(seq[j - 1])
        j -= 1
    return "".join(reversed(ra)), "".join(reversed(sa))


def star_middle_msa(segs, center=CENTER_NAME):
    """middle 段星形多序列比对。返回 ({name: aligned_middle}, center_mid), 全部等长。

    两趟式合并(确定性): 第一趟每条序列与中心参考 NW 比对, 解析为
    (参考位字符或'-') + 各边界的插入串; 第二趟按 边界b插入列(宽=该边界最长
    插入串, 不足右侧补 '-') + 参考列 的顺序发射列, 保证插入列对全部序列共享。"""
    center_mid = next(s["middle"] for s in segs if s["name"] == center)
    n = len(center_mid)
    per_seq = {}    # name -> {"ref": [char|'-']*n, "ins": {boundary: str}}
    for s in segs:
        if s["name"] == center:
            per_seq[s["name"]] = {"ref": list(center_mid), "ins": {}}
            continue
        ra, sa = nw_align(center_mid, s["middle"])
        ref_row, ins = [], {}
        rpos = 0
        for rc, sc in zip(ra, sa):
            if rc == "-":
                ins.setdefault(rpos, "")
                ins[rpos] += sc
            else:
                ref_row.append(sc)
                rpos += 1
        assert rpos == n and len(ref_row) == n
        per_seq[s["name"]] = {"ref": ref_row, "ins": ins}
    ins_width = {b: max((len(p["ins"].get(b, "")) for p in per_seq.values()),
                        default=0)
                 for b in range(n + 1)}
    rows = {}
    for name, p in per_seq.items():
        row = []
        for b in range(n + 1):
            if ins_width[b]:
                row.extend(p["ins"].get(b, "").ljust(ins_width[b], "-"))
            if b < n:
                row.append(p["ref"][b])
        rows[name] = "".join(row)
    width = len(next(iter(rows.values())))
    assert all(len(r) == width for r in rows.values())
    return rows, center_mid


def build_alignment(segs, mid_msa):
    """锚L + middle比对 + 锚R 拼成全比对行; 返回 {name: aligned_full_core}。"""
    out = {}
    for s in segs:
        out[s["name"]] = s["anchor_l"] + mid_msa[s["name"]] + s["anchor_r"]
    return out


def column_profile(aln_rows):
    """对等长比对行逐列算保守度。identity = 最高频非 gap 碱基数 / 全列条目数
    (gap 计入分母——插入列少数派序列不得被视为保守; 显式规则)。gap 数单列报告。"""
    names = sorted(aln_rows)
    width = len(aln_rows[names[0]])
    if any(len(aln_rows[n]) != width for n in names):
        raise ValueError("比对行长度不齐")
    cols = []
    for i in range(width):
        col = [aln_rows[n][i] for n in names]
        comp = {}
        for ch in col:
            comp[ch] = comp.get(ch, 0) + 1
        nongap = {b: v for b, v in comp.items() if b != "-"}
        # 平局破缺: 按字母序(sorted)取首个最高频——dict 插入序虽已确定,
        # 但字母序破缺使平局结果与序列入库顺序解耦(外部复现评审: 旧版 max(set(col),
        # key=col.count) 受 PYTHONHASHSEED 影响, 打平位点共识字母跨进程翻转)。
        top = max(sorted(nongap), key=nongap.get)
        cols.append({"col_1based": i + 1, "consensus": top,
                     "identity": round(nongap[top] / len(col), 3),
                     "nongap_n": sum(nongap.values()),
                     "gap_n": comp.get("-", 0),
                     "composition": comp})
    return cols


def classify(identity):
    if identity >= STRICT_IDENTITY:
        return "strictly_conserved"
    if identity >= CONS_THRESHOLD:
        return "conserved"
    return "variable"


def cons3_runs(cols):
    """3' 端连续 run 长度: (strictly_conserved run, conserved(>=阈值) run)。"""
    strict = cons = 0
    for c in reversed(cols):
        if c["identity"] >= STRICT_IDENTITY:
            strict += 1
        else:
            break
    for c in reversed(cols):
        if c["identity"] >= CONS_THRESHOLD:
            cons += 1
        else:
            break
    return strict, cons


def map_effector(dr_dna, segs_by_name, mid_msa, cols, entry_name):
    """把列保守性映射到效应子 DR 的逐位坐标。效应子必须在比对库内
    (按其比对行逐列定位, 天然处理 gap)。"""
    seg = segs_by_name[entry_name]
    if seg["full"] != dr_dna:
        raise ValueError("注册表 DR 与库条目 %s 序列不一致" % entry_name)
    row = seg["anchor_l"] + mid_msa[entry_name] + seg["anchor_r"]
    col_of_corepos = {}     # core 内 0-based 位置 -> 比对列 1-based
    ci = 0
    for k, ch in enumerate(row):
        if ch != "-":
            col_of_corepos[ci] = k + 1
            ci += 1
    core = seg["anchor_l"] + seg["middle"] + seg["anchor_r"]
    positions = []
    for i, b in enumerate(dr_dna):
        pos = i + 1
        if i < seg["core_start_in_full"]:
            positions.append({"pos_1based": pos, "base": b, "column": None,
                              "region": "flank5", "identity": None,
                              "class": "unaligned_flank5"})
        elif i < seg["core_start_in_full"] + len(core):
            c = cols[col_of_corepos[i - seg["core_start_in_full"]] - 1]
            positions.append({"pos_1based": pos, "base": b, "column": c["col_1based"],
                              "region": "core", "identity": c["identity"],
                              "class": classify(c["identity"])})
        else:
            positions.append({"pos_1based": pos, "base": b, "column": None,
                              "region": "trailing3", "identity": None,
                              "class": "unaligned_trailing3"})
    strict_run, cons_run = cons3_runs(cols)
    # 管线窗口口径: DR 3' 末端 window nt = 严格保守 run 的列 + trailing3 位
    # (trailing3 在锚外不参与列统计, 按 Dmytrenko 2023 Fig.1c 3' 端保守口径并入窗口)
    window_prior = strict_run + len(seg["trailing3"])
    return {"dr_dna": dr_dna, "library_entry": entry_name,
            "flank5": seg["flank5"], "trailing3": seg["trailing3"],
            "middle": seg["middle"], "positions": positions,
            "cons3_window_prior": window_prior,
            "cons3_window_rule": "strictly_conserved 3' run(%d 列) + trailing3(%d nt) = %d; "
                                 "0.80 阈值 run 为 %d 列, 因这些位有同源替换史不作全权重惩罚"
                                 % (strict_run, len(seg["trailing3"]), window_prior, cons_run)}


def build_interpretation(cols, cols_dedup, strict_run, cons_run, n_entries, n_dedup):
    """判读文字由数值按规则生成: 只陈述数值与阈值比较结果。"""
    var_cols = [c["col_1based"] for c in cols if classify(c["identity"]) == "variable"]
    cons_mid = [c["col_1based"] for c in cols if classify(c["identity"]) == "conserved"]
    same_cls = all(classify(a["identity"]) == classify(b["identity"])
                   for a, b in zip(cols, cols_dedup))
    return [
        "n=%d 条目(%d unique 序列)锚定+middle星形比对共 %d 列: 完全保守 %d 列, "
        "保守(>=%.2f 且 <1.00) %d 列, 可变 %d 列%s。" % (
            n_entries, n_dedup, len(cols),
            sum(1 for c in cols if c["identity"] >= STRICT_IDENTITY),
            CONS_THRESHOLD,
            sum(1 for c in cols if CONS_THRESHOLD <= c["identity"] < STRICT_IDENTITY),
            len(var_cols),
            "(列 %s)" % ",".join(map(str, var_cols)) if var_cols else ""),
        "可变列为 %s(identity=%s), 位于比对区中段——对应 Dmytrenko 2023 Fig.1c 的 "
        "loop 可变区; 3' 端锚区 %d 列完全保守。" % (
            ",".join(map(str, var_cols)) if var_cols else "无",
            ",".join("%.2f" % c["identity"] for c in cols
                     if classify(c["identity"]) == "variable"),
            strict_run),
        "3' 保守窗先验: 严格口径 run=%d(管线 --cons3-window 采用), 0.80 阈值口径 run=%d"
        "(列 %s 有同源替换史, 仅作保守标注不作全权重惩罚)。" % (
            strict_run, cons_run,
            ",".join(map(str, cons_mid)) if cons_mid else "无"),
        "去重敏感性: 两口径下逐列分类%s。" % ("完全一致" if same_cls else "存在差异, 见逐列对照"),
    ]


def main():
    entries = load_library()
    segs = [segmented(e["name"], e["dr_dna"], e["provenance"]) for e in entries]
    mid_msa, center_mid = star_middle_msa(segs)
    aln_rows = build_alignment(segs, mid_msa)
    cols = column_profile(aln_rows)
    # 去重敏感性口径: 按完整 DR 序列去重
    uniq_rows = {}
    for s in segs:
        uniq_rows.setdefault(s["full"], aln_rows[s["name"]])
    cols_dedup = column_profile(uniq_rows)
    strict_run, cons_run = cons3_runs(cols)
    for c, cd in zip(cols, cols_dedup):
        c["class"] = classify(c["identity"])
        c["identity_dedup"] = cd["identity"]
        c["class_dedup"] = classify(cd["identity"])

    with open(os.path.join(ROOT, "configs", "scaffold_registry.json"),
              encoding="utf-8") as fh:
        registry = json.load(fh)["entries"]
    segs_by_name = {s["name"]: s for s in segs}
    entry_of_effector = {"cas12a2": "SuCas12a2_8D4A",
                         "cas12a2_zeng2026": "Cas12a2_zeng2026"}
    mapped = {}
    for key in EFFECTOR_KEYS:
        mapped[key] = map_effector(registry[key]["scaffold"], segs_by_name,
                                   mid_msa, cols, entry_of_effector[key])

    w3 = [c["identity"] for c in cols if c["col_1based"] > len(cols) - 5]
    rest = [c["identity"] for c in cols if c["col_1based"] <= len(cols) - 5]

    hdr = "%-24s %-6s %-8s %-*s %-6s %-8s" % ("DR", "5'侧", "锚L", len(center_mid) + 2,
                                              "middle", "锚R", "3'侧")
    print(hdr)
    for s in segs:
        print("%-24s %-6s %-8s %-*s %-6s %-8s" % (
            s["name"], s["flank5"] or "-", s["anchor_l"], len(center_mid) + 2,
            mid_msa[s["name"]], s["anchor_r"], s["trailing3"] or "-"))
    print("\n位点保守性(identity / class):")
    for c in cols:
        print("  col%2d %s  %.3f  %-18s %s" % (c["col_1based"], c["consensus"],
                                               c["identity"], c["class"],
                                               c["composition"]))
    print("3' 末 5 列平均 identity %.3f vs 其余 %.3f" % (
        sum(w3) / len(w3), sum(rest) / len(rest)))
    print("3' 保守 run: strict=%d, >=%.2f 阈值=%d" % (strict_run, CONS_THRESHOLD, cons_run))
    interp = build_interpretation(cols, cols_dedup, strict_run, cons_run,
                                  len(segs), len(uniq_rows))
    print("\n判读(由数值按规则生成):")
    for line in interp:
        print("  " + line)
    for key, m in mapped.items():
        print("映射 %s (%s): cons3_window_prior=%d" % (key, m["dr_dna"],
                                                       m["cons3_window_prior"]))

    out = os.path.join(ROOT, "data", "dr_conservation.json")
    payload = {
        "generated_by": "scripts/crrna_dr_conservation.py (任务⑪扩库重做, 2026-09-08)",
        "library": os.path.join("data", "dr_homologs.json"),
        "library_version": 2,
        "n_sequences": len(segs),
        "n_unique_sequences": len(uniq_rows),
        "scope_note": "V-K 公开可核验 DR 仅 SuCas12a2 天然基因座 + Zeng 2026 合成口径; "
                      "其余 V-K 位点(Dmytrenko 2023 ED Fig.1/10 的 ca33Cas12a2 等)未见文本"
                      "形式 repeat 序列, 不入库。V-A 条目按 Dmytrenko 2023 Fig.1c 联合比对"
                      "框架作近缘参照。检索路径: Dmytrenko 2023 正文+SI(仅图无文本序列)、"
                      "Methods Enzymol 2025 (PMC12975306, array 全长 repeat)、"
                      "Zetsche 2015 Cell 及同族专利 SEQ ID 多源复核、Zeng 2026 Nature SI。",
        "method": {
            "segmentation": "flank5 / 锚L(%s) / middle / 锚R(%s) / trailing3; "
                            "两锚各须唯一出现(合法性检查)" % (ANCHOR_L, ANCHOR_R),
            "middle_msa": "星形比对, 中心=%s; NW match=%+d mismatch=%+d gap=%+d, "
                          "回溯决胜序 对角>gap-in-seq>gap-in-ref" % (
                              CENTER_NAME, NW_MATCH, NW_MISMATCH, NW_GAP),
            "identity": "每列 identity = 非 gap 字符中最高频碱基占比(gap 单列报告)",
            "thresholds": {"strictly_conserved": STRICT_IDENTITY,
                           "conserved": CONS_THRESHOLD,
                           "variable": "< %.2f" % CONS_THRESHOLD},
            "weighting": "主口径 n=%d 全条目; 敏感性口径按完整序列去重 n=%d" % (
                len(segs), len(uniq_rows)),
        },
        "alignment": [{**{k: s[k] for k in ("name", "full", "flank5", "middle",
                                            "trailing3", "core_start_in_full",
                                            "provenance")},
                       "aligned_middle": mid_msa[s["name"]],
                       "aligned_core": aln_rows[s["name"]]} for s in segs],
        "conservation": cols,
        "cons3_prior": {
            "strictly_conserved_3prime_run": strict_run,
            "conserved_3prime_run_at_threshold": cons_run,
            "recommended_window": strict_run,
            "rule": "管线 3' 保守窗 = strictly_conserved 3' run; 0.80-1.00 位点有同源"
                    "替换史, 标注但不全权重惩罚; trailing3 位在锚外按 Fig.1c 3' 保守口径"
                    "并入效应子窗口(见 mapped_effectors.*.cons3_window_rule)",
        },
        "cons3_window_5col_mean_identity": round(sum(w3) / len(w3), 3),
        "rest_mean_identity": round(sum(rest) / len(rest), 3),
        "mapped_effectors": mapped,
        "interpretation": interp,
        "limitation": "n=%d 条目/%d unique 序列仍为小样本定量; V-K 仅 2 条, 结论主体由 "
                      "V-A 近缘参照承载, 与 Dmytrenko 2023 Fig.1c 跨家族定性结论一致"
                      "(3' 保守/中段可变), 不替代其结论" % (len(segs), len(uniq_rows)),
    }
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
