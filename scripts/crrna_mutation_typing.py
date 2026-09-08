"""突变类型判定(最长 ORF 法)(策划案 V3 §4.4 模块一靶标解析 / 表2 APC 入选理由
对齐整改, 2026-09-07)。

输入突变与野生转录本序列, 分别在两者中定位最长 ORF, 通过比较突变对 WT 最长
ORF(参考 ORF)的影响判定突变类型: 错义/无义/移码(边界类别——同义/框内 indel/
终止子丢失/非编码——如实单列, 不塞进三主类)。判定规则、阈值、边界处理全部显式
写在本模块(全仓唯一定义点), 判读文字由数值按显式规则生成, 不手写结论方向。

ORF 规则(显式声明):
  1. 标准密码子表; 起始密码子 ATG(AUG), 终止子 TAA/TAG/TGA(UAA/UAG/UGA);
     只扫描正义链 3 个阅读框——输入为 mRNA 转录本, 不做反向互补链;
  2. ORF = 某 ATG 起, 到同阅读框第一个终止子(含终止子)为止; 若至序列末端仍无
     框内终止子, ORF 延伸至最后一个完整密码子, 末端不足 3 nt 的部分密码子
     丢弃, 记 terminated=False;
  3. 最长 ORF 选择(显式 tie-break): (a) length_nt 大者优先; (b) 并列时起始位置
     靠前者(0-based 小)优先; (c) 再并列时阅读框序号小者优先。序列中每条 ATG
     独立评估(同框靠后的 ATG 产生更短 ORF, 自然落选);
  4. 突变事件模型(与 crrna_design_agent.find_mutation 同口径): 仅支持单事件——
     等长=碱基替换(多处差异时取落在参考 ORF 内的第一处判定, 全部差异如实列于
     输出); 不等长=共同前后缀夹逼出的单段纯插入/缺失; 替换与 indel 混合等复杂
     事件拒绝判定并如实报错, 不猜;
  5. 类型判定(突变对参考 ORF = WT 最长 ORF 的影响; MUT 侧报告同一 ATG 起点的
     对应 ORF(mut_orf_same_start, 坐标系 = MUT 转录本)与独立找到的最长 ORF
     (mut_longest_orf)双份证据):
     - 差异全部落在参考 ORF 之外 -> non_coding(UTR/非编码, 如实报告);
     - 单碱基替换在 ORF 内: 密码子翻译不变 -> synonymous(同义);
       sense->stop -> nonsense(无义, MUT 对应 ORF 在该提前终止子处结束,
       ORF 变短); sense->sense 单氨基酸替换且 ORF 起止/长度不变 -> missense
       (错义); stop->sense -> stop_lost(如实单列);
     - 替换摧毁起始密码子本身 -> start_lost(如实单列);
     - 单段插入/缺失在 ORF 内: |Δlen|%3!=0 -> frameshift(移码, 阅读框自突变
       点起位移, MUT 对应 ORF 终止位置改变且非单点替换); |Δlen|%3==0 ->
       inframe_indel(框内插入缺失, 如实单列, 非三主类)。

APC 序列核查(2026-09-07, 本整改附带发现, 如实记录): 旧 data/agent/apc_q1312x
序列对的 C>T 编辑(转录本 3946 位)由 crrna_panel_targets.apc_design 以
seq.find("ATG") 锚定——NM_000038.6 在 5'UTR 存在上游 ATG(13-15 位), 真正
CDS 起始于 60 位(最长 ORF 0-based 59..8591, 2843 aa), 47 nt 偏移使旧编辑落在
真实阅读框密码子 1296(GCA>GTA, A1296V 错义), 并非设计的无义。2026-09-08
用户拍板正式采纳核查结论: APC 靶点整体切换为真实无义 Q1328*(按同一设计
规则 MCR 1300-1450 首个 CAG->TAG, 最长 ORF 锚点, 密码子 1328, c.3982C>T,
转录本 4041 位), panel 键名更名 apc_q1328x, 序列对与下游(湿实验 CSV /
演示 / Chai 增量 / 选型特征)全部重生成; 旧键名 apc_q1312x 仅存在于 git
历史, 本模块对当前序列对的分型输出为无义 Q1328*(提前终止于密码子 1328)。

用法(库): import crrna_mutation_typing as mt; mt.classify(wt_seq, mut_seq)
自检: python scripts/crrna_mutation_typing.py <wt.fa> <mut.fa>
"""
import json
import sys

# 标准密码子表(DNA 口径; 与 scripts/crrna_panel_targets.py 同一构造)
CODON_TABLE = {}
_b = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON_TABLE[_x + _y + _z] = _aas[_i]
            _i += 1

START = "ATG"
STOPS = ("TAA", "TAG", "TGA")

TYPE_LABELS = {
    "missense": "错义",
    "nonsense": "无义",
    "frameshift": "移码",
    "synonymous": "同义",
    "inframe_indel": "框内插入缺失",
    "stop_lost": "终止子丢失",
    "start_lost": "起始密码子丢失",
    "non_coding": "非编码区(参考ORF外)",
    "complex": "复杂多事件(不支持)",
}


def norm(seq):
    """归一为 DNA 大写串(接受 RNA/混合输入; 去空白)。"""
    return "".join(str(seq).split()).upper().replace("U", "T")


def find_orfs(seq):
    """枚举正义链全部 ATG 起始的 ORF(规则见模块 docstring 第 1-2 条)。

    返回 list[dict], 坐标一律 0-based: start0 = ATG 首位; stop0 = 框内第一个
    终止子首位(无则 None); end0 = 排他末端(含终止子或最后完整密码子);
    length_nt = end0 - start0。"""
    seq = norm(seq)
    n = len(seq)
    orfs = []
    for frame in range(3):
        i = frame
        while i + 3 <= n:
            if seq[i:i + 3] == START:
                j = i
                stop0 = None
                while j + 3 <= n:
                    if seq[j:j + 3] in STOPS:
                        stop0 = j
                        break
                    j += 3
                if stop0 is not None:
                    end0 = stop0 + 3
                else:
                    end0 = i + (n - i) // 3 * 3   # 最后完整密码子, 部分密码子丢弃
                orfs.append({"start0": i, "frame": frame, "stop0": stop0,
                             "end0": end0, "length_nt": end0 - i,
                             "terminated": stop0 is not None})
            i += 3
    return orfs


def longest_orf(seq):
    """最长 ORF(tie-break: 长度降序 -> 起点升序 -> 阅读框序号升序); 无 ATG 返回 None。"""
    orfs = find_orfs(seq)
    if not orfs:
        return None
    return min(orfs, key=lambda o: (-o["length_nt"], o["start0"], o["frame"]))


def orf_from(seq, start0):
    """从指定起点(0-based)取 ORF; 该处非 ATG 时返回 None(如实)。"""
    seq = norm(seq)
    if start0 < 0 or seq[start0:start0 + 3] != START:
        return None
    for o in find_orfs(seq):
        if o["start0"] == start0:
            return o
    return None


def translate(seq, start0, end0):
    """翻译 [start0, end0) 区段(end0 对齐密码子边界); 终止子译为 '*'。"""
    seq = norm(seq)
    return "".join(CODON_TABLE.get(seq[i:i + 3], "?")
                   for i in range(start0, end0, 3))


def _orf_public(o):
    """ORF 的 JSON 口径(1-based; protein_aa 不含终止子)。"""
    if o is None:
        return None
    return {"start_1based": o["start0"] + 1, "frame": o["frame"],
            "stop_1based": o["stop0"] + 1 if o["stop0"] is not None else None,
            "end_1based": o["end0"], "length_nt": o["length_nt"],
            "protein_aa": o["length_nt"] // 3 - (1 if o["terminated"] else 0),
            "terminated": o["terminated"]}


def locate_event(wt, mut):
    """单事件突变定位(与 crrna_design_agent.find_mutation 同口径, 扩到单段 indel)。

    等长: kind='substitution', diffs0 = 全部差异位(0-based, 升序);
    不等长: 共同前后缀夹逼, 恰为纯插入/缺失(一侧中段为空)时返回
    kind='insertion'/'deletion' + pos0(首个分歧位)与两段序列; 其余(混合复杂
    事件)抛 ValueError, 如实拒绝, 不猜。"""
    wt, mut = norm(wt), norm(mut)
    if len(wt) == len(mut):
        diffs = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
        if not diffs:
            raise ValueError("WT 与 MUT 序列完全一致, 无突变可判定")
        return {"kind": "substitution", "diffs0": diffs, "delta": 0,
                "wt_seg": "".join(wt[i] for i in diffs),
                "mut_seg": "".join(mut[i] for i in diffs)}
    p = 0
    while p < min(len(wt), len(mut)) and wt[p] == mut[p]:
        p += 1
    s = 0
    while s < min(len(wt), len(mut)) - p and wt[len(wt) - 1 - s] == mut[len(mut) - 1 - s]:
        s += 1
    wt_mid = wt[p:len(wt) - s]
    mut_mid = mut[p:len(mut) - s]
    if wt_mid and mut_mid:
        raise ValueError(
            "复杂突变(替换+indel 混合或多段): WT 中段 %s... vs MUT 中段 %s..., "
            "本模块只支持单事件, 拒绝判定" % (wt_mid[:12], mut_mid[:12]))
    kind = "insertion" if not wt_mid else "deletion"
    return {"kind": kind, "pos0": p, "delta": len(mut) - len(wt),
            "wt_seg": wt_mid, "mut_seg": mut_mid}


def classify(wt, mut):
    """最长 ORF 法突变类型判定主入口(规则见模块 docstring 第 5 条)。

    返回 dict: mutation_type(+_label) / event / 双序列 ORF 证据 / 提前终止子
    (nonsense)或移码位移(frameshift)数值 / reading(由数值按显式规则生成)。"""
    wt, mut = norm(wt), norm(mut)
    ref = longest_orf(wt)
    if ref is None:
        raise ValueError("WT 转录本中未找到任何 ATG 起始的 ORF, 无法判定")
    ev = locate_event(wt, mut)
    mlong = longest_orf(mut)
    out = {"rule": ("最长 ORF 法: 参考 ORF=WT 最长 ORF(正义链 3 框, ATG 起, 首个框内"
                    "终止子止, tie-break 长度/起点/框序); 类型由突变对参考 ORF 的"
                    "影响按显式规则判定(scripts/crrna_mutation_typing.py docstring)"),
           "event": {k: v for k, v in ev.items()},
           "wt_orf": _orf_public(ref),
           "mut_longest_orf": _orf_public(mlong),
           "premature_stop": None, "aa_change": None,
           "codon_no": None, "wt_codon": None, "mut_codon": None,
           "wt_aa": None, "mut_aa": None}

    def _finish(mtype, mref, reading):
        out["mutation_type"] = mtype
        out["mutation_type_label"] = TYPE_LABELS[mtype]
        out["mut_orf_same_start"] = _orf_public(mref)
        wlen = out["wt_orf"]["protein_aa"]
        mlen = out["mut_orf_same_start"]["protein_aa"] if mref else None
        out["orf_length_change_aa"] = (mlen - wlen) if mlen is not None else None
        out["reading"] = reading
        return out

    if ev["kind"] == "substitution":
        diffs = ev["diffs0"]
        in_orf = [d for d in diffs if ref["start0"] <= d < ref["end0"]]
        out["event"]["diffs_1based"] = [d + 1 for d in diffs]
        out["event"]["diffs_in_ref_orf_1based"] = [d + 1 for d in in_orf]
        if not in_orf:
            return _finish("non_coding", orf_from(mut, ref["start0"]),
                           "全部 %d 处差异(%s)落在 WT 最长 ORF(%d..%d)之外, "
                           "不影响参考 ORF 翻译 -> 非编码区事件"
                           % (len(diffs), [d + 1 for d in diffs][:6],
                              ref["start0"] + 1, ref["end0"]))
        d = in_orf[0]
        if len(diffs) > 1:
            out["event"]["note"] = ("检出 %d 处差异, 取参考 ORF 内第一处(1-based %d)"
                                    "判定, 其余如实列出" % (len(diffs), d + 1))
        codon_no = (d - ref["start0"]) // 3 + 1
        ci = ref["start0"] + (codon_no - 1) * 3
        wtc, mutc = wt[ci:ci + 3], mut[ci:ci + 3]
        waa, maa = CODON_TABLE[wtc], CODON_TABLE[mutc]
        out.update({"codon_no": codon_no, "wt_codon": wtc, "mut_codon": mutc,
                    "wt_aa": waa, "mut_aa": maa,
                    "aa_change": "%s%d%s" % (waa, codon_no, maa)})
        mref = orf_from(mut, ref["start0"])
        if d < ref["start0"] + 3:
            return _finish("start_lost", mref,
                           "差异(1-based %d)落在参考 ORF 起始密码子(ATG, %d..%d)内: "
                           "%s>%s, MUT 该处 %s, 起始密码子%s -> 起始密码子丢失"
                           % (d + 1, ref["start0"] + 1, ref["start0"] + 3,
                              wt[d], mut[d], mutc,
                              "被摧毁" if mref is None else "仍存在(同义)"))
        if waa == maa:
            return _finish("synonymous", mref,
                           "参考 ORF %d..%d(%d aa)在 MUT 中起止与长度不变; 密码子 %d "
                           "%s>%s 翻译不变(%s) -> 同义"
                           % (ref["start0"] + 1, ref["end0"],
                              out["wt_orf"]["protein_aa"], codon_no, wtc, mutc, waa))
        if maa == "*":
            mlen = mref["length_nt"] // 3 - (1 if mref["terminated"] else 0)
            wlen = out["wt_orf"]["protein_aa"]
            out["premature_stop"] = {"codon_no": codon_no,
                                     "transcript_pos_1based": ci + 1,
                                     "wt_codon": wtc, "mut_codon": mutc}
            out["truncated_fraction"] = round(1.0 - mlen / wlen, 4)
            return _finish("nonsense", mref,
                           "密码子 %d %s>%s 变为终止子: MUT 对应 ORF 在该处提前终止"
                           "(蛋白 %d aa -> %d aa, 截短 %.1f%%; WT ORF %d..%d vs MUT "
                           "%d..%d) -> 无义"
                           % (codon_no, wtc, mutc, wlen, mlen,
                              (1.0 - mlen / wlen) * 100,
                              ref["start0"] + 1, ref["end0"],
                              mref["start0"] + 1, mref["end0"]))
        if waa == "*":
            return _finish("stop_lost", mref,
                           "密码子 %d %s>%s: 终止子变为有义密码子(%s), ORF 延伸 -> "
                           "终止子丢失" % (codon_no, wtc, mutc, maa))
        return _finish("missense", mref,
                       "参考 ORF %d..%d(%d aa)在 MUT 中起止与长度不变; 密码子 %d "
                       "%s>%s(%s>%s), 单氨基酸替换 %s -> 错义"
                       % (ref["start0"] + 1, ref["end0"], out["wt_orf"]["protein_aa"],
                          codon_no, wtc, mutc, waa, maa, out["aa_change"]))

    # 单段插入/缺失
    p, delta = ev["pos0"], ev["delta"]
    out["event"]["pos_1based"] = p + 1
    if not (ref["start0"] <= p < ref["end0"]):
        mstart = ref["start0"] + (delta if p <= ref["start0"] else 0)
        return _finish("non_coding", orf_from(mut, mstart),
                       "%s(%+d nt, 首个分歧位 1-based %d)落在 WT 最长 ORF(%d..%d)"
                       "之外, 不影响参考 ORF 翻译 -> 非编码区事件"
                       % ({"insertion": "插入", "deletion": "缺失"}[ev["kind"]],
                          delta, p + 1, ref["start0"] + 1, ref["end0"]))
    mref = orf_from(mut, ref["start0"])
    kind_zh = {"insertion": "插入", "deletion": "缺失"}[ev["kind"]]
    if delta % 3 != 0:
        mstop = mref["stop0"] + 1 if (mref and mref["stop0"] is not None) else None
        wstop = ref["stop0"] + 1 if ref["stop0"] is not None else None
        return _finish("frameshift", mref,
                       "%s %+d nt(|Δ|%%3!=0)于参考 ORF 内(首个分歧位 1-based %d): "
                       "阅读框自该处起位移, MUT 对应 ORF 终止位置改变(WT 终止子 "
                       "1-based %s -> MUT %s; 非单点替换) -> 移码"
                       % (kind_zh, delta, p + 1, wstop, mstop))
    aa_delta = delta // 3
    return _finish("inframe_indel", mref,
                   "%s %+d nt(|Δ|%%3==0)于参考 ORF 内(首个分歧位 1-based %d): "
                   "阅读框不变, 蛋白长度 %+d aa -> 框内插入缺失(非错义/无义/移码"
                   "三主类, 如实单列)" % (kind_zh, delta, p + 1, aa_delta))


def feature_entry(typing):
    """模块一靶标解析特征向量的突变类型维(紧凑口径; 完整 ORF 证据在
    design.json 的 mutation.typing 块)。"""
    t = typing
    ps = t.get("premature_stop")
    return {"mutation_type": t["mutation_type"],
            "mutation_type_label": t["mutation_type_label"],
            "aa_change": t.get("aa_change"),
            "premature_stop_codon": ps["codon_no"] if ps else None,
            "orf_wt_protein_aa": t["wt_orf"]["protein_aa"],
            "orf_mut_protein_aa": (t["mut_orf_same_start"] or {}).get("protein_aa"),
            "orf_length_change_aa": t.get("orf_length_change_aa")}


def _read_fasta_single(path):
    with open(path, encoding="utf-8") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def main():
    if len(sys.argv) != 3:
        raise SystemExit("用法: python crrna_mutation_typing.py <wt.fa> <mut.fa>")
    t = classify(_read_fasta_single(sys.argv[1]), _read_fasta_single(sys.argv[2]))
    print(json.dumps(t, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
