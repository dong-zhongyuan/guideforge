# -*- coding: utf-8 -*-
"""Chai-1 骨架-靶标组合界面特征差量提取(策划案 V3 §4.1 末段 / 表1「选型特征」
对齐整改, 2026-09-09, 任务⑫)。

策划案原文要求: 「由于本项目骨架分型依赖靶标特征, 结构层特别计算每个
骨架-靶标组合的界面特征差量, 作为选型分类器的关键输入」; 表1 蛋白预测层
产出「ipTM、界面PAE、接触残基互补度」作为「选型特征」。本模块把该要求落到
现有 Chai 矩阵 JSON 上:

数据源(优先序, 显式声明于产物 source_matrix):
  1. data/chai_matrix_v3aligned.json(V3 对齐 4 靶矩阵, 由
     scripts/crrna_chai_increment_merge.py 在增量折叠回收后生成, 尚不存在);
  2. data/chai_matrix_4t.json(现行 4 靶 x 8 骨架 = 32 组合矩阵, 当前数据源)。

字段口径(以矩阵 JSON 实际存在的字段为准):
  iptm_mean    aggregate ipTM 5 模型均值(整体复合物)
  prot_crRNA   蛋白-crRNA 链对 ipTM 均值(界面主口径, 见 README Chai 节:
               crRNA-靶 RNA 双链预测稳定性存疑, 仅作参考)
  crRNA_target crRNA-靶RNA 链对 ipTM 均值(参考口径; 本矩阵靶 RNA 为合成
               完全互补构造, 全矩阵 0.45-0.49, 与旧矩阵原生靶 0.03-0.2
               是靶构造口径差, 两表分开引用)
  clash_frac   5 模型中带链间冲突的比例(接触冲突指标)
  iptm_sd      aggregate ipTM 5 模型 sd(模型间噪声代理; prot-crRNA 链对的
               per-model sd 未存盘, 2sd null 判定以此为代理, 与
               crrna_chai_matrix_summary.py 同口径)
  表1 的「界面PAE」与「接触残基互补度」未随矩阵存盘(矩阵未保留 PAE 与
  逐残基接触指标)——本特征表以实际存盘字段为准, 接触相关维度由两条链对
  ipTM 与 clash_frac 承载, 该映射关系如实声明于此, 不补造。

差量定义(显式): 对非 WT 骨架 s 与靶标 t,
  delta_<f>(s,t) = <f>(s,t) - <f>(WT,t)  (同靶 WT 组合为参照)
WT 组合值收入 wt_reference 块; 差量 4 维 = DELTA_FEATURES, 为选型契约列
(crrna_train_selector.INTERFACE_DELTA_FEATURES 经 import 共用本定义)。

判读规则(全部判读文字由数值按显式规则生成, 阈值常量 import 自
crrna_chai_matrix_summary, 不各自复制字面量; 方向词由数值符号决定):
  逐组合(prot-crRNA 主口径): |Δ| <= NEUTRAL_ABS(0.03) -> 界面中性;
    |Δ| >= COLLAPSE_ABS(0.25) -> 界面塌陷/大幅上升(方向由符号定);
    其间 -> 界面中等下降/上升; |Δ| <= NULL_SD_K(2) x iptm_sd -> 噪声内(判 null)
  逐骨架跨靶: 复用 crrna_chai_matrix_summary.scaffold_reading(与
    data/chai_matrix_4t_summary.json 的 reading 同源同函数, 逐字一致)。

方法学局限(如实声明, 随产物落盘):
  1. 自模板复述局限: 矩阵全部分数来自 100% 一致自模板(8D4A链A)注入——
     模板=查询自身的 WT 蛋白结构, 高 ipTM 是模板合规性的平凡结果(sanity
     check), 界面差量只描述「自模板口径下模型对骨架扰动的数值表型」,
     不能作为界面可预测性证据; template-free / scrambled-template 对照
     (AF3/Protenix 独立引擎)补齐前, 不得引为分型假说的结构证据;
  2. WT prot-crRNA ipTM 四靶均值 ~0.49 低于可用界面预测区间(~0.5-0.6),
     即使撇开自模板问题, 绝对水平也不足以支撑界面结论;
  3. B_break_compensate 行为已退役茎破坏混杂臂的表型(2026-09-04 起湿实验
     产物已换去混杂新臂, 重折叠待 Chai 机时), 不得引为当前分子的界面结论;
  4. 接入定位: 界面差量作为结构描述特征(选型特征维)进入选型链路——
     文献训练行(Han/Fig1f/Teng)无结构数据, 文献决策树在这 4 维上不能
     分裂, 故当前经显式附加层消费(与 target_abundance_tier 先例同:
     crrna_design_agent.target_context 的 interface_deltas 块 + webapp
     panel/design 页同源展示), IVT 8x4 矩阵训练集 schema 将含这 4 列,
     实测活性回填后才进入选型分类器训练; 不捏造训练标签。

运行: python scripts/crrna_chai_interface_features.py
输出: data/chai_interface_deltas.json
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

from crrna_chai_matrix_summary import (  # noqa: E402
    COLLAPSE_ABS, NEUTRAL_ABS, NULL_SD_K, scaffold_reading)

MATRIX_CANDIDATES = ["chai_matrix_v3aligned.json", "chai_matrix_4t.json"]
OUT_JSON = os.path.join(DATA, "chai_interface_deltas.json")

# 矩阵行实际存盘字段; 前 4 个参与差量, iptm_sd 为噪声代理, n_models 溯源
RAW_FIELDS = ("iptm_mean", "iptm_sd", "prot_crRNA", "crRNA_target",
              "clash_frac", "n_models")
DELTA_FIELDS = ("iptm_mean", "prot_crRNA", "crRNA_target", "clash_frac")
DELTA_FEATURES = ["delta_" + f for f in DELTA_FIELDS]

# 智能体 panel 键 -> 矩阵靶标列名(显式映射, 不做大小写猜测:
# 矩阵 APC 列为 "APC_Q1328x", 其余全大写)。KRAS_G12C 在现行矩阵中但不在
# V3 表2 panel(干实验附加列); APC_Q1328x 在 panel 但缺席现行矩阵(增量
# 折叠回填后由 v3aligned 矩阵覆盖)。APC 靶点 2026-09-08 经阅读框核查由旧
# Q1312x(实为 A1296V 错义)切换为真实无义 Q1328*, panel 键与矩阵列同步更名。
PANEL_TO_MATRIX_TARGET = {"tp53_r248q": "TP53_R248Q",
                          "kras_g12d": "KRAS_G12D",
                          "tp53_r273h": "TP53_R273H",
                          "apc_q1328x": "APC_Q1328x"}

LIMITATION = (
    "界面差量为 Chai-1 100% 一致自模板(8D4A链A)口径的模板合规性数值: "
    "模板=查询自身的 WT 结构, 高分为自模板复述的平凡结果(sanity check); "
    "本表只作为结构描述特征(选型特征维)记录与展示, 不构成界面可预测性证据, "
    "template-free/scrambled-template 对照(AF3/Protenix 独立引擎)补齐前不得"
    "引为分型假说的结构证据。WT prot-crRNA ipTM≈0.49 亦低于可用界面预测"
    "区间(~0.5-0.6); B_break_compensate 行为已退役混杂臂表型(见方法学注释)。")

SELECTOR_INTEGRATION = (
    "接入定位: 结构描述特征(选型特征维)经显式附加层进入选型链路——文献训练行"
    "(Han/Fig1f/Teng)无结构数据, 文献决策树在这 4 维上不能分裂(与 "
    "target_abundance_tier 先例同), 不捏造训练标签; 当前消费点 = "
    "crrna_design_agent.target_context 的 interface_deltas 块(随 design.json "
    "归档)与 webapp panel/design 页同源展示; IVT 8x4 矩阵训练集 schema 将含 "
    "这 4 列, 实测活性回填后进入选型分类器训练。")


def matrix_scaffold_of(library_desc):
    """orientation_library「+」分隔名 -> Chai 矩阵紧凑名(A8C+U15G -> A8CU15G)。"""
    return library_desc.replace("+", "")


def _classify(delta):
    """prot-crRNA 差量分级(方向词由数值符号决定, 不手写)。"""
    a = abs(delta)
    if a <= NEUTRAL_ABS:
        return "界面中性"
    if a >= COLLAPSE_ABS:
        return "界面塌陷" if delta < 0 else "界面大幅上升"
    return "界面中等下降" if delta < 0 else "界面中等上升"


def combo_reading(delta, iptm_sd):
    """逐组合判读: 全部由 delta 数值与显式阈值规则生成。"""
    d_pc = delta["delta_prot_crRNA"]
    noise = ("|Δ| <= %dsd 噪声代理(判 null)" % NULL_SD_K
             if abs(d_pc) <= NULL_SD_K * iptm_sd else
             "超 %dsd 噪声代理" % NULL_SD_K)
    return ("[主口径] Δprot-crRNA %+.4f → %s(%s); [参考] Δagg ipTM %+.4f / "
            "ΔcrRNA-靶 %+.4f / Δclash %+.4f"
            % (d_pc, _classify(d_pc), noise, delta["delta_iptm_mean"],
               delta["delta_crRNA_target"], delta["delta_clash_frac"]))


def _pick_source():
    for name in MATRIX_CANDIDATES:
        path = os.path.join(DATA, name)
        if os.path.isfile(path):
            return path
    raise SystemExit("缺 Chai 矩阵(候选: %s)" % ", ".join(MATRIX_CANDIDATES))


def build(matrix_path=None):
    """从 Chai 矩阵 JSON 构建界面特征差量表(全部数值源自矩阵, 不手写)。"""
    matrix_path = matrix_path or _pick_source()
    mat = json.load(open(matrix_path, encoding="utf-8"))
    rows = mat["rows"]
    targets = sorted({r["target"] for r in rows})
    scaffolds = sorted({r["scaffold"] for r in rows})
    by = {(r["scaffold"], r["target"]): r for r in rows}

    # 网格完整性: 期望 骨架 x 靶标 全网格, 缺失组合显式列出
    grid_missing = [(s, t) for s in scaffolds for t in targets
                    if (s, t) not in by]
    wt_reference = {}
    for t in targets:
        wt = by.get(("WT", t))
        if wt is None:
            raise SystemExit("矩阵缺 WT 参照组合(WT, %s), 差量无锚点" % t)
        wt_reference[t] = {f: wt[f] for f in RAW_FIELDS}

    combos = []
    for s in scaffolds:
        if s == "WT":
            continue
        for t in targets:
            if (s, t) not in by:
                continue
            r = by[(s, t)]
            raw = {f: r[f] for f in RAW_FIELDS}
            delta = {"delta_" + f: round(r[f] - wt_reference[t][f], 4)
                     for f in DELTA_FIELDS}
            combos.append({"target": t, "scaffold": s, "raw": raw,
                           "delta": delta,
                           "reading": combo_reading(delta, r["iptm_sd"])})

    # 逐骨架跨靶判读: 与 chai_matrix_4t_summary.json 同源同函数(逐字一致)
    scaf_reading = {}
    for s in scaffolds:
        if s == "WT":
            continue
        deltas = {t: round(by[(s, t)]["prot_crRNA"] - wt_reference[t]["prot_crRNA"], 4)
                  for t in targets if (s, t) in by}
        sds = {t: by[(s, t)]["iptm_sd"] for t in deltas}
        if len(deltas) == len(targets):
            scaf_reading[s] = scaffold_reading(deltas, sds)

    # panel 覆盖: 智能体四靶标键 -> 矩阵列, 缺席键显式标注原因
    panel_cov = {}
    for key, mt in PANEL_TO_MATRIX_TARGET.items():
        if mt in targets:
            missing_s = [s for s in scaffolds if s != "WT"
                         and (s, mt) not in by]
            panel_cov[key] = {"matrix_target": mt, "covered": True,
                              "n_combos": len(scaffolds) - 1 - len(missing_s),
                              "missing_scaffolds": missing_s}
        else:
            panel_cov[key] = {
                "matrix_target": mt, "covered": False,
                "n_combos": 0, "missing_scaffolds": sorted(scaffolds),
                "reason": ("矩阵 %s 无靶标列 %s(现有列: %s); 该列随增量折叠"
                           "回填(crrna_chai_increment_merge.py -> "
                           "chai_matrix_v3aligned.json, 待 Chai 机时)后重跑"
                           "本脚本即自动覆盖"
                           % (os.path.basename(matrix_path), mt,
                              "/".join(targets)))}

    return {
        "generated_by": "scripts/crrna_chai_interface_features.py",
        "task": ("策划案 V3 §4.1 末段/表1 对齐(任务⑫): 骨架-靶标组合界面特征"
                 "差量作为选型分类器关键输入(结构描述特征口径)"),
        "source_matrix": {
            "file": os.path.relpath(matrix_path, ROOT).replace(os.sep, "/"),
            "engine": mat.get("engine"),
            "n_rows": len(rows), "targets": targets, "scaffolds": scaffolds,
            "preference": "chai_matrix_v3aligned.json 存在时优先, 否则 "
                          "chai_matrix_4t.json(候选序见 MATRIX_CANDIDATES)"},
        "delta_features": DELTA_FEATURES,
        "delta_definition": ("delta_<f>(骨架s, 靶标t) = <f>(s,t) - <f>(WT,t), "
                             "同靶 WT 组合为参照; 4 维 = %s" % DELTA_FEATURES),
        "feature_field_notes": {
            "iptm_mean": "aggregate ipTM 5 模型均值(整体复合物)",
            "prot_crRNA": "蛋白-crRNA 链对 ipTM(界面主口径)",
            "crRNA_target": "crRNA-靶RNA 链对 ipTM(参考口径; 本矩阵为合成完全"
                            "互补靶 0.45-0.49, 与旧矩阵原生靶分开引用)",
            "clash_frac": "5 模型带链间冲突比例(接触冲突指标)",
            "iptm_sd": "aggregate ipTM 5 模型 sd, 模型间噪声代理(prot-crRNA "
                       "per-model sd 未存盘)",
            "unpersisted": "表1 的界面PAE 与接触残基互补度未随矩阵存盘——以实际"
                           "存盘字段为准, 接触相关维度由链对 ipTM 与 clash_frac "
                           "承载, 不补造"},
        "methodology": {
            "self_template": ("自模板复述局限: 矩阵全部分数来自 100% 一致自模板"
                              "(8D4A链A)注入, 模板=查询自身的 WT 结构, 高 ipTM 为"
                              "模板合规性平凡结果(sanity check); 差量只描述自模板"
                              "口径下模型对骨架扰动的数值表型, 不能作为界面可预测"
                              "性证据"),
            "usability_band": ("WT prot-crRNA ipTM 四靶均值 %.4f 低于可用界面预测"
                               "区间(~0.5-0.6)" % (
                                   sum(wt_reference[t]["prot_crRNA"]
                                       for t in targets) / len(targets))),
            "retired_arm": ("B_break_compensate 行折叠对象为已退役茎破坏混杂臂"
                            "(2026-09-04 起湿实验侧已换去混杂新臂, 重折叠待 Chai "
                            "机时), 其'界面塌陷'是该混杂的可预测表型, 不得引为"
                            "当前分子的界面结论"),
            "controls_pending": ["template-free 运行", "scrambled-template 运行",
                                 "AF3/Protenix 独立引擎对照(data/af3_results/, "
                                 "data/protenix_results/ 待回填)"],
            "selector_integration": SELECTOR_INTEGRATION},
        "reading_rules": {
            "neutral_abs": NEUTRAL_ABS, "collapse_abs": COLLAPSE_ABS,
            "null_sd_k": NULL_SD_K,
            "note": "阈值常量 import 自 crrna_chai_matrix_summary(单一定义); "
                    "逐骨架跨靶判读与 chai_matrix_4t_summary.json 的 reading "
                    "同源同函数(scaffold_reading), 逐字一致"},
        "wt_reference": wt_reference,
        "combos": combos,
        "scaffold_reading": scaf_reading,
        "coverage": {
            "matrix_grid": {"expected": len(scaffolds) * len(targets),
                            "present": len(rows), "missing": grid_missing},
            "panel": panel_cov},
    }


def load_deltas(path=OUT_JSON):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def target_block(doc, panel_key):
    """模块三附加层消费点: panel 靶标的界面差量特征块(agent/webapp 同源)。

    返回 dict: available=True 时携带 WT 参照 + 全部非 WT 骨架差量行;
    available=False 时携带显式缺失原因(如 APC 列待增量折叠回填)。"""
    cov = (doc or {}).get("coverage", {}).get("panel", {}).get(panel_key)
    base = {"panel_key": panel_key,
            "source": (doc or {}).get("source_matrix", {}).get("file"),
            "delta_features": DELTA_FEATURES,
            "delta_definition": (doc or {}).get("delta_definition"),
            "limitation": LIMITATION}
    if not doc or not cov:
        return dict(base, available=False,
                    missing_reason="界面差量特征表无 panel 键 %s 的覆盖记录"
                                   % panel_key)
    if not cov.get("covered"):
        return dict(base, available=False, matrix_target=cov["matrix_target"],
                    missing_reason=cov["reason"])
    mt = cov["matrix_target"]
    rows = [c for c in doc["combos"] if c["target"] == mt]
    rows.sort(key=lambda c: c["delta"]["delta_prot_crRNA"])
    return dict(base, available=True, matrix_target=mt,
                wt_reference=doc["wt_reference"][mt],
                scaffolds=[{"scaffold": c["scaffold"], "delta": c["delta"],
                            "reading": c["reading"]} for c in rows],
                note="Δ 按主口径(prot-crRNA)升序; 自模板口径结构描述特征, "
                     "非界面可预测性证据(见 limitation)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--matrix", default=None,
                    help="Chai 矩阵路径(默认按候选序自动选择: v3aligned 优先)")
    ap.add_argument("--out", default=OUT_JSON)
    args = ap.parse_args()

    payload = build(args.matrix)
    json.dump(payload, open(args.out, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("输出 -> %s" % args.out)
    print("数据源: %s (%d 行)" % (payload["source_matrix"]["file"],
                                payload["source_matrix"]["n_rows"]))
    print("差量组合 %d 个; 网格缺失 %s" % (
        len(payload["combos"]), payload["coverage"]["matrix_grid"]["missing"]))
    for key, cov in payload["coverage"]["panel"].items():
        print("  [panel] %-12s -> %-12s %s" % (
            key, cov["matrix_target"],
            "覆盖 %d 组合" % cov["n_combos"] if cov["covered"]
            else "缺失: " + cov["reason"]))
    for s, txt in payload["scaffold_reading"].items():
        print("[%s] %s" % (s, txt))


if __name__ == "__main__":
    main()
