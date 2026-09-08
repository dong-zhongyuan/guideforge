"""选型分类器训练: Han 2025 数据 + 本仓库候选库(2026-09-02, Fig1g 尺度版)。

两段式:
  第一段(文献预训练): Han 2025 工具箱变体 × 结构特征 × Fig1g 实测活性
    (归一化 GFP 表达, 低=抑制强; 与筛选图同尺度, 优于旧 RBS0/33 端点)
    -> 决策树回归(可解释): "什么 DR 结构特征区间活性最优"
  第二段(应用): 把第一段模型应用到本仓库 zengdr 口径全部通过变体
    -> 按文献模型预测活性排序候选库
  附加(Sanger 耐受集体检): 54 条真实 Sanger 变体特征经模型预测,
    与 204 候选分布对比 —— 检查模型是否把真实耐受变体排到离谱位置。

诚实边界: 干净监督对 = 工具箱 7 条 + 主文 Fig.1f 视觉转录 9 条
  (S3/S5/S7/S10/S14/S16/S20/S21/FL4, 2026-09-03 入库, 置信度分级见
  han2025_dataset.json fig1f_pairs; high=双读一致+Sanger 交叉校验,
  medium=单完整读+交叉校验)。145 个编号的完整映射仍未随论文发表
  (数值源数据已穷尽查证; 完整映射待作者回复, 见 docs/ 邮件草稿);
  LbCas12a CRISPRi ≠ Cas12a2 杀伤, 迁移假设已在 README 声明;
  决策树 = 答辩可解释。fig1f 留出验证(训练 7 条工具箱, 检验 9 条
  fig1f)作为转录配对的独立一致性体检随输出报告。
运行: python scripts/crrna_train_selector.py
输出: data/selector_model.json + data/selector_ranked_candidates.csv
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
from scipy.stats import rankdata
from sklearn.tree import DecisionTreeRegressor, export_text

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

DATA = os.path.join(ROOT, "data")

FEATURES = ["ddG_dr", "bp_dist", "cross_nt", "p_fold", "spacer_up"]

# 模块三靶标上下文特征(策划案 V3 §4.4 模块一 -> 模块三契约列):
# target_abundance_tier 由模块一靶标解析供给(data/target_abundance_tiers.json,
# 显式分档规则见 scripts/crrna_target_abundance.py); mutation_type 由模块一
# 最长 ORF 法突变类型判定供给(错义/无义/移码, 规则与边界处理显式声明于
# scripts/crrna_mutation_typing.py docstring, 2026-09-07 任务⑩)。
# 诚实边界: 文献训练行(Han/Fig1f/Teng)无靶标上下文, 决策树在这两维上都不能
# 分裂(mutation_type 为类目维, 文献行更没有突变类型标注); 丰度档维当前经显式
# 阈值门控(crrna_target_abundance.activation_gate)进入推荐链路, 突变类型维
# 当前经模块一判读(最长 ORF 法)进入特征向量与输出展示, IVT 8x4 矩阵训练集
# schema 将含这两列。webapp 端与本文档共用同一 FEATURES/TARGET_FEATURES 定义
# (crrna_agent_webapp.py 从本模块导入, 两端不得各自复制字面量)。
TARGET_FEATURES = ["target_abundance_tier", "mutation_type"]

# 结构层界面特征差量契约列(策划案 V3 §4.1 末段 / 表1「选型特征」, 2026-09-09
# 任务⑫): 每个骨架-靶标组合的 Chai-1 界面特征(aggregate ipTM / prot-crRNA
# 链对 ipTM / crRNA-靶链对 ipTM / 链间冲突比例)相对同靶 WT 组合的差量, 由
# scripts/crrna_chai_interface_features.py 从 Chai 矩阵 JSON 提取落盘
# (data/chai_interface_deltas.json, 差量定义与方法学注释见该文件)。
# 诚实边界(与上方 target_abundance_tier 先例同): 文献训练行(Han/Fig1f/Teng)
# 无结构数据, 决策树在这 4 维上不能分裂, 不捏造训练标签; 当前经显式附加层
# 进入推荐链路(crrna_design_agent.target_context 的 interface_deltas 块 +
# webapp panel/design 页同源展示), IVT 8x4 矩阵训练集 schema 将含这 4 列。
# 单一定义在提取模块(DELTA_FEATURES), 此处 import 共用, 两端不得各自复制
# 字面量; 数值为 100% 自模板口径的模板合规性表型, 仅作结构描述特征,
# 不构成界面可预测性证据(round-2/round-3 评审降级口径)。
from crrna_chai_interface_features import (  # noqa: E402
    DELTA_FEATURES as INTERFACE_DELTA_FEATURES)


def spearman(x, y):
    """tie-aware Spearman(midranks), 2026-09 round-3 R2 修复: 弃用手搓
    argsort-of-argsort(并列不取平均秩), 统一走 scipy.stats.spearmanr。
    与 scripts/crrna_creutzburg_reverse.py:spearman 同写法。"""
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic)


def load_han():
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    X, y, tags, confs = [], [], [], []
    for r in d["toolbox_features"]:
        if "fig1g" not in r:
            continue
        X.append([float(r[f]) for f in FEATURES])
        y.append(float(r["fig1g"]))
        tags.append(r["tag"])
        confs.append("toolbox")
    for r in d.get("fig1f_pairs") or []:
        if r.get("fig1g") is None:
            continue
        X.append([float(r[f]) for f in FEATURES])
        y.append(float(r["fig1g"]))
        tags.append(r["tag"])
        confs.append("fig1f:" + r["confidence"])
    return np.array(X), np.array(y), tags, d, confs


DEWEIRDT_SCAN = os.path.join(DATA, "raw", "deweirdt2020_dr_scan.json")


def load_deweirdt2020():
    """DeWeirdt 2020 alt-DR 扫描全量特征表(由 crrna_han2025_features.py 生成)。"""
    if not os.path.exists(DEWEIRDT_SCAN):
        return None
    return json.load(open(DEWEIRDT_SCAN, encoding="utf-8"))


ENDPOINTS = {"fig1g": "fig1g", "cis_end": "cis_end", "trans_end": "trans_end",
             "rbs0": "rbs0", "rbs33": "rbs33"}
# 方向统一: True = 数值越大越优(fig1g/rbs 为抑制读数取负号), 供池化模型共用语义
ENDPOINT_HIGHER_BETTER = {"fig1g": False, "cis_end": True, "trans_end": True,
                          "rbs0": False, "rbs33": False}


def load_pooled():
    """同源池化训练集(等湿实验数据的同源先验版, 2026-09-03):

    全部可用同源活性观测, 逐终点 z 标准化并统一"越大越优"方向:
      - Han 2025 工具箱 5 终点 x 7 条(fig1g 抑制 / cis / trans 切割 / RBS0/RBS33 两档表达)
      - Han 2025 Fig.1f 视觉转录 9 条(并入 fig1g 终点组)
      - Teng 2019 4n96 序数对(WT=-1, 4n96=+1; 跨物种 Fn/Cas12a 口径, 序数级)
    Tian 2025 RRS 96 条不入池(伪结外特征盲端, 见函数体注释), 作外部验证。
    已尝试并放弃的终点(SPR KD / SupFig8-10 原始荧光): 传感器预处理不可校验 /
    增长混杂方向与 Fig1g 相反——拒绝理由记录于 han2025_dataset.json。
    返回 (X, y, meta), meta 每行 = {dataset, tag, endpoint}。
    """
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    X, y, meta = [], [], []
    for ep, key in ENDPOINTS.items():
        vals = [(r, "han2025") for r in d["toolbox_features"] if key in r]
        if ep == "fig1g":  # Fig.1f 视觉转录 9 条并入 fig1g 终点组(2026-09-03)
            vals += [(r, "han2025_fig1f") for r in d.get("fig1f_pairs") or []
                     if r.get(key) is not None]
        if len(vals) < 5:
            continue
        ys = np.array([float(r[key]) for r, _ in vals])
        z = (ys - ys.mean()) / max(ys.std(ddof=0), 1e-9)
        if not ENDPOINT_HIGHER_BETTER[ep]:
            z = -z
        for (r, src), zi in zip(vals, z):
            X.append([float(r[f]) for f in FEATURES])
            y.append(float(zi))
            meta.append({"dataset": src, "tag": r["tag"],
                         "endpoint": ep})
    # 注(2026-09-03): Tian 2025 RRS 单点扫描 96 条不入池——RRS 区(DR 5' 端 4nt)
    # 活性由假结配对承载, 伪结外特征按构造不可表示(与 Creutzburg Sp8 盲区同源),
    # 混入训练实测把 Han 各终点 LOEO 全部拖差; 该集改作池化模型的独立外部验证,
    # 见 main() 的 tian2025 外部验证段。
    # Teng 2019 4n96 序数对(自身数据集内 z=±1; 特征以 canonical spacer 上下文计算)
    try:
        import crrna_scaffold_design as core
        import RNA
        tj = json.load(open(os.path.join(DATA, "teng4n96_validation.json"),
                            encoding="utf-8"))
        wt_top = json.load(open(os.path.join(
            DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
        spacer = core.to_rna(wt_top["spacer_fixed_dna"])
        # Teng 为 18nt PDB 口径 DR; 以其自身 WT 为基准算特征
        drs = {"Teng_WT": core.to_rna(tj["provenance"]["wt_dr_dna"]),
               "Teng_4n96": core.to_rna(tj["provenance"]["v4n96_dr_dna"])}
        base = drs["Teng_WT"]
        base_ss, _ = core.fold(base + spacer)
        _, base_mfe = core.fold(base)
        base_stem = core.stem_pairs_of(core.fold(base)[0])
        for tag, dr in drs.items():
            full = dr + spacer
            ss, _ = core.fold(full)
            _, dr_mfe = core.fold(dr)
            cross_nt, _ = core.cross_pairs(full, len(dr))
            pf = core.stem_intact_prob(full, base_stem) if base_stem else 1.0
            _, sp_up, _ = core.pf_stats(full, len(dr), 7)
            X.append([round(dr_mfe - base_mfe, 2),
                      RNA.bp_distance(base_ss, ss), cross_nt,
                      round(pf, 5), round(sp_up, 3)])
            y.append(-1.0 if tag == "Teng_WT" else 1.0)
            meta.append({"dataset": "teng2019", "tag": tag,
                         "endpoint": "ordinal_4n96"})
    except Exception as e:
        print("Teng 序数对跳过: %s" % e)
    # 注(2026-09-03): DeWeirdt 2020 alt-DR 35,883 条已尝试子样入池并回退——
    # 120 条 LFC 秩分层子样(占池 72%)实测把其余 5 终点 LOEO 全部拖负
    # (rbs33 +0.96->-0.75, rbs0 +0.68->-0.64, trans +0.25->-0.71,
    # fig1g +0.04->-0.70; 仅 cis -0.86->+0.32), 单终点淹没多终点共识语义,
    # 且 DeWeirdt 终点自身 LOEO 亦为负(-0.36)。拒绝入池, 改作 main() 的
    # 独立外部验证(与 Tian 2025 同模式)。
    return np.array(X), np.array(y), meta


def tolerance_neighborhood(row_feats):
    """54 条真实 Sanger 耐受变体的特征分布邻域度(半监督门控)。

    对候选的每个特征取其在耐受集分布中的百分位 p, 邻域度
    = 1 - mean(|p-50|)/50; 1 = 每维都落在真实耐受变体分布中心。
    """
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    tol = d.get("sanger_tolerance") or []
    if len(tol) < 10:
        return None
    devs = []
    for f, v in zip(FEATURES, row_feats):
        col = np.array([float(r[f]) for r in tol])
        p = float((col < v).mean() * 100.0)
        devs.append(abs(p - 50.0) / 50.0)
    return round(1.0 - float(np.mean(devs)), 3)


def load_han_endpoints():
    """同源多终点训练集: fig1g(抑制) + cis_end + trans_end(切割, 与旁切杀伤最近缘)"""
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    out = {}
    for ep, key in ENDPOINTS.items():
        X, y, tags = [], [], []
        src = [(r, "toolbox") for r in d["toolbox_features"] if key in r]
        if ep == "fig1g":
            src += [(r, "fig1f") for r in d.get("fig1f_pairs") or []
                    if r.get(key) is not None]
        for r, _tag in src:
            X.append([float(r[f]) for f in FEATURES])
            y.append(float(r[key]))
            tags.append(r["tag"])
        if len(y) >= 5:
            out[ep] = (np.array(X), np.array(y), tags)
    return out


def rank_family(models):
    """对 8 员跨型候选族(orientation_library, 含 compensatory)逐终点预测排序"""
    import crrna_scaffold_design as core
    import RNA
    lib = json.load(open(os.path.join(DATA, "orientation_library.json"),
                         encoding="utf-8"))
    spacer = core.to_rna(lib["spacer"])
    dr_wt = core.to_rna(lib["orientations"][0]["wt_construct"][:19]) if \
        lib["orientations"][0].get("wt_construct") else None
    wt_top = json.load(open(os.path.join(
        DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
    wt_dr = core.to_rna(wt_top["top"][0]["dr_seq"])
    wt_full_ss, _ = core.fold(wt_dr + spacer)
    _, wt_dr_mfe = core.fold(wt_dr)
    wt_stem = core.stem_pairs_of(core.fold(wt_dr)[0])
    fam = [{"desc": "WT", "dr": wt_dr}]
    seen = {"WT"}
    for p in lib["orientations"]:
        b = p["best"]
        if b["desc"] in seen:
            continue
        seen.add(b["desc"])
        fam.append({"desc": b["desc"], "dr": core.to_rna(b["dr_seq"])})
    rows = []
    for m in fam:
        dr = m["dr"]
        full = dr + spacer
        ss, _ = core.fold(full)
        _, dr_mfe = core.fold(dr)
        cross_nt, _ = core.cross_pairs(full, len(dr))
        p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
        _, sp_up, _ = core.pf_stats(full, len(dr), 7)
        feats = [round(dr_mfe - wt_dr_mfe, 2),
                 RNA.bp_distance(wt_full_ss, ss), cross_nt,
                 round(p_fold, 5), round(sp_up, 3)]
        row = {"scaffold": m["desc"],
               "features": dict(zip(FEATURES, feats))}
        for ep, dt in models.items():
            row["pred_" + ep] = round(float(dt.predict([feats])[0]), 4)
        rows.append(row)
    return rows


def load_our_candidates():
    import crrna_scaffold_design as core
    import RNA
    path = os.path.join(DATA, "tp53_r248q_zengdr.v6.variants.csv")
    wt_top = json.load(open(
        os.path.join(DATA, "tp53_r248q_zengdr.v6.top.json"), encoding="utf-8"))
    wt_dr = core.to_rna(wt_top["top"][0]["dr_seq"])
    spacer = core.to_rna(wt_top["spacer_fixed_dna"])
    wt_full_ss, _ = core.fold(wt_dr + spacer)
    wt_dr_ss, wt_dr_mfe = core.fold(wt_dr)
    wt_stem = core.stem_pairs_of(wt_dr_ss)
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["passed"] not in ("True", "true", "1"):
                continue
            dr = core.to_rna(r["dr_seq"])
            full = dr + spacer
            ss, mfe = core.fold(full)
            dr_ss, dr_mfe = core.fold(dr)
            cross_nt, _ = core.cross_pairs(full, len(dr))
            p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
            _, sp_up, _ = core.pf_stats(full, len(dr), 7)
            bp_d = RNA.bp_distance(wt_full_ss, ss)
            rows.append({
                "desc": r["desc"], "score": float(r["score"]),
                "ddG_dr": round(dr_mfe - wt_dr_mfe, 2),
                "bp_dist": bp_d, "cross_nt": cross_nt,
                "p_fold": round(p_fold, 5),
                "spacer_up": round(sp_up, 3),
                "construct_dna": r["construct_dna"]})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    args = ap.parse_args()

    # === 第一段: 文献预训练(Fig1g 尺度) ===
    X, y, tags, han, confs = load_han()
    print("=== 训练集(Fig1g 归一化 GFP, 低=抑制强; 工具箱 7 + Fig1f 转录 9) ===")
    print("%-6s %-16s" % ("var", "source"), " ".join("%10s" % f for f in FEATURES),
          "  Fig1g")
    for i, t in enumerate(tags):
        print("%-6s %-16s" % (t, confs[i]),
              " ".join("%10.3f" % v for v in X[i]), "  %6.4f" % y[i])

    # 决策树(小样本: max_leaf=4, min_samples_leaf=1)
    dt = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                               random_state=42)
    dt.fit(X, y)

    print("\n=== 决策树规则(可解释) ===")
    print(export_text(dt, feature_names=FEATURES, decimals=3))

    # 训练集拟合度(tie-aware: 树预测在本库仅个位数 distinct 值, 并列取平均秩)
    pred_train = dt.predict(X)
    rho_train = spearman(pred_train, y)
    print("训练集 Spearman rho = %.3f (tie-aware)" % rho_train)

    # fig1f 留出验证: 只用 7 条工具箱训练, 检验 9 条 Fig1f 转录配对
    fig1f_holdout = None
    idx_tb = [i for i, c in enumerate(confs) if c == "toolbox"]
    idx_f1f = [i for i, c in enumerate(confs) if c.startswith("fig1f")]
    if len(idx_tb) >= 5 and len(idx_f1f) >= 5:
        dt_tb = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                                      random_state=42)
        dt_tb.fit(X[idx_tb], y[idx_tb])
        pr = dt_tb.predict(X[idx_f1f])
        # tie-aware midranks(scipy rankdata, 与 spearmanr 同值): 4 叶树预测大量
        # 并列, 旧 argsort-of-argsort 对并列不取平均秩(2026-09-04 round-3 统一)
        rk = rankdata(pr)
        ry = rankdata(y[idx_f1f])
        rho_ho = float(np.corrcoef(rk, ry)[0, 1])
        mae_ho = float(np.mean(np.abs(pr - y[idx_f1f])))
        y_tr = y[idx_tb]
        # 2026-09-04 round-3 修复(泄漏): 旧口径 thr=(y.max()-y.min())/3 用含 9 条
        # 留出测试标签的全量 y 计算(y.max()=S20=0.928 本身是测试点), 是测试集的
        # 函数, "预登记"名不副实; 现只用训练折 y 值域。登记见 docs/preregistration.md §C。
        thr = float((y_tr.max() - y_tr.min()) / 3)  # 训练折可用阈: 训练 y 值域 1/3
        # 外推归因 + 逐点误差(2026-09-03: MAE 超阈的结构化诊断, 不改阈值不改数据)
        per_point = [{"tag": tags[k], "confidence": confs[k],
                      "y_true": round(float(y[k]), 4),
                      "y_pred": round(float(pr[j]), 4),
                      "abs_err": round(float(abs(pr[j] - y[k])), 4),
                      "label_in_train_range": bool(y[k] <= y_tr.max())}
                     for j, k in enumerate(idx_f1f)]
        n_out = sum(1 for p in per_point if not p["label_in_train_range"])
        # bootstrap 95% CI(固定种子, 对检验点重抽样 10000 次)
        rng = np.random.default_rng(42)
        n_ho = len(idx_f1f)
        boot_rho, boot_mae = [], []
        boot_nan = 0
        for _ in range(10000):
            b = rng.integers(0, n_ho, n_ho)
            if len(set(b.tolist())) < 3:
                continue
            c = float(np.corrcoef(rk[b], ry[b])[0, 1])
            if np.isnan(c):
                # 并列秩重抽样零方差 → ρ 未定义, 剔除并计数(tie-aware 连带修正)
                boot_nan += 1
                continue
            boot_rho.append(c)
            boot_mae.append(float(np.mean(np.abs(pr[b] - y[idx_f1f][b]))))
        rho_ci = [round(float(np.percentile(boot_rho, q)), 3)
                  for q in (2.5, 97.5)]
        mae_ci = [round(float(np.percentile(boot_mae, q)), 4)
                  for q in (2.5, 97.5)]
        fig1f_holdout = {
            "design": "train toolbox n=%d, test fig1f n=%d" % (
                len(idx_tb), len(idx_f1f)),
            "spearman": round(rho_ho, 3), "mae": round(mae_ho, 4),
            "spearman_method": "scipy rankdata midranks(与 scipy.stats.spearmanr "
                               "同值, tie-aware); 2026-09-04 round-3 统一, 旧口径 "
                               "argsort-of-argsort 对并列不取平均秩",
            "spearman_boot_ci95": rho_ci, "mae_boot_ci95": mae_ci,
            "spearman_boot_excluded_zero_variance": boot_nan,
            "mae_threshold_trainfold": round(thr, 4),
            "mae_threshold_note": "2026-09-04 round-3 修复: 旧字段 "
                "mae_threshold_prereg=(y.max()-y.min())/3=0.1617 用含 9 条留出测试"
                "标签的全量 y 计算(y.max()=S20=0.928 本身是测试点), 泄漏测试集信息, "
                "'预登记'名不副实; 现口径仅用 7 条训练折 y 值域的 1/3。该阈仍为数据"
                "导出(训练折)相对阈, 非外部绝对标准; 判据登记见 "
                "docs/preregistration.md §C",
            "train_y_range": [round(float(y_tr.min()), 4),
                              round(float(y_tr.max()), 4)],
            "test_labels_above_train_range": n_out,
            "per_point": per_point,
            "extrapolation_note": "决策树预测上界=训练 y 最大值; 全部 %d 条 fig1f "
                "标签高于训练上界 %.3f → MAE 超阈主因是 Fig1g 主图与 Fig1f 附图两图版"
                "间的水平偏移叠加树模型不可外推, 非排序能力失败; 本用途(先验排序)"
                "的有效判据是秩一致性(Spearman 及其 CI)" % (n_out, y_tr.max()),
            "note": "Fig1f 视觉转录配对对工具箱模型的独立一致性体检; "
                    "Spearman>0 且 MAE<训练折 y 值域的 1/3 视为转录配对可用; "
                    "MAE 超阈后果: fig1f 配对仅用于秩级体检, 绝对水平不采用"}
        print("\n=== fig1f 留出验证(训练%d工具箱 -> 检验%d转录对) ===" % (
            len(idx_tb), len(idx_f1f)))
        print("Spearman = %+.3f [CI95 %+.3f, %+.3f], MAE = %.4f [CI95 %.4f, %.4f] "
              "(训练折阈值 %.4f, 训练折 y 范围 %.3f-%.3f)" % (
                  rho_ho, rho_ci[0], rho_ci[1], mae_ho, mae_ci[0], mae_ci[1],
                  thr, y_tr.min(), y_tr.max()))
        print("外推归因: %d/%d 条检验标签高于训练上界 %.3f" % (
            n_out, n_ho, y_tr.max()))

    # 特征重要性
    print("\n=== 特征重要性 ===")
    for f, imp in sorted(zip(FEATURES, dt.feature_importances_),
                         key=lambda x: -x[1]):
        print("  %-12s %.3f" % (f, imp))

    # === 第二段: 应用到候选库 ===
    cands = load_our_candidates()
    X_ours = np.array([[float(r[f]) for f in FEATURES] for r in cands])
    preds = dt.predict(X_ours)

    ranked = sorted(zip([r["desc"] for r in cands], preds,
                        [float(r["ddG_dr"]) for r in cands],
                        [float(r["score"]) for r in cands]),
                    key=lambda x: x[1])  # 低=好

    print("\n=== 文献模型对本仓库 %d 条通过变体的排序(前 20) ===" % len(cands))
    print("%-16s %10s %10s %10s" % ("变体", "文献预测", "DDR-ΔΔG", "管线得分"))
    for name, p, ddg, score in ranked[:20]:
        print("%-16s %10.1f %10.2f %10.4f" % (name, p, float(ddg), float(score)))

    # 与管线综合分的 Spearman(tie-aware, 2026-09 round-3 R3 修复)
    # 旧口径 np.argsort(np.argsort()) 对并列不取平均秩且默认排序不稳定:
    # 决策树在本候选库上仅输出个位数 distinct 预测值(204 条大量并列),
    # tie-blind 统计量随 numpy 版本/输入顺序摆动(同数据实测: 本机文件序 +0.587,
    # 逆序 -0.493, 快照 8c33354 的 artifact 为 -0.872), 属纯伪影, 不构成
    # "两排序系统反向"的证据。改用 scipy midrank 并附敏感度记录。
    from scipy.stats import spearmanr as _spearmanr
    our_scores = np.array([float(s) for _, _, _, s in ranked])
    lit_preds = np.array([p for _, p, _, _ in ranked])
    _sp = _spearmanr(lit_preds, our_scores)
    rho, rho_p = float(_sp.statistic), float(_sp.pvalue)
    rho_tb_fwd = float(np.corrcoef(np.argsort(np.argsort(lit_preds)),
                                   np.argsort(np.argsort(our_scores)))[0, 1])
    rho_tb_rev = float(np.corrcoef(np.argsort(np.argsort(lit_preds[::-1])),
                                   np.argsort(np.argsort(our_scores[::-1])))[0, 1])
    n_leaf_vals = len(set(np.round(lit_preds, 6)))
    # 顶部集合重合(K=12): 方向口径 lit_pred 低=好, score 高=优
    sel_top = {nm for nm, *_ in ranked[:12]}
    pipe_top = {r["desc"] for r in
                sorted(cands, key=lambda r: -float(r["score"]))[:12]}
    top12_overlap = {"k": 12,
                     "intersection": len(sel_top & pipe_top),
                     "jaccard": round(len(sel_top & pipe_top)
                                      / len(sel_top | pipe_top), 3),
                     "selector_top12": sorted(sel_top),
                     "pipeline_top12": sorted(pipe_top)}
    print("\n文献模型 vs 管线综合分 Spearman(tie-aware) = %+.3f (p=%.4f)"
          % (rho, rho_p))
    print("  方向口径: lit_pred 低=预测活性强, score 高=优 -> 排序一致性读 "
          "goodness 序相关 = %+.3f" % -rho)
    print("  tie-blind 旧口径敏感度: 正序 %+.3f / 逆序 %+.3f (distinct 预测值 %d 个,"
          "伪影; 快照 8c33354 artifact -0.872 同源)" % (rho_tb_fwd, rho_tb_rev,
                                                       n_leaf_vals))
    print("  TOP-12 重合: %d/12 (Jaccard %.3f)" % (top12_overlap["intersection"],
                                                  top12_overlap["jaccard"]))

    # === Sanger 耐受集体检: 真实变体在模型预测分布中的位置 ===
    tol_stats = None
    if han.get("sanger_tolerance"):
        tol = [r for r in han["sanger_tolerance"]
               if r["group"] in ("flank", "loop")]  # 与 Cas12a2 DR 口径近缘
        X_tol = np.array([[float(r[f]) for f in FEATURES] for r in tol])
        tol_pred = dt.predict(X_tol)
        cand_pred = dt.predict(X_ours)
        pct = float(np.mean([(cand_pred < p).mean() for p in tol_pred]))
        print("\n=== Sanger 耐受集体检(n=%d, flank+loop) ===" % len(tol))
        print("真实耐受变体预测活性中位数 = %.4f (工具箱锚点范围 %.4f-%.4f)"
              % (np.median(tol_pred), y.min(), y.max()))
        print("耐受变体在 204 候选预测分布中的平均分位 = %.1f%%" % (pct * 100))
        tol_stats = {
            "n": len(tol),
            "pred_median": round(float(np.median(tol_pred)), 4),
            "pred_min": round(float(tol_pred.min()), 4),
            "pred_max": round(float(tol_pred.max()), 4),
            "avg_percentile_in_candidates": round(pct, 3),
            "by_group": {g: round(float(np.median(
                [p for r, p in zip(tol, tol_pred) if r["group"] == g])), 4)
                for g in sorted({r["group"] for r in tol})}}

    # === 同源多终点: cis/trans 切割终点模型 + 8 员家族逐终点排序 ===
    multi = {ep: (X_e, y_e) for ep, (X_e, y_e, _) in
             load_han_endpoints().items()}
    models = {"fig1g": dt}
    print("\n=== 同源多终点模型(cis/trans = LbCas12a 切割, trans 与旁切杀伤最近缘) ===")
    for ep, (X_e, y_e) in multi.items():
        m = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                                  random_state=42)
        m.fit(X_e, y_e)
        models[ep] = m
        print("  %-9s n=%d  训练范围 %.4f-%.4f" % (ep, len(y_e),
                                                  y_e.min(), y_e.max()))
    fam_rows = rank_family(models)
    # 池化模型 + 耐受集邻域度(同源数据全量训练, 等湿实验回流前的先验版)
    Xp, yp, pmeta = load_pooled()
    pooled = DecisionTreeRegressor(max_leaf_nodes=6, min_samples_leaf=2,
                                   random_state=42)
    pooled.fit(Xp, yp)
    print("\n=== 同源池化模型(n=%d 观测 = Han 5 终点 x7 + Fig1f 转录 x9(fig1g 组) + Teng 序数对) ===" % len(yp))
    print(export_text(pooled, feature_names=FEATURES, decimals=3))
    # 留一终点交叉验证(方向一致性)
    loeo = {}
    for ep in sorted({m["endpoint"] for m in pmeta}):
        idx_tr = [i for i, m in enumerate(pmeta) if m["endpoint"] != ep]
        idx_te = [i for i, m in enumerate(pmeta) if m["endpoint"] == ep]
        if len(idx_te) < 5:
            continue
        m2 = DecisionTreeRegressor(max_leaf_nodes=6, min_samples_leaf=2,
                                   random_state=42)
        m2.fit(Xp[idx_tr], yp[idx_tr])
        pr = m2.predict(Xp[idx_te])
        loeo[ep] = round(spearman(pr, yp[idx_te]), 3)
        print("  LOEO %-9s (n=%d): Spearman(tie-aware)=%+.3f" % (
            ep, len(idx_te), loeo[ep]))
    for r in fam_rows:
        fv = [r["features"][f] for f in FEATURES]
        r["pred_pooled_z"] = round(float(pooled.predict([fv])[0]), 3)
        r["tolerance_neighborhood"] = tolerance_neighborhood(fv)
    print("\n=== 8 员族: 池化模型 z 预测 + 真实耐受集邻域度 ===")
    print("%-22s %10s %12s" % ("scaffold", "pooled_z", "tol_neighbor"))
    for r in sorted(fam_rows, key=lambda x: -x["pred_pooled_z"]):
        print("%-22s %10.3f %12s" % (r["scaffold"], r["pred_pooled_z"],
                                     r["tolerance_neighborhood"]))
    print("\n=== 8 员跨型候选族逐终点预测(fig1g/cis 低=抑制强; trans 高=切割强) ===")
    print("%-22s %10s %10s %10s" % ("scaffold", "fig1g", "cis_end", "trans_end"))
    for r in fam_rows:
        print("%-22s %10.4f %10.3f %10.3f" % (
            r["scaffold"], r.get("pred_fig1g", float("nan")),
            r.get("pred_cis_end", float("nan")),
            r.get("pred_trans_end", float("nan"))))

    # === Tian 2025 RRS 外部验证(2026-09-03): 池化模型(训练集不含 Tian)
    #     预测 96 条独立 RRS 单点扫描, 量化伪结外特征对 RRS 区的盲端 ===
    rrs_check = None
    rrs = [r for r in han.get("tian2025_rrs_pairs") or []
           if r.get("rrs_rel") is not None]
    if len(rrs) >= 20:
        Xr = np.array([[float(r[f]) for f in FEATURES] for r in rrs])
        yr = np.array([float(r["rrs_rel"]) for r in rrs])
        pr = pooled.predict(Xr)
        rho_all = spearman(pr, yr)
        per_crna, n_undef = [], 0
        for nn in sorted({int(r["crRNA"]) for r in rrs}):
            idx = [i for i, r in enumerate(rrs) if int(r["crRNA"]) == nn]
            if len(idx) >= 5:
                # tie-aware 且在子集内重取秩(旧口径沿用全局秩切片, 并列处理不一致);
                # 子集内预测全并列(树模型在该 crRNA 上无区分)时 ρ 未定义, 剔除并计数
                rho_c = spearman(pr[idx], yr[idx])
                if np.isnan(rho_c):
                    n_undef += 1
                else:
                    per_crna.append(rho_c)
        rrs_check = {
            "n": len(rrs),
            "spearman_pooled_vs_rrs_rel": round(rho_all, 3),
            "per_crna_spearman_median": round(float(np.median(per_crna)), 3)
            if per_crna else None,
            "per_crna_n_valid": len(per_crna),
            "per_crna_n_undefined_all_tied": n_undef,
            "method": "scipy.stats.spearmanr(tie-aware midranks), 逐 crRNA 在子集"
                      "内重取秩; 2026-09-04 round-3 统一",
            "design": "池化模型训练集不含任何 Tian 2025 数据(独立外部验证)",
            "interpretation": "RRS 区(DR 5' 端 4nt)活性由假结配对承载(Tian 2025 "
                              "Fig.1: U+3/U+4 与茎环互作), 伪结外特征按构造不可表示"
                              "——Spearman≈0 为预期盲端, 与 Creutzburg Sp8 盲区同源; "
                              "结论: 选型器排序不适用于 DR 5' 端 RRS 位点变体, "
                              "该位点区候选须按位置规则另行排除"}
        print("\n=== Tian2025 RRS 外部验证(n=%d) ===" % len(rrs))
        print("整体 Spearman = %+.3f; 逐 crRNA 中位 = %+.3f (有效 %d, 全并列剔除 %d)"
              % (rho_all, rrs_check["per_crna_spearman_median"] or float("nan"),
                 len(per_crna), n_undef))
        print("解读: 伪结外特征对 RRS 区盲端(预期内), 选型器不适于 RRS 位点变体")

    # === DeWeirdt 2020 alt-DR 外部验证(2026-09-03): 池化模型(训练集不含
    #     DeWeirdt)预测全量 35,883 条变体两方向, 检验同源迁移对大规模 DR
    #     功能扫描的外推(茎/环区在变体表示范围内, 与 RRS 盲端互补的对照) ===
    dw_check = None
    scan = load_deweirdt2020()
    if scan:
        dw_check = {"design": "池化模型训练集不含任何 DeWeirdt 2020 数据"
                              "(独立外部验证, 全量无剔除)",
                    "method": "scipy.stats.spearmanr(tie-aware midranks); "
                              "2026-09-04 round-3 统一",
                    "direction": "LFC 低=活性强, 模型 z 高=预测优 -> 报告 "
                                 "Spearman(pred, -LFC), 正=方向正确"}
        for ori, fkey, lkey in (("127", "features_127", "lfc127"),
                                ("128", "features_128", "lfc128")):
            idx = [i for i, r in enumerate(scan["rows"])
                   if r.get(fkey) and r.get(lkey) is not None]
            Xd = np.array([[float(scan["rows"][i][fkey][f]) for f in FEATURES]
                           for i in idx])
            yd = np.array([-float(scan["rows"][i][lkey]) for i in idx])
            pr = pooled.predict(Xd)
            rho_dw = spearman(pr, yd)
            # 同数据集 5 折 CV(归因诊断: 区分"特征在该体系内无信息"
            # 与"跨数据集迁移失败"; 树超参与池化模型一致)
            from sklearn.model_selection import KFold
            pr_cv = np.zeros(len(idx))
            for tr, te in KFold(n_splits=5, shuffle=True,
                                random_state=42).split(Xd):
                m3 = DecisionTreeRegressor(max_leaf_nodes=6,
                                           min_samples_leaf=2,
                                           random_state=42)
                m3.fit(Xd[tr], yd[tr])
                pr_cv[te] = m3.predict(Xd[te])
            rho_cv = spearman(pr_cv, yd)
            dw_check["ori" + ori] = {
                "n": len(idx),
                "spearman_pred_vs_neglfc": round(rho_dw, 3),
                "within_dataset_cv5_spearman": round(rho_cv, 3)}
            print("\n=== DeWeirdt2020 alt-DR 外部验证(pRDA_%s, n=%d) ==="
                  % (ori, len(idx)))
            print("Spearman(池化预测 vs -LFC) = %+.3f; 同数据集 5 折 CV = %+.3f"
                  % (rho_dw, rho_cv))
    else:
        # 2026-09-04: 本机缺 data/raw/deweirdt2020_dr_scan.json 时沿用既有
        # homolog_training.json 的全量数据运行结果并如实标注口径, 防止静默丢失
        # (旧行为: dw_check=None 直接把上一轮的 DeWeirdt 外部验证结果覆盖为 null)
        prev_path = os.path.join(DATA, "homolog_training.json")
        prev_dw = None
        if os.path.exists(prev_path):
            try:
                prev_dw = json.load(open(prev_path, encoding="utf-8")).get(
                    "deweirdt2020_external_check")
            except Exception:
                prev_dw = None
        if prev_dw:
            dw_check = dict(prev_dw)
            dw_check["preserved_from_previous_run"] = (
                "本机缺 data/raw/deweirdt2020_dr_scan.json, 本次未重算; 保留自上一"
                "全量数据运行(旧 tie-blind Spearman 口径), 待数据齐备后按 "
                "tie-aware 口径重算")
            print("\n=== DeWeirdt2020 外部验证: 本机缺 data/raw, 沿用上一全量"
                  "运行结果并标注(见 JSON preserved_from_previous_run) ===")

    # === 界面特征差量(任务⑫): 结构层骨架-靶标组合特征, 附加层接入声明 ===
    if_path = os.path.join(DATA, "chai_interface_deltas.json")
    if_block = {"status": "missing",
                "note": "缺 %s —— 先运行 scripts/crrna_chai_interface_features.py"
                        % if_path}
    if os.path.isfile(if_path):
        ifdoc = json.load(open(if_path, encoding="utf-8"))
        cov = ifdoc.get("coverage", {}).get("panel", {})
        if_block = {
            "status": "loaded",
            "source": "data/chai_interface_deltas.json",
            "source_matrix": ifdoc["source_matrix"]["file"],
            "n_combos": len(ifdoc.get("combos", [])),
            "panel_covered": sorted(k for k, v in cov.items()
                                    if v.get("covered")),
            "panel_missing": sorted(k for k, v in cov.items()
                                    if not v.get("covered"))}

    # 输出
    out = {"model": "DecisionTreeRegressor(max_leaf=4)",
           "usage_scope": "文献先验排序(相对次序)专用; 绝对预测值不可引用——"
                          "fig1f 留出 MAE 超训练折可用阈(归因见 fig1f_holdout_check; "
                          "旧'预登记'口径含测试标签泄漏, 2026-09-04 修复说明见 "
                          "mae_threshold_note)",
           "training_data": "Han 2025 (s41467-025-64010-z) LbCas12a "
                            "toolbox n=7 + Fig.1f 视觉转录 n=9 "
                            "(置信度分级见 han2025_dataset.json fig1f_pairs)",
           "fig1f_holdout_check": fig1f_holdout,
           "features": FEATURES,
           "target_context_features": TARGET_FEATURES,
           "target_context_note": "靶标上下文契约列(模块一供给)不经文献决策树"
               "消费——文献训练行无靶标上下文, 决策树在这两维上都不能分裂: "
               "靶RNA丰度档位经显式阈值门控 crrna_target_abundance.activation_gate "
               "进入推荐链路(数据源 data/target_abundance_tiers.json); 突变类型"
               "(mutation_type, 类目维: 错义/无义/移码)经最长 ORF 法判定进入模块一"
               "特征向量与输出展示(scripts/crrna_mutation_typing.py, 完整 ORF 证据"
               "在 data/agent/*.design.json 的 mutation.typing 块); IVT 8x4 矩阵"
               "训练集 schema 将含这两列",
           "interface_delta_features": INTERFACE_DELTA_FEATURES,
           "interface_delta_note": "结构层界面特征差量(策划案 V3 §4.1 末段/表1 "
               "选型特征, 任务⑫): 每个骨架-靶标组合的 Chai-1 界面特征相对同靶 "
               "WT 组合的差量(4 维契约列 = crrna_chai_interface_features."
               "DELTA_FEATURES, import 共用单一定义)。诚实边界: 文献训练行"
               "(Han/Fig1f/Teng)无结构数据, 文献决策树在这 4 维上不能分裂——"
               "与 target_abundance_tier 先例相同, 当前经显式附加层进入推荐链路"
               "(crrna_design_agent.target_context 的 interface_deltas 块与 "
               "webapp panel/design 页同源展示), 不捏造训练标签; IVT 8x4 矩阵"
               "训练集 schema 将含这 4 列, 实测活性回填后进入选型分类器训练。"
               "数值为 100% 一致自模板口径的模板合规性表型(局限声明见 "
               "data/chai_interface_deltas.json methodology), 仅作结构描述"
               "特征, 不构成界面可预测性证据(round-2/round-3 评审降级口径)",
           "interface_delta_data": if_block,
           "feature_importance": {f: float(imp) for f, imp in
                                  zip(FEATURES, dt.feature_importances_)},
           "tree_rules": export_text(dt, feature_names=FEATURES),
           "train_spearman": round(float(rho_train), 3),
           "applied_to": "tp53_r248q_zengdr.v6 (%d passed)" % len(cands),
           "vs_pipeline_score_spearman": round(float(rho), 3),
           "vs_pipeline_score_spearman_pval": round(float(rho_p), 4),
           "vs_pipeline_score_method": "scipy.stats.spearmanr(tie-aware midranks); "
               "2026-09 round-3 R3 修复: 旧口径 argsort(argsort()) 对并列不取平均秩,"
               "且决策树预测仅个位数 distinct 值, 结果随 numpy 版本/输入顺序摆动",
           "vs_pipeline_score_direction": "lit_pred 低=预测活性强(Fig1g 抑制读数), "
               "pipeline score 高=优; 两排序系统的一致性应读 goodness 序相关 = "
               "-vs_pipeline_score_spearman",
           "vs_pipeline_score_tieblind_artifact": {
               "forward_order": round(rho_tb_fwd, 3),
               "reversed_order": round(rho_tb_rev, 3),
               "snapshot_8c33354_artifact": -0.872,
               "distinct_pred_values": n_leaf_vals,
               "note": "tie-blind 统计量在少 distinct 值预测上随输入顺序/numpy 版本"
                       "摆动(+0.587/-0.493/-0.872 同源伪影), 不构成两排序系统"
                       "真实反向分歧的证据; 以 tie-aware 值为准"},
           "top12_overlap_pipeline_vs_selector": top12_overlap,
           "ranking_authority": "pipeline_score",
           "ranking_authority_reason":
               "合成候选排序以管线打分(score_variant, 同一 Cas12a2 zeng2026 构建的 "
               "ViennaRNA 结构口径)为准。(1) 本选型器是 n=16 LbCas12a CRISPRi 终点"
               "训练的文献先验, fig1f 留出 MAE 超训练折可用阈、LOEO 在 cis 端点反预测"
               "(%.3f, 本运行计算值), 按 usage_scope 仅供先验排序参考; (2) 决策树在本候选库仅"
               "输出个位数 distinct 预测值, 对多数保守变体无区分力; (3) 管线打分经 "
               "filter+random 基线检验(data/selection_baseline.cas12a2_zeng2026.json)"
               "确认相对'硬过滤+随机'有大效应排序信息量。" % loeo["cis_end"],
           "divergence_note":
               "修正后两排序并非反向: tie-aware rho(lit_pred, score)=%+.3f "
               "(p=%.4f), 按方向口径折算 goodness 序相关=%+.3f, 仅为弱残余分歧; "
               "TOP-12 重合 %d 条(Jaccard %.3f), 顶部集合大体一致。残余分歧的可能"
               "解释: 选型器偏好 DR 单独强稳定化(ddG_dr<=-1.15 叶节点), 而管线在"
               "全长口径下 ddG 仅 0.1 罚/0.3 赏; 且训练终点为 LbCas12a CRISPRi GFP "
               "抑制, 与 Cas12a2 旁切杀伤存在体系迁移差距。快照 8c33354 artifact 的 "
               "-0.872 为 tie-blind 伪影(见 vs_pipeline_score_tieblind_artifact)。"
               % (rho, rho_p, -rho, top12_overlap["intersection"],
                  top12_overlap["jaccard"]),
           "sanger_tolerance_check": tol_stats,
           "multi_endpoint": {"endpoints": sorted(models),
                              "note": "cis/trans 终点比 = 同源(LbCas12a)切割先验, "
                                      "trans 与 Cas12a2 旁切杀伤最近缘; "
                                      "fig1g n=16(含 Fig1f 转录 9), cis/trans n=7"},
           "family_ranking_8": fam_rows,
           "top10_lit_model": [{"desc": n, "lit_pred": round(p, 4),
                                "ddG_dr": float(d), "score": float(s)}
                               for n, p, d, s in ranked[:10]]}
    json.dump(out, open(os.path.join(DATA, "selector_model.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump({"family": fam_rows, "models": sorted(models)},
              open(os.path.join(DATA, "selector_family_ranking.json"), "w",
                   encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump({
        "status": "同源先验全量训练(等 IVT 实测替换); 池化目标 = 逐终点 z "
                  "标准化并统一'越大越优'(fig1g/rbs 取负, cis/trans 原向)",
        "n_observations": len(yp),
        "composition": {"han2025_5endpoints_x7": 35,
                        "han2025_fig1f_fig1g": sum(
                            1 for m in pmeta if m["dataset"] == "han2025_fig1f"),
                        "teng2019_ordinal": sum(
                            1 for m in pmeta if m["dataset"] == "teng2019")},
        "tian2025_external_check": rrs_check,
        "deweirdt2020_external_check": dw_check,
        "fig1f_note": "fig1g 终点组并入 Fig.1f 视觉转录 9 条(meta.dataset="
                      "han2025_fig1f 可追溯; 置信度 high 2/medium 7, 序列均经 "
                      "Sanger 独立源交叉校验)",
        "pooled_tree_rules": export_text(pooled, feature_names=FEATURES),
        "pooled_feature_importance": {f: float(imp) for f, imp in
                                      zip(FEATURES, pooled.feature_importances_)},
        "loeo_spearman": loeo,
        "spearman_method": "scipy.stats.spearmanr(tie-aware midranks); "
                           "2026-09-04 round-3 统一: LOEO/Tian2025/DeWeirdt 各 ρ "
                           "弃用手搓 argsort-of-argsort(并列不取平均秩), "
                           "与旧 tie-blind 口径的数值差异为并列秩处理所致",
        "rejected_endpoints": {
            "spr_kd": "传感图相位预处理无法对论文定性序(canonical 最强)校验, 拟合不收敛",
            "supfig8_10_raw_fluor": "0mM 诱导下面板间仍 3 倍差、10mM 方向与 "
                                    "Fig1g 相反——增长混杂未归一, 拒绝入库",
            "deweirdt2020_pooling": "120 条 LFC 秩分层子样入池实测(2026-09-03): "
                                    "单终点占池 72%, 其余 5 终点 LOEO 全拖负"
                                    "(rbs33 +0.96→−0.75, rbs0 +0.68→−0.64, "
                                    "trans +0.25→−0.71, fig1g +0.04→−0.70; "
                                    "deweirdt 自身 −0.36; 仅 cis 转正), "
                                    "拒绝入池, 全量改作独立外部验证"},
        "family_ranking_pooled": fam_rows},
        open(os.path.join(DATA, "homolog_training.json"), "w",
             encoding="utf-8"), ensure_ascii=False, indent=1)

    with open(os.path.join(DATA, "selector_ranked_candidates.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank_lit", "desc", "lit_pred_activity",
                    "ddG_dr", "pipeline_score"] +
                   [f"feat_{x}" for x in FEATURES])
        for i, (n, p, d, s) in enumerate(ranked, 1):
            row_idx = [j for j, r in enumerate(cands) if r["desc"] == n]
            feats = [cands[row_idx[0]][f] for f in FEATURES] if row_idx else []
            w.writerow([i, n, round(p, 4), float(d), float(s)] + feats)

    print("\n输出 -> data/selector_model.json / selector_ranked_candidates.csv")


if __name__ == "__main__":
    main()
