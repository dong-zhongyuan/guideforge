"""第 2 层(同族实验层)证据记录生成器。

文献自报数值为带出处的常量(引用 Creutzburg 2020 原文, 属事实引用而非本项目计算);
本项目复算统计自动读取 data/creutzburg_reverse_validation.json(脚本产出),
两者并列写入——杜绝"只引论文显著性"或"手改 JSON"两种状态。
运行: python scripts/crrna_layer2_evidence.py
输出: data/creutzburg_layer2_validation.json(整文件再生, 勿手改)
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

PAPER_MEASURES = [
    {"name": "stem-forming potential",
     "definition": "DR 茎 5 对碱基配对概率乘积(ViennaRNA bpp)",
     "pearson_rho": 0.31, "p_value": "<0.01",
     "guideforge_analog": "dp_fold(stem_intact_prob, 同为 bpp 茎配对概率乘积)"},
    {"name": "loop-accessibility",
     "definition": "DR loop 区平均未配对概率(ViennaRNA)",
     "pearson_rho": 0.29, "p_value": "<0.01",
     "guideforge_analog": "spacer_unpaired(全长平均未配对概率)"},
    {"name": "spacer opening energy",
     "definition": "spacer 前 19nt 打开能量(RNAplfold)",
     "pearson_rho": -0.38, "p_value": "<0.01",
     "guideforge_analog": "spacer_unpaired + d_seed(负向对应)"},
]


def main():
    rev_path = os.path.join(ROOT, "data", "creutzburg_reverse_validation.json")
    rev = json.load(open(rev_path, encoding="utf-8"))
    stats = rev["contexts"]["mature19"]["stats"]

    payload = {
        "date": "2026-09-01",
        "layer": "2 同族实验层(Creutzburg 2020 反向验证)",
        "source": "Creutzburg et al. 2020, NAR 48(6):3228-3243, PMC7102956, FnCas12a+AsCas12a",
        "finding": {
            "description": ("并列口径(引用时两句必须同时出现, 不得只引前半)。"
                            "(1) 论文自报: 该文计算并报告三个折叠特征与 FnCas12a 活性显著相关"
                            "(pre-crRNA 口径, n约25, 茎含手工定义假结配对)——支持[特征类型]有同族实验依据; "
                            "(2) 本项目复算(脚本 scripts/crrna_creutzburg_reverse.py, 数据集 "
                            "data/creutzburg_2020_dataset.json): 同文 Fig 1D 子集 n=14, 活性为读图近似值, "
                            "功效不足, 特征均未达显著且方向混合。复算不构成对论文结论的否定, 也不构成独立验证。"),
            "measures": PAPER_MEASURES,
            "n_spacers_paper": "~25 (图 1D/2A/2C/4B/6C 五组构建, 论文口径)",
            "recomputation_n14": {
                "mature19": stats,
                "note": "seed5_up 方向正确不显著; spacer_up 方向反转; stem_prob 各口径不一致; "
                        "Sp8(唯一 0% 活性)未被任何现有特征抓住: 活性态含 5' 假结在伪结外空间不可表示, 侵占型嵌套中间态好坏 guide 均大量出现(Sp8 93/318 vs Sp12 146/391)故伪结外交叉指标按构造无区分力——盲区开放"},
        },
        "conclusion": (
            "并列结论: (1) 特征[类型]有同族实验支持(论文自报 n约25 显著, p<0.01, rho 0.29-0.38); "
            "(2) 本项目子集复算(n=14, 活性图读)功效不足未达显著(方向混合), 且暴露 Sp8 盲区。"
            "两者并列引用; 该层支持的是特征类型与活性相关, 不支持具体排序预测。"),
        "machine_generated_by": "scripts/crrna_layer2_evidence.py (勿手改本文件)",
    }
    out = os.path.join(ROOT, "data", "creutzburg_layer2_validation.json")
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
