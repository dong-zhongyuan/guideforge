"""MD 预备数据聚合器: prep.json 汇总 + 分析输入派生(茎配对由 ViennaRNA 官方折叠计算)。

产出(均脚本生成, 勿手改):
  data/md_prep_manifests.json          四体系 prep 清单汇总
  data/md/<tag>/analysis_input.json    DR 茎配对(0基)等分析输入, 供
                                       crrna_md_analysis.py 读取(替代脚本内硬编码)
运行: python scripts/crrna_md_manifest.py
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import RNA  # noqa: E402
import crrna_scaffold_design as core  # noqa: E402

TAGS = ("wt", "A7C_U14G", "A7G_U14C", "U6C_A15G")


def stem_pairs_from_dr(dr_rna):
    """DR 单独 MFE 折叠(官方 RNA.fold)解析茎配对, 0 基 DR 内坐标。"""
    db, _mfe = core.fold(dr_rna)
    stack, pairs = [], []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.append((stack.pop(), i))
    return pairs, db


def main():
    plan = json.load(open(os.path.join(ROOT, "data", "md_analysis_plan.json"),
                          encoding="utf-8"))
    manifests = {}
    for tag in TAGS:
        pj = os.path.join(ROOT, "data", "md", tag, "prep.json")
        prep = json.load(open(pj, encoding="utf-8"))
        manifests[tag] = prep
        dr_rna = core.to_rna(prep["wt_dr_rna"])
        pairs, db = stem_pairs_from_dr(dr_rna)
        out = {"tag": tag, "dr_rna": dr_rna, "dr_len": len(dr_rna),
               "dr_mfe_db": db, "dr_stem_pairs_0based": pairs,
               "contact_json_hint": "data/8D4A_dr_contacts_*.json (5A 口径)",
               "criteria": {"preregistered": plan["preregistered"],
                            "M1_pass_delta_A": 0.5,   # 预注册: 候选茎RMSD中位数 <= WT+0.5A
                            "M2_pass_delta_pp": 10.0,  # 预注册: 接触存活率 >= WT-10个百分点
                            "main_criteria": plan["main_criteria"],
                            "source_plan": "data/md_analysis_plan.json"}}
        path = os.path.join(ROOT, "data", "md", tag, "analysis_input.json")
        json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("%s: DR 茎配对 %s" % (tag, pairs))
    manifests["runtime"] = ("GROMACS 2026.3 (容器 ~/dzy/envs/gmx, nompi_cuda, sm_86); "
                            "OpenMM(amber14 全套)参数化经 ParmEd 导出; "
                            "EM(仅 nb GPU)/NVT/NPT/MD(nb+pme+bonded GPU); 产出 10 ns")
    manifests["generated_by"] = "scripts/crrna_md_manifest.py (勿手改)"
    out = os.path.join(ROOT, "data", "md_prep_manifests.json")
    json.dump(manifests, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % out)


if __name__ == "__main__":
    main()
