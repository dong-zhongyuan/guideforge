"""AlphaFold Server 就绪共折叠检查输入生成(策划案 V3 表1 蛋白预测层, 2026-09-04)。

定位: V3 表1/第1周计划写 AF3; 本机无 AF3 权重(需 DeepMind 审批)与 GPU 机时,
故产出 AlphaFold Server (alphafoldserver.com) 可直接上传的 JSON 任务文件
(dialect=alphafoldserver, version=1, 官方 server/README.md 口径); 人工提交后
把结果 zip 放回 data/af3_results/, 用 scripts/crrna_af3_collect.py 回填对比。
**格式勘误(2026-09-05)**: 旧版误用 AF3 开源代码库的 dialect=alphafold3(带链 id),
服务器只认 alphafoldserver 方言(顶层为任务列表, 实体键 proteinChain/rnaSequence,
链 id 由服务器自动分配); 且服务器**默认启用 PDB 模板**——8D4A/8D49 在库, 不设
useStructureTemplate=false 则预测**并非** template-free, 旧版"恰好 template-free"
的说法不成立。本版已强制蛋白链 useStructureTemplate=false, 才是真正的
template-free 独立对照层: 若 AF3 template-free 下 prot-crRNA 界面同样恢复,
"界面可预测性"才具备首个独立证据。

输入: WT + TOP-12 候选(data/tp53_r248q_zengdr.v6.top.json, 13 套, V3 docx 口径);
蛋白链 = SuCas12a2(1207aa, 取自 data/chai_input_WT.fasta, 与 Chai 矩阵同源);
靶 RNA = R248Q protospacer 窗口 24nt + PFS 5nt(CAGAG 突变体语境, 与 Chai 4t 同构)。

运行: python scripts/crrna_af3_input.py
输出: data/af3_inputs/GF_<desc>.json x 13 + GF_all13_alphafoldserver.json(合并单文件,
服务器支持一次上传多任务) + manifest.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

COMP = str.maketrans("ACGT", "TGCA")


def revcomp_rna(spacer_dna):
    return spacer_dna.translate(COMP)[::-1]


def read_chai_fasta(path):
    """读 Chai 头部约定的三链 FASTA, 返回 (蛋白, crRNA, 靶) 序列。"""
    seqs, name, chunks = {}, None, []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:], []
        elif line:
            chunks.append(line)
    if name:
        seqs[name] = "".join(chunks)
    prot = next(v for k, v in seqs.items() if k.startswith("protein|"))
    crna = next(v for k, v in seqs.items() if k == "rna|name=crRNA")
    return prot, crna


def main():
    prot, _wt_crna = read_chai_fasta(os.path.join(DATA, "chai_input_WT.fasta"))
    top = json.load(open(os.path.join(DATA, "tp53_r248q_zengdr.v6.top.json"),
                         encoding="utf-8"))["top"]
    spacer = json.load(open(os.path.join(DATA, "tp53_r248q_scan.summary.json"),
                            encoding="utf-8"))["spacer_dna"]
    pfs = json.load(open(os.path.join(DATA, "tp53_r248q_scan.summary.json"),
                         encoding="utf-8"))["sites"][0]["pfs"]
    target_rna = (revcomp_rna(spacer) + pfs).replace("T", "U")  # 窗口 24nt + 3' 下游 PFS 5nt(R248Q 突变体 CAGAG)

    out_dir = os.path.join(DATA, "af3_inputs")
    os.makedirs(out_dir, exist_ok=True)
    files = []
    jobs = []
    for e in top[:13]:
        desc = "WT" if e.get("desc") in (None, "WT") else e["desc"]
        crna_rna = e["construct_dna"].replace("T", "U")
        name = "GF_" + str(desc).replace("+", "_")
        job = {
            "name": name,
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": prot, "count": 1,
                                  "useStructureTemplate": False}},
                {"rnaSequence": {"sequence": crna_rna, "count": 1}},
                {"rnaSequence": {"sequence": target_rna, "count": 1}},
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
        fp = os.path.join(out_dir, name + ".json")
        json.dump([job], open(fp, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        jobs.append(job)
        files.append({"file": os.path.basename(fp), "desc": desc,
                      "crRNA_nt": len(crna_rna)})

    combined = os.path.join(out_dir, "GF_all13_alphafoldserver.json")
    json.dump(jobs, open(combined, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    json.dump({
        "generated_by": "scripts/crrna_af3_input.py (2026-09-05, alphafoldserver 方言重发版)",
        "purpose": ("策划案 V3 表1 AF3 复合物定性终检输入; AlphaFold Server 上传即用; "
                    "蛋白链 useStructureTemplate=false 强制 template-free"
                    "(服务器默认启用 PDB 模板, 8D4A/8D49 在库, 不设此项则非 template-free)"
                    " = Chai-1 自模板口径的独立对照层"),
        "dialect_note": ("dialect=alphafoldserver(官方 server/README.md 口径, 顶层任务列表, "
                         "链 id 服务器自动分配 A/B/C); 旧版误用 alphafold3 方言已废弃"),
        "chains": {"A": "SuCas12a2 蛋白 1207aa(data/chai_input_WT.fasta, 与 Chai 矩阵同源; 模板禁用)",
                   "B": "crRNA 构建 = DR 变体 19nt + R248Q spacer 24nt",
                   "C": "靶 RNA = R248Q protospacer 窗口 24nt + PFS CAGAG 5nt(与 Chai 4t 同构)"},
        "submission": ("https://alphafoldserver.com/ 逐文件上传, 或一次上传合并文件 "
                       "GF_all13_alphafoldserver.json(含全部 13 任务); "
                       "步骤与限额见 docs/af3_server_runbook.md; 结果 zip 放回 data/af3_results/"),
        "jobs": files,
    }, open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8"),
        ensure_ascii=False, indent=1)
    print("AF3 任务 %d 套(alphafoldserver 方言, template-free) -> %s (蛋白 %daa, 靶 RNA %dnt); 合并文件 %s" % (
        len(files), out_dir, len(prot), len(target_rna), os.path.basename(combined)))


if __name__ == "__main__":
    main()
