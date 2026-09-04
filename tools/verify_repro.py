# -*- coding: utf-8 -*-
"""GuideForge 交付包一键复现验证(Tier 0, 纯 CPU、无网络)。

逐项输出 PASS/FAIL:
  a. 依赖导入(RNA/numpy/Bio/scipy/sklearn/openpyxl/matplotlib/flask)
  b. pytest 全量测试
  c. 核心管线 smoke(全部输出到临时目录, 不污染包内 data/):
     1. crrna_scaffold_design.py 最小运行(仅 ViennaRNA)
     2. crrna_specificity_scan.py 小包内 fasta 扫描
     3. crrna_ivt_template.py 重生成模板+订单表, 与包内 data/ 同名文件逐字节 diff
     4. crrna_chai_matrix_summary.py 重跑, 与包内 data/chai_matrix_4t_summary.json
        逐字节 diff(确定性产物, 无时间戳/随机数)

用法: python tools/verify_repro.py
退出码: 0=全绿, 1=有 FAIL
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SPACER = "GTTCATGCCGCCCATGCAGGAACT"

IMPORTS = [("RNA", "viennarna"), ("numpy", "numpy"), ("Bio", "biopython"),
           ("scipy", "scipy"), ("sklearn", "scikit-learn"),
           ("openpyxl", "openpyxl"), ("matplotlib", "matplotlib"),
           ("flask", "flask")]

results = []


def report(name, ok, detail=""):
    results.append((name, ok, detail))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         (" — " + detail) if detail else ""), flush=True)


def run(cmd, **kw):
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8",
               PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", **kw)


def check_imports():
    for mod, pkg in IMPORTS:
        try:
            __import__(mod)
            report("依赖导入 %s (%s)" % (mod, pkg), True)
        except Exception as e:  # noqa: BLE001
            report("依赖导入 %s (%s)" % (mod, pkg), False, str(e))


def check_pytest():
    r = run([sys.executable, "-m", "pytest", "tests/", "-q",
             "-p", "no:cacheprovider"])
    tail = (r.stdout.strip().splitlines() or [""])[-1]
    report("pytest tests/ -q", r.returncode == 0, tail)


def file_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def smoke_scaffold_design(tmp):
    prefix = os.path.join(tmp, "run1")
    r = run([sys.executable, os.path.join("scripts", "crrna_scaffold_design.py"),
             "--effector", "cas12a2_zeng2026", "--spacer", SPACER,
             "--topk", "12", "--out-prefix", prefix])
    if r.returncode != 0:
        report("smoke: crrna_scaffold_design 最小运行", False,
               r.stderr.strip()[-300:])
        return
    outs = [prefix + s for s in (".variants.csv", ".top.json", ".top.fasta")]
    missing = [p for p in outs if not os.path.isfile(p)]
    if missing:
        report("smoke: crrna_scaffold_design 最小运行", False,
               "缺产物 %s" % missing)
        return
    top = json.load(open(prefix + ".top.json", encoding="utf-8"))["top"]
    ok = len(top) >= 2 and top[0]["desc"] == "WT"
    report("smoke: crrna_scaffold_design 最小运行", ok,
           "TOP 候选 %d 条, 首条 %s" % (len(top), top[0]["desc"]))


def smoke_specificity_scan(tmp):
    fa = os.path.join(tmp, "_wt_mut.fa")
    with open(fa, "wb") as out:
        for name in ("tp53_mrna_NM000546.fa", "tp53_r248q_mrna_NM000546.fa"):
            out.write(file_bytes(os.path.join(ROOT, "data", name)))
    prefix = os.path.join(tmp, "scan1")
    r = run([sys.executable, os.path.join("scripts", "crrna_specificity_scan.py"),
             "--spacer", SPACER, "--effector", "cas12a2_zeng2026",
             "--fasta", fa, "--out", prefix])
    if r.returncode != 0:
        report("smoke: crrna_specificity_scan 小 fasta", False,
               r.stderr.strip()[-300:])
        return
    outs = [prefix + s for s in (".sites.csv", ".summary.json")]
    missing = [p for p in outs if not os.path.isfile(p)]
    if missing:
        report("smoke: crrna_specificity_scan 小 fasta", False,
               "缺产物 %s" % missing)
        return
    summ = json.load(open(prefix + ".summary.json", encoding="utf-8"))
    report("smoke: crrna_specificity_scan 小 fasta", True,
           "summary 键 %d 个" % len(summ))


def byte_diff(name, regenerated, reference):
    if file_bytes(regenerated) == file_bytes(reference):
        return True, "逐字节一致"
    a, b = file_bytes(regenerated), file_bytes(reference)
    return False, ("不一致: 重生成 %d 字节 vs 包内 %d 字节" % (len(a), len(b)))


def smoke_ivt_template(tmp):
    t = os.path.join(tmp, "ivt_round1_template.csv")
    o = os.path.join(tmp, "ivt_round1_order_sheet.csv")
    r = run([sys.executable, os.path.join("scripts", "crrna_ivt_template.py"),
             "--out", t, "--out-order-sheet", o])
    if r.returncode != 0:
        report("smoke: crrna_ivt_template 重生成", False,
               r.stderr.strip()[-300:])
        return
    for path, ref_name in ((t, "ivt_round1_template.csv"),
                           (o, "ivt_round1_order_sheet.csv")):
        ok, detail = byte_diff(ref_name, path,
                               os.path.join(ROOT, "data", ref_name))
        report("smoke: ivt 重生成一致 (%s)" % ref_name, ok, detail)


def smoke_chai_summary(tmp):
    out = os.path.join(tmp, "chai_matrix_4t_summary.json")
    r = run([sys.executable, os.path.join("scripts", "crrna_chai_matrix_summary.py"),
             "--out", out])
    if r.returncode != 0:
        report("smoke: crrna_chai_matrix_summary 重跑", False,
               r.stderr.strip()[-300:])
        return
    ok, detail = byte_diff("chai_matrix_4t_summary.json", out,
                           os.path.join(ROOT, "data",
                                        "chai_matrix_4t_summary.json"))
    report("smoke: chai 矩阵汇总重生成一致", ok, detail)


def main():
    print("GuideForge 复现验证 — 包根: %s" % ROOT)
    print("Python: %s" % sys.executable)
    check_imports()
    check_pytest()
    tmp = tempfile.mkdtemp(prefix="gf_verify_")
    print("smoke 输出临时目录: %s" % tmp)
    smoke_scaffold_design(tmp)
    smoke_specificity_scan(tmp)
    smoke_ivt_template(tmp)
    smoke_chai_summary(tmp)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print("=" * 60)
    print("总计 %d 项, PASS %d, FAIL %d" % (len(results),
                                            len(results) - n_fail, n_fail))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
