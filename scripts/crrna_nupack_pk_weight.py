"""NUPACK 含假结配分竞争检验(Sp8 判别第四模型, 2026-09-02)。

工具: NUPACK 3.2.2 pfunc/-pseudo(Dirks-Pierce 2003 假结配分; 手册确认
pairs/mfe 不支持 -pseudo, 故唯一原语为配分函数)。
指标(单一 NUPACK 标尺):
  W_pk = (Z_pk - Z_nopk)/Z_pk = 假结构象的系综权重占比
依据: FnCas12a 体系活性态含 canonical pseudoknot(Creutzburg 自述),
Sp8 的失活=spacer 侵占(嵌套可表示); 若侵占压过假结活性态, W_pk 应最低。
判别检验(先定后跑): Sp8 的 W_pk 在 pre 口径 14 条中应排最末 2 位;
否则该模型亦 FAIL, out_of_model_domain 判定(四模型)成立。

编排: jobs 模式写任务文件(容器侧 /dawn); 本机管道执行
nupack_pk_runner.py 输出到宿主 /tmp/nupack_out.json; assemble 模式装配。
运行: python scripts/crrna_nupack_pk_weight.py jobs|assemble
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import crrna_scaffold_design as core  # noqa: E402

DAWN = "/public/home/mengxl/openclaw/dawn"
OUT_TMP = "/tmp/nupack_out.json"
DATA = os.path.join(ROOT, "data")


def build_jobs():
    ds = json.load(open(os.path.join(DATA, "creutzburg_2020_dataset.json"),
                        encoding="utf-8"))
    dr18 = core.to_rna(ds["dr18_dna"])
    pre5 = core.to_rna(ds["pre5_flank_dna"])
    jobs = {}
    for nm, sp in ds["spacers_21nt_dna"].items():
        jobs["%s_pre" % nm] = pre5 + dr18 + core.to_rna(sp)
    d = json.load(open(os.path.join(DATA, "tp53_r248q_pdbdr.v2.top.json"),
                       encoding="utf-8"))
    jobs["gf_WT"] = core.to_rna(d["top"][0]["construct_dna"])
    for e in d["top"]:
        if e.get("desc") not in (None, "WT") and e.get("rank", 99) <= 3:
            jobs["gf_%s" % e["desc"].replace("+", "_")] = core.to_rna(e["construct_dna"])
    json.dump(jobs, open(os.path.join(DAWN, "pk_jobs.json"), "w"))
    print("任务 %d 条 -> /dawn/pk_jobs.json" % len(jobs))


def assemble():
    res = json.load(open(OUT_TMP))
    ds = json.load(open(os.path.join(DATA, "creutzburg_2020_dataset.json"),
                        encoding="utf-8"))
    act = ds["activity_fig1d_pct_approx"]
    rows = []
    for nm in sorted(act, key=lambda s: int(s[2:])):
        r = res["%s_pre" % nm]
        try:
            z_pk, z_nopk = r["pk"]["Z"], r["nopk"]["Z"]
            w_pk = (z_pk - z_nopk) / z_pk
        except Exception:
            w_pk, z_pk, z_nopk = None, None, None
        rows.append({"name": nm, "activity": act[nm],
                     "Z_pk": z_pk, "Z_nopk": z_nopk,
                     "W_pk": round(w_pk, 4) if w_pk is not None else None})
    rows_ok = [r for r in rows if r["W_pk"] is not None]
    rows_ok.sort(key=lambda r: r["W_pk"])
    print("%-5s %4s %12s %9s" % ("name", "act", "Z_pk", "W_pk"))
    for r in rows_ok:
        print("%-5s %4d %12.3e %9.4f" % (r["name"], r["activity"],
                                         r["Z_pk"], r["W_pk"]))
    rank = [r["name"] for r in rows_ok].index("Sp8") + 1
    verdict = "PASS" if rank <= 2 else "FAIL"
    print("\n判别检验: Sp8 W_pk 升序名次 %d/%d -> %s" % (rank, len(rows_ok), verdict))

    import numpy as np
    acts = np.array([r["activity"] for r in rows_ok], float)
    wpk = np.array([r["W_pk"] for r in rows_ok], float)
    rho = float(np.corrcoef(np.argsort(np.argsort(wpk)),
                            np.argsort(np.argsort(acts)))[0, 1])
    print("W_pk vs activity Spearman = %+.3f" % rho)

    gf = {k: {"W_pk": round((v["pk"]["Z"] - v["nopk"]["Z"]) / v["pk"]["Z"], 4)
              if "Z" in v["pk"] and "Z" in v["nopk"] else None,
              "G_pk": v["pk"].get("G_kcal"), "G_nopk": v["nopk"].get("G_kcal")}
          for k, v in res.items() if k.startswith("gf_")}
    payload = {"engine": "NUPACK 3.2.2 pfunc/-pseudo (Dirks-Pierce 假结配分)",
               "metric": "W_pk=(Z_pk-Z_nopk)/Z_pk, 假结构象权重占比",
               "preregistered": "Sp8 W_pk 应排末 2 位(pre 口径); 2026-09-02 先定后跑",
               "rows": rows, "sp8_rank_asc": rank, "verdict": verdict,
               "spearman_Wpk_activity": round(rho, 3), "guideforge": gf,
               "boundary": "pairs/mfe 不支持 -pseudo(手册+实测崩溃), 故无 pk 配对级信息"}
    path = os.path.join(DATA, "nupack_pk_weight.json")
    json.dump(payload, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("jobs", "assemble"))
    args = ap.parse_args()
    (build_jobs if args.mode == "jobs" else assemble)()


if __name__ == "__main__":
    main()
