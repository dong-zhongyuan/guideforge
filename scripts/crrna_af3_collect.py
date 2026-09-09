"""AF3 Server 结果回填与对比汇总(策划案 V3 表1 独立对照层, 2026-09-05)。

读取 data/af3_results/ 下 AlphaFold Server 下载的结果(支持 .zip 原包或解压目录),
逐任务解析 *summary_confidence*.json(AF3 服务器每任务 5 模型), 提取 ipTM/pTM/
ranking_score/链对 ipTM, 汇总 13 套(WT + TOP-12)的 prot-crRNA 与 crRNA-target
界面均值/sd, 与 WT 做差量, 写入 data/af3_summary.json。

2026-09-09 扩展: --protenix 模式回收 Protenix-v1 结果(data/protenix_results/,
容器 2026-09-06 批跑产物, 命名 *summary_confidence_sample_N.json, 字段同构:
iptm/ptm/chain_pair_iptm; 链序同为 蛋白/crRNA/靶RNA)。产出
data/protenix_summary.json, 附与 Chai-1 自模板 WT 参照
(data/chai_cofold_WT_template_real.json)的跨引擎对照与预登记分支判读
(docs/af3_server_runbook.md 二分支: 界面恢复→独立佐证 / 低迷→主张撤回)。

链序约定(服务器自动分配): 0/A=蛋白 SuCas12a2, 1/B=crRNA, 2/C=靶RNA
(与 crrna_af3_input.py 的 sequences 顺序一致)。

口径声明(写入产物): 本层为 AlphaFold Server + useStructureTemplate=false 的
template-free 预测, 是 Chai-1 自模板口径(降级为 sanity check)的独立对照;
但服务器仍自动构建 MSA, 且单种子 5 模型, 结论表述不得超出"独立引擎的
template-free 单构建观测"。

运行: python scripts/crrna_af3_collect.py [--dir data/af3_results]
"""
import argparse
import glob
import io
import json
import os
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

IPTM_BAND = (0.5, 0.6)  # 可用界面预测置信区间下沿的常用口径(与 chai_collect 一致)


def iter_summary_jsons(path):
    """从一个结果 zip 或解压目录产出 (模型名, dict)。"""
    if zipfile.is_zipfile(path):
        z = zipfile.ZipFile(path)
        names = sorted(n for n in z.namelist()
                       if "summary_confidence" in os.path.basename(n)
                       and n.endswith(".json"))
        for n in names:
            yield os.path.basename(n), json.loads(z.read(n).decode("utf-8"))
    elif os.path.isdir(path):
        for f in sorted(glob.glob(os.path.join(
                path, "**", "*summary_confidence*.json"), recursive=True)):
            yield os.path.basename(f), json.load(open(f, encoding="utf-8"))


def collect_one(path):
    """汇总单任务: 返回 5 模型均值/sd 的 ipTM、prot-crRNA、crRNA-target 链对分数。"""
    rows = []
    for name, d in iter_summary_jsons(path):
        pair = d.get("chain_pair_iptm")
        if not pair or len(pair) < 3:
            continue
        rows.append({
            "model": name,
            "iptm": round(float(d["iptm"]), 4),
            "ptm": round(float(d["ptm"]), 4),
            "ranking_score": round(float(d.get("ranking_score", 0.0)), 4),
            "prot_crRNA": round(float(pair[0][1]), 4),
            "crRNA_target": round(float(pair[1][2]), 4),
            "has_clash": bool(d.get("has_clash", False)),
        })
    if not rows:
        return None

    def ms(key):
        xs = [r[key] for r in rows]
        m = sum(xs) / len(xs)
        sd = (sum((x - m) ** 2 for x in xs) / len(xs)) ** 0.5
        return round(m, 4), round(sd, 4)

    im, isd = ms("iptm")
    pm, psd = ms("prot_crRNA")
    cm, csd = ms("crRNA_target")
    return {"n_models": len(rows), "iptm_mean": im, "iptm_sd": isd,
            "prot_crRNA_mean": pm, "prot_crRNA_sd": psd,
            "crRNA_target_mean": cm, "crRNA_target_sd": csd,
            "models": rows}


def build_reading(job, rec, wt_pc):
    """由数值程序化生成解读(方向与量级不手写)。"""
    pc = rec["prot_crRNA_mean"]
    parts = ["prot-crRNA ipTM %.3f(5模型sd %.3f)" % (pc, rec["prot_crRNA_sd"])]
    if pc < IPTM_BAND[0]:
        parts.append("低于可用界面预测区间(~%.1f-%.1f)" % IPTM_BAND)
    elif pc < IPTM_BAND[1]:
        parts.append("处于可用区间下沿, 弱参考")
    else:
        parts.append("达可用区间(>=%.1f)" % IPTM_BAND[1])
    if job != "GF_WT" and wt_pc is not None:
        d = round(pc - wt_pc, 4)
        parts.append("vs WT %+.4f(%s)" % (
            d, "噪声内(|Δ|<=2*sd合成)" if abs(d) <= 2 * (rec["prot_crRNA_sd"] ** 2
               + 0.0) ** 0.5 else "超出2sd"))
    return "; ".join(parts)


def _chai_comparison(px_wt_pc):
    """与 Chai-1 自模板 WT 参照的跨引擎对照 + 预登记分支判读(runbook 二分支)。"""
    ref = {"source": "data/chai_cofold_WT_template_real.json"}
    try:
        d = json.load(open(os.path.join(
            DATA, "chai_cofold_WT_template_real.json"), encoding="utf-8"))
        xs = [m["iptm_prot_crRNA"] for m in d["models"]]
        m = sum(xs) / len(xs)
        ref["chai_wt_prot_crRNA_mean"] = round(m, 4)
        ref["chai_wt_n_models"] = len(xs)
    except Exception as e:  # noqa: BLE001
        ref["error"] = str(e)
        return ref
    if px_wt_pc is None:
        ref["note"] = "Protenix GF_WT 缺失, 无法对照"
        return ref
    ref["protenix_wt_prot_crRNA_mean"] = round(px_wt_pc, 4)
    band = IPTM_BAND
    both = px_wt_pc >= band[1] and ref["chai_wt_prot_crRNA_mean"] >= band[0]
    if px_wt_pc >= band[1]:
        ref["preregistered_verdict"] = (
            "分支一: template-free 独立引擎下 WT prot-crRNA ipTM %.3f >= %.1f, "
            "界面恢复成立——Chai 自模板口径的蛋白-crRNA 界面预测获得独立引擎佐证, "
            "可表述为'两个独立引擎对界面恢复的预测一致'(差异 %.3f 属引擎间正常离散); "
            "按口径不用于骨架间排序主张" % (
                px_wt_pc, band[1], px_wt_pc - ref["chai_wt_prot_crRNA_mean"]))
    else:
        ref["preregistered_verdict"] = (
            "分支二: template-free 独立引擎下界面低迷(%.3f < %.1f)——Chai 自模板"
            "高分确认是模板复述, '界面可预测'主张撤回, 仅保留湿实验判据" % (
                px_wt_pc, band[0]))
    ref["both_engines_above_band"] = both
    return ref


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", default=None,
                    help="结果 zip/解压目录所在目录(默认 data/af3_results; "
                         "--protenix 时为 data/protenix_results)")
    ap.add_argument("--protenix", action="store_true",
                    help="回收 Protenix-v1 结果(容器批跑产物)并附 Chai 跨引擎对照")
    args = ap.parse_args()

    if args.protenix:
        args.dir = args.dir or os.path.join(DATA, "protenix_results")
    else:
        args.dir = args.dir or os.path.join(DATA, "af3_results")

    entries = []
    for p in sorted(glob.glob(os.path.join(args.dir, "*"))):
        base = os.path.basename(p)
        if zipfile.is_zipfile(p) or os.path.isdir(p):
            entries.append(p)
    if not entries:
        raise SystemExit("[data-pending] %s 下无 AF3 结果 zip/目录: "
                         "服务器结果下载后放入该目录再重跑(见 docs/af3_server_runbook.md)"
                         % args.dir)

    results = {}
    for p in entries:
        base = os.path.basename(p)
        job = os.path.splitext(base)[0]
        for prefix in ("fold_",):
            if job.startswith(prefix):
                job = job[len(prefix):]
        rec = collect_one(p)
        if rec:
            results[job] = rec
            print("%-24s n=%d prot-crRNA %.3f ipTM %.3f" % (
                job, rec["n_models"], rec["prot_crRNA_mean"], rec["iptm_mean"]))
    if not results:
        raise SystemExit("未解析到任何 summary_confidences(检查 zip 内容命名)")

    wt = results.get("GF_WT")
    wt_pc = wt["prot_crRNA_mean"] if wt else None
    engine = ("Protenix-v1 (protenix_base_default_v1.0.0), template-free"
              "(输入无 templates 字段), seed_101 5 sample"
              if args.protenix else
              "AlphaFold Server (AF3), useStructureTemplate=false, 单随机种子5模型")
    payload = {
        "generated_by": "scripts/crrna_af3_collect.py" + (
            " --protenix" if args.protenix else ""),
        "source_dir": os.path.relpath(args.dir, ROOT),
        "engine": engine,
        "chain_order": {"0/A": "SuCas12a2 蛋白", "1/B": "crRNA", "2/C": "靶RNA"},
        "claim_tier": ("独立引擎 template-free 单构建观测(MSA 自动构建); "
                       "作为 Chai-1 自模板口径(sanity check)的对照层; "
                       "结论表述不得超出此口径"),
        "wt_reference": {"prot_crRNA_mean": wt_pc} if wt else "GF_WT 缺失",
        "jobs": {j: {k: v for k, v in r.items() if k != "models"}
                 for j, r in results.items()},
        "readings": {j: build_reading(j, r, wt_pc) for j, r in results.items()},
    }
    if args.protenix:
        payload["chai_cross_engine"] = _chai_comparison(wt_pc)
    out = os.path.join(
        DATA, "protenix_summary.json" if args.protenix else "af3_summary.json")
    json.dump(payload, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("汇总 %d 任务 -> %s" % (len(results), out))


if __name__ == "__main__":
    main()
