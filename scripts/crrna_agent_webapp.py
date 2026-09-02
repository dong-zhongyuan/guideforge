"""crRNA 设计智能体 Web 演示(策划案 V3 §4.4, 2026-09-02)。

单文件 Flask 应用:
  GET  /                 演示页(输入 spacer 或选 panel 靶标)
  POST /api/design       {spacer: "..."} -> 骨架分型选择(实算, 秒级)
  GET  /api/panel        -> 四靶标 panel 预计算设计(data/agent/*.design.json)
  POST /api/mut          {mut_fasta, wt_fasta} -> tilling 端到端(分钟级)

启动(服务器): nohup python scripts/crrna_agent_webapp.py --port 8899 &
访问(本机):   ssh -L 8899:localhost:8899 srv  然后打开 http://localhost:8899
边界: 结构口径分型建议 + 文献先验活性排序(Han 2025 Fig1g 尺度, n=7 外部锚点,
  非活性保证; 与 crrna_train_selector.py 同一模型同一定义)。
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from flask import Flask, jsonify, request  # noqa: E402
from sklearn.tree import DecisionTreeRegressor  # noqa: E402

import crrna_scaffold_design as core  # noqa: E402
import crrna_design_agent as agent  # noqa: E402

app = Flask(__name__)

DATA = os.path.join(ROOT, "data")
DR = core.to_rna(agent.get_scaffold("cas12a2_zeng2026"))
MODEL = agent.load_type_models(
    os.path.join(DATA, "context_typing", "typing.clusters.json"),
    os.path.join(DATA, "context_typing", "spacers.txt"),
    DR)

# ---- 模块三: 文献先验选择器(与 crrna_train_selector.py 同一配方) ----
FEATURES = ["ddG_dr", "bp_dist", "cross_nt", "p_fold", "spacer_up"]


def _train_selector():
    d = json.load(open(os.path.join(DATA, "han2025_dataset.json"),
                       encoding="utf-8"))
    X, y = [], []
    for r in d["toolbox_features"]:
        if "fig1g" not in r:
            continue
        X.append([float(r[f]) for f in FEATURES])
        y.append(float(r["fig1g"]))
    dt = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                               random_state=42)
    dt.fit(np.array(X), np.array(y))
    return dt


SELECTOR = _train_selector()


def _load_library_drs():
    """orientation_library.fasta 的 7 成员(含 WT) DR 序列"""
    drs = []
    for line in open(os.path.join(DATA, "orientation_library.fasta")):
        if line.startswith(">"):
            label = line[1:].strip().split("|")
            drs.append([label[0] if label[0] == "WT" else label[1], ""])
        elif line.strip() and drs:
            full = core.to_rna(line.strip())
            drs[-1][1] = full[:len(DR)]
    return drs


LIBRARY = _load_library_drs()


def lit_rank_scaffolds(spacer_rna):
    """对 7 成员骨架逐个算特征并给文献模型预测(Fig1g 尺度, 低=抑制强)"""
    import RNA
    wt_full_ss, _ = core.fold(DR + spacer_rna)
    _, wt_dr_mfe = core.fold(DR)
    wt_stem = core.stem_pairs_of(core.fold(DR)[0])
    out = []
    for desc, dr in LIBRARY:
        full = dr + spacer_rna
        ss, _ = core.fold(full)
        _, dr_mfe = core.fold(dr)
        cross_nt, _ = core.cross_pairs(full, len(dr))
        p_fold = core.stem_intact_prob(full, wt_stem) if wt_stem else 1.0
        _, sp_up, _ = core.pf_stats(full, len(dr), 7)
        feats = [round(dr_mfe - wt_dr_mfe, 2),
                 RNA.bp_distance(wt_full_ss, ss), cross_nt,
                 round(p_fold, 5), round(sp_up, 3)]
        pred = float(SELECTOR.predict(np.array([feats]))[0])
        out.append({"scaffold": desc, "lit_pred_fig1g": round(pred, 4),
                    "features": dict(zip(FEATURES, feats))})
    out.sort(key=lambda r: r["lit_pred_fig1g"])
    for i, r in enumerate(out, 1):
        r["rank"] = i
    return out

PANEL = {"TP53-R248Q": "tp53_r248q", "KRAS-G12C": "kras_g12c",
         "KRAS-G12D": "kras_g12d", "TP53-R273H": "tp53_r273h",
         "APC-Q1312x": "apc_q1338x"}

PAGE = """<!doctype html><html><head><meta charset="utf-8">
<title>GuideForge crRNA 设计智能体</title>
<style>
body{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em;
background:#f7f9fb;color:#20303c}
h1{font-size:1.4em} .card{background:#fff;border:1px solid #d8e2ea;border-radius:8px;
padding:1em 1.2em;margin:1em 0}
textarea,input{width:100pc;box-sizing:border-box;padding:.5em;font-family:monospace}
button{background:#0e6ba8;color:#fff;border:0;border-radius:6px;padding:.55em 1.4em;
cursor:pointer;margin-top:.5em}
pre{background:#f0f4f8;padding:.8em;border-radius:6px;overflow-x:auto;white-space:pre-wrap}
.small{color:#62788a;font-size:.85em}
</style></head><body>
<h1>GuideForge — Cas12a2 crRNA 设计智能体(演示)</h1>
<div class="small">输出 = 结构口径分型建议 + 文献先验活性排序(Han 2025 Fig1g 尺度, n=7 外部锚点, 非活性保证)。</div>
<div class="card"><b>① spacer 模式</b>(项目本体口径, 秒级)<br>
spacer (17-25nt ACGT):<br>
<input id="sp" value="GTTCATGCCGCCCATGCAGGAACT">
<button onclick="go()">设计</button><pre id="o1"></pre></div>
<div class="card"><b>② panel 模式</b>(预计算端到端案例)<br>
<select id="pt">@@OPTIONS@@</select>
<button onclick="panel()">查看</button><pre id="o2"></pre></div>
<script>
async function go(){const r=await fetch('/api/design',{method:'POST',
headers:{'Content-Type':'application/json'},
body:JSON.stringify({spacer:document.getElementById('sp').value})});
document.getElementById('o1').textContent=JSON.stringify(await r.json(),null,1);}
async function panel(){const t=document.getElementById('pt').value;
const r=await fetch('/api/panel?target='+encodeURIComponent(t));
document.getElementById('o2').textContent=JSON.stringify(await r.json(),null,1);}
</script></body></html>"""


@app.get("/")
def index():
    opts = "".join("<option>%s</option>" % k for k in PANEL)
    return PAGE.replace("@@OPTIONS@@", opts)


@app.post("/api/design")
def design():
    sp = (request.json or {}).get("spacer", "")
    sp = sp.upper().replace("U", "T")
    if not 17 <= len(sp) <= 25 or set(sp) - set("ACGT"):
        return jsonify({"error": "spacer 须 17-25nt ACGT"}), 400
    try:
        sel = agent.select_scaffold_type(core.to_rna(sp), DR, MODEL)
        rep = MODEL["type_dr"].get(sel["nearest_type"], {})
        return jsonify({"input_spacer": sp,
                        "scaffold_type": sel["nearest_type"],
                        "confidence": sel["confidence"],
                        "features": sel["features"],
                        "type_representative": rep.get("representative_desc"),
                        "scaffold_ranking_lit": lit_rank_scaffolds(
                            core.to_rna(sp)),
                        "note": "结构分型建议 + 文献先验活性排序"
                                "(Han 2025 Fig1g 尺度, n=7 锚点, 低=抑制强,"
                                "非活性保证)"})
    except Exception as e:
        return jsonify({"error": str(e)[:200]}), 500


@app.get("/api/panel")
def panel():
    key = request.args.get("target", "")
    f = os.path.join(DATA, "agent", "%s.design.json" % PANEL.get(key, ""))
    if not os.path.isfile(f):
        return jsonify({"error": "未知靶标"}), 404
    d = json.load(open(f, encoding="utf-8"))
    top_spacer = d["designs"][0]["spacer_dna"] if d.get("designs") else None
    extra = {}
    if top_spacer:
        extra = {"scaffold_ranking_lit": lit_rank_scaffolds(
            core.to_rna(top_spacer))}
    # 模块四: 细胞层预测(实测标定虚拟细胞, Scholz 2026 剂量曲线口径)
    vc_path = os.path.join(DATA, "virtual_cell_prior.json")
    if os.path.isfile(vc_path):
        vc = json.load(open(vc_path, encoding="utf-8"))
        gene = "TP53" if key.startswith("TP53") else "KRAS"
        tag = "tp53" if gene == "TP53" else "kras"
        lines = [{"cell_line": r["cell_line"], "rpkm": r["rpkm"].get(gene),
                  "pred_survival_mid_pct": r.get("surv_%s_mid" % tag)}
                 for r in vc.get("named_cell_lines", [])]
        extra["virtual_cell"] = {
            "calibration": vc.get("calibration"),
            "status": vc.get("status"),
            "named_lines": lines,
            "note": "预测存活率基于 Scholz 2026 实测剂量曲线(RNP 体系移植边界"
                    "见 status); 细胞系选系前须复核 DepMap 突变状态"}
    return jsonify({"target": key, "design": d, **extra})


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--port", type=int, default=8899)
    args = ap.parse_args()
    app.run(host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()
