"""crRNA 设计智能体 Web 演示(策划案 V3 §4.4, 2026-09-02)。

单文件 Flask 应用:
  GET  /                 演示页(输入 spacer 或选 panel 靶标)
  POST /api/design       {spacer: "..."} -> 骨架分型选择(实算, 秒级)
  GET  /api/panel        -> 四靶标 panel 预计算设计(data/agent/*.design.json)
  POST /api/mut          {mut_fasta, wt_fasta} -> tilling 端到端(分钟级)

启动(服务器): nohup python scripts/crrna_agent_webapp.py --port 8899 &
访问(本机):   ssh -L 8899:localhost:8899 srv  然后打开 http://localhost:8899
边界: 输出为结构口径排序与分型建议, 不含活性预测(README 路径A校准)。
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from flask import Flask, jsonify, request  # noqa: E402

import crrna_scaffold_design as core  # noqa: E402
import crrna_design_agent as agent  # noqa: E402

app = Flask(__name__)

DATA = os.path.join(ROOT, "data")
DR = core.to_rna(agent.get_scaffold("cas12a2_zeng2026"))
MODEL = agent.load_type_models(
    os.path.join(DATA, "context_typing", "typing.clusters.json"),
    os.path.join(DATA, "context_typing", "spacers.txt"),
    DR)

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
<div class="small">输出为结构口径的骨架分型建议与排序, 不含活性预测(见仓库 README 校准声明)。</div>
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
                        "note": "结构口径分型建议; 非活性预测"})
    except Exception as e:
        return jsonify({"error": str(e)[:200]}), 500


@app.get("/api/panel")
def panel():
    key = request.args.get("target", "")
    f = os.path.join(DATA, "agent", "%s.design.json" % PANEL.get(key, ""))
    if not os.path.isfile(f):
        return jsonify({"error": "未知靶标"}), 404
    d = json.load(open(f, encoding="utf-8"))
    return jsonify({"target": key, "design": d})


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--port", type=int, default=8899)
    args = ap.parse_args()
    app.run(host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()
