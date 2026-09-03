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

PAGE = """<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GuideForge — Cas12a2 crRNA 设计智能体</title>
<style>
:root{
--canvas:#FAFAF8;--ink:#1C1B19;--muted:#7A7468;--line:#E7E4DD;
--surface:#FFFFFF;--accent:#0F5F52;--accent-soft:#E3EFEB;
--warn-bg:#FBF3DB;--warn-tx:#8A5D00;--err-bg:#FDEBEC;--err-tx:#9F2F2D;
--mono:"JetBrains Mono","SF Mono",Consolas,"Courier New",monospace;
--sans:"Helvetica Neue","Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;
}
*{box-sizing:border-box}
body{margin:0;background:var(--canvas);color:var(--ink);font-family:var(--sans);
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:56px 24px 40px}
.eyebrow{font-family:var(--mono);font-size:12px;letter-spacing:.12em;
text-transform:uppercase;color:var(--accent);margin-bottom:10px}
h1{font-size:30px;font-weight:600;letter-spacing:-.02em;line-height:1.2;
margin:0 0 8px;text-wrap:balance}
.sub{color:var(--muted);margin:0 0 14px;max-width:65ch}
.disclaimer{display:inline-block;background:var(--warn-bg);color:var(--warn-tx);
font-size:12.5px;border-radius:6px;padding:6px 12px;margin-bottom:36px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:start}
@media(max-width:820px){.grid{grid-template-columns:1fr}}
.card{background:var(--surface);border:1px solid var(--line);border-radius:12px;
padding:24px 28px 28px}
.card-head{display:flex;align-items:baseline;gap:10px;margin-bottom:18px}
.card-head h2{font-size:17px;font-weight:600;margin:0;flex:1}
.tag{font-family:var(--mono);font-size:11px;letter-spacing:.08em;color:var(--accent);
background:var(--accent-soft);border-radius:4px;padding:2px 7px}
.meta{font-size:12px;color:var(--muted)}
label{display:block;font-size:13px;color:var(--muted);margin-bottom:6px}
input,select{width:100%;box-sizing:border-box;padding:9px 11px;font-family:var(--mono);
font-size:13.5px;border:1px solid var(--line);border-radius:6px;background:#FDFDFC;
color:var(--ink);outline:none;transition:border-color .2s,box-shadow .2s}
input:focus,select:focus{border-color:var(--accent);
box-shadow:0 0 0 3px rgba(15,95,82,.12)}
.row{display:flex;align-items:center;gap:12px;margin-top:14px}
button{background:var(--ink);color:#fff;border:0;border-radius:6px;
padding:9px 20px;font-size:13.5px;font-family:var(--sans);cursor:pointer;
transition:background .2s,transform .1s}
button:hover{background:#3A3833}
button:active{transform:scale(.98)}
button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.hint{font-size:12.5px;color:var(--muted)}
.result{margin-top:20px}
.kv{display:grid;grid-template-columns:auto 1fr;gap:4px 16px;font-size:13.5px;
margin-bottom:16px}
.kv dt{color:var(--muted)}.kv dd{margin:0}
.seq{font-family:var(--mono);font-size:13px;word-break:break-all}
.badge{display:inline-block;font-family:var(--mono);font-size:12px;
background:var(--accent-soft);color:var(--accent);border-radius:4px;
padding:2px 8px}
table{width:100%;border-collapse:collapse;font-size:13px;margin:6px 0 4px}
th{font-family:var(--mono);font-size:11px;letter-spacing:.06em;color:var(--muted);
text-align:left;font-weight:500;padding:6px 8px;border-bottom:1px solid var(--line)}
td{padding:7px 8px;border-bottom:1px solid var(--line);
font-variant-numeric:tabular-nums}
tr:last-child td{border-bottom:0}
td.num,th.num{text-align:right;font-family:var(--mono);font-size:12.5px}
.bar{height:4px;border-radius:2px;background:var(--accent);margin-top:4px}
.bar-track{width:100%;height:4px;border-radius:2px;background:var(--line);
margin-top:4px}
.note{font-size:12px;color:var(--muted);margin-top:12px}
details{margin-top:14px}
summary{font-size:12.5px;color:var(--muted);cursor:pointer}
pre{background:#F4F3F0;border:1px solid var(--line);border-radius:6px;
padding:12px;font-size:12px;overflow-x:auto;white-space:pre-wrap;
font-family:var(--mono)}
.err{background:var(--err-bg);color:var(--err-tx);border-radius:6px;
padding:10px 14px;font-size:13px}
.sk{border-radius:6px;background:linear-gradient(90deg,#F0EEE9 25%,#FAFAF8 50%,#F0EEE9 75%);
background-size:200% 100%;animation:sh 1.2s infinite;height:14px;margin:10px 0}
@keyframes sh{to{background-position:-200% 0}}
footer{margin-top:48px;padding-top:18px;border-top:1px solid var(--line);
font-size:12px;color:var(--muted)}
@media(prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style></head><body>
<div class="wrap">
<header>
<div class="eyebrow">GuideForge · crRNA Design Agent</div>
<h1>Cas12a2 crRNA 设计智能体</h1>
<p class="sub">输入 spacer 序列或选择 panel 靶标，输出骨架分型建议与文献先验活性排序。</p>
<div class="disclaimer">口径说明：结构口径分型建议 + 文献先验活性排序（Han 2025 Fig1g 尺度，n=7 外部锚点，非活性保证）</div>
</header>
<main class="grid">
<section class="card">
<div class="card-head"><span class="tag">模式一</span><h2>spacer 设计</h2><span class="meta">实算 · 秒级</span></div>
<label for="sp">spacer 序列（17–25 nt，ACGT）</label>
<input id="sp" value="GTTCATGCCGCCCATGCAGGAACT" spellcheck="false" autocomplete="off">
<div class="row"><button onclick="go()">运行设计</button><span class="hint" id="hint1"></span></div>
<div class="result" id="r1"></div>
</section>
<section class="card">
<div class="card-head"><span class="tag">模式二</span><h2>panel 靶标</h2><span class="meta">预计算 · 端到端</span></div>
<label for="pt">靶标</label>
<select id="pt">@@OPTIONS@@</select>
<div class="row"><button onclick="panel()">查看设计</button><span class="hint" id="h2"></span></div>
<div class="result" id="r2"></div>
</section>
</main>
<footer>GuideForge · Cas12a2 crRNA 骨架分型 · 演示环境</footer>
</div>
<script>
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
function loading(el,h){el.innerHTML='<div class="sk"></div><div class="sk" style="width:70%"></div><div class="sk" style="width:45%"></div>';if(h)h.textContent='计算中…';}
function fail(el,h,msg){el.innerHTML='<div class="err">'+esc(msg)+'</div>';if(h)h.textContent='';}
function rankTable(rows){
if(!rows||!rows.length)return'';
const vals=rows.map(r=>r.lit_pred_fig1g),lo=Math.min.apply(null,vals),hi=Math.max.apply(null,vals);
let h='<table><thead><tr><th>#</th><th>骨架</th><th class="num">Fig1g 预测</th><th class="num">ddG</th><th class="num">bp_dist</th><th class="num">cross_nt</th><th class="num">p_fold</th><th class="num">spacer_up</th></tr></thead><tbody>';
for(const r of rows){const f=r.features||{};const w=hi>lo?((r.lit_pred_fig1g-lo)/(hi-lo)*100):50;
h+='<tr><td class="num">'+r.rank+'</td><td>'+esc(r.scaffold)+'</td>'+
'<td class="num">'+r.lit_pred_fig1g.toFixed(4)+'<div class="bar-track"><div class="bar" style="width:'+w.toFixed(1)+'%"></div></div></td>'+
'<td class="num">'+(f.ddG_dr!=null?f.ddG_dr:'')+'</td><td class="num">'+(f.bp_dist!=null?f.bp_dist:'')+'</td>'+
'<td class="num">'+(f.cross_nt!=null?f.cross_nt:'')+'</td><td class="num">'+(f.p_fold!=null?f.p_fold:'')+'</td>'+
'<td class="num">'+(f.spacer_up!=null?f.spacer_up:'')+'</td></tr>';}
return h+'</tbody></table><div class="note">Fig1g 尺度：数值越低 = 预测抑制越强（Han 2025 文献先验，n=7 锚点）。</div>';}
function raw(d){return '<details><summary>查看原始 JSON</summary><pre>'+esc(JSON.stringify(d,null,1))+'</pre></details>';}
async function go(){loading($('r1'),$('hint1'));
const sp=$('sp').value.trim();
let r,d;
try{r=await fetch('/api/design',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({spacer:sp})});d=await r.json();}
catch(e){fail($('r1'),$('hint1'),'网络请求失败，请检查服务是否在线。');return;}
if(!r.ok||d.error){fail($('r1'),$('hint1'),d.error||('请求失败（HTTP '+r.status+'）'));return;}
$('hint1').textContent='';
let h='<dl class="kv">'+
'<dt>输入 spacer</dt><dd class="seq">'+esc(d.input_spacer)+'</dd>'+
'<dt>骨架分型</dt><dd><span class="badge">'+esc(d.scaffold_type)+'</span> &nbsp;置信度 '+(d.confidence!=null?d.confidence:'—')+'</dd>'+
(d.type_representative?'<dt>型代表</dt><dd class="seq">'+esc(d.type_representative)+'</dd>':'')+
'</dl>';
h+='<label>骨架文献先验活性排序</label>'+rankTable(d.scaffold_ranking_lit);
if(d.note)h+='<div class="note">'+esc(d.note)+'</div>';
h+=raw(d);$('r1').innerHTML=h;}
async function panel(){loading($('r2'),null);
const t=$('pt').value;let r,d;
try{r=await fetch('/api/panel?target='+encodeURIComponent(t));d=await r.json();}
catch(e){fail($('r2'),null,'网络请求失败，请检查服务是否在线。');return;}
if(!r.ok||d.error){fail($('r2'),null,d.error||('请求失败（HTTP '+r.status+'）'));return;}
let h='<dl class="kv"><dt>靶标</dt><dd><span class="badge">'+esc(d.target)+'</span></dd>';
const dz=(d.design&&d.design.designs)||[];
if(dz.length){h+='<dt>候选 spacer</dt><dd>'+dz.length+' 条，首位：</dd>';h+='<dt></dt><dd class="seq">'+esc(dz[0].spacer_dna)+'</dd>';}
h+='</dl>';
if(d.scaffold_ranking_lit){h+='<label>骨架文献先验活性排序</label>'+rankTable(d.scaffold_ranking_lit);}
const vc=d.virtual_cell;
if(vc){const cal=vc.calibration||{};
const calTxt=(cal.ec50_fpkm!=null)?('EC50 '+cal.ec50_fpkm+' FPKM（95% CI '+cal.ec50_ci95[0]+'–'+cal.ec50_ci95[1]+'），Hill '+cal.hill+'，R²='+cal.r2+'，n='+cal.n_points+'，'+esc(cal.source||'')):esc(vc.calibration||'—');
h+='<label style="margin-top:14px">细胞层预测（虚拟细胞）</label><dl class="kv">'+
'<dt>标定</dt><dd>'+calTxt+'</dd><dt>状态</dt><dd>'+esc(vc.status||'—')+'</dd></dl>';
const nl=vc.named_lines||[];
if(nl.length){h+='<table><thead><tr><th>细胞系</th><th class="num">RPKM</th><th class="num">预测存活 %</th></tr></thead><tbody>';
for(const x of nl){h+='<tr><td>'+esc(x.cell_line)+'</td><td class="num">'+(x.rpkm!=null?x.rpkm:'—')+'</td><td class="num">'+(x.pred_survival_mid_pct!=null?x.pred_survival_mid_pct:'—')+'</td></tr>';}
h+='</tbody></table>';}
if(vc.note)h+='<div class="note">'+esc(vc.note)+'</div>';}
h+=raw(d);$('r2').innerHTML=h;}
window.addEventListener('DOMContentLoaded',()=>{
if(location.hash==='#demo-spacer')go();
if(location.hash==='#demo-panel')panel();});
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
