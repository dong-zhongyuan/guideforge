"""crRNA 设计智能体 Web 演示(策划案 V3 §4.4, 2026-09-02)。

单文件 Flask 应用:
  GET  /                 演示页(输入 spacer 或选 panel 靶标)
  POST /api/design       {spacer: "...", target: 可选panel键} -> 骨架分型选择(实算, 秒级);
                         提供 target 时附靶RNA丰度档位(模块一特征维)+门控判读
                         + 界面特征差量(结构层选型特征维附加层, Chai-1 自模板
                           口径, data/chai_interface_deltas.json, 与离线 agent
                           同源; 可用时文献先验排序行附 Δprot-crRNA 列);
                         响应均携带 degradation 块(阴性退化判据状态, 见下)
  GET  /api/panel        -> 四靶标 panel 预计算设计(data/agent/*.design.json)
                         + 丰度档位(CCLE 18q3 实测, data/target_abundance_tiers.json)
                         + 突变类型判定(最长 ORF 法, 模块一特征维, 与离线
                           crrna_mutation_typing 同源; 证据取自 design.json
                           mutation.typing 块)
                         + 界面特征差量(任务⑫, 与离线同源; 矩阵未覆盖的靶标
                           如实返回 available=False + 缺失原因)
                         + degradation 块(与离线智能体同源的判据状态)
  POST /api/mut          {mut_fasta, wt_fasta} -> tilling 端到端(分钟级)

阴性退化模式(策划案 V3 §5.3 末段/表5 风险1, 2026-09-06): 启动时按预登记
go/no-go 判据(docs/preregistration.md §A)从分型产物重算"骨架最优解与靶标
特征是否显著相关"(agent.interaction_assessment); 判据未通过(阴性)时
/api/design 推荐从分型退化为通用型单骨架(agent.universal_scaffold 全场
最优, 与离线智能体同源同口径), 页面以醒目横幅如实展示退化状态、判据数值
与"科学结论依然完整可交付"的交付口径, 不掩饰。--clusters/--spacers 可指向
其它分型产物复算判据(阴性场景演示)。

启动(服务器): nohup python scripts/crrna_agent_webapp.py --port 8899 &
访问(本机):   ssh -L 8899:localhost:8899 srv  然后打开 http://localhost:8899
边界: 结构口径分型建议 + 文献先验活性排序(Han 2025 Fig1g 尺度, n=7 外部锚点,
  非活性保证; 特征表与 crrna_train_selector.py 同一定义, webapp 侧为 n=7
  轻量重训, 正式先验排序以 data/selector_model.json 为准)。丰度分档阈值显式
  声明于 scripts/crrna_target_abundance.py(Scholz 2026 EC50 95%CI 锚点)。
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

from flask import Flask, jsonify, request  # noqa: E402
from sklearn.tree import DecisionTreeRegressor  # noqa: E402

import crrna_scaffold_design as core  # noqa: E402
import crrna_design_agent as agent  # noqa: E402
import crrna_target_abundance as tabund  # noqa: E402
import crrna_mutation_typing as mtyping  # noqa: E402
import crrna_chai_interface_features as chai_if  # noqa: E402
from crrna_train_selector import (FEATURES, TARGET_FEATURES,  # noqa: E402
                                  INTERFACE_DELTA_FEATURES)

app = Flask(__name__)

DATA = os.path.join(ROOT, "data")
DR = core.to_rna(agent.get_scaffold("cas12a2_zeng2026"))
CLUSTERS_JSON = os.path.join(DATA, "context_typing", "typing.clusters.json")
SPACERS_TXT = os.path.join(DATA, "context_typing", "spacers.txt")
MODEL = agent.load_type_models(CLUSTERS_JSON, SPACERS_TXT, DR)

# 阴性退化判据(策划案 V3 §5.3 末段/表5 风险1, 2026-09-06): 启动时按预登记
# §A go/no-go 判据从分型产物重算一次, 全站共用; mode=degraded_universal 时
# 推荐从分型退化为通用型单骨架(agent.universal_scaffold, 与离线智能体同源
# 同口径), 页面如实展示退化状态与量化证据, 不掩饰。
DEGRADATION = agent.interaction_assessment(CLUSTERS_JSON)

# 模块一靶标解析特征向量的丰度维(panel 靶标 -> 丰度档位 + 门控判读)
ABUNDANCE = tabund.load_abundance(os.path.join(DATA, "target_abundance_tiers.json"))

# 结构层界面特征差量(任务⑫, 选型特征维附加层): 与离线 agent
# (crrna_design_agent.target_context)读同一 JSON、走同一 target_block 函数,
# 两端不得各自复制口径; 矩阵无该靶标列的 panel 键(如 APC 待增量折叠回填)
# 如实返回 available=False + 缺失原因
CHAI_IF = chai_if.load_deltas(os.path.join(DATA, "chai_interface_deltas.json"))

# ---- 模块三: 文献先验选择器(FEATURES/TARGET_FEATURES 经 import 与
#      crrna_train_selector.py 共用同一定义; webapp 侧为 n=7 工具箱锚点
#      轻量重训, 仅作演示, 正式先验排序以 data/selector_model.json 为准) ----


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
                    "features": dict(zip(FEATURES, feats)), "dr": dr})
    out.sort(key=lambda r: r["lit_pred_fig1g"])
    for i, r in enumerate(out, 1):
        r["rank"] = i
    return out

PANEL = {"TP53-R248Q": "tp53_r248q", "KRAS-G12D": "kras_g12d",
         "TP53-R273H": "tp53_r273h", "APC-Q1328x": "apc_q1328x"}

PAGE = """<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="description" content="GuideForge：Cas12a2 crRNA 骨架分型设计智能体。输入 spacer 或选择 panel 靶标，输出骨架分型建议与文献先验活性排序。">
<link rel="icon" href="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'><path d='M11 27V13a5 5 0 0 1 10 0v14' stroke='%230F5F52' stroke-width='2.6' stroke-linecap='round' fill='none'/><path d='M11 21.5h10M11 16.5h10' stroke='%230F5F52' stroke-width='2' stroke-linecap='round' opacity='.55'/></svg>">
<script src="/static/gsap.min.js"></script>
<script src="/static/MotionPathPlugin.min.js"></script>
<title>GuideForge · Cas12a2 crRNA 设计智能体</title>
<style>
:root{
--canvas:#F6F4EE;--ink:#1C1B19;--muted:#6E685C;--line:#E2DED4;
--surface:#FFFFFF;--surface-2:#FBF9F5;--accent:#0F5F52;--accent-deep:#0B4A40;
--accent-soft:#E3EEEA;
--warn-bg:#FBF3DB;--warn-tx:#7A5200;--warn-edge:#D9A93B;
--err-bg:#FDEBEC;--err-tx:#9F2F2D;
--mono:"JetBrains Mono","SF Mono",Consolas,"Courier New",monospace;
--sans:"Helvetica Neue","Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;
--serif:Georgia,"Times New Roman","Songti SC","STSong",SimSun,serif;
}
*{box-sizing:border-box}
[hidden]{display:none!important}
::selection{background:var(--accent-soft);color:var(--accent-deep)}
::-webkit-scrollbar{width:10px}
::-webkit-scrollbar-track{background:transparent}
::-webkit-scrollbar-thumb{background:#CFC8B9;border-radius:5px;
border:2px solid var(--canvas)}
.frame{position:fixed;inset:14px;border:1px solid var(--line);
pointer-events:none;z-index:4}
html{scroll-behavior:smooth}
body{margin:0;background:var(--canvas);color:var(--ink);font-family:var(--sans);
line-height:1.65;-webkit-font-smoothing:antialiased}
body::before{content:"";position:fixed;inset:0;pointer-events:none;
background:
radial-gradient(1100px 480px at 88% -12%,rgba(15,95,82,.055),transparent 62%),
radial-gradient(880px 420px at -12% 112%,rgba(28,27,25,.045),transparent 58%);}
body::after{content:"";position:fixed;inset:0;pointer-events:none;opacity:.5;
background-image:url("data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' width='160' height='160'><filter id='n'><feTurbulence type='fractalNoise' baseFrequency='0.9' numOctaves='2'/><feColorMatrix type='saturate' values='0'/><feComponentTransfer><feFuncA type='linear' slope='0.035'/></feComponentTransfer></filter><rect width='100%25' height='100%25' filter='url(%23n)'/></svg>");}
.wrap{position:relative;max-width:1120px;margin:0 auto;padding:60px 28px 48px}
.skip{position:absolute;left:-9999px;top:0}
.skip:focus{left:16px;top:12px;z-index:10;background:var(--ink);color:#fff;
padding:8px 14px;border-radius:6px;font-size:13px;text-decoration:none}
header{margin-bottom:42px}
.brand{display:flex;align-items:center;gap:11px;margin-bottom:24px}
.brand svg{width:30px;height:30px;flex:none}
.wordmark{font-size:15.5px;font-weight:600;letter-spacing:-.01em}
.wordmark em{font-style:normal;font-weight:450;color:var(--muted)}
h1{font-family:var(--serif);font-size:clamp(28px,3.6vw,38px);font-weight:600;
letter-spacing:-.015em;line-height:1.18;margin:0 0 12px;text-wrap:balance}
.sub{color:var(--muted);margin:0 0 20px;max-width:60ch;font-size:15px}
.disclaimer{display:block;max-width:660px;background:var(--warn-bg);
color:var(--warn-tx);font-size:12.5px;line-height:1.55;
border-left:3px solid var(--warn-edge);border-radius:0 8px 8px 0;
padding:9px 14px}
.layers{display:grid;grid-template-columns:repeat(4,1fr);margin-top:38px;
border-top:1px solid var(--line);padding-top:20px}
.layer{padding-right:18px}
.layer+.layer{border-left:1px solid var(--line);padding-left:18px}
.layer b{display:block;font-size:13.5px;font-weight:600;letter-spacing:-.01em;
margin-bottom:3px}
.layer span{font-size:12px;color:var(--muted);line-height:1.5}
@media(max-width:860px){.layers{grid-template-columns:1fr 1fr;gap:18px 0}
.layer:nth-child(3){border-left:0;padding-left:0}}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:start}
@media(max-width:860px){.grid{grid-template-columns:1fr}}
.card{background:var(--surface);border:1px solid var(--line);border-radius:14px;
padding:26px 28px 30px;
box-shadow:0 1px 0 rgba(255,255,255,.7) inset,0 18px 40px -24px rgba(28,27,25,.22);
transition:border-color .25s;animation:rise .5s ease both}
.card:nth-child(2){animation-delay:.09s}
@keyframes rise{from{opacity:0;transform:translateY(10px)}}
.card:hover{border-color:#CBC5B7}
.card-head{display:flex;align-items:baseline;gap:10px;margin-bottom:20px;
padding-bottom:14px;border-bottom:1px solid var(--line)}
.card-head h2{font-size:17px;font-weight:600;letter-spacing:-.01em;margin:0;flex:1}
.tag{font-family:var(--mono);font-size:11px;letter-spacing:.06em;
color:var(--accent-deep);background:var(--accent-soft);border-radius:4px;
padding:2px 7px}
.meta{font-size:12px;color:var(--muted)}
label{display:block;font-size:13px;font-weight:500;color:var(--muted);
margin-bottom:7px}
input,select{width:100%;box-sizing:border-box;padding:10px 12px;
font-family:var(--mono);font-size:13.5px;border:1px solid var(--line);
border-radius:8px;background:var(--surface-2);color:var(--ink);outline:none;
transition:border-color .2s,box-shadow .2s,background .2s}
input:focus,select:focus{border-color:var(--accent);background:#fff;
box-shadow:0 0 0 3px rgba(15,95,82,.14)}
select{appearance:none;-webkit-appearance:none;padding-right:34px;
background-image:url("data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' width='10' height='6'><path d='M1 1l4 4 4-4' stroke='%236E685C' stroke-width='1.5' fill='none' stroke-linecap='round'/></svg>");
background-repeat:no-repeat;background-position:right 13px center}
input.invalid{border-color:var(--err-tx);
box-shadow:0 0 0 3px rgba(159,47,45,.12)}
.row{display:flex;align-items:center;gap:12px;margin-top:16px}
button{background:var(--ink);color:#fff;border:0;border-radius:8px;
padding:10px 22px;font-size:13.5px;font-weight:600;font-family:var(--sans);
cursor:pointer;transition:background .2s,transform .1s;
box-shadow:0 1px 0 rgba(255,255,255,.09) inset,0 8px 18px -10px rgba(28,27,25,.4)}
button:hover{background:var(--accent)}
button:active{transform:scale(.98)}
button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.hint{font-size:12.5px;color:var(--muted)}
.result{margin-top:22px}
.result:empty::before{content:"运行后在此显示设计结果";display:block;
border:1px dashed var(--line);border-radius:8px;padding:20px 18px;
font-size:12.5px;color:var(--muted);text-align:center}
.result>*{animation:rise .35s ease both}
.kv{display:grid;grid-template-columns:auto 1fr;gap:0 16px;font-size:13.5px;
margin-bottom:18px;border-top:1px solid var(--line)}
.kv dt,.kv dd{padding:7px 0;border-bottom:1px solid var(--line)}
.kv dt{color:var(--muted)}.kv dd{margin:0}
.seq{font-family:var(--mono);font-size:13px;letter-spacing:.02em;
word-break:break-all}
.seq .dr{color:var(--muted)}
.seq .sp{color:var(--accent-deep);font-weight:600;background:var(--accent-soft);
border-radius:3px;padding:0 2px}
.badge{display:inline-block;font-family:var(--mono);font-size:12px;
background:var(--accent-soft);color:var(--accent-deep);border-radius:4px;
border:1px solid rgba(15,95,82,.22);padding:2px 8px}
table{width:100%;border-collapse:collapse;font-size:13px;margin:8px 0 6px}
th{font-family:var(--mono);font-size:11px;letter-spacing:.05em;color:var(--muted);
text-align:left;font-weight:500;padding:6px 8px;white-space:nowrap;
border-bottom:1.5px solid #D7D2C5}
td{padding:7px 8px;border-bottom:1px solid var(--line);
font-variant-numeric:tabular-nums}
tbody tr{transition:background .15s}
tbody tr:hover{background:var(--surface-2)}
tbody tr:first-child td{font-weight:600}
tbody tr:first-child td.num:first-child{color:var(--accent)}
tr:last-child td{border-bottom:0}
td.num,th.num{text-align:right;font-family:var(--mono);font-size:12.5px}
.bar{height:3px;border-radius:2px;background:var(--accent);margin-top:4px;
max-width:110px;transition:width .6s ease}
.note{font-size:12px;color:var(--muted);margin-top:12px}
details{margin-top:16px}
summary{font-size:12.5px;color:var(--muted);cursor:pointer;user-select:none}
summary::before{content:"›";display:inline-block;margin-right:6px;
transition:transform .2s}
details[open] summary::before{transform:rotate(90deg)}
summary:hover{color:var(--ink)}
summary:focus-visible{outline:2px solid var(--accent);outline-offset:2px;
border-radius:4px}
pre{background:#F3F1EB;border:1px solid var(--line);border-radius:8px;
padding:12px 14px;font-size:12px;overflow-x:auto;white-space:pre-wrap;
font-family:var(--mono)}
.err{background:var(--err-bg);color:var(--err-tx);border-radius:8px;
padding:10px 14px;font-size:13px}
.sk{border-radius:6px;background:linear-gradient(90deg,#EDE9E0 25%,#F8F5EF 50%,#EDE9E0 75%);
background-size:200% 100%;animation:sh 1.2s infinite;height:14px;margin:10px 0}
@keyframes sh{to{background-position:-200% 0}}
footer{margin-top:56px;padding-top:20px;border-top:1px solid var(--line);
font-size:12px;color:var(--muted);display:flex;justify-content:space-between;
gap:12px 24px;flex-wrap:wrap}
.view{position:relative;z-index:1}
#entry{min-height:100dvh;display:flex;overflow:hidden}
.splash{position:absolute;inset:0;z-index:2;display:flex;flex-direction:column;
align-items:center;justify-content:center;gap:24px;cursor:pointer;
outline:none;background:var(--canvas);border-radius:12px}
.splash:focus-visible{outline:2px solid var(--accent);outline-offset:-6px}
.scene{width:min(92vw,1060px)}
.cellpng,.celldead{user-select:none}
.bead{fill:url(#gTeal)}
.bead-r{fill:url(#gRed)}
.bp{stroke:rgba(20,40,35,.18);stroke-width:.6}
.bp0{fill:#E85454}.bp1{fill:#4D7DE8}.bp2{fill:#4FBB5A}.bp3{fill:#E8C84D}
.bp4{fill:#B44FD6}
.bp-mut{fill:#F5DE4C;stroke:#B8931B;stroke-width:.9}
.splash-cap{text-align:center}
.splash-cap b{display:block;font-family:var(--serif);font-size:21px;
font-weight:600;letter-spacing:.01em}
.splash-cap span{font-size:12.5px;color:var(--muted)}
body.gsap .entry-inner>.brand,body.gsap .entry-mid,body.gsap .masthead,
body.gsap .entry-art,body.gsap #app .wrap,body.gsap #app .card,
body.gsap .result>*{animation:none}
.display .line{display:block;overflow:hidden;padding-bottom:.08em}
.display .ch{display:inline-block;will-change:transform}
.entry-inner{position:relative;z-index:1;width:100%;max-width:1120px;
margin:0 auto;padding:44px 28px 40px;display:flex;flex-direction:column;flex:1}
.entry-inner>.brand{animation:rise .6s ease both}
.entry-mid{margin:auto 0;padding:7vh 0;animation:rise .6s ease .1s both}
.masthead{animation:rise .6s ease .22s both}
.entry-art{position:absolute;right:0;top:44%;transform:translateY(-50%);
width:min(33vw,400px);max-height:58vh;color:var(--accent);pointer-events:none;
animation:fadein .9s ease .22s both}
@keyframes fadein{from{opacity:0}}
.display{font-family:var(--serif);font-size:clamp(40px,5.8vw,68px);font-weight:600;
letter-spacing:-.015em;line-height:1.12;margin:0 0 10px;max-width:11em;
text-wrap:balance}
.display .cjk{font-size:.84em;letter-spacing:.005em}
body.gsap .display .line:last-child{font-size:.84em}
.sub-en{font-family:var(--serif);font-style:italic;font-size:15px;
color:var(--muted);letter-spacing:.015em;line-height:1.4;margin:0 0 26px}
.lede{color:var(--muted);font-size:16px;max-width:52ch;margin:0 0 32px}
.enter{font-size:14.5px;padding:12px 32px}
.masthead{display:grid;grid-template-columns:repeat(3,1fr);gap:24px;
border-top:1px solid var(--line);padding-top:22px}
.stat b{display:block;font-family:var(--serif);font-size:30px;font-weight:600;
letter-spacing:-.01em;color:var(--ink)}
.stat span{display:block;font-size:12.5px;color:var(--muted);line-height:1.5;
margin-top:3px}
.stat+.stat{border-left:1px solid var(--line);padding-left:24px}
.back{display:inline-block;font-size:12.5px;color:var(--muted);
text-decoration:none;padding:6px 14px;border:1px solid var(--line);
border-radius:8px;background:var(--surface);transition:border-color .2s,color .2s}
.back:hover{color:var(--accent);border-color:var(--accent)}
.back:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.topbar{display:flex;align-items:center;justify-content:space-between;
margin-bottom:44px}
.topbar .brand{margin-bottom:0}
.topbar .brand svg{width:24px;height:24px}
.topbar .wordmark{font-size:14px}
#app .wrap{animation:rise .45s ease both}
@media(max-width:860px){
.frame{inset:8px}
.masthead{grid-template-columns:1fr;gap:14px}
.stat+.stat{border-left:0;padding-left:0}
.entry-art{right:-30%;opacity:.14}
.entry-mid{padding:5vh 0}}
body.dark{
--canvas:#151310;--ink:#ECE8DE;--muted:#98927F;--line:#332F26;
--surface:#1E1B16;--surface-2:#25221A;--accent:#3DA88B;--accent-deep:#8FD8BF;
--accent-soft:#1B332A;
--warn-bg:#33290F;--warn-tx:#E0C47C;--warn-edge:#7A6420;
--err-bg:#3A1B1A;--err-tx:#DE8E8B;
}
body.dark::before{background:
radial-gradient(1100px 480px at 88% -12%,rgba(61,168,139,.12),transparent 62%),
radial-gradient(880px 420px at -12% 112%,rgba(236,232,222,.045),transparent 58%);}
body.dark .card{box-shadow:0 1px 0 rgba(255,255,255,.05) inset,
0 18px 40px -24px rgba(0,0,0,.55)}
body.dark .card:hover{border-color:#4A4436}
body.dark button{background:#E9E4D8;color:#171410;
box-shadow:0 1px 0 rgba(255,255,255,.25) inset,0 8px 18px -10px rgba(0,0,0,.6)}
body.dark button:hover{background:var(--accent);color:#0C0F0D}
body.dark th{border-bottom-color:#3B362B}
body.dark pre{background:#211E17}
body.dark .sk{background:linear-gradient(90deg,#26221B 25%,#2F2A20 50%,#26221B 75%);
background-size:200% 100%}
body.dark .badge{border-color:rgba(61,168,139,.4)}
body.dark ::-webkit-scrollbar-thumb{background:#3C372C;border-color:var(--canvas)}
body.dark .brand svg path{stroke:#3DA88B}
body.dark .entry-art>path{stroke:#4E4A40}
body.dark .entry-art g[opacity] path{stroke:#CFC9BB}
body.dark .entry-art g[fill="#FFFFFF"] circle{fill:#24211A;stroke:#CFC9BB}
body.dark .entry-art g[fill="#1C1B19"] text{fill:#ECE8DE}
body.dark .entry-art g[fill="#E3EEEA"] circle{fill:#0E332A;stroke:#3DA88B}
body.dark .entry-art g[fill="#0B4A40"] text{fill:#9BDEC9}
body.dark .entry-art>text{fill:#98927F}
.themetog{position:fixed;top:26px;right:30px;z-index:6;background:transparent;
color:var(--muted);border:1px solid var(--line);border-radius:8px;
padding:6px 14px;font-size:12px;font-weight:500;box-shadow:none}
.themetog:hover{background:transparent;color:var(--accent);
border-color:var(--accent)}
.themetog:active{transform:scale(.98)}
@media(prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style></head><body>
<div class="frame" aria-hidden="true"></div>
<button class="themetog" id="themetog" type="button">深色演示</button>
<div class="view" id="entry">
<div class="splash" id="splash" role="button" tabindex="0" aria-label="Cas12a2 靶向杀伤机制动画，点击跳过">
<svg class="scene" viewBox="0 0 1200 675" role="img" aria-label="RNP 复合物切割突变转录本，靶细胞死亡示意">
<defs>
<radialGradient id="gTeal" cx="35%" cy="28%" r="80%">
<stop offset="0%" stop-color="#CDEFE6"/><stop offset="48%" stop-color="#5EBBA6"/><stop offset="100%" stop-color="#2C7A68"/>
</radialGradient>
<radialGradient id="gRed" cx="35%" cy="28%" r="80%">
<stop offset="0%" stop-color="#F6B1A6"/><stop offset="48%" stop-color="#D95A4B"/><stop offset="100%" stop-color="#9E2F24"/>
</radialGradient>
</defs>
<g class="cell">
<image class="cellpng" href="/static/cell_alive.png" x="0" y="0" width="1200" height="675"/>
<image class="celldead" href="/static/cell_dead.png" x="0" y="0" width="1200" height="675" opacity="0"/>
<path id="swim" d="M-260 330 C 60 250 260 150 420 185 C 490 200 560 212 600 228" fill="none" style="visibility:hidden"/>
<path class="crack" d="M600 333 l-34 44 24 36 -38 33 20 32" fill="none" style="stroke:var(--err-tx)" stroke-width="2.5" opacity="0"/>
</g>
<g class="strand-l"></g>
<g class="strand-r"></g>
<circle class="mutsite" cx="600" cy="331" r="5.5" style="fill:var(--err-tx)" opacity=".9"/>
<circle class="flash" cx="600" cy="333" r="8" fill="none" style="stroke:var(--accent)" stroke-width="3" opacity="0"/>
<circle class="frag" cx="594" cy="324" r="4" style="fill:var(--warn-edge)" opacity="0"/>
<circle class="frag" cx="608" cy="326" r="3.4" style="fill:var(--warn-edge)" opacity="0"/>
<g class="pts">
<circle class="pt" cx="600" cy="330" r="3.5" style="fill:var(--accent)" opacity="0"/><circle class="pt" cx="600" cy="330" r="3" style="fill:var(--accent)" opacity="0"/><circle class="pt" cx="600" cy="330" r="4" style="fill:var(--accent)" opacity="0"/><circle class="pt" cx="600" cy="330" r="3" style="fill:var(--accent)" opacity="0"/><circle class="pt" cx="600" cy="330" r="3.5" style="fill:var(--accent)" opacity="0"/><circle class="pt" cx="600" cy="330" r="3" style="fill:var(--accent)" opacity="0"/>
</g>
<g class="rnp"><g transform="translate(600 252)">
<image class="protimg" href="/static/cas12a2_rnp.png" x="-71" y="-143" width="270" height="228"/>
</g></g>
</svg>
<div class="splash-cap"><b>Cas12a2 靶向杀伤机制</b><span id="capdetail"></span></div>
</div>
<div class="entry-inner" hidden>
<svg class="entry-art" viewBox="0 0 480 745" role="img" aria-label="crRNA 结构示意：DR 骨架茎环与 spacer 尾">
<path d="M150 572 L172 522 L190 468 L205 410 L205 150 L179 115 L216 84 L264 84 L300 115 L275 150 L275 410 L292 466 L310 518 L322 572 L300 622 L268 656 L228 676 L188 684 L150 676 L118 656 L96 624" stroke="#B9B3A6" stroke-width="1.6" fill="none"/>
<path d="M96 624 C78 646 58 658 34 664" stroke="#B9B3A6" stroke-width="1.6" stroke-dasharray="1 7" stroke-linecap="round" fill="none"/>
<g stroke="#1C1B19" stroke-width="1.6" opacity=".28">
<path d="M205 150 L275 150"/><path d="M205 215 L275 215"/><path d="M205 280 L275 280"/><path d="M205 345 L275 345"/><path d="M205 410 L275 410"/>
</g>
<path d="M291 280 L350 280" stroke="#B9B3A6" stroke-width="1.2"/>
<text x="358" y="276" font-size="12.5" fill="#6E685C" font-family="Helvetica Neue,Segoe UI,PingFang SC,Microsoft YaHei,sans-serif">DR 骨架 · 19 nt</text>
<text x="358" y="293" font-size="11" fill="#A39D8F" font-family="Helvetica Neue,Segoe UI,PingFang SC,Microsoft YaHei,sans-serif">茎 5 bp · 环 4 nt</text>
<text x="240" y="728" font-size="12.5" fill="#6E685C" text-anchor="middle" font-family="Helvetica Neue,Segoe UI,PingFang SC,Microsoft YaHei,sans-serif">spacer · 24 nt（示意为前 8 nt）</text>
<text x="120" y="592" font-size="11" fill="#A39D8F" font-family="JetBrains Mono,Consolas,monospace">5'</text>
<g font-family="Helvetica Neue,Segoe UI,PingFang SC,Microsoft YaHei,sans-serif" font-size="13" font-weight="600" text-anchor="middle">
<g fill="#FFFFFF" stroke="#1C1B19" stroke-width="1.4">
<circle cx="150" cy="572" r="16"/><circle cx="172" cy="522" r="16"/><circle cx="190" cy="468" r="16"/>
<circle cx="205" cy="410" r="16"/><circle cx="205" cy="345" r="16"/><circle cx="205" cy="280" r="16"/><circle cx="205" cy="215" r="16"/><circle cx="205" cy="150" r="16"/>
<circle cx="179" cy="115" r="16"/><circle cx="216" cy="84" r="16"/><circle cx="264" cy="84" r="16"/><circle cx="300" cy="115" r="16"/>
<circle cx="275" cy="150" r="16"/><circle cx="275" cy="215" r="16"/><circle cx="275" cy="280" r="16"/><circle cx="275" cy="345" r="16"/><circle cx="275" cy="410" r="16"/>
<circle cx="292" cy="466" r="16"/><circle cx="310" cy="518" r="16"/>
</g>
<g fill="#1C1B19" stroke="none">
<text x="150" y="576.5">A</text><text x="172" y="526.5">A</text><text x="190" y="472.5">T</text>
<text x="205" y="414.5">T</text><text x="205" y="349.5">T</text><text x="205" y="284.5">C</text><text x="205" y="219.5">T</text><text x="205" y="154.5">A</text>
<text x="179" y="119.5">C</text><text x="216" y="88.5">T</text><text x="264" y="88.5">G</text><text x="300" y="119.5">T</text>
<text x="275" y="154.5">G</text><text x="275" y="219.5">T</text><text x="275" y="284.5">G</text><text x="275" y="349.5">A</text><text x="275" y="414.5">G</text>
<text x="292" y="470.5">A</text><text x="310" y="522.5">T</text>
</g>
<g fill="#E3EEEA" stroke="#0F5F52" stroke-width="1.4">
<circle cx="322" cy="572" r="16"/><circle cx="300" cy="622" r="16"/><circle cx="268" cy="656" r="16"/><circle cx="228" cy="676" r="16"/>
<circle cx="188" cy="684" r="16"/><circle cx="150" cy="676" r="16"/><circle cx="118" cy="656" r="16"/><circle cx="96" cy="624" r="16"/>
</g>
<g fill="#0B4A40" stroke="none">
<text x="322" y="576.5">G</text><text x="300" y="626.5">T</text><text x="268" y="660.5">T</text><text x="228" y="680.5">C</text>
<text x="188" y="688.5">A</text><text x="150" y="680.5">T</text><text x="118" y="660.5">G</text><text x="96" y="628.5">C</text>
</g>
</g>
</svg>
<div class="brand">
<svg viewBox="0 0 32 32" fill="none" aria-hidden="true"><path d="M11 27V13a5 5 0 0 1 10 0v14" stroke="#0F5F52" stroke-width="2.6" stroke-linecap="round"/><path d="M11 21.5h10M11 16.5h10" stroke="#0F5F52" stroke-width="2" stroke-linecap="round" opacity=".55"/></svg>
<span class="wordmark">GuideForge <em>· Cas12a2 crRNA Design Agent</em></span>
</div>
<div class="entry-mid">
<h1 class="display">Cas12a2 crRNA<br><span class="cjk">骨架分型设计智能体</span></h1>
<p class="sub-en">Target-adaptive scaffold typing for Cas12a2 crRNA design</p>
<p class="lede">骨架最优解随靶标特征变化。输入目标突变序列，输出 spacer 设计、骨架选型与完整 crRNA 构建。</p>
<button class="enter" onclick="enterApp()">进入演示</button>
</div>
<div class="masthead">
<div class="stat"><b>32</b><span>8 骨架 × 4 靶标 IVT 验证矩阵</span></div>
<div class="stat"><b>0.874</b><span>Chai-1 共折叠 ipTM（8D4A 自模板注入）</span></div>
<div class="stat"><b>0.867</b><span>虚拟细胞剂量标定 R²（Scholz 2026）</span></div>
</div>
</div>
</div>
<div class="view" id="app" hidden>
<div class="wrap">
<a class="skip" href="#main">跳到主内容</a>
<div class="topbar">
<div class="brand">
<svg viewBox="0 0 32 32" fill="none" aria-hidden="true"><path d="M11 27V13a5 5 0 0 1 10 0v14" stroke="#0F5F52" stroke-width="2.6" stroke-linecap="round"/><path d="M11 21.5h10M11 16.5h10" stroke="#0F5F52" stroke-width="2" stroke-linecap="round" opacity=".55"/></svg>
<span class="wordmark">GuideForge <em>· Cas12a2 crRNA Design Agent</em></span>
</div>
<a class="back" href="#" onclick="backHome();return false;">返回首页</a>
</div>
<header>
<h1>Cas12a2 crRNA 设计智能体</h1>
<p class="sub">输入 spacer 序列或选择 panel 靶标，输出骨架分型建议与文献先验活性排序。</p>
<div class="disclaimer">口径说明：结构口径分型建议 + 文献先验活性排序（Han 2025 Fig1g 尺度，n=7 外部锚点，非活性保证）；含阴性退化模式（策划案 §5.3/表5 风险1）：分型 go/no-go 判据未通过时退化为通用型单骨架推荐，页面如实标注并给量化证据</div>
<div class="layers">
<div class="layer"><b>蛋白预测层</b><span>Chai-1 三元共折叠，ipTM 0.874</span></div>
<div class="layer"><b>AI 改造层</b><span>受约束生成 × 贝叶斯优化闭环</span></div>
<div class="layer"><b>虚拟细胞层</b><span>剂量标定 EC50 11.3 FPKM，R² 0.867</span></div>
<div class="layer"><b>智能体层</b><span>本页双模式端到端演示</span></div>
</div>
</header>
<main class="grid" id="main">
<section class="card">
<div class="card-head"><span class="tag">模式一</span><h2>spacer 设计</h2><span class="meta">实算 · 秒级</span></div>
<label for="sp">spacer 序列（17-25 nt，ACGT）</label>
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
<footer><span>GuideForge · Cas12a2 crRNA 骨架分型（演示环境）</span><span>8 骨架 × 4 靶标 IVT 矩阵口径 · 数据为脚本实算 / 预计算</span></footer>
</div>
</div>
<noscript><style>#entry{display:none}#app{display:block}</style></noscript>
<script>
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
function loading(el,h){el.innerHTML='<div class="sk"></div><div class="sk" style="width:70%"></div><div class="sk" style="width:45%"></div>';if(h)h.textContent='计算中…';}
function fail(el,h,msg){el.innerHTML='<div class="err">'+esc(msg)+'</div>';if(h)h.textContent='';}
function rankTable(rows){
if(!rows||!rows.length)return'';
const vals=rows.map(r=>r.lit_pred_fig1g),lo=Math.min.apply(null,vals),hi=Math.max.apply(null,vals);
const hasIf=rows.some(r=>r.interface_delta&&r.interface_delta.delta_prot_crRNA!=null);
let h='<table><thead><tr><th>#</th><th>骨架</th><th class="num">Fig1g 预测</th><th class="num">ddG</th><th class="num">bp_dist</th><th class="num">cross_nt</th><th class="num">p_fold</th><th class="num">spacer_up</th>'+(hasIf?'<th class="num">Δprot-crRNA</th>':'')+'</tr></thead><tbody>';
for(const r of rows){const f=r.features||{};const w=hi>lo?((r.lit_pred_fig1g-lo)/(hi-lo)*100):50;
const dpc=r.interface_delta?r.interface_delta.delta_prot_crRNA:null;
h+='<tr><td class="num">'+r.rank+'</td><td>'+esc(r.scaffold)+'</td>'+
'<td class="num">'+r.lit_pred_fig1g.toFixed(4)+'<div class="bar" style="width:'+w.toFixed(1)+'%"></div></td>'+
'<td class="num">'+(f.ddG_dr!=null?f.ddG_dr:'')+'</td><td class="num">'+(f.bp_dist!=null?f.bp_dist:'')+'</td>'+
'<td class="num">'+(f.cross_nt!=null?f.cross_nt:'')+'</td><td class="num">'+(f.p_fold!=null?f.p_fold:'')+'</td>'+
'<td class="num">'+(f.spacer_up!=null?f.spacer_up:'')+'</td>'+
(hasIf?'<td class="num">'+(dpc!=null?((dpc>0?'+':'')+dpc.toFixed(4)):'—')+'</td>':'')+'</tr>';}
return h+'</tbody></table><div class="note">Fig1g 尺度：数值越低 = 预测抑制越强（Han 2025 文献先验，n=7 锚点）。'+(hasIf?' Δprot-crRNA 列为 Chai-1 界面差量（自模板口径结构描述特征，主口径链对 ipTM 相对同靶 WT 差量，非界面可预测性证据）。':'')+'</div>';}
function raw(d){return '<details><summary>查看原始 JSON</summary><pre>'+esc(JSON.stringify(d,null,1))+'</pre></details>';}
function abBlock(ta){
if(!ta||!ta.feature)return'';
const f=ta.feature,g=ta.gate||{};
let h='<label style="margin-top:14px">靶RNA丰度档位（模块一特征维 · CCLE 18q3 实测）</label><dl class="kv">'+
'<dt>基因@细胞系</dt><dd>'+esc(f.gene)+'@'+esc(f.cell_line)+'</dd>'+
'<dt>实测丰度</dt><dd>'+f.total_rpkm+' RPKM（杂合突变场景 '+f.mut_rpkm_het50+' RPKM）</dd>'+
'<dt>丰度档位</dt><dd><span class="badge">tier '+f.target_abundance_tier+' · '+esc(f.target_abundance_tier_label)+'</span></dd></dl>';
if(g.note)h+='<div class="note">'+esc(g.note)+'</div>';
return h;}
function mtBlock(mt){
if(!mt)return'';
const w=mt.wt_orf||{},m=mt.mut_orf_same_start||{},ps=mt.premature_stop;
let h='<label style="margin-top:14px">突变类型判定（模块一特征维 · 最长 ORF 法）</label><dl class="kv">'+
'<dt>类型</dt><dd><span class="badge">'+esc(mt.mutation_type)+' · '+esc(mt.mutation_type_label)+'</span>'+(mt.aa_change?' <span class="badge">'+esc(mt.aa_change)+'</span>':'')+'</dd>'+
'<dt>WT 最长 ORF</dt><dd>转录本 '+w.start_1based+'..'+w.end_1based+'（'+w.protein_aa+' aa）</dd>';
if(m)h+='<dt>MUT 同起点 ORF</dt><dd>转录本 '+m.start_1based+'..'+m.end_1based+'（'+m.protein_aa+' aa）'+(mt.orf_length_change_aa?'，变化 '+mt.orf_length_change_aa+' aa':'')+'</dd>';
if(ps)h+='<dt>提前终止子</dt><dd>密码子 '+ps.codon_no+'（'+esc(ps.wt_codon)+'→'+esc(ps.mut_codon)+'，转录本 '+ps.transcript_pos_1based+' 位）'+(mt.truncated_fraction!=null?'，截短 '+(mt.truncated_fraction*100).toFixed(1)+'%':'')+'</dd>';
h+='</dl>';
if(mt.reading)h+='<div class="note">'+esc(mt.reading)+'</div>';
return h;}
function ifBlock(ib){
if(!ib)return'';
if(ib.available===false)return '<label style="margin-top:14px">界面特征差量（结构层选型特征维 · Chai-1 骨架-靶标组合）</label><div class="note">该靶标暂无 Chai 矩阵覆盖：'+esc(ib.missing_reason||'')+'</div>';
const wt=ib.wt_reference||{};
const sg=v=>(v==null?'':((v>0?'+':'')+Number(v).toFixed(4)));
let h='<label style="margin-top:14px">界面特征差量 Δ vs 同靶 WT（结构层选型特征维 · Chai-1 自模板口径）</label><dl class="kv">'+
'<dt>矩阵靶标</dt><dd>'+esc(ib.matrix_target||'')+'</dd>'+
'<dt>WT 参照</dt><dd>agg ipTM '+wt.iptm_mean+' · prot-crRNA '+wt.prot_crRNA+' · crRNA-靶 '+wt.crRNA_target+' · clash '+(wt.clash_frac!=null?wt.clash_frac:'')+'</dd></dl>';
h+='<table><thead><tr><th>骨架</th><th class="num">Δagg ipTM</th><th class="num">Δprot-crRNA</th><th class="num">ΔcrRNA-靶</th><th class="num">Δclash</th><th>判读</th></tr></thead><tbody>';
for(const r of (ib.scaffolds||[])){const d=r.delta||{};
h+='<tr><td>'+esc(r.scaffold)+'</td><td class="num">'+sg(d.delta_iptm_mean)+'</td><td class="num">'+sg(d.delta_prot_crRNA)+'</td><td class="num">'+sg(d.delta_crRNA_target)+'</td><td class="num">'+sg(d.delta_clash_frac)+'</td><td>'+esc(r.reading||'')+'</td></tr>';}
h+='</tbody></table>';
if(ib.limitation)h+='<div class="note">'+esc(ib.limitation)+'</div>';
return h;}
function dgrBlock(dg){
if(!dg)return'';
const rel=dg.n_common<=dg.threshold?'≤':'>';
const nums='判据数值（预登记 §A go/no-go）：全型公共 TOP-'+dg.topk+' 子集 '+dg.n_common+' 条（'+rel+' 阈 '+dg.threshold+'），型间两两 Jaccard 均值 '+(dg.pairwise_jaccard_mean!=null?dg.pairwise_jaccard_mean:'n/a');
if(dg.mode==='degraded_universal'){
let h='<div class="disclaimer" style="margin-top:14px"><b>退化模式（阴性结果，如实标注）</b>：'+esc(dg.reading)+'<br>'+esc(nums);
if(dg.robustness_evidence)h+='<br>'+esc(dg.robustness_evidence.reading);
if(dg.deliverable_note)h+='<br>'+esc(dg.deliverable_note);
return h+'</div>';}
return '<div class="note">选型模式：分型推荐。'+esc(nums)+'</div>';}
async function go(){
const inp=$('sp');const sp=inp.value.trim().toUpperCase();
if(!/^[ACGTU]{17,25}$/.test(sp)){inp.classList.add('invalid');
fail($('r1'),$('hint1'),'spacer 须为 17-25 nt 的 ACGT(U) 序列，请检查输入。');return;}
inp.classList.remove('invalid');inp.value=sp;
loading($('r1'),$('hint1'));
let r,d;
try{r=await fetch('/api/design',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({spacer:sp})});d=await r.json();}
catch(e){fail($('r1'),$('hint1'),'网络请求失败，请检查服务是否在线。');return;}
if(!r.ok||d.error){fail($('r1'),$('hint1'),d.error||('请求失败（HTTP '+r.status+'）'));return;}
$('hint1').textContent='';
const degraded=d.degradation&&d.degradation.mode==='degraded_universal';
let h='<dl class="kv">'+
'<dt>输入 spacer</dt><dd class="seq">'+esc(d.input_spacer)+'</dd>'+
(degraded?
'<dt>骨架选型</dt><dd><span class="badge" style="background:var(--warn-bg);color:var(--warn-tx);border-color:var(--warn-edge)">退化模式 · 通用型单骨架</span>（分型诊断：最近型 '+esc(String(d.typing_nearest_type))+'，置信度 '+(d.confidence!=null?d.confidence:'n/a')+'，不参与推荐）</dd>'
:'<dt>骨架分型</dt><dd><span class="badge">'+esc(d.scaffold_type)+'</span> &nbsp;置信度 '+(d.confidence!=null?d.confidence:'n/a')+'</dd>')+
(!degraded&&d.type_representative?'<dt>型代表</dt><dd class="seq">'+esc(d.type_representative)+'</dd>':'')+
'</dl>';
h+=dgrBlock(d.degradation);
if(d.recommended){const c=d.recommended;
h+='<label>完整 crRNA 构建（'+(degraded?'通用型单骨架 · 退化模式':'推荐骨架')+' '+esc(String(c.scaffold))+'）</label><dl class="kv">'+
'<dt>crRNA 5′→3′</dt><dd class="seq"><span class="dr">'+esc(c.dr_rna)+'</span><span class="sp">'+esc(c.spacer_rna)+'</span></dd>'+
'<dt>合成订购（DNA）</dt><dd class="seq"><span class="dr">'+esc(c.dr_rna.replace(/U/g,'T'))+'</span><span class="sp">'+esc(c.spacer_rna.replace(/U/g,'T'))+'</span></dd>'+
'<dt>组成</dt><dd>DR 骨架 '+c.dr_len+' nt + spacer '+c.spacer_len+' nt ＝ '+(c.dr_len+c.spacer_len)+' nt</dd></dl>';}
h+='<label>骨架文献先验活性排序</label>'+rankTable(d.scaffold_ranking_lit);
if(d.target_abundance)h+=abBlock(d.target_abundance);
h+=ifBlock(d.interface_deltas);
if(d.note)h+='<div class="note">'+esc(d.note)+'</div>';
h+=raw(d);$('r1').innerHTML=h;animateResult($('r1'));}
async function panel(){loading($('r2'),null);
const t=$('pt').value;let r,d;
try{r=await fetch('/api/panel?target='+encodeURIComponent(t));d=await r.json();}
catch(e){fail($('r2'),null,'网络请求失败，请检查服务是否在线。');return;}
if(!r.ok||d.error){fail($('r2'),null,d.error||('请求失败（HTTP '+r.status+'）'));return;}
let h='<dl class="kv"><dt>靶标</dt><dd><span class="badge">'+esc(d.target)+'</span></dd>';
const dz=(d.design&&d.design.designs)||[];
if(dz.length){h+='<dt>候选 spacer</dt><dd>'+dz.length+' 条，首位：</dd>';h+='<dt></dt><dd class="seq">'+esc(dz[0].spacer_dna)+'</dd>';
if(dz[0].construct_dna){h+='<dt>完整 crRNA（DNA，'+esc(dz[0].dr_desc||'')+' 骨架）</dt><dd class="seq"><span class="dr">'+esc(dz[0].dr_dna||'')+'</span><span class="sp">'+esc(dz[0].spacer_dna)+'</span></dd>';}}
h+='</dl>';
h+=dgrBlock(d.degradation||(d.design&&d.design.degradation));
h+=mtBlock(d.mutation_typing);
h+=abBlock(d.target_abundance);
h+=ifBlock(d.interface_deltas);
if(d.scaffold_ranking_lit){h+='<label>骨架文献先验活性排序</label>'+rankTable(d.scaffold_ranking_lit);}
const vc=d.virtual_cell;
if(vc){const cal=vc.calibration||{};
const calTxt=(cal.ec50_fpkm!=null)?('EC50 '+cal.ec50_fpkm+' FPKM（95% CI '+cal.ec50_ci95[0]+'-'+cal.ec50_ci95[1]+'），Hill '+cal.hill+'，R²='+cal.r2+'，n='+cal.n_points+'，'+esc(cal.source||'')):esc(vc.calibration||'n/a');
h+='<label style="margin-top:14px">细胞层预测（虚拟细胞）</label><dl class="kv">'+
'<dt>标定</dt><dd>'+calTxt+'</dd><dt>状态</dt><dd>'+esc(vc.status||'n/a')+'</dd></dl>';
const nl=vc.named_lines||[];
if(nl.length){h+='<table><thead><tr><th>细胞系</th><th class="num">RPKM</th><th class="num">预测存活 %</th></tr></thead><tbody>';
for(const x of nl){h+='<tr><td>'+esc(x.cell_line)+'</td><td class="num">'+(x.rpkm!=null?x.rpkm:'n/a')+'</td><td class="num">'+(x.pred_survival_mid_pct!=null?x.pred_survival_mid_pct:'n/a')+'</td></tr>';}
h+='</tbody></table>';}
if(vc.note)h+='<div class="note">'+esc(vc.note)+'</div>';}
h+=raw(d);$('r2').innerHTML=h;animateResult($('r2'));}
let motionOK=false;
if(typeof window.gsap!=='undefined'){const mm=gsap.matchMedia();
mm.add('(prefers-reduced-motion: no-preference)',()=>{motionOK=true;
document.body.classList.add('gsap');
return()=>{motionOK=false;document.body.classList.remove('gsap');};});
if(/[?&]fast/.test(location.search)){gsap.ticker.lagSmoothing(0);}}
function splitTitle(){const h=document.querySelector('.display');
if(!h||h.dataset.split)return;h.dataset.split='1';
h.setAttribute('aria-label',h.textContent.trim());
const lines=[];let cur=[];
h.childNodes.forEach(n=>{
if(n.nodeName==='BR'){lines.push(cur);cur=[];}
else{cur.push.apply(cur,n.textContent.split(''));}});
lines.push(cur);
h.innerHTML=lines.map(cs=>'<span class="line" aria-hidden="true">'+
cs.map(c=>'<span class="ch">'+(c===' '?'&nbsp;':esc(c))+'</span>').join('')+
'</span>').join('');}
function countUp(){document.querySelectorAll('.masthead .stat b').forEach(el=>{
const t=el.textContent.trim();const target=parseFloat(t);
if(isNaN(target))return;
const dec=(t.split('.')[1]||'').length;const o={v:0};
gsap.to(o,{v:target,duration:1.2,ease:'power2.out',delay:1.05,
onUpdate:()=>{el.textContent=o.v.toFixed(dec);},
onComplete:()=>{el.textContent=t;}});});}
function magnet(){const b=document.querySelector('.enter');
if(!b||b.dataset.mag)return;b.dataset.mag='1';
const xTo=gsap.quickTo(b,'x',{duration:.35,ease:'power3.out'});
const yTo=gsap.quickTo(b,'y',{duration:.35,ease:'power3.out'});
b.addEventListener('mousemove',e=>{const r=b.getBoundingClientRect();
xTo((e.clientX-r.left-r.width/2)*.16);yTo((e.clientY-r.top-r.height/2)*.28);});
b.addEventListener('mouseleave',()=>{xTo(0);yTo(0);});}
function animateResult(el){
if(!motionOK||!el)return;
const kv=el.querySelectorAll('.kv dt,.kv dd');
if(kv.length)gsap.from(kv,{autoAlpha:0,y:8,stagger:.03,duration:.35,clearProps:'all'});
const rows=el.querySelectorAll('tbody tr');
if(rows.length)gsap.from(rows,{autoAlpha:0,y:10,stagger:.05,duration:.4,clearProps:'all'});
const bars=el.querySelectorAll('.bar');
if(bars.length)gsap.from(bars,{scaleX:0,transformOrigin:'left center',
duration:.7,ease:'power3.out',stagger:.06,clearProps:'transform'});}
function animateEntry(){
if(!motionOK)return;
splitTitle();magnet();
const tl=gsap.timeline({defaults:{ease:'power2.out'}});
tl.from('.entry-inner .brand',{autoAlpha:0,y:14,duration:.5,clearProps:'all'},0)
.from('.display .ch',{yPercent:112,autoAlpha:0,duration:.65,ease:'power3.out',
stagger:.024,clearProps:'all'},.12)
.from('.sub-en',{autoAlpha:0,y:14,duration:.5,clearProps:'all'},.62)
.from('.lede',{autoAlpha:0,y:18,duration:.55,clearProps:'all'},.7)
.from('.enter',{autoAlpha:0,y:14,duration:.45,clearProps:'all'},.84)
.from('.masthead .stat',{autoAlpha:0,y:16,stagger:.09,duration:.5,clearProps:'all'},.98);
countUp();
const fig=document.querySelector('.entry-art');
if(!fig)return;
const bb=fig.querySelector('path');const L=bb.getTotalLength();
const stub=fig.querySelector('path[stroke-dasharray]');
const ftl=gsap.timeline({delay:.15});
ftl.fromTo(bb,{strokeDasharray:L,strokeDashoffset:L},
{strokeDashoffset:0,duration:1.3,ease:'power2.inOut',
clearProps:'strokeDasharray,strokeDashoffset'},0)
.from(fig.querySelectorAll('g[opacity] path'),{autoAlpha:0,duration:.3,stagger:.06},.55)
.from(fig.querySelectorAll('circle'),{scale:0,transformOrigin:'50% 50%',
duration:.35,ease:'back.out(1.6)',stagger:.04,clearProps:'transform'},.5)
.from(fig.querySelectorAll('g text'),{autoAlpha:0,duration:.25,stagger:.04,clearProps:'all'},.62)
.from(fig.querySelectorAll(':scope>text'),{autoAlpha:0,y:6,duration:.4,stagger:.12,clearProps:'all'},1.6);
if(stub)ftl.from(stub,{autoAlpha:0,duration:.5,clearProps:'all'},1.7);}
function animateApp(){
if(!motionOK)return;
gsap.from('#app .topbar',{autoAlpha:0,y:-10,duration:.4,ease:'power2.out',clearProps:'all'});
gsap.from('#app header>*',{autoAlpha:0,y:16,stagger:.07,duration:.5,ease:'power2.out',clearProps:'all'});
gsap.from('#app .card',{autoAlpha:0,y:18,stagger:.1,duration:.55,ease:'power2.out',clearProps:'all'});
gsap.from('#app footer',{autoAlpha:0,duration:.5,delay:.35,clearProps:'all'});}
let opened=false;
function isDark(){return document.body.classList.contains('dark');}
function startScene(){
const sp=$('splash'),inner=document.querySelector('.entry-inner');
if(!motionOK){opened=true;sp.hidden=true;inner.hidden=false;return;}
const cap=$('capdetail');
const cr=document.querySelector('.crack');
const cl=cr.getTotalLength();
gsap.set(cr,{strokeDasharray:cl,strokeDashoffset:cl});
function bead(x,y,r,cls){return '<circle class="'+cls+'" cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="'+r+'"/>';}
function cap2(x,y1,y2,cls){const t=Math.min(y1,y2),h=Math.abs(y2-y1);
return '<rect class="bp '+cls+'" x="'+(x-2.6).toFixed(1)+'" y="'+t.toFixed(1)+'" width="5.2" height="'+Math.max(h,1.5).toFixed(1)+'" rx="2.6"/>';}
function helix(sel,x0,x1,y,mutX,c0,c1){const g=document.querySelector(sel);if(!g)return;
let bk='',rg='',cr='',mu='';
for(let x=x0;x<=x1+.1;x+=6.5){const s=Math.sin((x-330)/46*Math.PI*2);
bk+=bead(x,y+s*15.5,5.4,'bead')+bead(x,y-s*15.5,5.4,'bead');}
let i=0;
for(let x=x0+11.5;x<x1-3;x+=23){const s=Math.sin((x-330)/46*Math.PI*2);
rg+=cap2(x,y+s*15.5,y,'bp'+(i%5))+cap2(x,y,y-s*15.5,'bp'+((i+2)%5));i++;}
if(mutX)mu=cap2(mutX,y-15.5,y+15.5,'bp-mut');
if(c0<c1)for(let x=c0;x<=c1;x+=7.2){cr+=bead(x,y,4.3,'bead-r');}
g.innerHTML=rg+mu+cr+bk;}
helix('.strand-l',330,596,332,585,506,596);
helix('.strand-r',604,868,332,0,604,694);
const hasPath=typeof MotionPathPlugin!=='undefined';
if(hasPath)gsap.registerPlugin(MotionPathPlugin);
const wd=setTimeout(()=>{if(opened)return;opened=true;
if(window._scene){window._scene.kill();window._scene=null;}
sp.hidden=true;inner.hidden=false;animateEntry();
setTimeout(()=>{try{gsap.globalTimeline.getChildren(true,true,true)
.forEach(t=>t.progress(1));}catch(e){}},3500);},6000);
const tq=new URLSearchParams(location.search).get('t');
const tl=gsap.timeline();
window._scene=tl;
tl.add(()=>{cap.textContent='RNP 复合物识别突变转录本';},.05);
if(hasPath){
tl.to('.rnp',{duration:1.5,ease:'power2.out',
motionPath:{path:'#swim',align:'#swim',alignOrigin:[.5,.5]}},.1);
}else{
tl.fromTo('.rnp',{x:-540,y:-90,rotation:-9},
{x:0,y:0,rotation:0,duration:1.35,ease:'power2.out'},.1);}
tl.to('.rnp',{y:-8,duration:.3,ease:'sine.out'},1.6)
.fromTo('.mutsite',{scale:1,transformOrigin:'50% 50%'},
{scale:1.9,duration:.28,yoyo:true,repeat:1,ease:'power2.inOut'},1.65)
.to('.rnp',{y:22,scaleY:.92,scaleX:1.05,transformOrigin:'50% 100%',
duration:.22,ease:'power3.in'},1.95)
.to('.rnp',{scaleY:1,scaleX:1,duration:.38,ease:'elastic.out(1,.45)'},2.17)
.fromTo('.flash',{attr:{r:8},autoAlpha:.9},
{attr:{r:58},autoAlpha:0,duration:.55,ease:'power2.out'},2.1)
.add(()=>{cap.textContent='Cas12a2 激活 · 旁切效应开启';},2.12)
.to('.mutsite',{autoAlpha:0,duration:.15},2.12)
.to('.strand-l',{x:-78,y:9,rotation:-6,svgOrigin:'596 333',opacity:.55,
duration:.3,ease:'power4.in'},2.12)
.to('.strand-r',{x:78,y:9,rotation:6,svgOrigin:'604 333',opacity:.55,
duration:.3,ease:'power4.in'},2.12)
.to('.strand-l',{x:-66,rotation:-4.5,duration:.4,ease:'power2.out'},2.42)
.to('.strand-r',{x:66,rotation:4.5,duration:.4,ease:'power2.out'},2.42)
.fromTo('.frag',{autoAlpha:1},{x:i=>i?52:-48,y:i=>i?-40:-26,autoAlpha:0,
duration:.7,ease:'power2.out'},2.15)
.to('.pt',{opacity:1,duration:.1,stagger:.04},2.3)
.to('.pt',{x:()=>gsap.utils.random(-170,170),y:()=>gsap.utils.random(-110,90),
autoAlpha:0,duration:.9,ease:'power2.out',stagger:.05},2.35)
.add(()=>{cap.textContent='靶细胞死亡';},3.1)
.to(cr,{strokeDashoffset:0,opacity:1,duration:.45,ease:'power2.out'},3.15)
.to('.celldead',{opacity:1,duration:.9,ease:'sine.inOut'},3.15)
.to('.cell',{scale:.955,y:8,svgOrigin:'600 347',duration:.9,ease:'sine.inOut'},3.15)
.to('.rnp',{opacity:.72,duration:.9,ease:'sine.inOut'},3.15);
tl.to(sp,{autoAlpha:0,y:-10,duration:.5,ease:'power2.in'},4.3)
.add(()=>{opened=true;window._scene=null;
sp.hidden=true;gsap.set(sp,{clearProps:'all'});inner.hidden=false;animateEntry();},4.55);
if(tq!==null){clearTimeout(wd);tl.pause();tl.seek(parseFloat(tq));}}
function openRNP(){
if(opened)return;opened=true;
if(window._scene){window._scene.kill();window._scene=null;}
const sp=$('splash'),inner=document.querySelector('.entry-inner');
if(!motionOK){sp.hidden=true;inner.hidden=false;return;}
gsap.to(sp,{autoAlpha:0,duration:.4,ease:'power2.in',
onComplete:()=>{sp.hidden=true;gsap.set(sp,{clearProps:'all'});}});
setTimeout(()=>{sp.hidden=true;},700);
inner.hidden=false;animateEntry();}
(function(){const sp=$('splash');if(!sp)return;
sp.addEventListener('click',openRNP);
sp.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){
e.preventDefault();openRNP();}});})();
(function(){const t=$('themetog');if(!t)return;
if(/[?&]theme=dark/.test(location.search))document.body.classList.add('dark');
const sync=()=>{t.textContent=isDark()?'浅色':'深色演示';};
t.addEventListener('click',()=>{document.body.classList.toggle('dark');sync();});
sync();})();
function show(v){const e=$('entry'),a=$('app');
if(v==='app'){
if(motionOK&&!e.hidden){
gsap.to(e,{autoAlpha:0,y:-14,duration:.35,ease:'power2.in',
onComplete:()=>{e.hidden=true;gsap.set(e,{clearProps:'all'});}});
setTimeout(()=>{e.hidden=true;},600);}
else{e.hidden=true;}
a.hidden=false;window.scrollTo(0,0);animateApp();}
else{a.hidden=true;e.hidden=false;window.scrollTo(0,0);
if(opened)animateEntry();}}
function enterApp(){history.replaceState(null,'','#app');show('app');}
function backHome(){history.replaceState(null,'','#');show('entry');}
window.addEventListener('DOMContentLoaded',()=>{
const h=location.hash;
if(h==='#demo-spacer'){show('app');go();}
else if(h==='#demo-panel'){show('app');panel();}
else if(h==='#app'){show('app');}
else{startScene();}});
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
    # 可选靶标上下文: 提供时附模块一丰度档维 + 模块三门控判读
    # + 结构层界面特征差量附加层(任务⑫, 与离线 agent.target_context 同源)
    tkey = (request.json or {}).get("target")
    tctx = None
    iblock = None
    if tkey:
        feat = tabund.target_feature(ABUNDANCE, tkey)
        if feat is None:
            return jsonify({"error": "未知 panel 靶标 %s" % tkey}), 400
        tctx = {"feature": feat, "gate": tabund.activation_gate(feat)}
        iblock = chai_if.target_block(CHAI_IF, tkey)
    try:
        sel = agent.select_scaffold_type(core.to_rna(sp), DR, MODEL)
        rep = MODEL["type_dr"].get(sel["nearest_type"], {})
        sp_rna = core.to_rna(sp)
        ranking = lit_rank_scaffolds(sp_rna)
        if iblock and iblock.get("available"):
            # 选型特征维附加入排序行: 库骨架名(+分隔)映射矩阵紧凑名,
            # 矩阵未覆盖的骨架如实不带该字段(不补造)
            dmap = {r["scaffold"]: r["delta"] for r in iblock["scaffolds"]}
            for r in ranking:
                dlt = dmap.get(chai_if.matrix_scaffold_of(r["scaffold"]))
                if dlt is not None:
                    r["interface_delta"] = dlt
        degraded = DEGRADATION["mode"] == "degraded_universal"
        if degraded:
            # 退化模式: 推荐 = 通用型单骨架(分型产物全场最优, 与离线 agent 同源);
            # 分型结果仅作诊断(typing_nearest_type), 文献先验排序仍附作证据
            uni = agent.universal_scaffold(MODEL)
            best = {"scaffold": uni["desc"], "dr": core.to_rna(uni["dr_dna"])}
            payload_type = {"scaffold_type": None,
                            "typing_nearest_type": sel["nearest_type"],
                            "type_representative": None}
        else:
            best = ranking[0]
            payload_type = {"scaffold_type": sel["nearest_type"],
                            "type_representative": rep.get("representative_desc")}
        cr_rna = best["dr"] + sp_rna
        return jsonify({"input_spacer": sp,
                        "confidence": sel["confidence"],
                        "features": sel["features"],
                        "target_abundance": tctx,
                        "interface_deltas": iblock,
                        "degradation": DEGRADATION,
                        **payload_type,
                        "recommended": {"scaffold": best["scaffold"],
                                        "dr_rna": best["dr"],
                                        "spacer_rna": sp_rna,
                                        "crrna_rna": cr_rna,
                                        "crrna_dna": cr_rna.replace("U", "T"),
                                        "dr_len": len(best["dr"]),
                                        "spacer_len": len(sp_rna)},
                        "scaffold_ranking_lit": ranking,
                        "note": ("退化模式(阴性结果): 通用型单骨架推荐 + 文献先验"
                                 "活性排序附作证据" if degraded else
                                 "结构分型建议 + 文献先验活性排序"
                                 "(Han 2025 Fig1g 尺度, n=7 锚点, 低=抑制强,"
                                 "非活性保证)")})
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
    # 模块一丰度档维(CCLE 18q3 实测 + 显式分档, 与离线 agent 同源同口径)
    feat = tabund.target_feature(ABUNDANCE, PANEL.get(key, ""))
    if feat:
        extra["target_abundance"] = {"feature": feat,
                                     "gate": tabund.activation_gate(feat)}
    # 结构层界面特征差量(任务⑫, 选型特征维附加层): 与离线 agent
    # target_context 读同一 JSON、走同一 target_block; 矩阵未覆盖的靶标
    # (APC 列待增量折叠回填)如实返回 available=False + 缺失原因
    extra["interface_deltas"] = chai_if.target_block(CHAI_IF, PANEL.get(key, ""))
    # 模块一突变类型维(最长 ORF 法, 2026-09-07 任务⑩): 判定与 ORF 证据在
    # 预计算阶段由 crrna_mutation_typing 生成并随 design.json 归档, webapp
    # 侧只读不重构(与离线同源的证据载体 = design.json mutation.typing 块)
    mtp = (d.get("mutation") or {}).get("typing")
    if mtp:
        extra["mutation_typing"] = mtp
    # 模块四: 细胞层预测(实测标定虚拟细胞, Scholz 2026 剂量曲线口径)
    vc_path = os.path.join(DATA, "virtual_cell_prior.json")
    if os.path.isfile(vc_path):
        vc = json.load(open(vc_path, encoding="utf-8"))
        # 虚拟细胞只标定了 TP53/KRAS 两条存活曲线; APC 靶标无对应预测,
        # 如实不展示细胞系存活表(丰度信息见 target_abundance 块)
        gene = feat["gene"] if feat else (
            "TP53" if key.startswith("TP53") else "KRAS")
        if gene in ("TP53", "KRAS"):
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
    # 阴性退化状态(与离线智能体同源的启动时评估): 无论 panel 预计算设计
    # 是否含 degradation 块, 响应均携带当前判据状态, 页面如实展示
    extra["degradation"] = DEGRADATION
    return jsonify({"target": key, "design": d, **extra})


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--port", type=int, default=8899)
    ap.add_argument("--clusters", default=None,
                    help="分型产物 clusters.json(默认 data/context_typing/"
                         "typing.clusters.json); 指向其它产物可复算退化判据"
                         "(阴性场景演示)")
    ap.add_argument("--spacers", default=None,
                    help="分型 spacer 表(默认 data/context_typing/spacers.txt)")
    args = ap.parse_args()
    if args.clusters or args.spacers:
        global MODEL, DEGRADATION
        cj = args.clusters or CLUSTERS_JSON
        MODEL = agent.load_type_models(cj, args.spacers or SPACERS_TXT, DR)
        DEGRADATION = agent.interaction_assessment(cj)
    app.run(host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()
