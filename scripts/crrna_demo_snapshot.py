"""答辩离线快照演示页生成器(2026-09-02)。

Web 演示页需要 ssh -L 端口转发, 答辩现场网络/设备不可控; 本脚本把
可演示的核心产物(五靶标 panel 设计/四取向候选族/补偿三臂/预注册信息)
固化进一个零依赖单文件 HTML(可 file:// 直接打开, 无服务器无网络)。
运行: python scripts/crrna_demo_snapshot.py -> demo_snapshot.html
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")

PANEL = [("TP53-R248Q", "tp53_r248q"), ("KRAS-G12C", "kras_g12c"),
         ("KRAS-G12D", "kras_g12d"), ("TP53-R273H", "tp53_r273h"),
         ("APC-Q1312x", "apc_q1338x")]

TPL = """<!doctype html><html><head><meta charset="utf-8">
<title>GuideForge 离线演示快照</title><style>
body{font-family:sans-serif;max-width:1000px;margin:2em auto;padding:0 1em;color:#20303c}
h1{font-size:1.3em}h2{font-size:1.05em;border-left:4px solid #0e6ba8;padding-left:.5em;margin-top:1.6em}
table{border-collapse:collapse;width:100%;font-size:.9em}
td,th{border:1px solid #d8e2ea;padding:.4em .6em;text-align:left;vertical-align:top}
th{background:#eef4f9}code{background:#f0f4f8;padding:0 .3em}
.small{color:#62788a;font-size:.85em}
</style></head><body>
<h1>GuideForge — Cas12a2 crRNA 骨架分型设计与 AI 智能体(离线快照)</h1>
<div class="small">生成时间 @@DATE@@ | 输出均为结构口径排序/分型建议, 不含活性预测
(游离态指标经三模型反向验证不构成活性预测器, 见仓库 README 校准声明)。</div>
<h2>① 五靶标 panel 端到端设计(tilling + 骨架分型, 智能体预计算)</h2>
@@PANEL@@
<h2>② 四取向跨型候选骨架族(V3 §4.2, 供体外矩阵)</h2>
@@ORIENT@@
<h2>③ MYCg1 补偿突变三臂(折叠竞争因果检验, 预注册预测已固化)</h2>
@@COMP@@
<h2>④ 预注册时间戳</h2>
<div>git tag <code>preregistration-20260901</code>(GitHub/Gitee 双远端): 合成清单+候选族 sha256、
打分权重出处、预测(A_break&gt;C_WT, B 臂回落; 活性主张一律不写)。</div>
<div class="small">在线演示: 服务器 <code>python scripts/crrna_agent_webapp.py</code> + 本机
<code>ssh -L 8899:localhost:8899 srv</code> -> http://localhost:8899 (spacer 实算分型)。</div>
</body></html>"""


def main():
    rows = []
    for label, key in PANEL:
        d = json.load(open(os.path.join(DATA, "agent", key + ".design.json"),
                           encoding="utf-8"))
        top = d["designs"][0] if d.get("designs") else {}
        rows.append("<tr><th>%s</th><td><code>%s</code></td><td>%s(%s)</td>"
                    "<td>%s</td><td><code>%s</code></td></tr>" % (
                        label, top.get("spacer_dna", top.get("spacer", "-")),
                        top.get("dr_desc", "-"), top.get("confidence", "-"),
                        top.get("tier", "-"), top.get("construct_dna", "-")))
    panel = ("<table><tr><th>靶标</th><th>spacer(A)</th><th>骨架型/置信</th>"
             "<th>tier</th><th>完整 crRNA 构建</th></tr>" + "".join(rows) + "</table>")

    lib = json.load(open(os.path.join(DATA, "orientation_library.json"),
                         encoding="utf-8"))
    rows = ["<tr><th>WT 参照</th><td colspan=3><code>AATTTCTACTGTTGTAGAT"
            "GTTCATGCCGCCCATGCAGGAACT</code></td></tr>"]
    for p in lib["orientations"]:
        b = p["best"]
        rows.append("<tr><th>%s</th><td>%s</td><td>DDR-alone ΔΔG %s</td>"
                    "<td><code>%s</code></td></tr>" % (
                        p["orientation"], b["desc"], b["ddG_dr"], b["construct_dna"]))
    orient = ("<table><tr><th>取向</th><th>变体</th><th>指标</th><th>构建</th></tr>"
              + "".join(rows) + "</table>")

    comp = json.load(open(os.path.join(DATA, "ivt_compensatory_panel.json"),
                          encoding="utf-8"))
    rows = []
    for a in comp["arms"]:
        rows.append("<tr><th>%s</th><td><code>%s</code></td><td>侵占螺旋 %d 对</td>"
                    "<td><code>%s</code></td></tr>" % (
                        a["arm"], a["construct_dna"][:19] + "…" + a["spacer_dna"][:8],
                        a["inv_max_run"], a["target_rna_dna"]))
    comps = ("<table><tr><th>臂</th><th>构建</th><th>结构核查</th><th>同源靶</th></tr>"
             + "".join(rows) + "</table><div class='small'>" + comp["prediction"] + "</div>")

    html = (TPL.replace("@@DATE@@", "2026-09-02")
               .replace("@@PANEL@@", panel)
               .replace("@@ORIENT@@", orient)
               .replace("@@COMP@@", comps))
    out = os.path.join(ROOT, "demo_snapshot.html")
    open(out, "w", encoding="utf-8").write(html)
    print("输出 -> %s (%d KB)" % (out, len(html) // 1024))


if __name__ == "__main__":
    main()
