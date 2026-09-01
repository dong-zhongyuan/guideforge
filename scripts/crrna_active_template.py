"""8D4A 结合态活性模板提取(联审方案: 假结竞争模型的模板, 2026-09-01)。

从 8D4A(SuCas12a2-crRNA-靶RNA 结合态) 链 B(crRNA 41nt, 前 18nt=DR) 的三维坐标
提取碱基对(Leontis-Westhof 简化距离判据):
  WC/摆动: N1(A/G)···N3(U/C)(或反向) < 4.0 A, 且 C1'-C1' < 12.5 A;
输出 DR 区全部配对, 并按"与茎配对交叉与否"分类 stem / pseudoknot;
另报 DR-spacer 交界配对(结合态的天然跨界, 与游离态对照)。
输出: data/cas12a2_active_template.json (供假结竞争模型作活性态模板)
运行(gpumd env): python scripts/crrna_active_template.py
"""
import json
import os
from itertools import combinations

from openmm.app import PDBxFile
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
CIF = os.path.join(ROOT, "data", "8D4A.cif")

N1N3 = {"A": "N1", "G": "N1", "U": "N3", "C": "N3"}
WOB_PARTNER = {"G": "O6", "U": "N3"}  # G-U 摆动的另一判据用 N1(G)···N3(U) 已覆盖


def main():
    f = PDBxFile(CIF)
    rna = max((c for c in f.topology.chains()
               if all(r.name in ("A", "U", "G", "C") for r in c.residues())),
              key=lambda c: len(list(c.residues())))
    res = list(rna.residues())
    n = len(res)
    # 每残基原子坐标
    atoms = {}
    pos = np.array([[a.x, a.y, a.z] for a in f.positions]) * 10.0  # nm->Å
    for a in f.topology.atoms():
        if a.residue in [r.index for r in []]:
            pass
    idx_of_res = {}
    for k, r in enumerate(res):
        idx_of_res[k] = []
    for a in f.topology.atoms():
        for k, r in enumerate(res):
            if a.residue == r:
                idx_of_res[k].append((a.name, pos[a.index]))
                break

    def get(k, name):
        for nm, p in idx_of_res[k]:
            if nm == name:
                return p
        return None

    pairs = []
    for i, j in combinations(range(n), 2):
        if j - i < 3:
            continue
        bi, bj = res[i].name, res[j].name
        ok = False
        # 简化 Watson-Crick/摆动: N1(嘌呤)···N3(嘧啶) < 4.0 A
        if bi in "AG" and bj in "UC":
            pi, pj = get(i, N1N3[bi]), get(j, N1N3[bj])
            ok = pi is not None and pj is not None and np.linalg.norm(pi - pj) < 4.0
        elif bj in "AG" and bi in "UC":
            pi, pj = get(j, N1N3[bj]), get(i, N1N3[bi])
            ok = pi is not None and pj is not None and np.linalg.norm(pi - pj) < 4.0
        if not ok:
            continue
        c1i, c1j = get(i, "C1'"), get(j, "C1'")
        if c1i is None or c1j is None or np.linalg.norm(c1i - c1j) > 12.5:
            continue
        pairs.append((i + 1, j + 1, bi, bj, round(float(np.linalg.norm(c1i - c1j)), 1)))

    # 茎配对(按 5bp 茎已知: 1基 5-9 / 14-18 对应 0基 4-8/13-17)与交叉分类
    dr_len = 18
    stem0 = [(4, 17), (5, 16), (6, 15), (7, 14), (8, 13)]
    stem1 = {(i + 1, j + 1) for i, j in stem0}

    def crosses(p, q):
        a, b = p
        c, d = q
        return (a < c < b < d) or (c < a < d < b)

    out = []
    for p1, p2, bi, bj, c1d in pairs:
        in_dr = p2 <= dr_len
        cross_pk = any(crosses((p1, p2), (s0 + 1, s1 + 1)) for s0, s1 in stem0
                       if (s0 + 1, s1 + 1) in stem1)
        region = "DR" if in_dr else ("DR-spacer" if p1 <= dr_len < p2 else "spacer")
        kind = "pseudoknot-crossing" if cross_pk else "nested"
        out.append({"i_1based": p1, "j_1based": p2, "bases": bi + "-" + bj,
                    "region": region, "class": kind, "C1p_dist_A": c1d})

    dr_pairs = [x for x in out if x["region"] == "DR"]
    print("链B(crRNA) %dnt, 几何判据提取配对 %d 对; DR 区内 %d 对:" % (n, len(out), len(dr_pairs)))
    for x in dr_pairs:
        print("  %2d-%2d %s %-6s %s  C1'=%.1f" % (
            x["i_1based"], x["j_1based"], x["bases"], x["class"], x["region"], x["C1p_dist_A"]))
    xs_pairs = [x for x in out if x["region"] == "DR-spacer"]
    print("DR-spacer 天然交界配对(结合态): %d 对: %s" % (
        len(xs_pairs), [(x["i_1based"], x["j_1based"]) for x in xs_pairs]))

    payload = {"source": "data/8D4A.cif 链B(crRNA) 三维几何提取",
               "criteria": "N1(嘌呤)-N3(嘧啶) <4.0A + C1'-C1' <12.5A (简化 LW 判据)",
               "dr_len": dr_len, "stem_pairs_1based": sorted(list(x) for x in stem1),
               "pairs": out,
               "note": "class=pseudoknot-crossing 指与 5bp 茎交叉的配对(假结元素候选); "
                       "供假结竞争模型作活性态模板"}
    path = os.path.join(ROOT, "data", "cas12a2_active_template.json")
    json.dump(payload, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("输出 -> %s" % path)


if __name__ == "__main__":
    main()
