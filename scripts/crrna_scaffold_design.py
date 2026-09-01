"""crRNA 骨架变体库生成与打分（Cas12a2 crRNA 骨架优化 · 核心计算模块）

定位:
  效应子 Cas12a2（RNA 靶向）；设计变量只有 crRNA 骨架（DR 茎环），spacer 向导序列固定不变。
  本模块做"候选库压缩"：生成结构保持型骨架突变体并按透明指标排序，输出 TOP-K 供湿实验
  初筛；不做活性预测（活性/特异性结论由体外生化与细胞实验闭环给出，见研究方案）。

输入:
  - 骨架/PFS 口径: configs/scaffold_registry.json（经 scaffold_registry.py 只读），默认 cas12a2
  - 固定 spacer: --spacer 传入（DNA/RNA 字母均可，T 自动转 U）
  - 蛋白接触约束（可选）: data/8D4A_dr_contacts.json；缺失时自动从 RCSB 下载
    8D4A.cif（Cas12a2 四元复合物, Bravo 2023）计算并缓存，可用 --no-contacts 关闭

策略（与申报书一致，--strategy 控制）:
  A. 结构保持局部突变扫描: DR 全部单点突变 + 抽样双点突变，ViennaRNA 折叠过滤
  A2. Struct2SeQ 逆折叠提案(--use-struct2seq): 以 WT DR 茎环为靶结构、WT 序列 upweight,
      子进程调 rnet 环境跑 toolbox/Struct2SeQ, 经 RNet-SS Jaccard 筛选
  A4. 茎区共变枚举(--use-covariation): 对 WT DR 折叠的每个茎配对做保持配对的双侧补偿
      突变(标准配对 GC/CG/AU/UA/GU/UG 互换), 单对全枚举 + 双对组合抽样
      (依据 Dmytrenko 2023 Fig.1c: 跨家族茎保守/loop 可变, 共变是同源 DR
      保持茎环拓扑的天然机制)
  筛选增强: --rnet-screen 对 WT+TOP-N 候选预测 RNet SHAPE 图谱并算与 WT 结构的一致性
      (Science 2026 aeg6829 的"先模拟实验再筛选"口径, eterna score 分量)
  B. 模拟退火: 以同一打分为目标，突变提案偏向 loop/非蛋白接触位点（--sa-steps）

硬过滤（默认值均可 CLI 调）:
  - 全长 construct 的 MFE 结构与 WT 的 base-pair 距离 <= --max-bp-dist（默认 4）
  - spacer 区平均未配对概率不低于 WT 值 - --spacer-unpaired-margin（默认 0.10，
    spacer 必须保持游离；因 spacer 固定，与 WT 相对比较即可排除 spacer 自身发卡的影响）
  - DR 变体不含 TTTT/GGGG（合成与转录纪律，与旧链一致）
  - 加工位点保护(--no-protect-processing 可关): Cas12a2 pre-crRNA 加工切点在 DR 3' 端
    紧下游(Dmytrenko 2023 ED Fig.3, 比 Cas12a 下游 1nt; 成熟边界由结构口径 18nt DR 锚定);
    DR 3' 末端 3nt + 交界后首位的配对状态须与 WT 一致, 防止成熟 crRNA 长度漂移

打分（透明启发式，仅用于候选排序，不是活性预测；权重 CLI 可调）:
  score = -w_bp*bp_dist - w_ddg*max(ddG,0) - w_contact*contact_sum
          - w_ens*max(d_ens,0) + w_hbond*hbond_net - w_cons3*cons3_frac
          + w_fold*dp_fold + w_seed*d_seed + w_stab*max(-ddG,0)
  其中 hbond_loss 为 8D4A 极性接触的碱基替换互补损失（几何不变近似）;
  hbond_net = hbond_gain - hbond_loss(有符号: 找回互补记正; gain 仅同名原子可算,
  新碱基独有极性原子无 WT 坐标, 系统性低估);
  dp_fold = P(茎完整)_variant - P(茎完整)_WT, P 用全长 bpp 矩阵茎配对概率乘积
  (独立性近似, 忽略配对间耦合; 机制依据 Bush 2023: 骨架正确折叠率→有效分子比例;
  代理指标, 未实验标定);
  d_seed = 种子区游离度变化(spacer 3' 端 seed_len=7nt 平均未配对概率之差;
  依据 Bravo 2023: 3' 端 7 碱基有序溶剂暴露为靶 RNA 结合种子; Liao 2018: DR
  下游结构缠住 spacer 抑制切割——游离态口径: 预测降低溶液态自缠结/加速种子暴露,
  非结合态活性保证);
  w_stab*max(-ddG,0) 为茎稳定化奖赏(共变强化轨道: 变体比 WT 更稳给正分;
  依据 Sudhakar 2023(修正后因果: 预折叠 crRNA 诱导蛋白构象变化驱动组装)的
  组装/均一性方向 + Teng 2019 4n96 存在性先例; 无"更稳->活性更高"直接实验证据,
  表述限定在组装效率/折叠均一性, 不写提高活性; 权重与 w_cons3 对称, 解决
  "纯罚分制下增强候选被保守窗+接触双重压制浮不出"的结构性问题);
  cons3_frac 为落在 DR 3' 保守窗的突变位点比例(Dmytrenko 2023 Fig.1c: 跨家族 3' 端
  高度保守、loop 可变, 故 3' 窗内突变给显式惩罚; 窗大小 --cons3-window 为项目设定,
  文献依据为定性结论; 功能佐证: Zhang 2025 DR 3' 端化学修饰可逆调控 Cas12a 活性);
  另报告 ddG、spacer_unpaired、ens_diversity、hbond_preserved 供人工挑选。

输出: <out>.variants.csv（全部打分行）/ <out>.top.json（TOP-K + WT 参考 + 参数血缘）/
      <out>.top.fasta（完整 construct，DNA 字母表，可直接送合成）。

用法:
  python crrna_scaffold_design.py --spacer <21nt> --topk 12 --out-prefix run1
"""
import argparse
import csv
import hashlib
import json
import os
import urllib.request

import numpy as np
import RNA

from scaffold_registry import get_entry

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
CONTACTS_JSON_LEGACY = os.path.join(ROOT, 'data', '8D4A_dr_contacts.json')


def contacts_json_path(dr_seq):
    """接触缓存按 DR 序列分文件（多 DR 并存时互不覆盖）。"""
    tag = hashlib.sha256(dr_seq.encode()).hexdigest()[:8]
    return os.path.join(ROOT, 'data', f'8D4A_dr_contacts_{tag}.json')
PDB_CIF = os.path.join(ROOT, 'data', '8D4A.cif')
PDB_URL = 'https://files.rcsb.org/download/8D4A.cif'
PDB_ID = '8D4A'
CONTACT_CUTOFF_A = 5.0

BASES = 'ACGU'
_fold_cache = {}

HBOND_CUTOFF_A = 3.5

# RNA 碱基极性原子表（氢键供体 d / 受体 a;2'OH 各碱基相同不计）
RNA_BASE_POLAR = {
    'A': {'N1': 'a', 'N3': 'a', 'N6': 'd', 'N7': 'a'},
    'U': {'N3': 'd', 'O2': 'a', 'O4': 'a'},
    'G': {'N1': 'd', 'N2': 'd', 'N3': 'a', 'O6': 'a', 'N7': 'a'},
    'C': {'N3': 'a', 'N4': 'd', 'O2': 'a'},
}

# 蛋白原子供受体的分类（b = 兼性, 如羟基/咪唑;主链 N=d, O=a 对所有残基生效）
PROT_SIDECHAIN_CLASS = {
    'ARG': {'NE': 'd', 'NH1': 'd', 'NH2': 'd'},
    'LYS': {'NZ': 'd'}, 'TRP': {'NE1': 'd'},
    'ASP': {'OD1': 'a', 'OD2': 'a'}, 'GLU': {'OE1': 'a', 'OE2': 'a'},
    'ASN': {'ND2': 'd', 'OD1': 'a'}, 'GLN': {'NE2': 'd', 'OE1': 'a'},
    'HIS': {'ND1': 'b', 'NE2': 'b'}, 'SER': {'OG': 'b'},
    'THR': {'OG1': 'b'}, 'TYR': {'OH': 'b'}, 'CYS': {'SG': 'b'},
}


def prot_atom_class(resname, atom_name):
    """蛋白原子供受体分类: 主链 N/O + 侧链表; 非极性返回 None。"""
    if atom_name == 'N':
        return 'd'
    if atom_name == 'O':
        return 'a'
    return PROT_SIDECHAIN_CLASS.get(resname, {}).get(atom_name)


def classes_complementary(prot_class, rna_class):
    if prot_class == 'b':
        return rna_class in ('d', 'a')
    return (prot_class == 'd' and rna_class == 'a') or \
           (prot_class == 'a' and rna_class == 'd')


def to_rna(seq):
    return str(seq).upper().replace('T', 'U')


def to_dna(seq):
    return str(seq).upper().replace('U', 'T')


def fold(seq):
    """MFE 折叠（带缓存）。返回 (dotbracket, mfe_kcal)。"""
    if seq not in _fold_cache:
        _fold_cache[seq] = RNA.fold(seq)
    return _fold_cache[seq]


def pf_stats(seq, dr_len, seed_len=7):
    """配分函数统计: (ensemble_diversity, spacer_mean_unpaired, seed_unpaired)。
    seed_unpaired = spacer 3' 端 seed_len nt 平均未配对概率
    (依据 Bravo 2023: Cas12a2 crRNA spacer 3' 端 7 碱基有序且溶剂暴露,
    为靶 RNA 结合的种子区——种子区被自缠结压住即损失可及性)。"""
    fc = RNA.fold_compound(seq)
    fc.pf()
    diversity = fc.mean_bp_distance()
    n = len(seq)
    P = np.array([list(r) for r in fc.bpp()])[:n + 1, :n + 1]
    paired = P.sum(axis=0) + P.sum(axis=1)
    unpaired = 1.0 - paired[1:]
    spacer_up = float(unpaired[dr_len:].mean()) if len(seq) > dr_len else 1.0
    seed_up = float(unpaired[max(dr_len, n - seed_len):].mean()) if len(seq) > dr_len else 1.0
    return float(diversity), spacer_up, seed_up


def wt_reference(dr, spacer, seed_len=7):
    """WT 全长参考: MFE 结构/能量、系综多样性、spacer 游离度、种子区游离度、
    DR 单独折叠、茎完整概率。"""
    full = dr + spacer
    ss, mfe = fold(full)
    diversity, spacer_up, seed_up = pf_stats(full, len(dr), seed_len)
    dr_ss, dr_mfe = fold(dr)
    p_fold = stem_intact_prob(full, stem_pairs_of(dr_ss))
    return {'full_len': len(full), 'mfe_struct': ss, 'mfe_kcal': round(mfe, 2),
            'ens_diversity': round(diversity, 2), 'spacer_mean_unpaired': round(spacer_up, 3),
            'seed_unpaired': round(seed_up, 3),
            'dr_only_struct': dr_ss, 'dr_only_mfe_kcal': round(dr_mfe, 2),
            'p_fold': p_fold}


def stem_positions(dr_struct, dr_len):
    """WT DR 折叠中成茎的位置集合（0-based）。"""
    return {i for i in range(min(dr_len, len(dr_struct))) if dr_struct[i] in '()'}


def compute_contacts(expected_dr, cif_path=PDB_CIF, cutoff=CONTACT_CUTOFF_A):
    """从 PDB mmCIF 提取 DR 每个核苷酸与蛋白链的接触（<=cutoff Å 重原子计数）。

    crRNA 链靠序列锚定: 5' 端必须等于 expected_dr（8D4A 里还有一条 28nt 靶 RNA,
    按长度选链会选错, 必须按序列）。DR = 该链 5' 端前 dr_len 个残基。
    返回 dict 并写缓存 JSON。
    """
    from Bio.PDB import MMCIFParser, NeighborSearch
    dr_len = len(expected_dr)
    if not os.path.isfile(cif_path):
        os.makedirs(os.path.dirname(cif_path), exist_ok=True)
        print(f'下载 {PDB_ID}.cif: {PDB_URL}')
        urllib.request.urlretrieve(PDB_URL, cif_path)
    structure = MMCIFParser(QUIET=True).get_structure(PDB_ID, cif_path)
    model = structure[0]
    aa_names = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
    protein_atoms, rna_chains = [], []
    for chain in model:
        residues = [r for r in chain if r.id[0] == ' ']
        if not residues:
            continue
        names = {r.resname.strip() for r in residues}
        if names <= aa_names:
            protein_atoms.extend([a for r in residues for a in r])
        elif names <= set('ACGU') and len(residues) >= dr_len:
            seq = ''.join(r.resname.strip() for r in residues)
            rna_chains.append((chain.id, seq, residues))
    if not protein_atoms:
        raise ValueError(f'{cif_path} 中未识别到蛋白链')
    # 定位 crRNA 链: 候选 = (错配数, expected 5' 悬垂, 链 5' 偏移, ...), 取最小者;
    # 覆盖 Zeng2026 DR 这类与结构 DR 相差 5' 延伸/个位替换的情况。
    # 悬垂/错配位点的极性接触不可转移, 置空或剔除。
    best = None
    for cid, s, rs in rna_chains:
        for eoff in (0, 1, 2):
            n = len(expected_dr) - eoff
            for coff in (0, 1, 2):
                if len(rs) - coff < n:
                    continue
                seg = s[coff:coff + n]
                mm = [k for k in range(n) if seg[k] != expected_dr[eoff + k]]
                if len(mm) <= 2 and (best is None or (len(mm), eoff, coff) < best[:3]):
                    best = (len(mm), eoff, coff, cid, rs[coff:coff + n], seg, mm)
    if best is None:
        preview = {cid: s[:25] for cid, s, _ in rna_chains}
        raise ValueError(
            f"{cif_path} 中无匹配(或可对齐)DR 的 crRNA 链; 候选链: {preview}")
    _, eoff, coff, chain_id, crna_residues, crna_seq, mm = best
    exact = (eoff == 0 and coff == 0 and not mm)
    transfer_info = None if exact else {
        'method': 'aligned_transfer',
        'expected_5p_overhang': eoff, 'chain_offset': coff,
        'mismatch_positions_1based': [eoff + k + 1 for k in mm],
        'note': "接触按序列对齐从结构 DR 转移; 错配位点与 5' 悬垂位点的极性接触已剔除/置空"}
    ns = NeighborSearch(protein_atoms)
    contacts, atom_contacts = {}, {}
    for pos0 in range(dr_len):
        i = pos0 - eoff
        res = crna_residues[i] if 0 <= i < len(crna_residues) else None
        if res is None:
            contacts[str(pos0 + 1)] = {'n_contact_residues': 0, 'residues': []}
            atom_contacts[str(pos0 + 1)] = []
            continue
        hits = set()
        polar = []
        base = res.resname.strip()
        for atom in res:
            for patom in ns.search(atom.coord, cutoff):
                pr = patom.get_parent()
                hits.add(f'{pr.get_parent().id}:{pr.resname.strip()}{pr.id[1]}')
                rna_class = RNA_BASE_POLAR.get(base, {}).get(atom.name.strip())
                prot_class = prot_atom_class(pr.resname.strip(), patom.name.strip())
                dist = float(np.linalg.norm(atom.coord - patom.coord))
                if rna_class and prot_class and dist <= HBOND_CUTOFF_A:
                    polar.append({'rna_atom': atom.name.strip(), 'rna_class': rna_class,
                                  'prot_atom': patom.name.strip(), 'prot_class': prot_class,
                                  'prot_res': f'{pr.get_parent().id}:{pr.resname.strip()}{pr.id[1]}',
                                  'dist': round(dist, 2)})
        contacts[str(pos0 + 1)] = {'n_contact_residues': len(hits), 'residues': sorted(hits)}
        atom_contacts[str(pos0 + 1)] = polar
    if transfer_info:
        for p in transfer_info['mismatch_positions_1based']:
            atom_contacts[str(p)] = []  # 碱基不同, 极性接触化学不可转移
    payload = {'pdb_id': PDB_ID, 'cutoff_A': cutoff, 'hbond_cutoff_A': HBOND_CUTOFF_A,
               'dr_len': dr_len, 'dr_seq': expected_dr, 'crna_chain': chain_id,
               'crna_chain_seq': crna_seq, 'transfer': transfer_info,
               'cif_sha256': hashlib.sha256(open(cif_path, 'rb').read()).hexdigest(),
               'note': 'Cas12a2 四元复合物(Bravo 2023); crRNA 链按 5\' 端 DR 序列锚定; '
                       'contacts 为 <=cutoff 内蛋白残基数(重原子), '
                       'atom_contacts 为 <=3.5A 的极性原子对(供受体分类, 用于碱基替换互补评估)',
               'contacts': contacts, 'atom_contacts': atom_contacts}
    out_path = contacts_json_path(expected_dr)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    return payload


def load_contacts(dr_seq, enabled=True):
    """加载/构建蛋白接触表；不可用时返回 None（接触项打分为 0 并在输出中声明）。"""
    if not enabled:
        return None
    payload = None
    out_path = contacts_json_path(dr_seq)
    if os.path.isfile(out_path):
        with open(out_path, encoding='utf-8') as fh:
            payload = json.load(fh)
    if payload is None or 'atom_contacts' not in payload or payload.get('dr_seq') != dr_seq:
        payload = compute_contacts(expected_dr=dr_seq)
    counts = {int(k): v['n_contact_residues'] for k, v in payload['contacts'].items()}
    return {'max': max(counts.values()) or 1, 'counts': counts,
            'atom_contacts': {int(k): v for k, v in payload['atom_contacts'].items()},
            'source': payload['pdb_id'], 'transfer': payload.get('transfer')}


def hbond_loss(dr, seq, mut_positions, contact):
    """碱基替换的氢键互补损失(几何不变近似: 只替换碱基类型, 原子几何不动)。

    对每个突变位点: WT 碱基在 8D4A 中的极性接触里, 变体碱基在同名原子上仍呈现
    互补供/受体的比例记为 preserved; 该位点损失 = 1 - preserved。无接触位点损失 0。
    返回值 (loss_sum, preserved_mean)。"""
    if contact is None or not mut_positions:
        return 0.0, 1.0
    loss, fracs = 0.0, []
    for p in mut_positions:
        wt_base, var_base = dr[p - 1], seq[p - 1]
        atom_contacts = contact['atom_contacts'].get(p, [])
        if not atom_contacts:
            continue
        kept = 0
        for ac in atom_contacts:
            var_class = RNA_BASE_POLAR.get(var_base, {}).get(ac['rna_atom'])
            if var_class and classes_complementary(ac['prot_class'], var_class):
                kept += 1
        frac = kept / len(atom_contacts)
        fracs.append(frac)
        loss += 1.0 - frac
    return loss, (sum(fracs) / len(fracs) if fracs else 1.0)


def hbond_gain(dr, seq, mut_positions, contact):
    """找回互补(正向项): 突变位点 WT 接触原子对中, WT 碱基不互补而变体碱基互补的对数。
    仅同名原子几何不变近似——变体碱基独有的极性原子无 WT 坐标, 系统性低估(保守方向)。"""
    if contact is None or not mut_positions:
        return 0
    gain = 0
    for p in mut_positions:
        wt_base, var_base = dr[p - 1], seq[p - 1]
        for ac in contact['atom_contacts'].get(p, []):
            wt_class = RNA_BASE_POLAR.get(wt_base, {}).get(ac['rna_atom'])
            var_class = RNA_BASE_POLAR.get(var_base, {}).get(ac['rna_atom'])
            if wt_class is None and var_class \
                    and classes_complementary(ac['prot_class'], var_class):
                gain += 1
    return gain


def stem_pairs_of(dr_ss):
    """DR 折叠点括号 → 茎配对列表 [(i,j)](0-based, DR 内坐标)。"""
    stack, pairs = [], []
    for i, ch in enumerate(dr_ss):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            pairs.append((stack.pop(), i))
    return pairs


def stem_intact_prob(full_seq, dr_pairs):
    """P(茎完整): 全长配分函数 bpp 矩阵上, WT DR 各茎配对概率的乘积(独立性近似,
    忽略配对间耦合——作变体间相对比较有效)。dr_pairs 为 DR 内坐标, 偏移至全长。"""
    fc = RNA.fold_compound(full_seq)
    fc.pf()
    bpp = fc.bpp()
    p = 1.0
    for i, j in dr_pairs:
        p *= bpp[i + 1][j + 1]
    return p


def contact_sum(mut_positions, contact):
    """变体的接触代价: 突变位点的蛋白接触计数归一化求和。"""
    if contact is None:
        return 0.0
    return sum(contact['counts'].get(p, 0) / contact['max'] for p in mut_positions)


def mut_desc(dr, seq):
    diff = [f'{dr[i]}{i + 1}{seq[i]}' for i in range(len(dr)) if dr[i] != seq[i]]
    return 'WT' if not diff else '+'.join(diff), [int(d[1:-1]) for d in diff]


def proc_window_ok(wt_ss, ss, dr_len, window):
    """加工位点保护(Dmytrenko 2023 ED Fig.3): Cas12a2 加工切点在 DR 3' 端紧下游,
    成熟边界由结构口径 18nt DR 锚定。要求 DR 3' 末端 window-1 个碱基 + 交界后首位
    在变体全长 MFE 中的配对状态与 WT 一致, 防止成熟 crRNA 长度漂移。"""
    idx = list(range(max(dr_len - (window - 1), 0), min(dr_len + 1, len(ss))))
    return all((wt_ss[i] != '.') == (ss[i] != '.') for i in idx)


STEM_PAIR_SETS = ('GC', 'CG', 'AU', 'UA', 'GU', 'UG')


def strategy_a4(dr, args, rng):
    """策略A4: 茎区共变枚举。对 WT DR 折叠的每个茎配对, 枚举保持配对的双侧补偿突变
    (标准配对互换), 单对全枚举 + 双对组合抽样。茎配对来源 = DR 单独折叠的点括号。"""
    ss = fold(dr)[0]
    pairs, stack = [], []
    for i, ch in enumerate(ss):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            pairs.append((stack.pop(), i))

    def apply_subs(subs):
        s = list(dr)
        for (i, j), (bi, bj) in subs.items():
            s[i], s[j] = bi, bj
        return ''.join(s)

    variants = []
    for i, j in pairs:  # 单对共变: 全枚举
        cur = dr[i] + dr[j]
        for p in STEM_PAIR_SETS:
            if p != cur:
                variants.append(apply_subs({(i, j): (p[0], p[1])}))
    combos = [(a, b) for x, a in enumerate(pairs) for b in pairs[x + 1:]]
    rng.shuffle(combos)
    for c in combos[:args.cov_double]:  # 双对共变: 抽样
        subs = {}
        for i, j in c:
            cur = dr[i] + dr[j]
            alt = [p for p in STEM_PAIR_SETS if p != cur]
            ch = alt[rng.integers(len(alt))]
            subs[(i, j)] = (ch[0], ch[1])
        variants.append(apply_subs(subs))
    return [v for v in variants if v != dr]


def score_variant(dr, seq, wt, spacer, args, contact, stem_pos):
    """折叠并打分一个 DR 变体。返回指标 dict（含 pass/score）。"""
    full = seq + spacer
    ss, mfe = fold(full)
    bp_dist = RNA.bp_distance(wt['mfe_struct'], ss)
    diversity, spacer_up, seed_up = pf_stats(full, len(dr), args.seed_len)
    ddg = mfe - wt['mfe_kcal']
    d_ens = diversity - wt['ens_diversity']
    desc, mut_pos = mut_desc(dr, seq)
    csum = contact_sum(mut_pos, contact)
    hloss, hfrac = hbond_loss(dr, seq, mut_pos, contact)
    hgain = hbond_gain(dr, seq, mut_pos, contact)
    hnet = hgain - hloss
    p_fold = stem_intact_prob(full, stem_pairs_of(wt['dr_only_struct']))
    dp_fold = p_fold - wt['p_fold']
    d_seed = seed_up - wt['seed_unpaired']
    dr_ss, dr_mfe = fold(seq)
    ddg_dr = dr_mfe - wt['dr_only_mfe_kcal']  # DR 单独折叠 ΔΔG(茎稳定化语义更干净的口径)
    ok = (bp_dist <= args.max_bp_dist
          and spacer_up >= wt['spacer_mean_unpaired'] - args.spacer_unpaired_margin
          and 'TTTT' not in to_dna(seq) and 'GGGG' not in to_dna(seq)
          and (args.no_protect_processing
               or proc_window_ok(wt['mfe_struct'], ss, len(dr), args.proc_window)))
    cons3_n = sum(1 for p in mut_pos if p > len(dr) - args.cons3_window)
    cons3_frac = cons3_n / max(args.cons3_window, 1)
    stab_dd = ddg_dr if args.stab_dr_only else ddg
    score = (-args.w_bp * bp_dist - args.w_ddg * max(ddg, 0.0)
             - args.w_contact * csum - args.w_ens * max(d_ens, 0.0)
             + args.w_hbond * hnet - args.w_cons3 * cons3_frac
             + args.w_fold * dp_fold + args.w_seed * d_seed
             + args.w_stab * max(-stab_dd, 0.0))
    return {'desc': desc, 'mut_positions': mut_pos, 'n_mut': len(mut_pos),
            'dr_seq': to_dna(seq), 'construct_dna': to_dna(full),
            'mfe_struct': ss, 'mfe_kcal': round(mfe, 2), 'ddG': round(ddg, 2),
            'dr_only_mfe_kcal': round(dr_mfe, 2), 'ddG_dr': round(ddg_dr, 2),
            'bp_dist': bp_dist, 'spacer_unpaired': round(spacer_up, 3),
            'ens_diversity': round(diversity, 2), 'd_ens': round(d_ens, 2),
            'contact_sum': round(csum, 3),
            'hbond_loss': round(hloss, 3), 'hbond_preserved': round(hfrac, 3),
            'hbond_gain': hgain, 'hbond_net': round(hnet, 3),
            'p_fold': round(p_fold, 5), 'dp_fold': round(dp_fold, 5),
            'seed_up': round(seed_up, 3), 'd_seed': round(d_seed, 3),
            'cons3_frac': round(cons3_frac, 3),
            'proc_ok': bool(args.no_protect_processing
                            or proc_window_ok(wt['mfe_struct'], ss, len(dr), args.proc_window)),
            'pareto': None, 'advantage': None,
            'shape_cons': None,
            'passed': bool(ok), 'score': round(score, 4)}


def pareto_flag(rows):
    """Pareto 最优标记(通过硬过滤的变体)。目标: DR 稳定(-ddG_dr)↑, spacer 游离↑,
    结构距离↓, 接触+氢键损失↓。总分排序会埋没"单维占优但有代价"的权衡型候选,
    Pareto 前沿保证它们在产出中可见(每条都有不可替代的优势维度)。"""
    passed = [r for r in rows if r['passed']]
    if not passed:
        return
    obj = np.array([[-r['ddG_dr'], r['spacer_unpaired'],
                     -r['bp_dist'], -(r['contact_sum'] + r['hbond_loss'])]
                    for r in passed])
    for k, r in enumerate(passed):
        dominated = np.any(np.all(obj >= obj[k], axis=1)
                           & np.any(obj > obj[k], axis=1))
        r['pareto'] = bool(not dominated)


def advantage_str(r, wt, seed_len=7):
    """候选优势论述(自动注释): 只陈述有文献/机制依据的定向改进与保持项。"""
    adv = []
    if r['d_seed'] > 0.02:
        adv.append(f"种子区(3'端{seed_len}nt)游离度 {wt['seed_unpaired']:.2f}→{r['seed_up']:.2f}")
    if r['spacer_unpaired'] - wt['spacer_mean_unpaired'] > 0.02:
        adv.append(f"spacer 游离度 {wt['spacer_mean_unpaired']:.2f}→{r['spacer_unpaired']:.2f}")
    if r['ddG_dr'] < -0.1:
        adv.append(f"DR 茎稳定化 {r['ddG_dr']:+.1f} kcal/mol")
    if r['dp_fold'] > 0.01:
        adv.append(f"茎完整概率 {r['dp_fold']:+.2f}")
    if r['hbond_net'] > 0:
        adv.append(f"氢键净找回 +{r['hbond_net']:.0f}")
    if r['bp_dist'] == 0 and r['hbond_loss'] == 0.0:
        adv.append("结构与蛋白极性接触保持")
    return '; '.join(adv) if adv else '与 WT 无显著差异'


def strategy_a(dr, args, rng):
    """结构保持局部突变扫描: 全部单点 + 分层抽样双点。"""
    variants = []
    for i in range(len(dr)):
        for b in BASES:
            if b != dr[i]:
                variants.append(dr[:i] + b + dr[i + 1:])
    positions = list(range(len(dr)))
    seen = set()
    while len(variants) < len(dr) * 3 + args.n_double:
        i, j = sorted(rng.choice(positions, size=2, replace=False).tolist())
        bi, bj = rng.choice(list(BASES)), rng.choice(list(BASES))
        if bi == dr[i] or bj == dr[j]:
            continue
        seq = dr[:i] + bi + dr[i + 1:j] + bj + dr[j + 1:]
        if seq in seen:
            continue
        seen.add(seq)
        variants.append(seq)
    return variants


def proposal_weights(dr, contact, stem_pos):
    """SA 提案权重: 茎位与高蛋白接触位点降低被扰动概率。"""
    w = np.ones(len(dr))
    for i in range(len(dr)):
        if i in stem_pos:
            w[i] *= 0.3
        if contact is not None:
            w[i] *= 1.0 - 0.5 * contact['counts'].get(i + 1, 0) / contact['max']
    return w / w.sum()


def shape_consistency(profile, target_db, threshold=0.5):
    """SHAPE-结构一致性（OpenKnot/eterna score 口径）: 靶结构的配对/非配对状态与
    预测反应性（反应性 < threshold 视为配对）一致的碱基比例 ×100。
    靶结构无假结, 故只取 eterna_score 分量, 不含 crossed-pair 分量。"""
    n = len(target_db)
    agree = sum(1 for c, r in zip(target_db, profile)
                if (c != '.') == (r < threshold))
    return 100.0 * agree / n


def rnet_shape_screen(entries, target_db, args):
    """对 [(key, construct_rna)] 预测 RNet SHAPE 并算与 target_db 的一致性。
    返回 {key: score}; 子进程失败返回 None（不阻断主流程, 输出中标注）。"""
    import subprocess
    workdir = os.path.join(ROOT, 'data', '_rnet_shape_run')
    os.makedirs(workdir, exist_ok=True)
    in_csv = os.path.join(workdir, 'in.csv')
    out_csv = os.path.join(workdir, 'out.csv')
    with open(in_csv, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['id', 'sequence'])
        for key, seq in entries:
            w.writerow([key, seq])
    probe = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'crrna_rnet_shape.py')
    try:
        proc = subprocess.run([args.s2s_python, probe, '--csv', in_csv,
                               '--out', out_csv],
                              capture_output=True, text=True, encoding='utf-8',
                              errors='replace', timeout=args.s2s_timeout)
        if proc.returncode != 0:
            print(f"[rnet-screen] 子进程失败, 跳过: {(proc.stderr or proc.stdout or '')[-300:]}")
            return None
        result = {}
        with open(out_csv, encoding='utf-8') as fh:
            for row in csv.DictReader(fh):
                profile = [float(x) for x in row['shape'].split(';')]
                result[row['id']] = round(shape_consistency(profile, target_db), 1)
        return result
    finally:
        for f in os.listdir(workdir):
            os.remove(os.path.join(workdir, f))
        os.rmdir(workdir)


def strategy_a2(dr, wt, args):
    """策略A2: Struct2SeQ 逆折叠提案（子进程调 rnet 环境）。
    以 WT DR 茎环结构为靶、WT 序列 upweight 生成候选，由官方 Jaccard/rescue 机制筛选。"""
    import subprocess
    workdir = os.path.join(ROOT, 'data', '_s2s_run')
    os.makedirs(workdir, exist_ok=True)
    target_csv = os.path.join(workdir, 'target.csv')
    with open(target_csv, 'w', encoding='utf-8') as fh:
        fh.write('Title,Dot-bracket,wild_type_sequence\n')
        fh.write(f"dr_wt,{wt['dr_only_struct']},{dr}\n")
    wrapper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'crrna_struct2seq_gen.py')
    proc = subprocess.run(
        [args.s2s_python, wrapper, '--target-csv', target_csv,
         '--out-dir', workdir, '--up-bias', str(args.s2s_up_bias)],
        capture_output=True, text=True, encoding='utf-8', errors='replace',
        timeout=args.s2s_timeout)
    if proc.returncode != 0:
        print(f"[struct2seq] 子进程失败, 跳过 A2(策略A/B 不受影响): "
              f"{(proc.stderr or proc.stdout or '')[-300:]}")
        return []
    seqs = set()
    with open(os.path.join(workdir, 'results.csv'), encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            s = row['sequence'].strip().upper().replace('T', 'U')
            if len(s) == len(dr) and set(s) <= set(BASES) and s != dr:
                seqs.add(s)
    print(f'[struct2seq] 逆折叠提案 {len(seqs)} 条唯一 DR 变体(经 RNet-SS Jaccard 筛选)')
    for f in os.listdir(workdir):
        os.remove(os.path.join(workdir, f))
    os.rmdir(workdir)
    return sorted(seqs)


def strategy_a3(dr, args):
    """策略A3: gRNAde 三维几何逆折叠提案（子进程调 rnet 环境，8D4A 结合态几何）。
    仅当注册表 DR 与 8D4A 结构 DR 完全一致时可用（几何按序列锚定）;
    其他 DR（如 cas12a2_zeng2026）无直接结合态几何，跳过并说明。"""
    pdb_dr = 'AUUUCUACUAUUGUAGAU'  # 8D4A B 链 5' 端（结构口径, 见 compute_contacts 锚定）
    if dr != pdb_dr:
        print(f'[grnade] 当前 DR ({dr}) 与 8D4A 结构 DR ({pdb_dr}) 不一致, '
              f'无直接结合态几何, 跳过 A3')
        return []
    import subprocess
    pdb_spacer = 'UGGAGCAACACCUGAAGAAGGCU'  # 8D4A 结构原生 spacer(几何上下文)
    target_db = fold(dr + pdb_spacer)[0]  # 结构链的 WT 二级结构(RNet-SS 复核一致)
    workdir = os.path.join(ROOT, 'data', '_grnade_run')
    os.makedirs(workdir, exist_ok=True)
    out_csv = os.path.join(workdir, 'grnade.csv')
    wrapper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'crrna_grnade_gen.py')
    pdb = os.path.join(ROOT, 'data', '8D4A_crna_chainB.pdb')
    proc = subprocess.run(
        [args.s2s_python, wrapper, '--pdb', pdb, '--dr-len', str(len(dr)),
         '--n-batches', str(args.grnade_batches), '--seed', str(args.seed),
         '--target-sec-struct', target_db,
         '--out', out_csv],
        capture_output=True, text=True, encoding='utf-8', errors='replace',
        timeout=args.s2s_timeout)
    if proc.returncode != 0:
        print(f"[grnade] 子进程失败, 跳过 A3(其他策略不受影响): "
              f"{(proc.stderr or proc.stdout or '')[-300:]}")
        return []
    seqs = set()
    with open(out_csv, encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            s = row['sequence'].strip().upper().replace('T', 'U')[:len(dr)]
            if len(s) == len(dr) and set(s) <= set(BASES) and s != dr:
                seqs.add(s)
    print(f'[grnade] 三维几何提案 {len(seqs)} 条唯一 DR 变体(8D4A 结合态)')
    for f in os.listdir(workdir):
        os.remove(os.path.join(workdir, f))
    os.rmdir(workdir)
    return sorted(seqs)


def strategy_b(dr, args, rng, wt, spacer, contact, stem_pos):
    """模拟退火: 以 score_variant 的 score 为目标, 记录访问过的唯一变体。"""
    weights = proposal_weights(dr, contact, stem_pos)
    current = dr
    cur = score_variant(dr, current, wt, spacer, args, contact, stem_pos)
    best, visited = cur['score'], {dr: cur}
    t0, t1 = 1.0, 0.05
    for step in range(args.sa_steps):
        temp = t0 * (t1 / t0) ** (step / max(args.sa_steps - 1, 1))
        i = rng.choice(len(dr), p=weights)
        b = BASES[rng.integers(4)]
        if b == current[i]:
            continue
        cand = current[:i] + b + current[i + 1:]
        row = visited.get(cand) or score_variant(dr, cand, wt, spacer, args, contact, stem_pos)
        visited[cand] = row
        delta = row['score'] - cur['score']
        if delta >= 0 or rng.random() < np.exp(delta / temp):
            current, cur = cand, row
            best = max(best, row['score'])
    return [s for s in visited if s != dr]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--spacer', required=True, help='固定 spacer 序列(17-25nt ACGT/U), 全程不变')
    ap.add_argument('--effector', default='cas12a2')
    ap.add_argument('--registry', default=None, help='注册表路径(默认 configs/scaffold_registry.json)')
    ap.add_argument('--strategy', choices=('A', 'B', 'both'), default='both')
    ap.add_argument('--n-double', type=int, default=200, help='策略A 双点突变抽样数')
    ap.add_argument('--sa-steps', type=int, default=300, help='策略B 模拟退火步数')
    ap.add_argument('--topk', type=int, default=12)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--max-bp-dist', type=float, default=4.0)
    ap.add_argument('--spacer-unpaired-margin', type=float, default=0.10)
    ap.add_argument('--w-bp', type=float, default=0.3)
    ap.add_argument('--w-ddg', type=float, default=0.1)
    ap.add_argument('--w-contact', type=float, default=1.0)
    ap.add_argument('--w-ens', type=float, default=0.05)
    ap.add_argument('--w-hbond', type=float, default=0.5,
                    help='氢键净变化权重(有符号: 找回互补记正, 丢失记负; 8D4A 极性接触, 几何不变近似)')
    ap.add_argument('--w-fold', type=float, default=0.2,
                    help='茎完整概率变化 dp_fold 权重(正向项; bpp 独立性近似; 置 0 关闭)')
    ap.add_argument('--w-seed', type=float, default=0.0,
                    help="种子区游离度变化 d_seed 权重(正向项, spacer 3' 端)。默认 0: 与 "
                         "Δspacer_up Spearman 共线 0.835 超过 0.8 预注册合并线, 机制列保留、"
                         "权重默认关闭(湿实验标定后决定独立价值); 且 d_seed 单独加权会被"
                         "极端多点变体套利(砸烂 DR→spacer 解放), 须与硬过滤耦合使用")
    ap.add_argument('--seed-len', type=int, default=7,
                    help="种子区长度 nt(Bravo 2023: 3' 端 7 碱基有序溶剂暴露)")
    ap.add_argument('--w-stab', type=float, default=0.3,
                    help='茎稳定化奖赏权重(max(-ddG,0) 正分; 组装/均一性口径, 非活性直证; 置 0 回到纯罚分制)')
    ap.add_argument('--w-cons3', type=float, default=0.3,
                    help="3'保守窗突变惩罚权重(Dmytrenko 2023 Fig.1c 跨家族 3' 端保守)")
    ap.add_argument('--cons3-window', type=int, default=5,
                    help="DR 3' 保守窗长度(nt)")
    ap.add_argument('--proc-window', type=int, default=4,
                    help='加工位点保护窗: DR 3\' 末端 window-1 nt + 交界后首位的配对状态须与 WT 一致')
    ap.add_argument('--no-protect-processing', action='store_true',
                    help='关闭加工位点保护过滤(Dmytrenko 2023 ED Fig.3 口径)')
    ap.add_argument('--use-covariation', action='store_true',
                    help='启用策略A4: 茎区共变枚举(保持配对的双侧补偿突变)')
    ap.add_argument('--cov-double', type=int, default=150,
                    help='策略A4 双对共变组合抽样数')
    ap.add_argument('--use-struct2seq', action='store_true',
                    help='启用策略A2: Struct2SeQ 逆折叠提案(子进程调 rnet 环境)')
    ap.add_argument('--s2s-python', default='/public/home/mengxl/dzy/envs/rnet/bin/python')
    ap.add_argument('--s2s-up-bias', type=float, default=0.85, help='Struct2SeQ 向 WT 偏倚强度')
    ap.add_argument('--s2s-timeout', type=float, default=3600)
    ap.add_argument('--rnet-screen', action='store_true',
                    help='对 WT+TOP 候选做 RNet SHAPE 一致性筛选(子进程调 rnet 环境)')
    ap.add_argument('--rnet-screen-n', type=int, default=30, help='SHAPE 筛选的候选数')
    ap.add_argument('--use-grnade', action='store_true',
                    help='启用策略A3: gRNAde 三维几何逆折叠提案(仅 PDB 结构口径 DR 可用)')
    ap.add_argument('--grnade-batches', type=int, default=4, help='gRNAde 采样批数(每批 64)')
    ap.add_argument('--no-contacts', action='store_true', help='关闭蛋白接触项(不推荐)')
    ap.add_argument('--stab-dr-only', action='store_true',
                    help='w_stab 改用 DR 单独折叠 ddG_dr(茎稳定化语义更干净; 默认否, 保持 v1.4 全长口径)')
    ap.add_argument('--out-prefix', default='crrna_scaffold_run')
    args = ap.parse_args()

    entry = get_entry(args.effector, args.registry)
    dr = to_rna(entry['scaffold'])
    spacer = to_rna(args.spacer)
    if not 17 <= len(spacer) <= 25 or set(spacer) - set(BASES):
        ap.error('--spacer 必须为 17-25nt ACGT/U')
    if entry.get('scaffold_side') != '5prime':
        ap.error(f"效应子 {args.effector} 骨架方位 {entry.get('scaffold_side')} 非 5prime, 暂不支持")
    rng = np.random.default_rng(args.seed)

    wt = wt_reference(dr, spacer, args.seed_len)
    print('=== WT 参考 ===')
    print(f"效应子: {entry['display_name']}")
    print(f"DR 骨架({len(dr)}nt): {to_dna(dr)}  spacer({len(spacer)}nt, 固定): {to_dna(spacer)}")
    print(f"全长 MFE: {wt['mfe_struct']}  {wt['mfe_kcal']} kcal/mol")
    print(f"DR 单独折叠: {wt['dr_only_struct']}  {wt['dr_only_mfe_kcal']} kcal/mol")
    print(f"系综多样性: {wt['ens_diversity']}  spacer 平均未配对: {wt['spacer_mean_unpaired']}")

    contact = None if args.no_contacts else load_contacts(dr, enabled=True)
    contact_msg = '关闭' if args.no_contacts else f"{contact['source']} (DR 位点接触计数已加载)"
    if contact is not None and contact.get('transfer'):
        t = contact['transfer']
        contact_msg += (f" [对齐转移: DR 5'悬垂 {t['expected_5p_overhang']}nt/链偏移 {t['chain_offset']}, "
                        f"错配位点 {t['mismatch_positions_1based']} 的极性接触已剔除]")
    print(f"蛋白接触约束: {contact_msg}")
    stem_pos = stem_positions(wt['dr_only_struct'], len(dr))

    pool, source = {}, {}

    def _add(seq, tag):
        pool.setdefault(seq, None)
        source.setdefault(seq, set()).add(tag)

    if args.strategy in ('A', 'both'):
        for seq in strategy_a(dr, args, rng):
            _add(seq, 'A')
    if args.strategy in ('B', 'both'):
        for seq in strategy_b(dr, args, rng, wt, spacer, contact, stem_pos):
            _add(seq, 'B')
    if args.use_struct2seq:
        for seq in strategy_a2(dr, wt, args):
            _add(seq, 'S2S')
    if args.use_covariation:
        for seq in strategy_a4(dr, args, rng):
            _add(seq, 'COV')
    if args.use_grnade:
        for seq in strategy_a3(dr, args):
            _add(seq, 'G3D')

    rows = []
    for seq in pool:
        row = score_variant(dr, seq, wt, spacer, args, contact, stem_pos)
        row['source'] = '+'.join(sorted(source[seq]))
        rows.append(row)
    pareto_flag(rows)
    for row in rows:
        row['advantage'] = advantage_str(row, wt, args.seed_len)
    rows.sort(key=lambda r: (not r['passed'], -r['score']))
    for rank, row in enumerate(rows, 1):
        row['rank'] = rank

    wt_row = {'desc': 'WT', 'mut_positions': [], 'n_mut': 0, 'dr_seq': to_dna(dr),
              'construct_dna': to_dna(dr + spacer), 'mfe_struct': wt['mfe_struct'],
              'mfe_kcal': wt['mfe_kcal'], 'ddG': 0.0,
              'dr_only_mfe_kcal': wt['dr_only_mfe_kcal'], 'ddG_dr': 0.0, 'bp_dist': 0,
              'spacer_unpaired': wt['spacer_mean_unpaired'],
              'ens_diversity': wt['ens_diversity'], 'd_ens': 0.0,
              'contact_sum': 0.0, 'hbond_loss': 0.0, 'hbond_preserved': 1.0,
              'hbond_gain': 0, 'hbond_net': 0.0,
              'p_fold': round(wt['p_fold'], 5), 'dp_fold': 0.0,
              'seed_up': wt['seed_unpaired'], 'd_seed': 0.0,
              'cons3_frac': 0.0, 'proc_ok': True,
              'pareto': None, 'advantage': '参考株',
              'shape_cons': None,
              'passed': True, 'score': 0.0, 'source': 'WT', 'rank': 0}

    n_pass = sum(r['passed'] for r in rows)
    n_pareto = sum(1 for r in rows if r['pareto'])
    print(f"\n=== 候选库 ===  生成 {len(rows)} 条唯一变体, 通过硬过滤 {n_pass} 条, "
          f"Pareto 最优 {n_pareto} 条")

    # 正向项区分力闸门: 非零比例 <5% 即声明该指标在本库无区分力(防虚假精度)
    for name, key, thresh in (('dp_fold', 'dp_fold', 0.01), ('hbond_gain', 'hbond_gain', 0.5),
                              ('d_seed', 'd_seed', 0.01)):
        nz = sum(1 for r in rows if abs(r.get(key) or 0) > thresh)
        frac = nz / len(rows) if rows else 0
        flag = ' [警告: 无区分力, 建议置零权重]' if frac < 0.05 else ''
        print(f"[闸门] {name} 非零(|>{thresh}|)比例 {frac:.1%}{flag}")

    if args.rnet_screen:
        screen_rows = [wt_row] + [r for r in rows if r['passed']][:args.rnet_screen_n]
        entries = [(('WT' if r is wt_row else f"r{r['rank']}"),
                    to_rna(r['construct_dna'])) for r in screen_rows]
        cons = rnet_shape_screen(entries, wt['mfe_struct'], args)
        if cons is not None:
            for r in screen_rows:
                r['shape_cons'] = cons.get('WT' if r is wt_row else f"r{r['rank']}")
            print(f"RNet SHAPE 一致性(对 WT 结构): WT = {cons.get('WT')}"
                  f"（变体列 shape_cons; < WT-10 者结构保持存疑）")
    print(f"{'rank':<5}{'desc':<14}{'ddG':>7}{'bp_d':>5}{'sp_up':>7}{'ens':>6}"
          f"{'cont':>6}{'hb':>6}{'sha':>7}{'score':>8}  pass/src")
    for row in rows[:args.topk]:
        sha = f"{row['shape_cons']:>7.1f}" if row['shape_cons'] is not None else f"{'-':>7}"
        print(f"{row['rank']:<5}{row['desc']:<14}{row['ddG']:>7.2f}{row['bp_dist']:>5}"
              f"{row['spacer_unpaired']:>7.3f}{row['ens_diversity']:>6.2f}"
              f"{row['contact_sum']:>6.2f}{row['hbond_loss']:>6.2f}"
              f"{sha}{row['score']:>8.3f}  {'Y' if row['passed'] else 'N'}/{row['source']}")

    out_csv = args.out_prefix + '.variants.csv'
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(wt_row.keys()))
        writer.writeheader()
        writer.writerow(wt_row)
        writer.writerows(rows)
    top = [wt_row] + [r for r in rows if r['passed']][:args.topk]
    out_json = args.out_prefix + '.top.json'
    payload = {
        'disclaimer': '候选库压缩结果, 非活性预测; 活性/特异性以体外生化与细胞实验为准; '
                      'dp_fold 与 hbond_net 为正向代理项(bpp 独立性近似/几何不变近似), '
                      '未经 Cas12a2 实验标定, 排序含义为先验更优而非已证提活',
        'effector': args.effector, 'registry_entry_version': entry.get('version'),
        'spacer_fixed_dna': to_dna(spacer), 'wt': wt,
        'params': {k: getattr(args, k) for k in
                   ('strategy', 'n_double', 'sa_steps', 'seed', 'max_bp_dist',
                    'spacer_unpaired_margin', 'w_bp', 'w_ddg', 'w_contact', 'w_ens',
                    'w_hbond', 'w_fold', 'w_seed', 'seed_len', 'w_stab', 'w_cons3', 'cons3_window', 'proc_window',
                    'no_protect_processing', 'use_covariation', 'cov_double',
                    'use_struct2seq', 's2s_up_bias', 'rnet_screen',
                    'rnet_screen_n', 'no_contacts', 'stab_dr_only', 'topk')},
        'n_variants': len(rows), 'n_passed': n_pass, 'n_pareto': n_pareto,
        'library_sha256': hashlib.sha256(
            json.dumps([r['dr_seq'] for r in rows], sort_keys=True).encode()).hexdigest(),
        'top': top,
    }
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    out_fasta = args.out_prefix + '.top.fasta'
    with open(out_fasta, 'w', encoding='utf-8') as fh:
        for row in top:
            fh.write(f">crRNA_{row['rank']}_{row['desc']}|{args.effector}|spacer_fixed\n"
                     f"{row['construct_dna']}\n")
    print(f"\n输出: {out_csv} / {out_json} / {out_fasta}")
    print("说明: top.json 含 WT(第 0 名) + TOP-K 通过硬过滤变体; spacer 全程固定, 变量只有 DR 骨架。")


if __name__ == '__main__':
    main()
