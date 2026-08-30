"""骨架注册脚本：新增/更新一个效应子的注册条目（2026-08-31）。

注册流程（docs/结构维度_平行设计框架_20260831.md §1.5）：
  1) 占位 spacer 拼接骨架（spcas9: spacer+scaffold；cas12a: handle+spacer）
  2) ViennaRNA 配分函数：MFE/质心/系综距离/逐位置 unpaired 概率（本 venv 内计算）
  3) RNet-SS 预测（经子进程调用独立 rnet 环境，本 venv 不引入 torch）
  4) 两路结构 F1 → 置信门控基线（confident / low_conf / unfolded）
  5) 条目（含参考折叠、基线指标、版本、输入 sha256）写入 configs/scaffold_registry.json

RNet 子进程路径：默认取本机部署（C:\\Users\\zhangdanjing\\rna-pseudoknot-ai），
可用环境变量 RNET_PYTHON / RNET_2D_SCRIPT 覆盖（开源/换机时必须设置）。
--skip-rnet 时仅 ViennaRNA 注册，门控强制 low_conf 并留注。

示例（cas12a 占位条目）:
  python scripts/register_scaffold.py --effector cas12a \
      --scaffold TAATTTCTACTAAGTGTAGAT --scaffold-side 5prime \
      --pam-rule TTTV --pam-side 5prime \
      --display-name "Cas12a (As 通用 handle, 占位)" \
      --bound-state-note "handle 假结为 Cas12a 蛋白诱导, 游离态不出现"
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

import RNA  # ViennaRNA 2.x（项目 venv 既有依赖）

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DEFAULT_REGISTRY = os.path.join(ROOT, 'configs', 'scaffold_registry.json')
DEFAULT_RNET_PYTHON = os.environ.get(
    'RNET_PYTHON', r'C:\Users\zhangdanjing\anaconda3\envs\rnet\python.exe')
DEFAULT_RNET_2D = os.environ.get(
    'RNET_2D_SCRIPT', os.path.join(ROOT, 'toolbox', 'rnet-inference', 'src', 'rnet_2d.py'))

OPEN2CLOSE = {'(': ')', '[': ']', '{': '}', '<': '>'}
CLOSE_SET = set(OPEN2CLOSE.values())


def dotbracket_to_bps(db):
    """点括号 → 碱基对集合 {(i,j)}（0-based，支持 () [] {} <> 四类，即可识别假结）。"""
    stacks, bps = {}, set()
    for i, ch in enumerate(db):
        if ch in OPEN2CLOSE:
            stacks.setdefault(ch, []).append(i)
        elif ch in CLOSE_SET:
            for op, st in stacks.items():
                if OPEN2CLOSE[op] == ch:
                    if not st:
                        raise ValueError(f'点括号不平衡（位置 {i} 多余 {ch}）')
                    bps.add((st.pop(), i))
                    break
    if any(stacks.values()):
        raise ValueError('点括号不平衡（存在未闭合括号）')
    return bps


def bp_f1(db1, db2):
    """两个点括号结构的碱基对集合 F1。"""
    a, b = dotbracket_to_bps(db1), dotbracket_to_bps(db2)
    if not a or not b:
        return 0.0
    tp = len(a & b)
    return 2.0 * tp / (len(a) + len(b))


def vienna_baseline(full_seq_rna, spacer_len, scaffold_side):
    """ViennaRNA 配分函数基线：MFE/质心/系综距离/spacer区平均unpaired概率。"""
    fc = RNA.fold_compound(full_seq_rna)
    mfe_struct, mfe = fc.mfe()
    fc.pf()
    centroid, ens_dist = fc.centroid()
    n = len(full_seq_rna)
    bpp = fc.bpp()
    unpaired = [1.0 - sum(bpp[i][j] for j in range(1, n + 1) if j != i)
                for i in range(1, n + 1)]
    if scaffold_side == '3prime':
        spacer_idx = range(0, spacer_len)
    else:
        spacer_idx = range(n - spacer_len, n)
    spacer_unp = sum(unpaired[i] for i in spacer_idx) / spacer_len
    return {
        'vienna_mfe_struct': mfe_struct,
        'vienna_mfe_kcal': round(float(mfe), 2),
        'vienna_centroid': centroid,
        'vienna_ens_diversity': round(float(ens_dist), 2),
        'spacer_mean_unpaired': round(spacer_unp, 3),
    }


def rnet_predict(full_seq_rna, rnet_python, rnet_2d, timeout=600):
    """经子进程调用 rnet 环境的 rnet_2d.py，解析 structure 与 confidence_f1。"""
    proc = subprocess.run(
        [rnet_python, rnet_2d, full_seq_rna],
        capture_output=True, text=True, encoding='utf-8', errors='replace',
        timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            f'RNet 子进程失败(rc={proc.returncode}): '
            f'stderr={ (proc.stderr or "")[-500:] } stdout={ (proc.stdout or "")[-300:] }')
    structure = conf = None
    for line in proc.stdout.splitlines():
        if line.startswith('structure:'):
            structure = line.split(':', 1)[1].strip()
        elif line.startswith('confidence_f1:'):
            conf = float(line.split(':', 1)[1].strip())
    if structure is None or conf is None:
        raise RuntimeError(f'RNet 输出解析失败: {proc.stdout[-500:]}')
    return structure, conf


def build_entry(args, baseline, rnet_result, gate, gate_f1):
    scaffold = args.scaffold.upper()
    spacer = args.placeholder_spacer.upper()
    full_dna = spacer + scaffold if args.scaffold_side == '3prime' else scaffold + spacer
    ref_db = baseline['vienna_centroid'] if rnet_result is None else rnet_result[0]
    inputs = '|'.join([args.effector, scaffold, spacer, args.scaffold_side,
                       args.pam_rule, args.pam_side, args.version])
    entry = {
        'effector': args.effector,
        'display_name': args.display_name,
        'scaffold': scaffold,
        'scaffold_len': len(scaffold),
        'scaffold_side': args.scaffold_side,
        'placeholder_spacer': spacer,
        'full_construct_dna': full_dna,
        'pam': {'rule': args.pam_rule.upper(), 'side': args.pam_side},
        'reference': {
            'full_len': len(full_dna),
            'dotbracket': ref_db,
            'source': ('RNet-SS × ViennaRNA pf 注册运行' if rnet_result
                       else '仅 ViennaRNA（--skip-rnet）'),
            'gate': gate,
            'baseline': {
                **{k: v for k, v in baseline.items()
                   if k not in ('vienna_mfe_struct', 'vienna_centroid')},
                'vienna_rnet_structure_f1': (round(gate_f1, 3)
                                             if gate_f1 is not None else None),
                'rnet_confidence_f1': (round(rnet_result[1], 3)
                                       if rnet_result else None),
            },
            'bound_state_note': args.bound_state_note,
        },
        'version': args.version,
        'registered_at': datetime.date.today().isoformat(),
        'inputs_sha256': hashlib.sha256(inputs.encode('utf-8')).hexdigest(),
    }
    return entry


def write_entry(registry_path, entry):
    if os.path.exists(registry_path):
        with open(registry_path, encoding='utf-8') as f:
            reg = json.load(f)
    else:
        reg = {'registry_version': 1, 'entries': {}}
    reg['entries'][entry['effector']] = entry
    os.makedirs(os.path.dirname(registry_path), exist_ok=True)
    with open(registry_path, 'w', encoding='utf-8') as f:
        json.dump(reg, f, ensure_ascii=False, indent=2)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--effector', required=True, help='效应子键名，如 cas12a / spcas9')
    p.add_argument('--scaffold', required=True, help='骨架序列（DNA 字母表，ACGT）')
    p.add_argument('--scaffold-side', choices=['3prime', '5prime'], default='3prime',
                   help="骨架相对 spacer 的位置：spcas9 为 3prime，cas12a 为 5prime")
    p.add_argument('--pam-rule', required=True, help='PAM 规则，如 NGG / TTTV')
    p.add_argument('--pam-side', choices=['3prime', '5prime'], default='3prime')
    p.add_argument('--display-name', default='', help='展示名（株系/来源说明）')
    p.add_argument('--placeholder-spacer', default='GACCCCCTCCACCCCGCCTC',
                   help='注册用占位 spacer（默认 20nt，两个先导实验验证游离）')
    p.add_argument('--version', default='v1')
    p.add_argument('--bound-state-note', default='',
                   help='结合态注释（如蛋白诱导假结；仅注释不作评分基准）')
    p.add_argument('--registry', default=DEFAULT_REGISTRY)
    p.add_argument('--rnet-python', default=DEFAULT_RNET_PYTHON)
    p.add_argument('--rnet-2d', default=DEFAULT_RNET_2D)
    p.add_argument('--skip-rnet', action='store_true',
                   help='跳过 RNet（仅 ViennaRNA 注册，门控强制 low_conf）')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    scaffold = args.scaffold.upper()
    if set(scaffold) - set('ACGT'):
        raise ValueError('--scaffold 含非 ACGT 字符')
    spacer = args.placeholder_spacer.upper()
    full_dna = spacer + scaffold if args.scaffold_side == '3prime' else scaffold + spacer
    full_rna = full_dna.replace('T', 'U')

    print(f'[1/4] 构建序列 {len(full_dna)}nt = spacer {len(spacer)} + scaffold {len(scaffold)} '
          f'({args.scaffold_side})')
    baseline = vienna_baseline(full_rna, len(spacer), args.scaffold_side)
    print(f"[2/4] ViennaRNA: MFE {baseline['vienna_mfe_kcal']} kcal/mol, "
          f"系综距离 {baseline['vienna_ens_diversity']}, "
          f"spacer unpaired {baseline['spacer_mean_unpaired']}")

    rnet_result, gate_f1 = None, None
    if args.skip_rnet:
        gate = 'low_conf'
        print('[3/4] --skip-rnet：仅 ViennaRNA 注册')
    else:
        print('[3/4] 调用 RNet-SS（子进程，可能需 1-2 分钟加载模型）...')
        rnet_result = rnet_predict(full_rna, args.rnet_python, args.rnet_2d)
        gate_f1 = bp_f1(baseline['vienna_centroid'], rnet_result[0])
        print(f'      RNet: {rnet_result[0]}  (自评 F1 {rnet_result[1]:.3f})')
        print(f'      两路结构 F1 = {gate_f1:.3f}')
        unfolded = baseline['spacer_mean_unpaired'] > 0.8 and \
            baseline['vienna_ens_diversity'] > 5.0
        if unfolded:
            gate = 'unfolded'
        else:
            gate = 'confident' if gate_f1 >= 0.8 else 'low_conf'

    entry = build_entry(args, baseline, rnet_result, gate, gate_f1)
    write_entry(args.registry, entry)
    print(f"[4/4] 已写入 {args.registry}: {args.effector} "
          f"(gate={gate}, version={args.version})")


if __name__ == '__main__':
    main()
