"""骨架注册表读取接口。

背景：平台不绑定单一固定骨架。每个效应子一条注册记录（序列/长度/PAM/参考结构/
置信基线/版本）。
注册记录由 scripts/register_scaffold.py 生成（RNet-SS × ViennaRNA 配分函数），
本模块只负责读取，不做任何计算，可在生产链安全 import（无重依赖）。
每个效应子一条注册记录，crRNA 模块（crrna_*）一律从注册表取骨架。
"""
from __future__ import annotations

import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DEFAULT_REGISTRY_PATH = os.path.join(ROOT, 'configs', 'scaffold_registry.json')

_REQUIRED_FIELDS = (
    'effector', 'scaffold', 'scaffold_len', 'scaffold_side',
    'pam', 'reference', 'version', 'registered_at',
)


def load_registry(path=None):
    """读取注册表 JSON，返回 dict。文件不存在时报错并提示先跑注册脚本。"""
    path = os.path.abspath(path or DEFAULT_REGISTRY_PATH)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f'骨架注册表不存在: {path}；请先运行 scripts/register_scaffold.py 注册条目')
    with open(path, encoding='utf-8') as f:
        reg = json.load(f)
    if 'entries' not in reg or not isinstance(reg['entries'], dict):
        raise ValueError(f'注册表格式错误（缺 entries）: {path}')
    return reg


def list_effectors(path=None):
    """返回已注册效应子名列表。"""
    return sorted(load_registry(path)['entries'].keys())


def get_entry(effector, path=None):
    """取单个效应子的完整注册记录，并做最小字段校验。"""
    reg = load_registry(path)
    try:
        entry = reg['entries'][effector]
    except KeyError:
        raise KeyError(
            f"效应子 '{effector}' 未注册；已注册: {sorted(reg['entries'])}") from None
    missing = [k for k in _REQUIRED_FIELDS if k not in entry]
    if missing:
        raise ValueError(f"注册条目 '{effector}' 缺字段: {missing}")
    return entry


def get_scaffold(effector, path=None):
    """返回骨架序列（DNA 字母表 T，项目约定；需 RNA 字母时调用方自行 T→U）。"""
    return get_entry(effector, path)['scaffold']


def get_reference_dotbracket(effector, path=None):
    """返回该效应子的游离态参考折叠（点括号，长度=占位spacer+骨架）。"""
    return get_entry(effector, path)['reference']['dotbracket']


def get_gate(effector, path=None):
    """返回注册时的置信门控基线（confident / low_conf / unfolded）。"""
    return get_entry(effector, path)['reference']['gate']


def get_pam(effector, path=None):
    """返回 PAM 规则 dict，含 rule 与 side 两个键。"""
    return get_entry(effector, path)['pam']


IUPAC_CODES = frozenset('ACGTURYSWKMBDHVN')
_PFS_REQUIRED = ('consensus', 'side', 'tolerant_mismatches')


def _validate_pfs(effector, pfs):
    """校验 pfs 字段可被扫描器直接消费（具体 ACGTU 共识 + 整数容忍度）。"""
    missing = [k for k in _PFS_REQUIRED if k not in pfs]
    if missing:
        raise ValueError(f"注册条目 '{effector}' 的 pfs 缺字段: {missing}")
    consensus = str(pfs['consensus']).upper().replace('U', 'T')
    if not consensus or set(consensus) - set('ACGT'):
        raise ValueError(
            f"注册条目 '{effector}' 的 pfs.consensus 必须为具体 ACGT/U 序列 "
            f"(扫描器逐字符比对, 不支持简并码/符号标签), 实为 {pfs['consensus']!r}")
    tol = pfs['tolerant_mismatches']
    if not isinstance(tol, int) or isinstance(tol, bool) or not 0 <= tol < len(consensus):
        raise ValueError(
            f"注册条目 '{effector}' 的 pfs.tolerant_mismatches 须为 [0, len(consensus)) 整数, "
            f"实为 {tol!r}")
    if pfs['side'] not in ('3prime', '5prime'):
        raise ValueError(f"注册条目 '{effector}' 的 pfs.side 非法: {pfs['side']!r}")


def get_pfs(effector, path=None):
    """返回 PFS（protospacer-flanking sequence）规则 dict 的副本。

    含 consensus（具体 ACGT 序列）、tolerant_mismatches（容忍错配数，int）、
    side，可选 note。PFS 是 RNA 靶向效应子（Cas12a2）靶 RNA 上的侧翼识别基序；
    DNA 靶向条目（cas12a/spcas9）无 pfs 字段，调用抛 KeyError。
    """
    entry = get_entry(effector, path)
    if 'pfs' not in entry:
        raise KeyError(
            f"效应子 '{effector}' 未注册 pfs 字段（仅 RNA 靶向条目有 PFS 概念）")
    pfs = dict(entry['pfs'])
    _validate_pfs(effector, pfs)
    return pfs


def pfs_mismatches(seq, consensus):
    """观察 PFS 与共识序列的汉明距离（U→T 归一）；长度不等返回 None（不可判）。"""
    seq = str(seq).upper().replace('U', 'T')
    ref = str(consensus).upper().replace('U', 'T')
    if len(seq) != len(ref):
        return None
    return sum(a != b for a, b in zip(seq, ref))


def resolve_pfs(effector, consensus=None, tolerant_mismatches=None, path=None):
    """解析生效 PFS 规则：默认注册表 pfs 字段，允许调用方（CLI）覆盖。

    返回 dict: consensus / tolerant_mismatches / side / note / source
    （source = 'registry:<effector>' 或 'cli-override' 组合）。
    """
    spec = get_pfs(effector, path)
    source = f'registry:{effector}'
    if consensus is not None:
        spec['consensus'] = str(consensus).upper().replace('U', 'T')
        source = 'cli-override'
    if tolerant_mismatches is not None:
        spec['tolerant_mismatches'] = int(tolerant_mismatches)
        if consensus is None:
            source = f'registry:{effector}+cli-tol'
        else:
            source = 'cli-override'
    _validate_pfs(effector, spec)
    spec['source'] = source
    return spec


if __name__ == '__main__':
    for name in list_effectors():
        e = get_entry(name)
        pfs = f", PFS {e['pfs']['consensus']}±{e['pfs']['tolerant_mismatches']}" \
            if 'pfs' in e else ''
        print(f"{name}: scaffold {e['scaffold_len']}nt ({e['scaffold_side']}), "
              f"PAM {e['pam']['rule']} ({e['pam']['side']}){pfs}, "
              f"gate={e['reference']['gate']}, version={e['version']}")
