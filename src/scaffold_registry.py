"""骨架注册表读取接口（2026-08-31 体系切换后引入）。

背景：平台不再绑定单一固定骨架。每个效应子一条注册记录（序列/长度/PAM/参考结构/
置信基线/版本），见 docs/结构维度_平行设计框架_20260831.md §1.5。
注册记录由 scripts/register_scaffold.py 生成（RNet-SS × ViennaRNA 配分函数），
本模块只负责读取，不做任何计算，可在生产链安全 import（无重依赖）。

与 configs/model_contract.json 的关系：model_contract.json 的 scaffold 字段是
SpCas9 锁定生产链的训练契约（83nt 校验见 contract.py），不动；注册表是
多效应子时代的骨架来源，crRNA 新模块（crrna_*）一律从注册表取骨架。
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


if __name__ == '__main__':
    for name in list_effectors():
        e = get_entry(name)
        print(f"{name}: scaffold {e['scaffold_len']}nt ({e['scaffold_side']}), "
              f"PAM {e['pam']['rule']} ({e['pam']['side']}), "
              f"gate={e['reference']['gate']}, version={e['version']}")
