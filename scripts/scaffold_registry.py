# -*- coding: utf-8 -*-
"""src/scaffold_registry.py 的 re-export shim(2026-09-02, 审计 P2 去重)。

真实实现只在 src/scaffold_registry.py(单一事实源); 本文件保留是因为
scripts/ 下多个脚本与 tests 以 'import scaffold_registry'(同目录)方式引用,
且独立运行(python scripts/xxx.py)时不设 PYTHONPATH=src。任何修改请改
src/scaffold_registry.py, 不要改本文件。
"""
import importlib.util
import os

_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'src', 'scaffold_registry.py')
_spec = importlib.util.spec_from_file_location('_src_scaffold_registry', _SRC)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
for _k in dir(_mod):
    if not _k.startswith('_'):
        globals()[_k] = getattr(_mod, _k)
del _spec, _mod, _k
