# -*- coding: utf-8 -*-
"""Chai 增量合并器契约测试(2026-09-06)。

锁定 scripts/crrna_chai_increment_merge.py 的合并口径:
  --selftest 必须全过(合成 11 折 npz -> 32 行 4靶x8骨架, G12C 列丢弃,
  退役 B 臂行被新臂替换, 复用 21 行与旧矩阵逐值一致);
  缺 --runs-dir 时必须以用法错误退出(防止空合并产出假矩阵)。
"""
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SCRIPT = os.path.join(ROOT, "scripts", "crrna_chai_increment_merge.py")


def run(args):
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8",
               PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, SCRIPT] + args, cwd=ROOT, env=env,
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace")


class TestChaiIncrementMerge(unittest.TestCase):

    def test_selftest_passes(self):
        r = run(["--selftest"])
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        self.assertIn("selftest", r.stdout)

    def test_missing_runs_dir_is_usage_error(self):
        r = run([])
        self.assertNotEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()
