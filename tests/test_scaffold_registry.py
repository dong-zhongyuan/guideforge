"""骨架注册表契约测试（2026-08-31）。

守护对象：configs/scaffold_registry.json + src/scaffold_registry.py 读取接口。
注册条目由 scripts/register_scaffold.py 生成；本测试不重新注册，只校验
现有条目的完整性与接口行为，防止注册表被手改坏。
"""
import os
import unittest

from scaffold_registry import (
    DEFAULT_REGISTRY_PATH, get_entry, get_gate, get_pam,
    get_reference_dotbracket, get_scaffold, list_effectors,
)

GATES = {'confident', 'low_conf', 'unfolded'}


class TestScaffoldRegistry(unittest.TestCase):
    def test_registry_file_exists(self):
        self.assertTrue(os.path.exists(DEFAULT_REGISTRY_PATH),
                        '注册表缺失，先跑 scripts/register_scaffold.py')

    def test_default_entries_present(self):
        effectors = list_effectors()
        for name in ('cas12a', 'spcas9'):
            self.assertIn(name, effectors, f'默认条目 {name} 未注册')

    def test_entry_fields_and_sequence(self):
        for name in list_effectors():
            e = get_entry(name)
            self.assertLessEqual(set(e['scaffold']), set('ACGT'),
                                 f'{name} 骨架含非 ACGT 字符')
            self.assertEqual(len(e['scaffold']), e['scaffold_len'],
                             f'{name} scaffold_len 与实际长度不符')
            self.assertIn(e['scaffold_side'], ('3prime', '5prime'))
            self.assertIn(e['pam']['side'], ('3prime', '5prime'))
            self.assertIn(e['reference']['gate'], GATES)
            self.assertTrue(e['inputs_sha256'])
            # 参考结构长度 = 占位 spacer + 骨架
            expect_len = len(e['placeholder_spacer']) + e['scaffold_len']
            self.assertEqual(len(get_reference_dotbracket(name)), expect_len,
                             f'{name} 参考结构长度与构建序列不符')

    def test_accessor_types(self):
        for name in list_effectors():
            self.assertIsInstance(get_scaffold(name), str)
            self.assertIsInstance(get_reference_dotbracket(name), str)
            self.assertIn(get_gate(name), GATES)
            self.assertIsInstance(get_pam(name), dict)

    def test_unregistered_effector_raises(self):
        with self.assertRaises(KeyError):
            get_entry('__no_such_effector__')

    def test_cas12a_layout(self):
        """cas12a 条目：骨架在 5' 端（handle-spacer），PAM 5' 侧。"""
        e = get_entry('cas12a')
        self.assertEqual(e['scaffold_side'], '5prime')
        self.assertEqual(e['pam']['side'], '5prime')
        self.assertTrue(e['full_construct_dna'].startswith(e['scaffold']))


if __name__ == '__main__':
    unittest.main()
