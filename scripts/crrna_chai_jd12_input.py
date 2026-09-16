# -*- coding: utf-8 -*-
"""JD12 定稿 6 构建的 Chai 共折叠输入准备(2026-09-16, 界面验证缺口补齐)。

缺口: U15G / A8C+U15G(DR 语境)在 JD12 spacer 上下文的蛋白-crRNA 界面
从未做共折叠验证(A1C 有 8D4A 零接触+双榜间接证据)。本脚本产出 6 条
Chai 输入 FASTA(蛋白 + crRNA + 靶 RNA 三链, 与 chai_matrix 同口径),
模板命中复用 dawn/chai_template_hits.m8(蛋白相同)。

输入源(单一事实源): data/wetlab_jd12_order.csv(下单表)。
输出: /public/home/mengxl/openclaw/dawn/chai_jd12/<name>.fasta
运行: PYTHONUTF8=1 python scripts/crrna_chai_jd12_input.py
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data')
OUTD = '/public/home/mengxl/openclaw/dawn/chai_jd12'

# 蛋白序列: 取自既有 chai_matrix 输入(单一事实源, 不重抄)
PROT_SRC = '/public/home/mengxl/openclaw/dawn/chai_input_WT.fasta'


def main():
    prot = None
    name = None
    for line in open(PROT_SRC, encoding='utf-8'):
        if line.startswith('>'):
            if prot:
                break
            name = line[1:].strip()
            prot = []
        elif prot is not None:
            prot.append(line.strip())
    prot = ''.join(prot)
    assert prot.startswith('MLHAFTNQ') and len(prot) > 1000, '蛋白序列异常'

    os.makedirs(OUTD, exist_ok=True)
    rows = list(csv.DictReader(open(
        os.path.join(DATA, 'wetlab_jd12_order.csv'), encoding='utf-8')))
    manifest = []
    for r in rows:
        crna = r['rna_crRNA_42nt']
        spacer = r['spacer_dna_23nt'].replace('T', 'U')
        # 靶 RNA = spacer 互补(与 alphafoldserver/af3 输入同口径: 反向互补)
        target = spacer.translate(str.maketrans('AUCG', 'UAGC'))[::-1]
        fname = r['oligo_name'].replace('JD12_', 'JD_')
        fa = os.path.join(OUTD, fname + '.fasta')
        with open(fa, 'w', encoding='utf-8') as f:
            f.write('>protein|name=SuCas12a2\n%s\n' % prot)
            f.write('>rna|name=crRNA\n%s\n' % crna)
            f.write('>rna|name=target\n%s\n' % target)
        manifest.append({'name': fname, 'scaffold': r['scaffold'],
                         'spacer_group': 'crRNA1' if 'crRNA1' in r['oligo_name'] else 'crRNA2',
                         'fasta': fa})
        print('%-24s %s %s' % (fname, r['scaffold'], crna[:19]))

    json.dump(manifest, open(os.path.join(DATA, 'chai_jd12_manifest.json'), 'w',
                             encoding='utf-8'), ensure_ascii=False, indent=1)
    print('-> %d 条输入 + manifest | 模板命中复用 dawn/chai_template_hits.m8'
          % len(manifest))


if __name__ == '__main__':
    main()
