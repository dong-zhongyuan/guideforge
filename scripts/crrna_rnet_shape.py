# -*- coding: utf-8 -*-
"""RNet SHAPE 化学图谱预测探针（crRNA 链；须在 rnet 环境运行）

由 crrna_scaffold_design.py --rnet-screen 经子进程调用。加载 toolbox 内
RibonanzaNet 主模型（经 Struct2SeQ Env 包装），对输入序列预测逐位 SHAPE
反应性（Science 2026 aeg6829 的"先模拟实验"环节所用口径）。

用法:
  python crrna_rnet_shape.py --csv in.csv --out out.csv
输入 in.csv 列: id,sequence（RNA, ACGU）; 输出 out.csv 列: id,shape（分号分隔浮点）
"""
import argparse
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S2S_DIR = os.path.join(ROOT, 'toolbox', 'Struct2SeQ')
RNET_W = os.path.join(ROOT, 'toolbox', 'rnet-inference', 'RibonanzaNet-Weights')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--batch', type=int, default=32)
    args = ap.parse_args()

    os.chdir(S2S_DIR)  # DQN_env 用相对路径读 test10_configs/pairwise.yaml
    sys.path.insert(0, S2S_DIR)

    # 同 crrna_struct2seq_gen.py: arnie 配置必须在任何 arnie import 前备好(绝对路径);
    # TMP 放包内 _arnie_tmp(服务器 /tmp 不一定可写, 且避免移交后路径失效)
    arnie_tmp = os.path.join(S2S_DIR, '_arnie_tmp')
    os.makedirs(arnie_tmp, exist_ok=True)
    arnie_cfg = os.path.join(S2S_DIR, 'arnie_file.txt')
    with open(arnie_cfg, 'w', encoding='utf-8') as fh:
        fh.write('linearpartition: . \nTMP: ' + arnie_tmp + '\n')
    os.environ['ARNIEFILE'] = arnie_cfg

    import numpy as np
    import torch
    from Env import DQN_env
    from Dataset import tokenize_sequence

    if not torch.cuda.is_available():
        torch.Tensor.cuda = lambda self, *a, **k: self  # type: ignore
        torch.nn.Module.cuda = lambda self, *a, **k: self  # type: ignore

    env = DQN_env(os.path.join(RNET_W, 'RibonanzaNet.pt'),
                  os.path.join(RNET_W, 'RibonanzaNet-SS.pt'),
                  use_gpu=torch.cuda.is_available())

    rows = list(csv.DictReader(open(args.csv, encoding='utf-8')))
    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['id', 'shape'])
        for start in range(0, len(rows), args.batch):
            chunk = rows[start:start + args.batch]
            src = torch.tensor([tokenize_sequence(r['sequence']) for r in chunk]).long()
            profiles = env.get_SHAPE(src)  # [B, L]
            for row, prof in zip(chunk, profiles):
                writer.writerow([row['id'],
                                 ';'.join(f'{float(x):.4f}' for x in np.ravel(prof))])


if __name__ == '__main__':
    main()
