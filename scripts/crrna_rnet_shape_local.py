# -*- coding: utf-8 -*-
"""RNet SHAPE 化学图谱预测探针（本地/Windows 版；纯 torch，不经 arnie）

与 crrna_rnet_shape.py 的差异: 不 import Struct2SeQ 的 Env/Dataset（其 arnie
配置文件解析在 Windows 绝对路径盘符冒号下报错），直接加载 RibonanzaNet 反应性
模型做前向推理。模型、权重、配置与原版完全一致（RibonanzaNet.pt +
test10_configs/pairwise.yaml），只需 RibonanzaNet.pt 一个权重（-SS.pt 仅供
结构预测，SHAPE 筛选用不到）。

用法:
  python crrna_rnet_shape_local.py --csv in.csv --out out.csv
输入 in.csv 列: id,sequence（RNA, ACGU）; 输出 out.csv 列: id,shape（分号分隔浮点）
"""
import argparse
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S2S_DIR = os.path.join(ROOT, 'toolbox', 'Struct2SeQ')
RNET_W = os.path.join(ROOT, 'toolbox', 'rnet-inference', 'RibonanzaNet-Weights')

NT2IDX = {nt: i for i, nt in enumerate('ACGU')}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--batch', type=int, default=32)
    args = ap.parse_args()

    sys.path.insert(0, S2S_DIR)
    import numpy as np
    import torch
    import yaml
    from Network_test10 import RibonanzaNet

    class Config:
        def __init__(self, **entries):
            self.__dict__.update(entries)

    with open(os.path.join(S2S_DIR, 'test10_configs', 'pairwise.yaml'),
              encoding='utf-8') as fh:
        config = Config(**yaml.safe_load(fh))
    model = RibonanzaNet(config)
    model.load_state_dict(torch.load(os.path.join(RNET_W, 'RibonanzaNet.pt'),
                                     map_location='cpu'))
    model.eval()

    rows = list(csv.DictReader(open(args.csv, encoding='utf-8')))
    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['id', 'shape'])
        for start in range(0, len(rows), args.batch):
            chunk = rows[start:start + args.batch]
            src = torch.tensor([[NT2IDX[c] for c in r['sequence'].upper()]
                                for r in chunk]).long()
            with torch.no_grad():
                profiles = model(src)[:, :, 0].cpu().numpy()
            for row, prof in zip(chunk, profiles):
                writer.writerow([row['id'],
                                 ';'.join(f'{float(x):.4f}' for x in np.ravel(prof))])


if __name__ == '__main__':
    main()
