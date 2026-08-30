# -*- coding: utf-8 -*-
"""gRNAde 三维几何逆折叠生成包装器（crRNA 链 · 策略A3 引擎；须在 rnet 环境运行）

由 crrna_scaffold_design.py --use-grnade 经子进程调用。与官方 design.py 的差异：
  - 官方 design.py 自带 Ribonanzanet 自洽评分（依赖未入包的 tools/ 目录），本包装器
    只做"几何约束生成 + perplexity 输出"，结构筛选统一由主管线（ViennaRNA/RNet）完成，
    口径单一、不重复造轮子；
  - 部分序列约束：spacer 位置按 8D4A 结构原生序列锁定（'DR 位点=_ 可设计'），
    仅 DR 区（5' 端 --dr-len 个位置）被设计。

输入几何：data/8D4A_crna_chainB.pdb（Cas12a2 结合态 crRNA 链，8D4A B 链）。
注意：模型按 backbone 几何生成序列；spacer 骨干来自 8D4A 结构，最终 construct 的
spacer 序列由主管线替换为项目固定 spacer（几何上下文不受影响）。

用法:
  python crrna_grnade_gen.py --pdb ../data/8D4A_crna_chainB.pdb --dr-len 18 \
      --n-batches 4 --batch-size 64 --out out.csv
输出 out.csv 列: sequence（全长）, perplexity, temperature
"""
import argparse
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRNADE_DIR = os.path.join(ROOT, 'toolbox', 'geometric-rna-design')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--dr-len', type=int, required=True)
    ap.add_argument('--n-batches', type=int, default=4)
    ap.add_argument('--batch-size', type=int, default=64)
    ap.add_argument('--temp-min', type=float, default=0.1)
    ap.add_argument('--temp-max', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--target-sec-struct', default='',
                    help='结构链的二级结构(点括号, 长度=链长); 传入后跳过从结构推算')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    sys.path.insert(0, GRNADE_DIR)
    os.chdir(GRNADE_DIR)  # src.constants.PROJECT_PATH 等相对路径解析
    # src.constants 从环境变量读路径(import 时生效); 本流程只需要 PROJECT_PATH
    # 指向 gRNAde 根目录, X3DNA/ETERNAFOLD 用不到(二级结构由 --target-sec-struct 传入)
    os.environ.setdefault('PROJECT_PATH', GRNADE_DIR)
    os.environ.setdefault('DATA_PATH', os.path.join(GRNADE_DIR, 'data'))
    os.environ.setdefault('X3DNA', '/nonexistent')
    os.environ.setdefault('ETERNAFOLD', '/nonexistent')

    import random
    import numpy as np
    import torch
    import torch.nn.functional as F
    import yaml

    from src.models import gRNAde
    from src.data.featurizer import RNAGraphFeaturizer
    from src.constants import NUM_TO_LETTER
    if args.target_sec_struct:
        # 二级结构由主管线传入(ViennaRNA 折叠, 与 RNet-SS 复核一致), 跳过 pdb 推算
        import src.data.data_utils as _du
        _du.pdb_to_sec_struct = lambda *a, **k: args.target_sec_struct

    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)

    with open(os.path.join(GRNADE_DIR, 'weights', 'config.yaml'), encoding='utf-8') as fh:
        cfg = {k: v['value'] for k, v in yaml.safe_load(fh).items()
               if isinstance(v, dict) and 'value' in v}

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    featurizer = RNAGraphFeaturizer(
        split='test', radius=cfg['radius'], top_k=cfg['top_k'], num_rbf=cfg['num_rbf'],
        num_posenc=cfg['num_posenc'], max_num_conformers=cfg['max_num_conformers'],
        noise_scale=cfg['noise_scale'], drop_prob_3d=cfg['drop_prob_3d'])
    model = gRNAde(
        node_in_dim=tuple(cfg['node_in_dim']), node_h_dim=tuple(cfg['node_h_dim']),
        edge_in_dim=tuple(cfg['edge_in_dim']), edge_h_dim=tuple(cfg['edge_h_dim']),
        num_layers=cfg['num_layers'], drop_rate=cfg['drop_rate'], out_dim=cfg['out_dim'])
    model.load_state_dict(torch.load(
        os.path.join(GRNADE_DIR, 'weights', 'gRNAde_drop3d@0.75_maxlen@500.h5'),
        map_location='cpu'))
    model = model.to(device).eval()

    _, raw = featurizer.featurize_from_pdb_file(args.pdb)
    native_seq = raw['sequence']
    if len(native_seq) <= args.dr_len:
        raise ValueError(f'结构序列 {len(native_seq)}nt 不大于 DR 长度 {args.dr_len}')
    partial = '_' * args.dr_len + native_seq[args.dr_len:]
    print(f'[grnade] 结构链 {len(native_seq)}nt; 设计 DR 前 {args.dr_len} 位, '
          f'spacer {len(native_seq) - args.dr_len} 位锁定为结构原生序列')

    letter_to_num = featurizer.letter_to_num
    idx = torch.tensor([letter_to_num[c] if c in letter_to_num else len(letter_to_num)
                        for c in partial], dtype=torch.long)
    logit_bias = F.one_hot(idx, num_classes=model.out_dim + 1).float()[:, :-1] * 100.0
    logit_bias = logit_bias.to(device)

    rows = []
    for _ in range(args.n_batches):
        temperature = float(np.random.uniform(args.temp_min, args.temp_max))
        data = featurizer(raw).to(device)
        samples, logits = model.sample(data, args.batch_size, temperature,
                                       logit_bias=logit_bias, return_logits=True)
        n_nodes = logits.shape[1]
        perplexity = torch.exp(F.cross_entropy(
            logits.view(args.batch_size * n_nodes, model.out_dim),
            samples.view(args.batch_size * n_nodes).long(), reduction='none',
        ).view(args.batch_size, n_nodes).mean(dim=1)).cpu().numpy()
        for s, ppx in zip(samples.cpu().numpy(), perplexity):
            rows.append((''.join(NUM_TO_LETTER[i] for i in s), float(ppx), temperature))

    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sequence', 'perplexity', 'temperature'])
        writer.writerows(rows)
    print(f'[grnade] 生成 {len(rows)} 条序列 -> {args.out}')


if __name__ == '__main__':
    main()
