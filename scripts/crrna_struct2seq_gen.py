# -*- coding: utf-8 -*-
"""Struct2SeQ 逆折叠生成包装器（crRNA 链 · 策略A2 引擎）

由 crrna_scaffold_design.py --use-struct2seq 经子进程调用（须在 rnet 环境运行，
见 toolbox/demo/run_struct2seq.py 的 CPU 兼容说明）。本包装器：
  1. 无 CUDA 时把 .cuda() 降级为 no-op（官方 generate.py 硬编码 .cuda()）
  2. 构造单条靶结构 CSV（WT DR 茎环点括号 + WT 序列），调用官方 generate.py
  3. 官方脚本按 WT-upweight + Jaccard（RNet-SS 预测结构 vs 靶结构）筛选输出

用法:
  python crrna_struct2seq_gen.py --target-csv t.csv --out-dir out/ --up-bias 0.85
输出: <out-dir>/results.csv（列: sequence, predicted_structure, jaccard, structure_match, ...）
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S2S_DIR = os.path.join(ROOT, 'toolbox', 'Struct2SeQ')
RNET_W = os.path.join(ROOT, 'toolbox', 'rnet-inference', 'RibonanzaNet-Weights')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target-csv', required=True, help='列为 Title,Dot-bracket,wild_type_sequence')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--up-bias', type=float, default=0.85, help='向 WT 序列偏倚的强度(0-1)')
    ap.add_argument('--weights', default=os.path.join(S2S_DIR, 'Struct2SeQ_SHAPE.pt'))
    args = ap.parse_args()

    os.chdir(S2S_DIR)  # generate.py 用相对路径读 test10_configs/
    sys.path.insert(0, S2S_DIR)

    # arnie 在 import 时(load_package_locations)就会找 arnie_file.txt(相对 CWD 或
    # ARNIEFILE), 先于 Env.py 的自建逻辑; 必须在任何 arnie import 前备好绝对路径配置。
    # TMP 放包内 _arnie_tmp(服务器 /tmp 不一定可写, 且避免移交后路径失效)
    arnie_tmp = os.path.join(S2S_DIR, '_arnie_tmp')
    os.makedirs(arnie_tmp, exist_ok=True)
    arnie_cfg = os.path.join(S2S_DIR, 'arnie_file.txt')
    with open(arnie_cfg, 'w', encoding='utf-8') as fh:
        fh.write('linearpartition: . \nTMP: ' + arnie_tmp + '\n')
    os.environ['ARNIEFILE'] = arnie_cfg

    import torch
    if not torch.cuda.is_available():
        torch.Tensor.cuda = lambda self, *a, **k: self  # type: ignore
        torch.nn.Module.cuda = lambda self, *a, **k: self  # type: ignore
        _orig_load = torch.load
        torch.load = lambda *a, **k: _orig_load(  # type: ignore
            *a, **dict(k, map_location=k.get('map_location', 'cpu')))

    os.makedirs(args.out_dir, exist_ok=True)
    sys.argv = [
        'generate.py',
        '--target_df', os.path.abspath(args.target_csv),
        '--output_csv', 'results.csv',
        '--out_folder', os.path.abspath(args.out_dir),
        '--gpu_id', '0', '--n_structures', '1',
        '--up_bias', str(args.up_bias),
        '--weights_path', args.weights,
        '--rnet_weights', os.path.join(RNET_W, 'RibonanzaNet.pt'),
        '--rnet_ss_weights', os.path.join(RNET_W, 'RibonanzaNet-SS.pt'),
    ]
    with open(os.path.join(S2S_DIR, 'generate.py'), encoding='utf-8') as fh:
        code = compile(fh.read(), 'generate.py', 'exec')
    exec(code, {'__name__': '__main__',
                '__file__': os.path.join(S2S_DIR, 'generate.py')})


if __name__ == '__main__':
    main()
