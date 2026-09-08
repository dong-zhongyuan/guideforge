#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复算审计: 不 import crrna_homolog_context_interaction, 独立实现重算
E1 rho(DeWeirdt test 类跨方向 Spearman)与 E3 Kendall W(Tian 8x12),
对照 data/homolog_context_interaction.json 落盘值。"""
import json
import numpy as np
from scipy.stats import spearmanr, rankdata

RAW = 'data/raw/deweirdt2020_dr_scan.json'
OUT = 'data/homolog_context_interaction.json'

# --- E1 独立复算 ---
d = json.load(open(RAW))
test = [r for r in d['rows'] if r['cls'] == 'test']
x = np.array([r['lfc127'] for r in test], dtype=float)
y = np.array([r['lfc128'] for r in test], dtype=float)
rho, p = spearmanr(x, y)
# 独立 bootstrap(不同 seed)
rng = np.random.default_rng(20260908 + 1)
n = len(x)
vals = []
for _ in range(200):
    idx = rng.integers(0, n, n)
    r, _ = spearmanr(x[idx], y[idx])
    if np.isfinite(r):
        vals.append(float(r))
ci = [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]

# --- E3 独立复算 ---
h = json.load(open('data/han2025_dataset.json'))
by_crna = {}
for r in h['tian2025_rrs_pairs']:
    by_crna.setdefault(r['crRNA'], {})[r['tag'].split('_crRNA')[0]] = r['rrs_rel']
mut_names = sorted(by_crna[1].keys())
mat = np.array([[by_crna[c][m] for m in mut_names] for c in sorted(by_crna)])
ranks = np.vstack([rankdata(row) for row in mat])
k = ranks.shape[0]
col_sums = ranks.sum(axis=0)
mean = k * (12 + 1) / 2.0
W = 12.0 * np.sum((col_sums - mean) ** 2) / (k * k * (12 ** 3 - 12))

# E3b 独立复算: 位置 1 三突变(1C/1G/1U) × 8 crRNA 的 Kendall W
p1 = [m for m in mut_names if m.startswith('1')]
mat1 = np.array([[by_crna[c][m] for m in p1] for c in sorted(by_crna)])
ranks1 = np.vstack([rankdata(row) for row in mat1])
k1, n1 = ranks1.shape
col_sums1 = ranks1.sum(axis=0)
mean1 = k1 * (n1 + 1) / 2.0
W1 = 12.0 * np.sum((col_sums1 - mean1) ** 2) / (k1 * k1 * (n1 ** 3 - n1))

out = json.load(open(OUT))
e1_arch = out['deweirdt']['e1']
tian_arch = out['tian']['kendall_w']
pos1_arch = out['tian']['pos1']['kendall_w']

print('E1 独立重算 rho=%.6f p=%.2e CI200=%s' % (rho, p, ci))
print('E1 归档      rho=%.6f      CI=%s' % (e1_arch['rho_cross'], e1_arch['bootstrap_ci95']))
print('E3 独立重算 W=%.6f' % W)
print('E3 归档      W=%.6f' % tian_arch['W'])
print('E3b pos1 独立重算 W=%.6f' % W1)
print('E3b pos1 归档      W=%.6f' % pos1_arch['W'])

ok1 = abs(rho - e1_arch['rho_cross']) < 1e-9
ok3 = abs(W - tian_arch['W']) < 1e-9
ok3b = abs(W1 - pos1_arch['W']) < 1e-9
# bootstrap CI 不同 seed 允许小幅波动(数量级一致即通过)
ok_ci = abs(ci[1] - e1_arch['bootstrap_ci95'][1]) < 0.02
print('AUDIT', 'PASS' if (ok1 and ok3 and ok3b and ok_ci) else 'FAIL',
      '| rho match:', ok1, '| W match:', ok3, '| pos1 W match:', ok3b,
      '| CI close:', ok_ci)
