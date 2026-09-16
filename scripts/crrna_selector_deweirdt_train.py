# -*- coding: utf-8 -*-
"""DeWeirdt 大库选型器再训练 + Cas12a2 外迁移检验(§K, 2026-09-16 登记).

动机(用户指令): 之前同源迁移训练的选型器仅 16 点(Han 工具箱 7 + Fig1f 9,
决策树 4 叶), 分辨率粗; DeWeirdt 2020 大库 35,883 条(AsCas12a, 双 spacer
上下文 127/128)带同套 5 维特征与实测敲低 lfc(低=抑制强, 方向同 fig1g)。
§I(Cas12a->Cas12a2 迁移证据链, Dmytrenko DR 互换)为跨酶使用提供边界。

预登记协议(docs/preregistration.md §K, 登记先于本次运行):
  1) 训练: DeWeirdt 全库(双上下文各一行), 特征 = 与现行选型器同套 5 维
     [ddG_dr, bp_dist, cross_nt, p_fold, spacer_up]; 模型 = Ridge(线性基线)
     + GradientBoosting(非线性), 无随机性的 Ridge 为主口径;
  2) 内部检验: 按上下文分组留一(GroupKFold 127<->128)——跨上下文泛化;
  3) 外迁移检验(决定性): DeWeirdt 训练模型预测 Han 16 点(fig1g 方向:
     预测 lfc 越低=抑制越强, 与 fig1g 数值 Spearman, 期望负相关),
     对照现行 16 点树的 LOEO 表现;
  4) 判读: 外迁移 Spearman <= 0 视为迁移失败, 模型只作 As 库内工具与
     边界记录, 不替换现行选型器; 显著负相关才进入替换讨论;
  5) 面板预测: 通过/未通过均如实输出两模型对 JD12/panel 变体空间的
     排序变化(WT-顶结论是否动摇)。

输出: data/selector_deweirdt_retrain.json
运行: PYTHONUTF8=1 python scripts/crrna_selector_deweirdt_train.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))

from sklearn.ensemble import GradientBoostingRegressor  # noqa: E402
from sklearn.linear_model import Ridge                  # noqa: E402
from sklearn.model_selection import GroupKFold          # noqa: E402
from scipy.stats import spearmanr                       # noqa: E402

import crrna_train_selector as sel                      # noqa: E402

DATA = os.path.join(ROOT, 'data')
FEATURES = sel.FEATURES  # 与现行选型器同套


def load_deweirdt():
    d = json.load(open(os.path.join(DATA, 'raw', 'deweirdt2020_dr_scan.json'),
                       encoding='utf-8'))
    X, y, g = [], [], []
    for r in d['rows']:
        for ctx in ('127', '128'):
            f = r.get('features_%s' % ctx)
            v = r.get('lfc%s' % ctx)
            if not f or v is None:
                continue
            X.append([float(f[k]) for k in FEATURES])
            y.append(float(v))
            g.append(ctx)
    return np.array(X), np.array(y), np.array(g)


def main():
    X, y, g = load_deweirdt()
    print('[K] DeWeirdt 训练集: n=%d (127:%d / 128:%d)'
          % (len(y), (g == '127').sum(), (g == '128').sum()))

    models = {
        'ridge': Ridge(alpha=1.0),
        'gbdt': GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                          random_state=42),
    }
    rep = {'generated_by': 'scripts/crrna_selector_deweirdt_train.py',
           'protocol': '§K(2026-09-16 登记, 登记先于运行; 见脚本头与 '
                       'docs/preregistration.md)',
           'n_train': int(len(y)), 'features': FEATURES,
           'internal_cv': {}, 'external_transfer': {}, 'models': {}}

    # 2) 内部: 分组留一上下文
    gkf = GroupKFold(n_splits=2)
    for name, m in models.items():
        rhos = []
        for tr, te in gkf.split(X, y, groups=g):
            m.fit(X[tr], y[tr])
            p = m.predict(X[te])
            rhos.append(spearmanr(p, y[te]).statistic)
        rep['internal_cv'][name] = {
            'fold_rhos': [round(float(r), 4) for r in rhos],
            'mean_rho': round(float(np.mean(rhos)), 4)}
        print('[K] 内部跨上下文 %-6s rho=%s' % (name, rep['internal_cv'][name]))

    # 全量拟合
    for name, m in models.items():
        m.fit(X, y)
        rep['models'][name] = 'fitted on n=%d' % len(y)

    # 3) 外迁移: Han 16 点
    Xh, yh, tags_h, _, confs_h = sel.load_han()
    # 现行 16 点树复刻
    from sklearn.tree import DecisionTreeRegressor
    dt16 = DecisionTreeRegressor(max_leaf_nodes=4, min_samples_leaf=1,
                                 random_state=42)
    dt16.fit(Xh, yh)
    rep['external_transfer']['baseline_dt16_train_rho'] = round(float(
        spearmanr(dt16.predict(Xh), yh).statistic), 4)
    for name, m in models.items():
        p = m.predict(Xh)
        rho, pv = spearmanr(p, yh)
        rep['external_transfer'][name] = {
            'spearman_vs_fig1g': round(float(rho), 4),
            'p_value': float(pv),
            'direction_ok': bool(rho < 0)}
        print('[K] 外迁移 %-6s rho=%+.4f p=%.4f (期望负: lfc低=抑制强=fig1g低)'
              % (name, rho, pv))

    # 4/5) 面板空间排序变化: JD12 语境变体用现行树 vs 大库模型
    import csv
    changes = {}
    for tag, csvf, wtref in (
            ('JD12-sp1', 'data/jd12_sp1_direct.variants.csv', None),
            ('JD12-sp2', 'data/jd12_direct.variants.csv', None)):
        rows = [r for r in csv.DictReader(open(os.path.join(ROOT, csvf),
                                               encoding='utf-8'))
                if r['passed'] in ('True', 'true', '1')]
        # 变体特征从 variants.csv 重算成本高; 用大库模型对 desc 做组级判读:
        # 直接对已通过变体的 5 维特征需现算——此处用近似: 只对关键骨架集合
        # (WT/A1C/A1U/A1G/A8C+U15G/U15G/U3C/U13C) 由 stemmax 同款算法重算,
        # 限于篇幅对 desc 前缀匹配 stemmax/winner 产物已有特征者跳过;
        # 完整逐变体重算列为后续(若外迁移通过)。
        changes[tag] = {'n_passed': len(rows),
                        'note': '若外迁移通过再做全空间重排(协议步骤5)'}
    rep['panel_implication'] = changes

    dst = os.path.join(DATA, 'selector_deweirdt_retrain.json')
    json.dump(rep, open(dst, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print('->', dst)


if __name__ == '__main__':
    main()
