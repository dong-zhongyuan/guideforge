# 湿实验分析链彩排报告(数据全部 SIMULATED, 不含任何真实实验数值)

- 生成: scripts/crrna_sim_rehearsal.py, 2026-09-16, seed=42
- 情景: H1 成立(A1C≈WT) + H2 阳性(A8C+U15G 每靶 +12pp), 噪声 sd=8pp, n=5
- 基线: HCT116/R248Q WT 70.4% / SW480/R273H 45.0%(Scholz 标定中位)

## anova (rc=0)

```
=== survival 两因素方差分析(a=3 骨架, b=2 靶标, n=3 重复) ===
  scaffold     SS=    247.99 df= 2 F=   3.277 p=0.07319 
  target       SS=   2941.45 df= 1 F=  77.739 p=1.371e-06 **
  interaction  SS=     36.39 df= 2 F=   0.481 p=0.6297 

结果 -> /public/home/mengxl/dzy/guideforge/data/ivt_matrix_anova.json
```

## twin (rc=0)

```
H4 校验(实测存活 vs Scholz 标定预测, WT 行逐靶):
  TP53-R248Q   WT 实测  69.2% vs 预测  70.4% (偏差  -1.2 pp; 标定移植边界 20pp)
  TP53-R273H   WT 实测  46.5% vs 预测  45.0% (偏差  +1.5 pp; 标定移植边界 20pp)
(--cell-csv 的 SI 对比属 IVT 口径, 细胞-only 下忽略)
结果 -> data/virtual_cell_check.json
```

## power (rc=0)

```
功效预分析 -> /public/home/mengxl/dzy/guideforge/data/power_analysis_cellpanel.json
  成对 MDD(80%功效, n=3 vs 3): 37.6 pp
  交互 F 功效(10/15/20/30pp): {'10pp': 0.159, '15pp': 0.321, '20pp': 0.539, '30pp': 0.897}
  重复数扫描 MDD: {'n=3': 37.6, 'n=4': 30.2, 'n=5': 26.1, 'n=6': 23.3, 'n=8': 19.7}
  结论: >=MDD 的骨架差异可下显著性结论; 更小差异按方向性/一致性叙述; 预测差 10-30pp 需按扫描表考虑加重复
```

## bayes (rc=0)

```
[同源先验池化(n=46 观测: {'fig1g': 16, 'cis_end': 7, 'trans_end': 7, 'rbs0': 7, 'rbs33': 7, 'ordinal_4n96': 2})] 提案 4 条(EI 降序):
  WT                 EI=126.1566  ddG_dr=0.0
  A1C                EI=126.1566  ddG_dr=0.0
  A1U                EI=126.1566  ddG_dr=0.0
  A1G                EI=126.1566  ddG_dr=0.0
输出 -> data/bayesopt_proposals_priorhan.json
```
