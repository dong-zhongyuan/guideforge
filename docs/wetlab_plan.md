# GuideForge 湿实验执行单（体外旁切初筛 + PFS 等位选择性鉴别）

对应赛道二评分项「实验验证与性能证据」（25 分）与管线闭环设计末环
（README 流程图：体外旁切初筛 → 细胞杀伤验证）。本文档为第一阶
段（体外）的完整执行方案。

## 1. 实验目标与评分要点映射

| 实验 | 验证主张 | 对应评分要点 |
|---|---|---|
| 实验一：体外旁切（trans-cleavage）初筛 | 10 条 panel 的三组假设（茎稳定化 / 种子区释放 / 综合 TOP）对切割活性的影响 | 核酸切割活性、剂量反应、重复实验 |
| 实验二：WT vs R248Q 双靶标鉴别 | 本项目核心机制主张——"突变产生 PFS（G>A）带来等位选择性" | 脱靶评估（等位级）、阳性/阴性对照 |
| 实验三（时间允许）：顺式切割凝胶验证 | 靶 RNA 的 cis 切割，旁切结果的独立佐证 | 核酸切割活性（第二读出） |

## 2. 合成 Panel（10 条，送单格式）

合成形式：**DNA 模板（HPLC 纯化），下游自做 T7 IVT 转 RNA**。
（直接合成 43nt RNA 贵 6–10 倍，初筛阶段不建议。）
模板 = T7 启动子 + DR(19nt) + spacer(24nt)。T7 启动子序列统一为
`TAATACGACTCACTATAGG`（具体拼接方式按所用 IVT 试剂盒说明，部分
试剂盒要求启动子后 G 起始，需核对）。

| 编号 | 名称 | 假设组 | DR + spacer（DNA 模板核心序列，5'→3'） |
|---|---|---|---|
| P01 | WT | 参考株（阳性基线） | AATTTCTACTGTTGTAGAT GTTCATGCCGCCCATGCAGGAACT |
| P02 | A8C+U15G | 茎稳定化 | AATTTCTCCTGTTGGAGAT GTTCATGCCGCCCATGCAGGAACT |
| P03 | A8G+U15C | 茎稳定化 | AATTTCTGCTGTTGCAGAT GTTCATGCCGCCCATGCAGGAACT |
| P04 | U7C+A16G | 茎稳定化 | AATTTCCACTGTTGTGGAT GTTCATGCCGCCCATGCAGGAACT |
| P05 | U13A+A16G | 种子区释放 | AATTTCTACTGTAGTGGAT GTTCATGCCGCCCATGCAGGAACT |
| P06 | U12G+G14A | 种子区释放 | AATTTCTACTGGTATAGAT GTTCATGCCGCCCATGCAGGAACT |
| P07 | A1C | 综合评分 TOP | CATTTCTACTGTTGTAGAT GTTCATGCCGCCCATGCAGGAACT |
| P08 | A1U | 综合评分 TOP | TATTTCTACTGTTGTAGAT GTTCATGCCGCCCATGCAGGAACT |
| P09 | A1G | 综合评分 TOP | GATTTCTACTGTTGTAGAT GTTCATGCCGCCCATGCAGGAACT |
| T01 | 靶标模板：TP53 WT 转录本片段 | — | 以 data/tp53_mrna_NM000546.fa 中靶点为中心取 ~300nt，加 T7 启动子（gBlock） |
| T02 | 靶标模板：TP53 R248Q 转录本片段 | — | 同上，以 data/tp53_r248q_mrna_NM000546.fa 为准（PFS 位 G>A） |

对照设计（随 panel 一并到位，不额外合成）：
- 阴性对照 1：无 crRNA（仅蛋白 + 靶标 + 探针）
- 阴性对照 2：无靶标（crRNA + 蛋白 + 探针）
- 阴性对照 3：WT crRNA × R248Q 错配 PFS 组合预期弱激活
- 阳性对照：P01（WT crRNA × R248Q 靶标，文献口径为有效向导）

## 3. 物料清单与预算

| 项目 | 规格建议 | 估算（RMB） |
|---|---|---|
| DNA 模板合成 ×10 | 43nt+T7 启动子，HPLC 纯化，各 2 OD | 600–1000 |
| 靶标 gBlock ×2 | ~300nt，含 T7 启动子 | 600–800 |
| T7 IVT 试剂盒 | HiScribe T7（NEB）或翌圣/近岸国产替代 | 1500–3000 |
| RNA 纯化 | RNeasy Mini 或磁珠法 | 300–600 |
| 荧光报告探针 | 5'-FAM-ssRNA(or ssDNA)-BHQ1-3'，5–8nt，1 OD | 400–600 |
| Cas12a2 蛋白 | 无货架产品，两条路线：①自行表达纯化（pET 构建 + Ni 柱，约 1500–3000 + 1–2 周人力）②合作实验室/公司渠道（约 2000–5000） | 1500–5000 |
| 缓冲体系与耗材 | 20 mM HEPES pH 7.5 / 150 mM KCl / 10 mM MgCl2 体系；无酶耗材、黑壁 96 孔板 | 500–800 |
| **合计** | | **约 5000–9000 元**（自表达蛋白路线） |

## 4. 实验一：体外旁切初筛（SOP 要点）

1. RNP 组装：Cas12a2 : crRNA = 1 : 1.2（摩尔比），室温 10 min。
2. 激活：加入靶标 RNA，剂量梯度 0 / 0.5 / 1 / 2 / 5 / 10 nM（6 个梯
   度，满足"剂量反应"评分要点）。
3. 旁切读出：加入 FQ 探针（终浓度 200 nM 过量），37°C 动力学读板，
   每 1–2 min 一个点，共 60–90 min。
4. 布局：10 条 panel × 6 剂量 × 3 重复 = 180 孔，加 4 组对照 ×3 重复
   ≈ 192 孔，一块半 96 孔板；建议分两块板，各板内自带全套对照（排
   除板间差异）。
5. 数据分析：线性段拟合 k_obs（初速），对每个 crRNA 做剂量-反应拟
   合；三重复给出均值 ± SD；组间比较用 ANOVA + 事后检验。统计方
   法预写入记录（满足"规范数据分析"要点）。

## 5. 实验二：WT vs R248Q 双靶标鉴别

1. 同一批 crRNA（P01–P09），分别以 T01（WT 转录本）与 T02（R248Q）
   为激活靶，旁切读出同实验一。
2. 关键读出：选择性窗口 = k_obs(R248Q) / k_obs(WT)。计算预测：PFS
   G>A 使 WT 转录本激活显著减弱（v6 三层特异性报告给出残余激活风
   险阈值，见 data/tp53_r248q_zengdr.v6.dr_specificity.json）。
3. 每条 crRNA 给出"激活效率 × 等位选择性"二维散点，作为候选排序的
   实测轴——这是 AI 设计主张的裁决图，也是区赛 PPT 的核心图。

## 6. 实验三（时间允许）：顺式切割凝胶验证

- 合成 5'-FAM 标记的短靶 RNA（~60nt），RNP + 靶 37°C 反应 0/5/15/
  30/60 min，尿素-PAGE 显影切割条带。作为旁切读出的独立佐证，堵
  评委"旁切不代表 cis 活性"的质疑。

## 6.5 第二阶段预案：细胞杀伤验证（区赛后启动）

启动判据：体外数据入库且区赛提交完成；仅带体外裁决出的最优 2–3 条
+ WT 对照（3–4 条规模），不全 panel 上细胞。

1. 体系：含 TP53-R248Q 杂合突变的肿瘤细胞，RNP（Cas12a2 + 候选
   crRNA）电转递送；野生型 TP53 细胞作对照。
2. 读出：
   - 细胞活力/杀伤曲线（CellTiter-Glo 或同类，剂量梯度 × 时间点）
   - 靶转录本切割验证（qPCR 或靶向测序，确认杀伤来自靶切割）
3. 对照：WT crRNA、无 crRNA、未处理；与体外对照体系一致。
4. 预算：细胞培养 + 电转/转染试剂 + 读出试剂，约 3000–8000 元；
   周期 3–4 周（含递送条件优化）。
5. 评分收益：体外 + 细胞两级证据齐备时，「实验验证与性能证据」
   （25 分）目标区间 18–22 分；细胞数据同时是「候选工具实质性改
   造」（20 分）的直接功能证据。

## 7. 里程碑倒排（以提交截止日 D 为锚）

| 时间 | 事项 |
|---|---|
| D-4 周 | 下单 DNA 模板 + gBlock；敲定 Cas12a2 蛋白来源（最长杆，最先做） |
| D-3 周 | 蛋白表达纯化启动 / 渠道到货；IVT 试剂与探针到货 |
| D-2 周 | IVT 制备 10 条 crRNA + 2 条靶标；纯化定量；预实验摸条件（P01 单点） |
| D-1 周 | 实验一 + 实验二正式上板；数据分析；裁决图产出 |
| D | 数据归档入库（raw 板读数 + 分析脚本），更新 README 闭环状态 |

## 8. 入库与可追踪要求（对应"实验记录可追踪"）

- 原始读板数据（导出 CSV）、IVT 质检（浓度/胶图）放入
  `data/wetlab/`，命名与干实验产出同规则（带日期与版本）。
- 分析脚本放 `scripts/`，复用主仓库的可复现惯例（固定随机性、输
  出含参数血缘）。
- 每个实验在本文档末尾补「执行记录」小节：日期、人员、批号、偏
  离项。

## 执行记录

（待填）
