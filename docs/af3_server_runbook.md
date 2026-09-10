# AlphaFold Server 提交操作手册（AF3 独立对照层，2026-09-05）

目的：把 `data/af3_inputs/` 的 13 套任务（WT + TOP-12 候选，SuCas12a2 + crRNA + 靶 RNA
三元复合物，蛋白链已强制 `useStructureTemplate=false`）提交到 AlphaFold Server，
结果回填后与 Chai-1 自模板口径对比，作为"界面可预测性"的首个独立证据层。

## 前置条件

- Google 账号；AlphaFold Server（https://alphafoldserver.com/ ）免费、仅限非商业用途。
- 任务额度：每账号每日有提交上限（历史上约 20 个/天，以页面实际提示为准）。
  本项目 13 个任务一天内可提完。
- 本层为 template-free 对照：输入 JSON 已禁用蛋白模板，**提交时不要再手动添加模板**。

## 提交步骤

方式 A（推荐，一次提完）：上传合并文件 `data/af3_inputs/GF_all13_alphafoldserver.json`
（顶层为 13 个任务的列表，服务器会拆成 13 个独立任务）。

方式 B（逐个提交）：在服务器页面选 "Upload JSON"，逐文件上传
`data/af3_inputs/GF_*.json`（13 个单任务文件，每个顶层为单元素列表）。

上传后服务器逐个排队运行（单任务通常数分钟到数十分钟，排队时间视负载）。
任务名 GF_xxx 会出现在 job history 表中。

## 结果回收

1. 每个任务完成后下载结果 zip（job history → Download）。
2. **不要解压改名**，直接把 13 个 zip 原样放进仓库 `data/af3_results/`目录
   （目录不存在就新建；该目录不入 git——结果文件大，属机时产物）。
3. 回填汇总：

```bash
export PYTHONUTF8=1   # Windows 必设
python scripts/crrna_af3_collect.py
```

产出 `data/af3_summary.json`：逐任务 5 模型的 ipTM / prot-crRNA / crRNA-靶 RNA
链对分数均值与 sd、vs WT 差量、程序化解读（阈值规则与 Chai 层一致：
prot-crRNA ipTM <0.5 低于可用区间，0.5–0.6 下沿弱参考，≥0.6 达可用区间）。

## 判读口径（预登记，防止事后叙事）

- **若 template-free AF3 下 WT 与候选的 prot-crRNA ipTM 普遍进入可用区间
  （≥0.5–0.6）**：Chai 矩阵的"界面恢复"获得首个独立引擎佐证，可表述为
  "两个独立引擎的蛋白-crRNA 界面预测一致"；
- **若 AF3 template-free 下界面同样低迷（<0.5，与 Chai 裸序列探针 0.02 接近）**：
  则 Chai 自模板高分确认是模板复述，"界面可预测"主张整体撤回，仅保留
  湿实验判据——这也是完整可交付的阴性结论；
- AF3 层分数只作**界面恢复与否**的定性对照，不用于骨架间排序主张
  （单构建观测、单随机种子、MSA 自动构建，噪声未标定）。

## 常见问题

- 上传报格式错误：确认上传的是 `data/af3_inputs/` 内 2026-09-05 之后生成的
  alphafoldserver 方言文件（顶层是 `[...]` 列表）。旧版 dialect=alphafold3
  文件（带链 id 字段）是给开源自托管 AF3 的，服务器不认，已废弃重发。
- 额度用完：分两天提交即可，任务相互独立。
- 结果 zip 内 `<job>_job_request.json` 记录了服务器实际使用的随机种子，
  需要复跑同一条件时以它为准。


## 本地 AF3 无 MSA 层（2026-09-10 补充，网络受阻期间的替代完成项）

alphafoldserver.com 谷歌托管网络不通（本机/A6000 均被阻断，无可用代理），
服务器提交暂缓。已用官方开源 alphafold3-3.0.4 + 官方 af3.bin.zst 权重
（storage.googleapis.com 直连公开可下，无需表单审批）在 A6000 容器完成
同一 13 套输入的本地推理：

- 口径：template-free 且**无 MSA**（蛋白/RNA 均单序列自比对输入，
  `--norun_data_pipeline`），seed 0 × 5 diffusion samples，triton 注意力，
  与 Protenix-v1 层同协议；
- 产物：`data/af3_local_summary.json`（`scripts/crrna_af3_collect.py --local`
  生成），结构快照 `data/af3_local_results/structures/`（每任务最佳样本 cif）；
- 结果：WT prot-crRNA ipTM 0.252（5 模型），全部 13 任务落在 0.23–0.33，
  低于可用区间（<0.5）——与 Protenix 分支二判读一致：无模板独立引擎下
  界面低迷，Chai 自模板高分属模板复述，维持"界面可预测"主张撤回、
  以湿实验判据为准的预登记结论；
- 服务器版（自动构建 MSA）与本地无 MSA 版是**不同口径**：前者仍是
  AlphaFold Server 恢复访问后的目标项，额度与流程见上文；
- 本地复跑：容器 `~/dzy/envs/af3/bin/python /dawn/af3_src/run_alphafold.py
  --json_path=/dawn/af3_jobs/<job>.json --model_dir=/dawn/af3_params
  --output_dir=<dir> --norun_data_pipeline`（输入 JSON 为 alphafold3 方言，
  含单序列自比对 MSA 字段；CCD pickle 已生成于该 env 的 site-packages，
  FetchContent 依赖经 ghfast.top git 代理拉取）。
