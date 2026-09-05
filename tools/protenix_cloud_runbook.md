# Protenix 云端运行手册(AF3 开源替代通道, 2026-09-05)

目的: 跑 13 套 template-free 共折叠(SuCas12a2 + crRNA 变体 + 靶 RNA),
作 Chai-1 自模板口径的独立对照层。引擎: Protenix-v1(`protenix_base_default_v1.0.0`,
368M 参数, 训练数据截止 2021-09-30 与 AF3 对齐, Apache-2.0 代码+权重,
bioRxiv 2026.02.05.703733)——DeepMind AF3 权重 gated 不可得时的合法 AF3 级替代。

## 1. 租机器

- AutoDL(或同类): 选 **RTX 4090 24GB** 起步(约 ¥2/小时), 稳妥选 A100 40GB
- 复合物规模 1279 残基(1207aa 蛋白 + 43nt + 29nt), 单折数分钟级; 13 套 < 2 小时
- 镜像选 PyTorch 2.x + CUDA 12.x 即可

## 2. 安装(云端)

```bash
pip install --upgrade protenix -i https://pypi.org/simple
# 权重: 首次 protenix pred 自动下载; 若失败见 GitHub README "Supported Models" 手动下载
```

## 3. 上传输入并批量跑

```bash
# 本地上传(在本机执行):
scp -r data/protenix_inputs root@<云机>:~/
# 云端执行:
mkdir -p ~/protenix_out
for f in ~/protenix_inputs/PX_*.json; do
  protenix pred -i "$f" -o ~/protenix_out -n protenix_base_default_v1.0.0
done
# MSA: 默认走内置检索(需联网); 若 MSA 服务不可达, 见官方 inference 文档的 MSA 开关。
# template-free 已由输入保证: JSON 不含 templates 字段。
```

## 4. 结果回传

```bash
# 云端打包: tar czf protenix_out.tgz -C ~ protenix_out
# 本机: scp 拉回后解到 data/protenix_results/(此目录当前为空, 待回填)
```

回填后通知 AI 会话: 写收集脚本聚合 13 套的 ipTM/PAE 到
`data/protenix_matrix_summary.json`, 与 `data/chai_matrix_4t_summary.json`
做跨引擎一致性对照, 并更新 README 结构层条目。

## 5. 记录义务(附件5 第三方平台条款)

提交材料的"第三方模型调用记录"里登记: Protenix 版本/权重版本/调用日期/
输入文件清单(见 data/protenix_inputs/manifest.json)/随机种子(输出目录自带)。
