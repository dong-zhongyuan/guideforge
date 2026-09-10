# 大文件清单(不进 git 历史)

以下文件为模型权重/基准数据, 存放在本仓库 Release `data-v1` 的资产中, 按相同相对路径放置即可。

| 路径 | 大小 | sha256 |
|---|---|---|
| `toolbox/OpenKnotAIDesignData/Data/OK7a_M2R_data.v4.5.1.csv` | 26534702 | 00d62f7d39b37e48507ca5a949ec8358121b306b3d680b2976c4369b3038a31b… † |
| `toolbox/OpenKnotAIDesignData/Data/OK7a_M2_data.v4.5.2.csv` | 30216667 | 33b353538c836ed28a366d3f5d5b2ddb10c378108505001ede3ba43781681a83… † |
| `toolbox/OpenKnotAIDesignData/Data/OpenKnotBench_data.v4.5.1.csv` | 139830358 | 22fa72df6d98f86fbf37a57e679f717628a5d5db3737df3f5519901bc9777901… † |
| `toolbox/Struct2SeQ/Struct2SeQ.pt` | 119736901 | e30294b2bce589dd6b1a5502d2c04fec96476e726f50ba6fa5807855ae791b80… † |
| `toolbox/Struct2SeQ/Struct2SeQ_SHAPE.pt` | 119736901 | 63f6bfd12d4d4f7ac5d81a3c31c3afb9dff533dbcd404064296836951b5294ed… † |
| `toolbox/geometric-rna-design/weights/gRNAde_drop3d@0.75_maxlen@500.h5` | 8726092 | d43454deaeec773644bfa52385b6fcc8264db256878fd9ca19ed8711ec62971c… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Deg.pt` | 45409198 | 8df072e9992cdff546f831095d8e518f7ba2421cbde371b0dad45a9f1476af37… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Drop.pt` | 45417124 | 8e423126b5686a48499d0a22d8b234e01092b757fe773ed4e6519dd8c1cf0bdb… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-SS.pt` | 45420400 | 626060952368affbf61b78b532d6166387094754b68bc0553da376f2d2b00d56… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet.pt` | 45406126 | c2aa45c14367863ece52d528d6c353ef40b66f7cb41539c19a042e87c7d3f215ce52d528d6c353ef40b66f7cb41539c19a042e87c7d3f215 ✅ |

✅ = 本地原件实测: 2026-09-04 用 hashlib 对本机文件重算, 完整 64 位哈希与原截断前缀一致, 可直接用于下载后完整性校验。

† = 本地缺原件, 校验值取自 Release `data-v1` 上传记录, 且原记录本身为**截断值**(仅前 16 位十六进制): 只能做前缀比对, 不足以做完整性校验; 完整 64 位哈希需从上传机重算后补录(切勿凭前缀推造)。

## 国内兜底(2026-09-02): 3 个超大文件已分片上传 Gitee Release

Gitee 单附件上限 100MB, 下列 3 件按 95MB 分片(part00/part01)附于 Gitee Release data-v1,
下载全部分片后拼回并按上表 sha256 校验(注意: 对应 3 行目前同为 † 截断前缀, 仅可前缀比对):

| 文件 | 分片 | 拼回命令 |
|---|---|---|
| OpenKnotBench_data.v4.5.1.csv | part00(95MB)+part01(38MB) | `cat OpenKnotBench_data.v4.5.1.csv.part* > OpenKnotBench_data.v4.5.1.csv` |
| Struct2SeQ.pt | part00(95MB)+part01(19MB) | `cat Struct2SeQ.pt.part* > Struct2SeQ.pt` |
| Struct2SeQ_SHAPE.pt | part00(95MB)+part01(19MB) | `cat Struct2SeQ_SHAPE.pt.part* > Struct2SeQ_SHAPE.pt` |

拼回后 sha256 已验证与原件一致(22fa72df…/e30294b2…/63f6bfd1…, 均为截断前缀, 见上表 † 说明)。
至此 10 件资产 GitHub(全量单件) 与 Gitee(7 单件 + 3 分片) 双源齐备。

## Release data-v2 (2026-09-10): Protenix-v1 13 折共折叠原始结果归档

| 路径 | 打包文件 | 大小(bytes) | sha256 |
|---|---|---|---|
| `data/protenix_results/` (13 折目录, 原 136MB) | `protenix_results_13fold.tar.gz` | 28950476 | 3b4eb0f7c6ea7b5c03639a3cf4c9b6ed55153cb6219a18651de120f60993ffba ✅ |

✅ = 完整 64 位哈希三重验证: (1)上传前在服务器对原件 tar.gz 实测; (2)与
GitHub Release 资产 API 返回的服务端 `digest` 字段逐字节一致; (3)从 Gitee
Release 下载回读再算一次, 逐字节一致(2026-09-10, Release data-v2 双源)。
下载后直接 `sha256sum` 对照即可。

内容: Protenix-v1 (protenix_base_default_v1.0.0) 13 折(WT + 12 变体)
共折叠原始输出, 每折含 `seed_101/`(5 样本结构+置信度) 与 `msa/` 目录;
汇总层为 `data/protenix_summary.json`(commit 48d7796, 已入 git)。
解包: `tar xzf protenix_results_13fold.tar.gz -C data/protenix_results/`。
生成链: 容器 2026-09-06 批跑 → 2026-09-09 回收汇总 → 2026-09-10 归档双源上传。
