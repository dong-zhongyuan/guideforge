# 大文件清单(不进 git 历史)

以下文件为模型权重/基准数据, 存放在本仓库 Release `data-v1` 的资产中, 按相同相对路径放置即可。

| 路径 | 大小 | sha256 |
|---|---|---|
| `toolbox/OpenKnotAIDesignData/Data/OK7a_M2R_data.v4.5.1.csv` | 26534702 | 00d62f7d39b37e48… † |
| `toolbox/OpenKnotAIDesignData/Data/OK7a_M2_data.v4.5.2.csv` | 30216667 | 33b353538c836ed2… † |
| `toolbox/OpenKnotAIDesignData/Data/OpenKnotBench_data.v4.5.1.csv` | 139830358 | 22fa72df6d98f86f… † |
| `toolbox/Struct2SeQ/Struct2SeQ.pt` | 119736901 | e30294b2bce589dd… † |
| `toolbox/Struct2SeQ/Struct2SeQ_SHAPE.pt` | 119736901 | 63f6bfd12d4d4f7a… † |
| `toolbox/geometric-rna-design/weights/gRNAde_drop3d@0.75_maxlen@500.h5` | 8726092 | d43454deaeec7736… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Deg.pt` | 45409198 | 8df072e9992cdff5… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Drop.pt` | 45417124 | 8e423126b5686a48… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-SS.pt` | 45420400 | 626060952368affb… † |
| `toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet.pt` | 45406126 | c2aa45c14367863ece52d528d6c353ef40b66f7cb41539c19a042e87c7d3f215 ✅ |

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
