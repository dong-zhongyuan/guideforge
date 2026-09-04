# -*- coding: utf-8 -*-
"""按 BIG_FILES.md 从 Gitee Release `data-v1` 下载 10 件权重/基准数据。

- 3 件超 100MB 的文件在 Gitee 侧为 part00/part01 分片, 本脚本自动下载并拼回;
- sha256 校验: 标 ✅ 的条目做全量 64 位比对; 标 † 的条目原记录只有前 16 位
  十六进制截断值, 只能做前缀比对——输出中会如实声明"仅前缀比对, 不构成完整
  完整性校验"。
- 已存在且校验通过的文件自动跳过; 校验失败的文件删除后重下。

用法:
  python tools/fetch_big_files.py            # 下载全部缺失项
  python tools/fetch_big_files.py --check    # 只校验本地已有文件, 不下载
  python tools/fetch_big_files.py --base-url <url>  # 换源(如 GitHub Release)
"""
import argparse
import hashlib
import os
import sys
import urllib.parse
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

DEFAULT_BASE = ("https://gitee.com/eastern-zhongyuan/guideforge/"
                "releases/download/data-v1/")

# (相对路径, 大小字节, sha256, 是否完整哈希(False=仅前16位截断前缀), 资产名或None=分片)
FILES = [
    ("toolbox/OpenKnotAIDesignData/Data/OK7a_M2R_data.v4.5.1.csv",
     26534702, "00d62f7d39b37e48", False,
     "toolbox__OpenKnotAIDesignData__Data__OK7a_M2R_data.v4.5.1.csv"),
    ("toolbox/OpenKnotAIDesignData/Data/OK7a_M2_data.v4.5.2.csv",
     30216667, "33b353538c836ed2", False,
     "toolbox__OpenKnotAIDesignData__Data__OK7a_M2_data.v4.5.2.csv"),
    ("toolbox/OpenKnotAIDesignData/Data/OpenKnotBench_data.v4.5.1.csv",
     139830358, "22fa72df6d98f86f", False, None),  # 分片
    ("toolbox/Struct2SeQ/Struct2SeQ.pt",
     119736901, "e30294b2bce589dd", False, None),  # 分片
    ("toolbox/Struct2SeQ/Struct2SeQ_SHAPE.pt",
     119736901, "63f6bfd12d4d4f7a", False, None),  # 分片
    ("toolbox/geometric-rna-design/weights/gRNAde_drop3d@0.75_maxlen@500.h5",
     8726092, "d43454deaeec7736", False,
     "gRNAde_drop3d@0.75_maxlen@500.h5"),
    ("toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Deg.pt",
     45409198, "8df072e9992cdff5", False,
     "toolbox__rnet-inference__RibonanzaNet-Weights__RibonanzaNet-Deg.pt"),
    ("toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-Drop.pt",
     45417124, "8e423126b5686a48", False,
     "toolbox__rnet-inference__RibonanzaNet-Weights__RibonanzaNet-Drop.pt"),
    ("toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet-SS.pt",
     45420400, "626060952368affb", False,
     "toolbox__rnet-inference__RibonanzaNet-Weights__RibonanzaNet-SS.pt"),
    ("toolbox/rnet-inference/RibonanzaNet-Weights/RibonanzaNet.pt",
     45406126,
     "c2aa45c14367863ece52d528d6c353ef40b66f7cb41539c19a042e87c7d3f215", True,
     "toolbox__rnet-inference__RibonanzaNet-Weights__RibonanzaNet.pt"),
]

PARTS = {  # 分片资产名(按拼接顺序)
    "toolbox/OpenKnotAIDesignData/Data/OpenKnotBench_data.v4.5.1.csv":
        ["OpenKnotBench_data.v4.5.1.csv.part00",
         "OpenKnotBench_data.v4.5.1.csv.part01"],
    "toolbox/Struct2SeQ/Struct2SeQ.pt":
        ["Struct2SeQ.pt.part00", "Struct2SeQ.pt.part01"],
    "toolbox/Struct2SeQ/Struct2SeQ_SHAPE.pt":
        ["Struct2SeQ_SHAPE.pt.part00", "Struct2SeQ_SHAPE.pt.part01"],
}


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(rel, size, digest, full):
    path = os.path.join(ROOT, rel)
    if not os.path.isfile(path):
        return "缺失"
    if os.path.getsize(path) != size:
        return "大小不符(%d != %d)" % (os.path.getsize(path), size)
    actual = sha256_of(path)
    if full:
        return "OK" if actual == digest else "sha256 不符"
    if actual.startswith(digest):
        return "OK(仅前16位前缀比对, 原记录为截断值, 不构成完整完整性校验)"
    return "sha256 前缀不符(前16位 %s != %s)" % (actual[:16], digest)


def download(url, dest):
    tmp = dest + ".part_download"
    print("  下载 %s" % url, flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "guideforge-fetch"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
        n = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
            print("\r    %.1f MB" % (n / 1e6), end="", flush=True)
    print()
    os.replace(tmp, dest)


def fetch_one(rel, size, digest, full, asset, base):
    path = os.path.join(ROOT, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if asset is not None:
        download(base + urllib.parse.quote(asset), path)
    else:
        part_paths = []
        for part in PARTS[rel]:
            pp = os.path.join(os.path.dirname(path), part)
            download(base + urllib.parse.quote(part), pp)
            part_paths.append(pp)
        with open(path + ".assembling", "wb") as out:
            for pp in part_paths:
                with open(pp, "rb") as f:
                    for chunk in iter(lambda: f.read(1 << 20), b""):
                        out.write(chunk)
        os.replace(path + ".assembling", path)
        for pp in part_paths:
            os.remove(pp)
        print("  已拼回 %d 个分片 -> %s" % (len(part_paths), rel))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="只校验, 不下载")
    ap.add_argument("--base-url", default=DEFAULT_BASE,
                    help="Release 资产下载基址(默认 Gitee data-v1)")
    args = ap.parse_args()

    n_fail = 0
    for rel, size, digest, full, asset in FILES:
        status = verify(rel, size, digest, full)
        tag = "✅全量" if full else "†前缀"
        if status.startswith("OK"):
            print("[跳过] %s (%s, 校验%s)" % (rel, tag, status))
            continue
        if args.check:
            print("[FAIL] %s: %s" % (rel, status))
            n_fail += 1
            continue
        print("[下载] %s" % rel)
        try:
            fetch_one(rel, size, digest, full, asset, args.base_url)
        except Exception as e:  # noqa: BLE001
            print("[FAIL] %s: 下载失败 %s" % (rel, e))
            n_fail += 1
            continue
        status = verify(rel, size, digest, full)
        if status.startswith("OK"):
            print("[OK] %s (%s, %s)" % (rel, tag, status))
        else:
            print("[FAIL] %s: 下载后校验失败: %s" % (rel, status))
            n_fail += 1
    if n_fail:
        print("失败 %d 项" % n_fail)
        sys.exit(1)
    print("全部 10 件就绪。注意: 标 † 的 9 件只做了 sha256 前 16 位前缀比对"
          "(原记录为截断值), 仅 RibonanzaNet.pt 为全量 64 位校验。")


if __name__ == "__main__":
    main()
