# -*- coding: utf-8 -*-
"""下载 GENCODE 转录组 FASTA(主管线特异性扫描依赖, 不进交付包)。

- GENCODE v47 transcripts -> data/raw/gencode.v47.transcripts.fa
- GENCODE v49 pc + lncRNA -> data/transcriptome/gencode.v49.{pc,lncRNA}_transcripts.fa
优先下载 .gz 压缩版再解压; 校验解压后文件大小(与原始下载实测值比对),
并把来源 URL 写入 data/download_log.json。

用法: python tools/fetch_transcriptome.py [--force]
"""
import argparse
import gzip
import json
import os
import shutil
import sys
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

EBI = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

# (目标相对路径, gz 来源 URL, 解压后期望大小字节(原始下载实测))
FILES = [
    ("data/raw/gencode.v47.transcripts.fa",
     EBI + "/release_47/gencode.v47.transcripts.fa.gz", 613241416),
    ("data/transcriptome/gencode.v49.pc_transcripts.fa",
     EBI + "/release_49/gencode.v49.pc_transcripts.fa.gz", 683175894),
    ("data/transcriptome/gencode.v49.lncRNA_transcripts.fa",
     EBI + "/release_49/gencode.v49.lncRNA_transcripts.fa.gz", 223740848),
]


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


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="已存在也重下")
    ap.add_argument("--keep-gz", action="store_true", help="保留下载的 .gz")
    args = ap.parse_args()

    log_path = os.path.join(ROOT, "data", "download_log.json")
    log = {"source": "tools/fetch_transcriptome.py", "files": []}
    n_fail = 0
    for rel, url, expect_size in FILES:
        dest = os.path.join(ROOT, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if os.path.isfile(dest) and not args.force:
            size = os.path.getsize(dest)
            status = "OK" if size == expect_size else "大小不符"
            print("[跳过] %s (%d 字节, %s)" % (rel, size, status))
            log["files"].append({"path": rel, "url": url, "size": size,
                                 "expected_size": expect_size,
                                 "status": "已存在-" + status})
            if status != "OK":
                n_fail += 1
            continue
        gz = dest + ".gz"
        try:
            download(url, gz)
            with gzip.open(gz, "rb") as fi, open(dest + ".inflating", "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
            os.replace(dest + ".inflating", dest)
            if not args.keep_gz:
                os.remove(gz)
        except Exception as e:  # noqa: BLE001
            print("[FAIL] %s: %s" % (rel, e))
            n_fail += 1
            continue
        size = os.path.getsize(dest)
        ok = size == expect_size
        print("[%s] %s (%d 字节, 期望 %d)" % ("OK" if ok else "FAIL",
                                             rel, size, expect_size))
        log["files"].append({"path": rel, "url": url, "size": size,
                             "expected_size": expect_size,
                             "status": "OK" if ok else "大小不符",
                             "fetched_at": time.strftime("%Y-%m-%d %H:%M:%S")})
        if not ok:
            n_fail += 1
    with open(log_path, "w", encoding="utf-8") as f:
        json.dump(log, f, ensure_ascii=False, indent=1)
    print("来源记录 -> %s" % log_path)
    if n_fail:
        print("失败 %d 项" % n_fail)
        sys.exit(1)


if __name__ == "__main__":
    main()
