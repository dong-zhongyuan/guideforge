#!/bin/bash
# MD 结果回传器(本机执行; 2026-09-02 正式化——此前 data/md/*/md_perf.log 由未入库的一次性脚本产出, 纪律违规)
# 用法(本机 Git Bash): bash scripts/crrna_md_retrieve.sh <tag...>
# 作用: 容器 ~/md_runs/<tag>/analysis.json 与 md.log 性能摘要 经本机双跳管道回传宿主入库。
#   轨迹文件留在容器(体积), 仅回传 KB 级结果。
set -u
SSH=C:/Windows/System32/OpenSSH/ssh.exe
for t in "$@"; do
  $SSH a6000-ct "cat ~/md_runs/$t/analysis.json" \
    | $SSH srv "cat > /public/home/mengxl/dzy/guideforge/data/md/$t/analysis.json"
  $SSH a6000-ct "grep -E 'Finished|Performance' ~/md_runs/$t/md.log | head -2" \
    | $SSH srv "cat > /public/home/mengxl/dzy/guideforge/data/md/$t/md_perf.log"
  echo "$t: analysis.json + md_perf.log 回传完成"
done
