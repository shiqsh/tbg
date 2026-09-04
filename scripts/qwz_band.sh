#!/bin/zsh

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

source "$ROOT/python/.venv/bin/activate"

(
    cd "$ROOT/fortran"
    fpm run --target qwz_band
)

python "$ROOT/plot/orbital-magnetization/qwz/band.py"

# $0              当前脚本
# dirname "$0"    脚本所在目录
# /..             上一级目录
# cd ...          进入上一级
# pwd             得到绝对路径
# $(...)          把命令输出保存下来
# ROOT=...        保存为 ROOT