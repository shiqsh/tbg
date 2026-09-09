#!/usr/bin/env python3
"""一条命令完成 BM 能带计算和绘图。

在项目根目录先执行 source python/.venv/bin/activate，随后运行：
    ./python/scripts/tbg/run_tbg_bands.py

两个任务按顺序执行：先生成最新能带数据，再读这份数据画图。
"""

import argparse
from pathlib import Path
import subprocess
import sys


# 当前脚本在 python/scripts/tbg/；绘图统一放到项目根目录的 plot/tbg/。
SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = PYTHON_ROOT.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-show', action='store_true', help='计算并保存图片，不弹出窗口')
    parser.add_argument('--data', type=Path,
                        default=PROJECT_ROOT / 'results' / 'tbg_bm_bands.npz',
                        help='本次计算保存的数据文件，也是绘图脚本的输入')
    parser.add_argument('--figure', type=Path,
                        default=PROJECT_ROOT / 'figure' / 'tbg' / 'bm_bands.png',
                        help='图片输出路径')
    args = parser.parse_args()

    # shebang 中的 env 从当前 PATH 寻找 python3。
    # 激活 .venv 后，PATH 中最先找到的是虚拟环境内的 Python。
    # sys.executable 是当前解释器路径，两个子脚本都使用它，保持环境一致。
    python = sys.executable
    print(f'Python 解释器：{python}', flush=True)

    # subprocess.run 会等待子脚本运行结束，而不是在后台启动后立刻返回。
    # check=True 表示出错就停止，因此计算失败时不会拿旧数据继续绘图。
    try:
        print('[1/2] 计算能带并保存数据……', flush=True)
        subprocess.run(
            [python, str(SCRIPT_DIR / 'tbg_band_structure.py'),
             '--output', str(args.data.resolve())],
            check=True,
        )

        print('[2/2] 读取数据并绘图……', flush=True)
        plot_command = [
            python, str(PROJECT_ROOT / 'plot' / 'tbg' / 'tbg_bands.py'),
            '--input', str(args.data.resolve()),
            '--output', str(args.figure.resolve()),
        ]
        if args.no_show:
            plot_command.append('--no-show')
        subprocess.run(plot_command, check=True)
    except subprocess.CalledProcessError as error:
        raise SystemExit(error.returncode) from None


if __name__ == '__main__':
    main()
