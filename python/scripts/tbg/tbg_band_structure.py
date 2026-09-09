#!/usr/bin/env python3
"""BM 能带计算的运行入口：设置参数 -> 调用模型 -> 保存原始数据。

从任意目录执行本文件均可；默认结果写到项目根目录 results/tbg_bm_bands.npz。
本脚本不画图，绘图入口是项目根目录下的 plot/tbg/tbg_bands.py。
"""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np


# 脚本位于 python/scripts/tbg/，可复用包位于 python/src/。
# 通过 __file__ 定位源码，不依赖终端当前目录，也不要求提前 pip install 项目。
# parents[0] 是 scripts/tbg，parents[1] 是 scripts，parents[2] 才是 python。
PYTHON_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = PYTHON_ROOT.parent
sys.path.insert(0, str(PYTHON_ROOT / 'src'))
from tbg import BMModel, BMParameters, make_k_path, solve_bands


# ===================== 修改计算参数的位置 =====================
PARAMETERS = BMParameters(
    theta_deg=1.09,       # 度；保留原 hbar*v_F 时的近魔角示例
    w0_mev=110.0,         # AA/BB 层间耦合，meV
    w1_mev=110.0,         # AB/BA 层间耦合，meV
    cutoff=4,            # 平面波截断；总矩阵维数 4*(2*cutoff+1)^2
)
POINTS_PER_SEGMENT = 60   # 三段路径，总点数 3*60+1=181
# =============================================================


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=PROJECT_ROOT / 'results' / 'tbg_bm_bands.npz',
                        help='结果文件路径；多组参数比较时可使用不同文件名')
    args = parser.parse_args()

    model = BMModel(PARAMETERS)
    k_points, distance, ticks = make_k_path(model, POINTS_PER_SEGMENT)
    energies = solve_bands(model, k_points)

    # 对单谷、单自旋 BM 模型，能量排序中间的两条带是中心能带。
    # ptp = max-min；这里打印的是沿此高对称路径的带宽，不是全 BZ 的严格极值。
    mid = model.dimension // 2
    widths = np.ptp(energies[:, mid - 1:mid + 1], axis=0)
    print(f'矩阵维数：{model.dimension}；动量点数：{len(k_points)}')
    print(f'中心两带沿路径带宽：{widths[0]:.6f}、{widths[1]:.6f} meV')

    # 将计算参数和单位一起保存，避免以后只有能带数组却忘记其对应模型。
    # JSON 字符串可以在 allow_pickle=False 的条件下读取，不保存 Python 对象。
    metadata = {
        'model': 'single-valley, single-spin BM',
        'parameters': asdict(PARAMETERS),
        'points_per_segment': POINTS_PER_SEGMENT,
        'energy_unit': 'meV',
        'momentum_unit': 'nm^-1',
        'energy_reference': 'raw eigenvalues, no energy shift',
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('wb') as handle:
        np.savez_compressed(
            handle, k_points_nm_inv=k_points, distance_nm_inv=distance,
            ticks_nm_inv=ticks, energies_mev=energies,
            metadata_json=json.dumps(metadata, ensure_ascii=False),
        )
    print(f'数据已保存：{output}')


# 只有直接运行本脚本时才启动计算；被其他程序 import 时不会扫描动量。
if __name__ == '__main__':
    main()
