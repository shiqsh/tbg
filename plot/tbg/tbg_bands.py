#!/usr/bin/env python3
"""读取已保存的 BM 能带数据并画图；不重新构造或对角化哈密顿量。

先运行 python/scripts/tbg/tbg_band_structure.py，再运行本脚本。
可以独立修改颜色、线宽和能量范围，然后重新运行本脚本。
"""

import argparse
import json
from pathlib import Path

import numpy as np


# 本脚本位于项目根目录的 plot/tbg/，不是 python/ 内部。
# 根据文件自身的位置定位 results/ 和 figure/，从任意工作目录执行都可找到它们。
# 绘图只依赖保存的数据，不需要导入 python/src/tbg 中的物理模型。
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def plot_bands(distance, ticks, energies_mev, metadata):
    """返回 Figure；输入来自计算脚本生成的 npz 文件。"""
    import matplotlib.pyplot as plt

    dimension = energies_mev.shape[1]
    mid = dimension // 2

    # make_k_path 的最后一点是 K_t。取该点中心两带的平均能量为 E_D，
    # 对所有动量、所有能带减去同一个数，只改变能量原点，不改变色散或带宽。
    dirac_energy = energies_mev[-1, mid - 1:mid + 1].mean()
    E = energies_mev - dirac_energy
    center = E[:, mid - 1:mid + 1]
    widths = np.ptp(center, axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax in axes:
        # 灰线展示中心附近的远端能带；蓝线、红线强调中心两带。
        for band in range(max(0, mid - 7), min(mid + 7, dimension)):
            ax.plot(distance, E[:, band], color='#9ba9b9', lw=1)
        ax.plot(distance, center[:, 0], color='#2166ac', lw=2)
        ax.plot(distance, center[:, 1], color='#d6604d', lw=2)
        for x in ticks:
            ax.axvline(x, color='0.85', lw=0.8, zorder=0)
        ax.axhline(0, color='0.55', ls='--', lw=0.7, zorder=0)
        ax.set(xticks=ticks, xticklabels=[r'$K_b$', r'$\Gamma$', r'$M$', r'$K_t$'],
               xlim=(distance[0], distance[-1]), ylabel=r'$E-E_D$ (meV)')
        ax.spines[['top', 'right']].set_visible(False)

    # 左图展示整体尺度，右图放大中心能带。放大坐标不会改变计算数据。
    axes[0].set(title='BM bands', ylim=(-120, 120))
    zoom = max(2.0, 1.4 * np.max(np.abs(center)))
    axes[1].set(title='Central two bands', ylim=(-zoom, zoom))
    axes[1].text(0.03, 0.97, f'Path bandwidths: {widths[0]:.3f}, {widths[1]:.3f} meV',
                 transform=axes[1].transAxes, va='top', fontsize=9)

    # 标题读取数据文件内的参数，防止改了脚本默认值却给旧数据标上新参数。
    p = metadata['parameters']
    fig.suptitle(rf"$\theta={p['theta_deg']:.3f}^\circ$, "
                 rf"$w_0={p['w0_mev']:g}$ meV, $w_1={p['w1_mev']:g}$ meV; "
                 rf"cutoff $={p['cutoff']}$, dimension $={dimension}$")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=PROJECT_ROOT / 'results' / 'tbg_bm_bands.npz')
    # 沿用项目当前已经存在的 figure/ 目录名称。
    parser.add_argument('--output', type=Path,
                        default=PROJECT_ROOT / 'figure' / 'tbg' / 'bm_bands.png')
    parser.add_argument('--no-show', action='store_true', help='仅保存图片，不弹出窗口')
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f'找不到数据文件 {args.input}；'
                     '请先运行 python/scripts/tbg/tbg_band_structure.py')
    if args.no_show:
        import matplotlib
        matplotlib.use('Agg')

    # np.load 只读取数组和 JSON 字符串；关闭 pickle 对象反序列化。
    with np.load(args.input, allow_pickle=False) as data:
        fig = plot_bands(
            data['distance_nm_inv'], data['ticks_nm_inv'], data['energies_mev'],
            json.loads(data['metadata_json'].item()),
        )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    print(f'图片已保存：{output}')
    if not args.no_show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == '__main__':
    main()
