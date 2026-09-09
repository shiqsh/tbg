# Python 数值计算

`src` 放可复用模块，相当于 Fortran 的 `module`。`scripts` 放具体计算任务，
相当于 Fortran 的 `program`。绘图单独读取已保存的数据，修改图形时不重新计算能带。

```text
tbg/                            # 项目根目录（这里只列 BM 相关文件）
├── python/
│   ├── src/tbg/
│   │   ├── __init__.py          # 对外导出模型和计算函数
│   │   └── bm_model.py          # BM 参数、H(k)、层位置算符、k 路径、本征值
│   └── scripts/tbg/
│       ├── run_tbg_bands.py     # 总入口：依次计算、绘图
│       └── tbg_band_structure.py # 设置参数、计算能带、保存结果
├── plot/tbg/
│   └── tbg_bands.py             # 只读取结果、绘制能带
├── results/tbg_bm_bands.npz     # 计算输出：原始能量、k 路径、参数、单位
└── figure/tbg/bm_bands.png      # 绘图输出
```

## 计算与绘图

在项目根目录 `/Users/shiqsh/Workspace/tbg` 执行：

```bash
source python/.venv/bin/activate
./python/scripts/tbg/run_tbg_bands.py
```

第一行激活虚拟环境，同一个终端会话通常只需执行一次。第二行是一条总命令：
先计算并保存数据，再读取这份数据画图。计算出错时会停止，不继续运行绘图。

三个入口脚本都有 `#!/usr/bin/env python3`，且已设置可执行权限。
shebang 是脚本第一行，其中 `#!` 必须是半角字符；它让系统从当前 `PATH`
寻找 `python3`。激活 `.venv` 后会优先找到该环境的 Python，所以不必每次手写
`python/.venv/bin/python`。总入口通过 `sys.executable` 让两个子脚本使用同一个解释器。

已经激活环境后，也可以分别运行：

```bash
./python/scripts/tbg/tbg_band_structure.py
./plot/tbg/tbg_bands.py
```

这两行在终端中是依次执行的。希望前一步成功后才执行下一步，可以写成：

```bash
./python/scripts/tbg/tbg_band_structure.py && ./plot/tbg/tbg_bands.py
```

图依赖计算结果，不应把这两个脚本并行启动。
总入口加 `--no-show` 可仅保存图片；`--data`、`--figure` 可指定数据和图片输出路径。

模型参数集中在 `python/scripts/tbg/tbg_band_structure.py` 的 `PARAMETERS` 中，
采样密度由 `POINTS_PER_SEGMENT` 控制。默认转角为 1.09 度、
`w0=w1=110 meV`、`cutoff=4`，总矩阵维数为 324。

计算结果写入 `results/tbg_bm_bands.npz`；图片写入 `figure/tbg/bm_bands.png`。
如需比较不同参数，用计算脚本的 `--output` 保留多个数据文件，
再通过绘图脚本的 `--input` 选择数据。绘图脚本支持 `--no-show` 仅保存图片。
同一路径再次运行会更新该结果文件。两个脚本都可以从其他工作目录用绝对路径运行。

绘图脚本集中在根目录的 `plot/`，BM 使用 `plot/tbg/`；
已有的 `plot/orbital-magnetization/` 保持独立、不受影响。
如果只调整图形，运行 `./plot/tbg/tbg_bands.py` 即可，不必重新运行总入口。

## 在其他计算中使用模型

现有目录采用源码直接导入方式，不要求安装成全局 Python 包。
从 `python/scripts/tbg/` 下的新脚本调用时（若直接放在 `python/scripts/`，应使用 `parents[1]`）：

```python
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from tbg import BMModel, BMParameters

model = BMModel(BMParameters(theta_deg=1.09, w0_mev=110, w1_mev=110))
H = model.hamiltonian([0.0, 0.0])   # 输入 k：nm^-1；输出 H：meV
energies, vectors = np.linalg.eigh(H)
z = model.layer_z_diagonal()        # z 算符对角线：nm
```

导入 `tbg` 不会扫描动量、写文件或画图。只有运行计算脚本才会开始计算。
核心模块只依赖 NumPy，绘图脚本另外依赖 Matplotlib；现有 `.venv` 已具备这些依赖。

## 阅读代码时的约定

- `hbar_vf_mev_nm` 表示 `hbar*v_F`，不是单独的速度，不能再乘一次 `hbar`。
- `cutoff` 对应 Fortran 的 `tr`，总态数为 `4*(2*cutoff+1)^2`。
- 基底按下层/上层、`n1`、`n2`、A/B 排列；详细索引公式写在核心模块中。
- `hamiltonian` 输出厄米矩阵；`eigh` 返回的本征矢按列排列。
- 这个版本令 Dirac 旋转角与莫尔几何转角一致，相当于原 `HTBG` 中 `te=theta`。
- `solve_bands` 返回原始能量；绘图才统一减去一个 Dirac 点参考能量。
- 打印的带宽是高对称路径上的带宽；判断整个 BZ 的平坦程度还需二维采样。
- 模型保存预计算结果；修改参数时创建新的 `BMParameters` 和 `BMModel`。
