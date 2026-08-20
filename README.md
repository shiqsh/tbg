# TBG — 扭转双层石墨烯连续模型数值研究

Numerical study of the continuum model for twisted bilayer graphene (TBG), including heterostrain effects and orbital magnetization.

## 1. 项目简介

本项目研究扭转双层石墨烯（Twisted Bilayer Graphene, TBG）的连续模型（continuum model，即 Bistritzer–MacDonald 模型）及其相关物理性质，采用 **Fortran + Python 双实现、相互验证** 的方式组织代码：

- **Fortran**（`fortran/`）：主要数值实现，由 `fpm` 管理，追求性能；
- **Python**（`python/`）：快速原型、交叉验证与绘图；
- 两套实现遵循**相同的物理约定**，结果写入共享的 `results/`，便于逐项比对。

计划研究的物理量：

- 能带结构（band structure）
- 态密度（density of states, DOS）
- Dirac 点（Dirac point）
- van Hove 奇异性（van Hove singularity, VHS）
- 费米面（Fermi surface）
- 轨道磁化（orbital magnetization）
- heterostrain 对上述物理量的影响

## 2. 当前状态

> 项目处于**早期搭建阶段**：Fortran 编译环境（stdlib / OpenMP / LAPACK）已通过自检，
> 物理计算模块**尚未实现**，`app/` 下多为占位文件。

| 部分 | 状态 | 说明 |
| --- | --- | --- |
| `fortran/` | ✅ 环境就绪 | `fpm build` / `fpm test` 可运行，`test/check.f90` 自检通过 |
| `fortran/src/` | 🚧 占位 | 仅有示例模块 `inplane_orbital_magnetization.f90`（Hello, world），物理模块待实现 |
| `fortran/app/` | 🚧 占位 | 程序骨架已建好，部分文件为空（0 字节），待填充 |
| `python/` | 🚧 环境就绪 | 仅 `requirements.txt`，`src/`、`scripts/`、`tests/` 均为空目录 |
| `data/` `results/` `figures/` `plot/` | ⏳ 空 | 目录已预留，尚未使用 |

## 3. 快速开始

### 3.1 Fortran

依赖：`gfortran`、`fpm`，以及 stdlib / OpenMP / BLAS / LAPACK（macOS 上可链接 Apple Accelerate）。

```bash
cd fortran
fpm build   # 编译项目（含 src/ 与 app/）
fpm test    # 运行环境自检：stdlib 精度类型、OpenMP 线程、LAPACK ZHEEV
fpm run     # 运行默认 app
```

### 3.2 Python

```bash
cd python
python3 -m venv .venv                 # 若尚未创建
source .venv/bin/activate
pip install -r requirements.txt
```

依赖（`requirements.txt`）：`numpy`、`scipy`、`matplotlib` 等。

## 4. 目录结构（现状）

```text
tbg/
├── README.md
├── .gitignore
├── .vscode/
│   └── settings.json          # Fortran 语言服务 / lint 配置
│
├── data/                      # (空) 长期保存的输入 / benchmark / reference 数据
├── results/                   # (空) 程序生成的数值结果
├── figures/                   # (空) 最终图片
├── plot/                      # (空) 绘图脚本目录（规划放置 Gnuplot 脚本）
│
├── fortran/                   # Fortran 数值实现（fpm 项目）
│   ├── fpm.toml               # 项目配置：stdlib / openmp / blas
│   ├── src/
│   │   └── inplane_orbital_magnetization.f90   # 占位模块（示例，待替换）
│   ├── app/
│   │   ├── tbg/                               # 无 heterostrain 的计算
│   │   │   └── band_structure.f90             # 占位程序
│   │   ├── heterostrain/                      # heterostrain 下的计算
│   │   │   ├── dirac.f90                      # (空) Dirac 点
│   │   │   ├── dos.f90                        # (空) 态密度
│   │   │   ├── fermi_surface.f90              # (空) 费米面
│   │   │   └── velocity.f90                   # (空) 速度 / 群速度
│   │   └── orbital-magnetization/
│   │       └── test_toy_model.f90             # (空) 轨道磁化 toy model
│   ├── test/
│   │   └── check.f90          # 环境自检：stdlib / OpenMP / LAPACK
│   └── build/                 # fpm 编译产物（不纳入版本控制）
│
└── python/                    # Python 数值实现
    ├── requirements.txt
    ├── .venv/                 # 虚拟环境（不纳入版本控制）
    ├── src/                   # (空) 可复用模块
    ├── scripts/               # (空) 计算脚本
    └── tests/                 # (空) 测试
```

## 5. 目标结构（规划）

以下为**规划中的设计**，模块落地后按此组织；`(空)` 目录即为此规划预留。

### 5.1 Fortran 模块划分（`fortran/src/`）

`src/` 只放可重复调用的 module，不放具体计算任务：

| 文件 | 职责 |
| --- | --- |
| `tbg_constants.f90` | 物理常数、Pauli 矩阵、晶格常数、Fermi velocity |
| `tbg_geometry.f90` | 石墨烯倒格矢、moiré 倒格矢、`q1 q2 q3`、转角、heterostrain 几何量 |
| `tbg_basis.f90` | plane-wave basis、`(layer, n1, n2)` ↔ 矩阵指标映射、cutoff 管理 |
| `tbg_hamiltonian.f90` | 单层 Dirac Hamiltonian、层间 tunneling、BM 连续模型、heterostrain 项 |
| `tbg_kpath.f90` | 高对称点、能带路径、倒空间采样 |
| `tbg_analytic.f90` | 8×8 / 20×20 等解析截断模型，作为 benchmark |

`fortran/app/` 中的 program 只负责：设置参数 → 调用 `src/` 模块 → 扫描动量等参数 → 写结果到 `results/`。物理模型与算法尽量放在 `src/`，不直接写进 program。

### 5.2 Python 镜像（`python/`）

```text
python/
├── src/tbg/
│   ├── constants.py    <->  tbg_constants.f90
│   ├── geometry.py     <->  tbg_geometry.f90
│   ├── basis.py        <->  tbg_basis.f90
│   ├── hamiltonian.py  <->  tbg_hamiltonian.f90
│   └── kpath.py        <->  tbg_kpath.f90
├── scripts/            # 计算脚本：band_structure / dos / dirac / vhs /
│                       #           fermi_surface / orbital_magnetization
├── plot/               # 绘图脚本：读取 results/，不重新计算 Hamiltonian
└── tests/
```

### 5.3 数据流

```text
fortran/app/  ─┐
               ├─→  results/  ──┬─→  python/plot/  ─┐
python/scripts/─┘               └─→  plot/ (gnuplot) ─┴─→  figures/
```

原则：数值计算一次、数据长期保存（`results/`）、绘图可反复修改（改线宽/颜色/坐标范围时**不需要**重算 Hamiltonian）。

## 6. 数据管理

| 目录 | 用途 | 是否入库 |
| --- | --- | --- |
| `data/input/` | 程序运行所需的外部输入（如 strain 参数、实验数据） | 视大小 |
| `data/benchmark/` | 已验证正确的少量基准数据，用于回归测试（如 `max\|E_python − E_fortran\|`） | ✅ 小文件入库 |
| `data/reference/` | 文献 / 实验 / 其他代码的外部参考数据 | 视大小 |
| `results/` | 程序生成的数值结果（可能很大） | ❌ 默认不入库 |
| `figures/` | 最终图片（.pdf / .png / .svg） | ❌ 默认不入库 |

区分：`benchmark/` 验证"我的代码有没有算错"，`reference/` 比较"我的结果与外部结果是否一致"。

## 7. Git 管理原则

入库：Fortran / Python 源代码、绘图脚本、`fpm.toml`、`requirements.txt`、README、`.gitignore`、VS Code 配置、小型 benchmark/reference 数据。

不入库：`fortran/build/`、`python/.venv/`、`__pycache__/`、大量 `results/`、自动生成的 `figures/`、临时文件（已由 `.gitignore` 覆盖）。

## 8. 物理与数值约定（待完善）

两套实现必须使用相同的物理定义，后续在本节逐步记录：

- valley / 上下层（top–bottom layer）convention
- sublattice basis
- 倒格矢与 moiré Brillouin zone convention
- `q1, q2, q3` convention
- heterostrain convention
- plane-wave cutoff convention
- 能量 / 动量 / 长度单位

## 9. 开发路线图

1. **实现 Fortran 连续模型核心模块**：constants → geometry → basis → hamiltonian → kpath；
2. **生成标准 benchmark**：用 Fortran 计算 pristine TBG 能带，存入 `data/benchmark/`；
3. **实现 Python 对应版本**，逐项比对能带等结果（`max|E_python − E_fortran|`）；
4. **heterostrain 系列计算**：能带 / Dirac 点 / VHS / DOS / Fermi surface；
5. **轨道磁化**（orbital magnetization）计算；
6. **绘图与论文级 figure**（`python/plot/` 与 `plot/`）。

---

## License

MIT（见 `fortran/fpm.toml`，Copyright 2026, Qian-Sheng Shi）。
