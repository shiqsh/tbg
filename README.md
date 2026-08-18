# TBG

本项目用于研究扭转双层石墨烯（Twisted Bilayer Graphene, TBG）的连续模型及相关物理性质。

目前主要包括：

- TBG 连续模型与能带结构；
- heterostrain 下的能带、Dirac 点和 van Hove singularity；
- density of states；
- Fermi surface；
- orbital magnetization；
- Fortran 与 Python 两套数值实现及相互验证。

---

## 1. 项目目录结构

```text
tbg/
│
├── README.md
├── .gitignore
│
├── fortran/
│   ├── fpm.toml
│   │
│   ├── src/
│   │   ├── tbg_constants.f90
│   │   ├── tbg_geometry.f90
│   │   ├── tbg_basis.f90
│   │   ├── tbg_hamiltonian.f90
│   │   ├── tbg_kpath.f90
│   │   └── tbg_analytic.f90
│   │
│   ├── app/
│   │   ├── pristine/
│   │   │   └── band_structure.f90
│   │   │
│   │   └── heterostrain/
│   │       ├── band_structure.f90
│   │       ├── dos.f90
│   │       ├── dirac.f90
│   │       ├── vhs.f90
│   │       ├── fermi_surface.f90
│   │       └── orbital_magnetization.f90
│   │
│   └── test/
│
├── python/
│   ├── .venv/
│   │
│   ├── src/
│   │   └── tbg/
│   │       ├── __init__.py
│   │       ├── constants.py
│   │       ├── geometry.py
│   │       ├── basis.py
│   │       ├── hamiltonian.py
│   │       └── kpath.py
│   │
│   ├── scripts/
│   │   ├── band_structure.py
│   │   ├── dos.py
│   │   ├── dirac.py
│   │   ├── vhs.py
│   │   ├── fermi_surface.py
│   │   └── orbital_magnetization.py
│   │
│   ├── plot/
│   │   ├── plot_band.py
│   │   ├── plot_dos.py
│   │   ├── plot_dirac.py
│   │   ├── plot_vhs.py
│   │   ├── plot_fermi_surface.py
│   │   └── plot_orbital_magnetization.py
│   │
│   └── tests/
│
├── gnuplot/
│   ├── band_structure.gp
│   ├── dos.gp
│   └── fermi_surface.gp
│
├── data/
│   ├── input/
│   ├── benchmark/
│   └── reference/
│
├── results/
│
└── figures/
```

---

## 2. Fortran

`fortran/` 用于保存基于 Fortran 的主要数值实现，并使用 `fpm` 管理项目。

### 2.1 `fortran/src/`

`src/` 中只放可重复调用的 Fortran module 和底层数值程序，不放具体的计算任务。

例如：

```text
fortran/src/
├── tbg_constants.f90
├── tbg_geometry.f90
├── tbg_basis.f90
├── tbg_hamiltonian.f90
├── tbg_kpath.f90
└── tbg_analytic.f90
```

各文件主要职责如下：

- `tbg_constants.f90`
  - 基本物理常数；
  - Pauli matrix；
  - graphene lattice constant；
  - Fermi velocity 等固定参数。

- `tbg_geometry.f90`
  - graphene reciprocal lattice；
  - moiré reciprocal lattice；
  - `q1, q2, q3`；
  - twist angle；
  - heterostrain 相关几何量。

- `tbg_basis.f90`
  - plane-wave basis；
  - `(layer, n1, n2)` 与 Hamiltonian matrix index 之间的映射；
  - basis cutoff 管理。

- `tbg_hamiltonian.f90`
  - 单层 Dirac Hamiltonian；
  - interlayer tunneling；
  - BM continuum Hamiltonian；
  - heterostrain Hamiltonian。

- `tbg_kpath.f90`
  - 高对称点；
  - band-structure 路径；
  - reciprocal-space sampling。

- `tbg_analytic.f90`
  - 8 × 8、20 × 20 等解析截断模型；
  - analytic benchmark；
  - 与完整连续模型进行比较。

---

### 2.2 `fortran/app/`

`app/` 中放可以直接运行的 Fortran program。

这些程序调用 `src/` 中已经实现好的 module 来完成具体计算任务。

例如：

```text
fortran/app/
├── pristine/
│   └── band_structure.f90
│
└── heterostrain/
    ├── band_structure.f90
    ├── dos.f90
    ├── dirac.f90
    ├── vhs.f90
    ├── fermi_surface.f90
    └── orbital_magnetization.f90
```

其中：

- `pristine/`
  - 无 heterostrain 的 TBG 计算；

- `heterostrain/`
  - heterostrain 下的各类计算。

`app/` 中的程序主要负责：

1. 设置计算参数；
2. 调用 `src/` 中的物理和数值模块；
3. 扫描 momentum 或其他参数；
4. 将最终数据写入 `results/`。

具体物理模型和算法应尽可能放在 `src/` 中，而不是直接写在 `app/` 的 program 中。

---

### 2.3 `fortran/test/`

`test/` 用于保存测试程序，包括：

- Hamiltonian Hermiticity；
- reciprocal lattice consistency；
- basis index consistency；
- numerical convergence；
- analytic benchmark；
- LAPACK；
- OpenMP；
- stdlib；
- Apple Accelerate。

测试程序不用于正式产生科研数据。

---

### 2.4 `fpm.toml`

Fortran 项目由 `fpm.toml` 管理。

当前主要依赖包括：

- `stdlib`；
- OpenMP；
- BLAS/LAPACK；
- macOS Apple Accelerate。

常用命令：

```bash
fpm build
```

```bash
fpm run
```

```bash
fpm test
```

---

## 3. Python

`python/` 用于保存 Python 版本的 TBG 数值程序以及数据处理和绘图程序。

Python 项目主要分成：

```text
python/
├── src/
├── scripts/
├── plot/
└── tests/
```

---

### 3.1 `python/src/tbg/`

这里放可以重复调用的 Python module。

例如：

```text
python/src/tbg/
├── constants.py
├── geometry.py
├── basis.py
├── hamiltonian.py
└── kpath.py
```

这些文件与 Fortran 中的 `src/` 基本对应。

例如：

```text
Fortran                         Python

tbg_constants.f90        <->   constants.py

tbg_geometry.f90         <->   geometry.py

tbg_basis.f90            <->   basis.py

tbg_hamiltonian.f90      <->   hamiltonian.py

tbg_kpath.f90            <->   kpath.py
```

Python 与 Fortran 两套实现应尽可能使用相同的物理 convention，以方便进行数值比较。

---

### 3.2 `python/scripts/`

`scripts/` 用于保存实际执行数值计算的 Python 脚本。

例如：

```text
python/scripts/
├── band_structure.py
├── dos.py
├── dirac.py
├── vhs.py
├── fermi_surface.py
└── orbital_magnetization.py
```

这些脚本主要负责：

1. 设置参数；
2. import `src/tbg/` 中的功能；
3. 执行计算；
4. 将结果写入 `results/`。

`scripts/` 中的文件主要用于“计算”，而不是用于定义底层 Hamiltonian 或绘图。

---

### 3.3 `python/plot/`

`plot/` 只保存 Python 绘图脚本。

例如：

```text
python/plot/
├── plot_band.py
├── plot_dos.py
├── plot_dirac.py
├── plot_vhs.py
├── plot_fermi_surface.py
└── plot_orbital_magnetization.py
```

这些程序读取 `results/` 中已经计算完成的数据，并生成最终图片。

基本流程为：

```text
python/scripts/
        |
        v
    results/
        |
        v
 python/plot/
        |
        v
    figures/
```

绘图脚本原则上不重新执行耗时的 TBG Hamiltonian 计算。

这样可以做到：

- 数值计算一次；
- 数据长期保存；
- 绘图可以反复修改；
- 修改线宽、字体、颜色、坐标范围时不需要重新计算 Hamiltonian。

---

### 3.4 Python 虚拟环境

Python 虚拟环境位于：

```text
python/.venv/
```

激活方式：

```bash
cd python
source .venv/bin/activate
```

`.venv/` 只属于当前计算机的本地 Python 环境，不上传 GitHub。

---

## 4. Gnuplot

`gnuplot/` 单独保存 Gnuplot 绘图脚本。

例如：

```text
gnuplot/
├── band_structure.gp
├── dos.gp
└── fermi_surface.gp
```

Gnuplot 与 Python plotting script 使用同一套 `results/` 数据。

基本流程为：

```text
results/
   |
   +------> python/plot/
   |
   +------> gnuplot/
                |
                v
             figures/
```

运行示例：

```bash
gnuplot gnuplot/band_structure.gp
```

---

## 5. Data

`data/` 保存需要长期保留的数据，而不是普通程序运行产生的大量临时结果。

目录结构：

```text
data/
├── input/
├── benchmark/
└── reference/
```

---

### 5.1 `data/input/`

保存程序运行所需要的外部输入数据。

例如：

```text
data/input/
├── strain_parameters.dat
├── experimental_data.dat
└── external_parameters.dat
```

对于可以直接由公式生成的 TBG 参数，没有必要专门保存为输入文件。

---

### 5.2 `data/benchmark/`

保存已经验证正确的少量 benchmark 数据。

这些数据主要用于验证后续代码修改是否改变了原有正确结果。

例如：

```text
data/benchmark/
├── pristine_band_theta_1p4_tr4.dat
├── hamiltonian_test.dat
└── eigenvalue_reference.dat
```

典型用途是比较：

```text
Fortran result
      |
      v
benchmark data
      ^
      |
Python result
```

例如检查：

```text
max |E_python - E_fortran|
```

是否足够小。

Benchmark 文件应满足：

- 文件较小；
- 结果已经确认正确；
- 可以长期使用；
- 应当上传 GitHub。

---

### 5.3 `data/reference/`

保存外部参考数据。

例如：

```text
data/reference/
├── literature_band.dat
├── published_parameters.dat
└── experimental_reference.csv
```

这些数据主要用于：

- 与文献结果比较；
- 与实验结果比较；
- 保存论文中的公开数据；
- 保存其他代码得到的参考结果。

区别可以简单理解为：

```text
benchmark/
    用来验证“我的代码有没有算错”

reference/
    用来比较“我的结果和外部结果是否一致”
```

---

## 6. Results

`results/` 保存程序运行生成的数值结果。

例如：

```text
results/
├── pristine/
│   └── band/
│
└── heterostrain/
    ├── band/
    ├── dos/
    ├── dirac/
    ├── vhs/
    ├── fermi_surface/
    └── orbital_magnetization/
```

典型数据包括：

- band structure；
- density of states；
- Dirac point；
- van Hove singularity；
- Fermi surface；
- orbital magnetization；
- strain-angle scan；
- chemical-potential scan。

这些结果原则上应该能够通过源代码重新生成。

由于数据量可能很大，`results/` 默认不上传 GitHub。

---

## 7. Figures

`figures/` 只保存最终生成的图片。

例如：

```text
figures/
├── band/
├── dos/
├── dirac/
├── vhs/
├── fermi_surface/
└── orbital_magnetization/
```

图片格式可能包括：

```text
.pdf
.png
.svg
```

绘图代码本身不放在 `figures/` 中。

绘图程序分别位于：

```text
python/plot/
```

和：

```text
gnuplot/
```

基本关系为：

```text
results/
    |
    v
plotting scripts
    |
    v
figures/
```

`figures/` 中自动生成的大量图片默认不上传 GitHub。

---

## 8. 项目数据流

整个项目的基本工作流程为：

```text
                   source code
                       |
           +-----------+-----------+
           |                       |
       fortran/                 python/
           |                       |
         app/                  scripts/
           |                       |
           +-----------+-----------+
                       |
                       v
                    results/
                       |
           +-----------+-----------+
           |                       |
     python/plot/               gnuplot/
           |                       |
           +-----------+-----------+
                       |
                       v
                    figures/
```

`data/` 则用于向计算程序提供长期保存的输入、benchmark 和 reference 数据：

```text
data/
  |
  +---- input
  |
  +---- benchmark
  |
  +---- reference
  |
  v
Fortran / Python calculations
```

---

## 9. Git 管理原则

建议上传 GitHub 的内容包括：

- Fortran 源代码；
- Python 源代码；
- Python 绘图程序；
- Gnuplot 绘图程序；
- `fpm.toml`；
- README；
- `.gitignore`；
- 小型 benchmark 数据；
- 小型 reference 数据；
- 项目相关 VS Code 配置。

原则上不上传：

- Fortran 编译文件；
- fpm `build/`；
- Python `.venv/`；
- Python cache；
- 大量 numerical results；
- 自动生成的大量 figures；
- 临时文件。

---

## 10. 物理与数值规范

后续应逐渐在本 README 中记录本项目采用的主要 convention，包括：

- valley convention；
- top / bottom layer convention；
- sublattice basis；
- reciprocal lattice convention；
- `q1, q2, q3` convention；
- moiré Brillouin zone convention；
- heterostrain convention；
- plane-wave cutoff convention；
- energy unit；
- momentum unit；
- length unit。

保证 Fortran 与 Python 两套实现始终使用相同的物理定义。

---

## 11. 当前开发目标

当前阶段主要目标为：

