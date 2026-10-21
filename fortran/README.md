# fortran — TBG Fortran 数值实现

TBG 连续模型的 Fortran 数值实现，由 `fpm` 管理。

## 常用命令

```bash
fpm build                       # 只编译项目，不运行
fpm test                        # 运行所有测试
fpm test --target check_linalg  # 运行指定测试
fpm run                         # 编译有改动的代码并运行默认 app
fpm run --target <name>         # 运行指定 app
fpm run --list                  # 查看所有可运行的 app
```

## 项目运行命令
```bash
cd /Users/shiqsh/Workspace/tbg/fortran
fpm run ...
```
需要进入到fortran所在的文件夹进行运行


## 日常运行流程

修改 `src/` 或 `app/` 中的代码后，通常直接运行：

```bash
fpm run
```

`fpm run` 会自动检查哪些源文件发生了变化，只重新编译需要更新的部分，然后运行程序。

因此通常不需要先执行：

```bash
fpm build
```

`fpm build` 主要用于只想检查项目能否成功编译、但暂时不运行程序的情况。

简单来说：

```text
修改代码后想运行      -> fpm run
只想检查能不能编译    -> fpm build
运行测试              -> fpm test
```

## build 目录

`build/` 中保存编译结果、`.mod`、`.o`、可执行文件，以及 `stdlib` 等依赖。

正常开发时不要删除 `build/`。

如果代码发生变化，直接：

```bash
fpm run
```

或：

```bash
fpm build
```

fpm 会自动进行增量编译。

如果执行：

```bash
rm -rf build
```

会把已经下载和编译的 `stdlib` 等依赖一起删除，下一次运行 fpm 时会重新下载并重新编译。

因此只有在编译缓存或依赖出现异常、需要彻底重新构建时，才使用：

```bash
rm -rf build
```

## 运行时间与 CPU 使用

普通运行：

```bash
fpm run
```

测量程序运行时间：

```bash
fpm run --runner "/usr/bin/time -l"
```

运行指定程序并测量时间：

```bash
fpm run --target <name> --runner "/usr/bin/time -l"
```

这对应于以前直接运行：

```bash
time ./executable
```

主要关注：

```text
real    实际运行时间
user    CPU 总计算时间
sys     系统调用时间
```

对于 OpenMP 并行程序，如果 `user` 明显大于 `real`，通常说明多个 CPU 核心同时参与了计算。

实时查看 CPU 占用：

```bash
top -o cpu
```

也可以使用 macOS 的 Activity Monitor（活动监视器）。

## 目录结构

- `src/` — 可复用物理模块（规划：`tbg_constants` / `tbg_geometry` / `tbg_basis` / `tbg_hamiltonian` / `tbg_kpath` / `tbg_analytic`）
- `app/` — 可执行程序（`tbg/` 无应变、`heterostrain/` 含应变、`orbital-magnetization/` 轨道磁化 toy model）
- `test/` — 测试与环境自检

当前为早期搭建阶段，`src/` 与 `app/` 中多为占位代码。详见根目录 [README](../README.md)。

## Heterostrain shifted Dirac-point energy

`src/heterostrain/dirac.f90` 提供 `FIND_DIRAC_POINTS`，`app/heterostrain/dirac_energy.f90` 负责两类参数扫描。应变 ε 输入扫描范围为 **0–0.2%**，传入物理模块前除以 100，即无量纲 **0–0.002**。默认 θ=1.4°；ε 扫描固定 φ=0°，φ 扫描固定 ε=0.1%。

在 `fortran/` 目录运行：

```bash
fpm run --target dirac_energy --profile release
```

默认同时生成两份数据，ε 包含 21 个点（步长 0.01%），φ 包含 91 个点（步长 2°），都包含起点和终点。默认输出目录为 `../data/tbg/fortran/dirac/`，程序自动创建目录。文件名记录扫描类别、θ、固定参数和点数（参数标签保留到小数点后六位）；相同配置重新运行会覆盖该配置的数据文件。默认文件名为 `epsilon_0-0p2pct_phi_0deg_theta_1p4deg_n21.dat` 和 `phi_0-180deg_epsilon_0p1pct_theta_1p4deg_n91.dat`。

位置参数依次为：

```text
dirac_energy [both|epsilon|phi] [phi_fixed_deg] [epsilon_fixed_percent] [theta_deg] [output_dir] [n_epsilon] [n_phi]
```

例如，仅计算 φ=30° 时的 ε 扫描，取 41 点；或仅计算 ε=0.15% 时的 φ 扫描，取 181 点：

```bash
fpm run --target dirac_energy --profile release -- epsilon 30 0.1 1.4 ../data/tbg/fortran/dirac 41 91
fpm run --target dirac_energy --profile release -- phi 0 0.15 1.4 ../data/tbg/fortran/dirac 21 181
fpm run --target dirac_energy -- --help
```

位置参数允许从末尾省略。只改后面的参数时，需要填入前面的默认值。角度参数单位为度，固定 ε 参数单位为百分比，因此 `0.1` 表示 0.1%，不是无量纲 0.1。

每个采样点输出两行，列顺序为：

| 列 | 含义 |
|---|---|
| `epsilon_percent`, `phi_deg` | 应变百分数与应变方向角（度） |
| `branch` | 1：从 top/origin 锥连续追踪；2：从 bottom/q1 锥连续追踪 |
| `kx/l`, `ky/l` | 动量除以 `MOIRE_BZ_SCALE(theta)`，与旧程序的归一化一致 |
| `E_D_meV` | 中央导带与价带能量的平均值，单位 meV |
| `gap_meV` | 同一动量的导带与价带能量差，单位 meV |
| `converged`, `valid` | 1 为真、0 为假；优化收敛及物理有效性标记 |
| `evaluations` | 本次锥搜索的哈密顿量/本征值评估次数 |

首个采样点使用模块的两个默认起点。之后，仅上一点标为 `valid` 的结果作为对应锥的连续追踪初值；无效锥重新从当前几何的名义位置搜索。两类扫描各自重新初始化。失败结果也保留在文件中，应结合 `valid`、`converged` 和残余 `gap_meV` 筛选数据。输出直接给出 `E_D`，不计算 `E_D/epsilon`，ε=0 不会发生除零。

文件头记录 θ、固定扫描参数、`tr`、哈密顿量维数、`w0/w1`、`hbar*vF`、应变参数和归一化动量尺度；模型参数继续使用 `constants.f90` 和 `strain.f90` 的定义。搜索成功与基底截断收敛属于不同问题；用于物理结论前仍需检查 `tr` 的收敛性。

模块只缓存同一组物理参数的 `H(0)/eV`，每次更新动量相关的 2×2 对角块；用 [LAPACK HEEVR](https://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f_source.html) 按索引求第 `NTBG`、`NTBG+1` 两个本征值，不求特征向量。局部二维 Nelder–Mead 搜索复用求解工作区，能量、能隙和坐标同步记录。失败或重复吸引域会尝试额外初值，再按与上一点的距离配对。这里搜索两条中心带的局部交点，不保证枚举全 BZ 内所有交点；branch 是连续轨迹标签，不代表始终属于某一层。

默认 `gap_tolerance=1e-8*eV`（J），`k_tolerance=1e-8`（以 `theta*|K_graphene|` 归一化），每个局部搜索最多 200 次迭代；达到动量收敛且残余能隙低于阈值才接受。两个交点间距小于 `max(10*k_tolerance,1e-8)*theta*|K_graphene|` 时按一个交点计数。默认值下实际去重尺度为 `1e-7*theta*|K_graphene|`。`n_found` 是有效且不同的点数；返回数组恒为两个诊断槽，不能用 `size(candidates)` 当有效点数。

模块入口（θ、φ 用弧度，ε 用无量纲应变；返回 k 用 m⁻¹，energy/gap 用 J）：

```fortran
call FIND_DIRAC_POINTS(theta, epsilon, phi, candidates, n_found)
! 可传 seeds(2,2)、gap_tolerance、k_tolerance、max_iterations。
call FIND_DIRAC_POINT(theta, epsilon, phi, candidate)
! 单点默认从原点搜索；可传 seed(2) 及上述收敛参数。
call CALCULATE_GAP(kx, ky, theta, epsilon, phi, gap, energy)
! energy 可省略，同一次求解计算 gap 与能量均值。
```

独立验证：

```bash
fpm test --target check_dirac
```

验证使用现有完整本征谱求解器交叉核对任意动量和优化交点的能量/能隙，并检查零应变、180°周期性、失败标记。

## Heterostrain DOS

`src/heterostrain/dos.f90` 保留原程序的线性三角形积分：在应变后的两个倒格矢张成的 `[0,1]²` 平行四边形上划分网格，每个小格分成两个三角形，并用三个顶点的本征值作线性插值。计算包含全部 `2*NTBG=324` 条能带，默认 `100×100` 个小格、每个方向 `101` 个顶点。DOS 归一化为单自旋、单谷的 `states/(eV·moire cell)`，一条完整能带积分为 1；没有乘自旋和谷的简并因子 4。三角形贡献本身已经有能量倒数的单位，因此删除旧程序结尾额外乘 `1/dE` 的操作。

输出采用 **bin-averaged analytic triangle DOS**：仍用同一线性三角形的 tent 分布，通过解析累计态数（CDF）的差计算每个能量区间的态数，再除以区间宽度。区间边界取相邻能量点的中点，首尾区间宽度为普通区间的一半，因此用输出数据作梯形积分正好对应区间态数之和。这样避免窄于能量步长的三角形在直接点采样时漏态，并正确处理三顶点等能或两端等能的情况。任意能量步长都能守恒计算能量窗口内的态数；完整能带被窗口包含时积分为 1。没有高斯展宽，也不对生成的曲线额外平滑。

默认参数为 `theta=1.4°`、`epsilon=0.1%`、`phi=90°`，能量范围 `[-50,50] meV`，2000 个区间、2001 个采样点。程序内部用弧度、无量纲应变和焦耳；`.dat` 两列分别为能量 eV 和 DOS `states/(eV·moire cell)`。数据头保存实际模型参数、网格、能量范围和归一化，绘图标题从这些数据头读取。

从项目根目录 `/Users/shiqsh/Workspace/tbg` 开始，先计算，再绘图：

```bash
cd fortran
OMP_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=1 fpm run --target strain_dos --profile release
cd ..
python/.venv/bin/python plot/tbg/strain_dos.py --no-show
```

默认数据为 `data/tbg/fortran/dos/dos_theta1p4_epsilon0p1_phi90_nk100_ne2000_eminm50_emax50.dat`；图片为 `figure/tbg/` 下同名 PNG 和 PDF。绘图显示完整计算能量范围和数值纵轴。省略 `--no-show` 可以打开图形窗口。绘图脚本按自身位置查找项目根目录，可以从任意工作目录运行。

计算主程序的可选位置参数为：

```text
strain_dos [theta_deg] [epsilon_percent] [phi_deg] [nk] [n_intervals] [E_min_meV] [E_max_meV] [output_dir]
```

其中 `nk` 表示每个倒格矢方向的小格数。`epsilon_percent=0.1` 表示 `0.1%`，内部应变为 `0.001`。例如，在 `fortran/` 中调整网格而保留默认物理参数：

```bash
fpm run --target strain_dos --profile release -- 1.4 0.1 90 150 2000 -50 50
```

输出文件名包含实际角度、应变、网格、能量区间数和能量端点。相同参数再次运行会替换对应数据文件。自定义数据的绘图命令，从项目根目录运行：

```bash
python/.venv/bin/python plot/tbg/strain_dos.py \
  --input data/tbg/fortran/dos/dos_theta1p4_epsilon0p1_phi90_nk150_ne2000_eminm50_emax50.dat \
  --output figure/tbg/strain_dos_nk150.png --no-show
```

`CAL_DOS_TRIANGLE` 的公开接口仍接收 `(theta, epsilon, phi, E_min, E_max, N_intervals, energy, density)`，角度为弧度、能量为 J、返回密度为 `1/J`；可选 `nkx`、`nky` 控制两个方向的小格数。默认值保持旧算法的 `100×100` 网格。

`OMP_NUM_THREADS=4` 并行处理动量点，`VECLIB_MAXIMUM_THREADS=1` 限制 macOS Accelerate 每次矩阵求解的线程数，避免并行层相互叠加。只求本征值，并复用缓存的 `H(0)` 和每线程 LAPACK 工作区。没有重写应变哈密顿量的物理公式。

回归验证：

```bash
fpm test --target check_dos
```

包括已知线性色散的 DOS、每条能带态数归一、能量分辨率、平带与重复顶点、跨越窗口边界的极窄三角形，以及完整本征谱对优化动量采样的交叉核对。
