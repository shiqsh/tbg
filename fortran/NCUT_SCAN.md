# 固定应变点的 ncut 收敛扫描

新主程序：`app/orbital-magnetization/tbg_strain/tbg_strain_chi_ncut.f90`。
原来的主程序和 module 不需要改动。

## 默认参数

- φ = 20°，ε = 0.05%（无量纲 `5.d-4`）。
- θ = 1.4°，μ = 0.5 meV，η = 0.25 meV，nk = 40。
- ncut = 20, 40, 60, 80, 100, 120, 140, 160。
- 当前 tr=4，N=324，允许的最大 ncut=N/2−1=161；160 合法。

上述参数集中在新主程序的参数区。所有 ncut 使用相同的 nk、物理参数与积分路径。

## 为什么内部仍有 epsilon 路径？

输出比较的是同一个固定物理点，但保留原程序的 reconstructed OM 需要：

`M_reconstructed(epsilon_fixed; ncut) = M_reference(0; ncut) + integral_0^epsilon_fixed chi(epsilon; ncut) d epsilon`。

因此 ncut 是外层循环；每个 ncut 内部计算 ε=0, 0.01%, …, 0.05% 六个路径点。
这与附件原程序 `linspace(0,0.002,21)` 的前六个节点对应。
每个路径点都会更新几何和 Hamiltonian；参考 OM 也按当前 ncut 重新计算。
不能把 χ 对 ncut 积分，χ 不是 dM/dncut。

这里 χ 仍是原程序 epsilon 分支对应的 dM/dε，不是 dM/dφ。
仅改变 ncut 不代表已经检验 nk、eta 或 epsilon 积分步长的收敛。

为忠实保留附件行为：direct OM 用 `CALCULATE_OM_CUTOFF_KET`，reference OM 用
`M_FROM_CHI_CUTOFF` 内的原 projector。没有悄悄切换为全 ket 参考，也不强制 M(0)=0。

## 运行

在 `/Users/shiqsh/Workspace/tbg/fortran` 目录，先只检查配置：

```sh
fpm run --target tbg_strain_chi_ncut -- --check-only
```

正式运行（外层 OpenMP 并行，限制 BLAS 线程，计算期间防止闲置休眠）：

```sh
OMP_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 caffeinate -i fpm run --target tbg_strain_chi_ncut
```

默认结果位于 `data/tbg/fortran/susceptibility/`，文件名包括物理参数、nk、ncut 范围、路径点数。
文件已存在时程序会拒绝覆盖；重跑可追加 `-- --output /完整路径/新文件名.dat`。

辅助选项 `--nk N --path-points N --ncut-last N` 可用于小规模试跑或改变网格，
不指定时按参数区运行。`--ncut-last` 必须在从20起、步长20的序列上。
例如 `--nk 2 --path-points 2 --ncut-last 40` **只适合检查流程，不是可信的物理精度**。

## 输出

主 `.dat` 文件：每完成一个 ncut 立即写入并 flush 一行，共15列：

`ncut epsilon_percent phi_deg Mx_OM My_OM Mx_from_chi My_from_chi chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y delta_Mx delta_My`

其中 `delta_M = M_from_chi - M_OM`。磁矩单位为 μB/cell；χ 为每单位无量纲应变的响应。

另一个 `.dat.path.dat` 文件逐点保存积分路径上的 FIX、FS、TOT，便于检查或事后重新积分。
此简化主程序不复用旧 `CHI_CHECKPOINT` 的 epsilon/phi 扫描格式，**不提供自动断点续算**；
中断后已完成的 ncut 行和 χ 路径行仍保留，但重新运行不会自动跳过这些点。
