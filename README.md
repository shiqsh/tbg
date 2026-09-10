# TBG — 扭转双层石墨烯数值计算

研究扭转双层石墨烯的 Bistritzer–MacDonald（BM）连续模型及其能带、应变效应和轨道磁化，并使用双层 QWZ 等简化模型检验计算方法。

**Fortran 负责正式数值计算；Python 用于 benchmark，并与 Fortran 的结果交叉验证。**
绘图与计算分开，绘图脚本统一放在根目录 `plot/`，便于读取已有结果、独立调整图形。

## 目录

```text
tbg/
├── fortran/              # 正式计算，使用 fpm 管理
│   ├── src/              # 可复用的物理模型与数值方法（module）
│   ├── app/              # 具体计算任务与参数扫描（program）
│   └── test/             # 数值方法测试与环境自检
├── python/               # benchmark 与交叉验证
│   ├── src/              # 对应模型的 Python 实现
│   └── scripts/          # 按模型组织的 benchmark 计算入口
├── plot/                 # 绘图脚本：tbg/、orbital-magnetization/ 等
├── scripts/              # 串联计算、绘图的 Shell 脚本
├── data/                 # 按课题保存的数据；现有 QWZ 计算结果在这里
├── results/              # 程序输出；现有 Python BM 能带数据在这里
└── figure/               # 输出图片
```

目前各算例的输出路径以对应程序为准，不假定所有数据都写入 `results/`。

## 运行

以下每组命令都从项目根目录执行。Fortran 需要 `gfortran` 和 `fpm`；
依赖配置见 [fpm.toml](fortran/fpm.toml) 和 [Python requirements](python/requirements.txt)。

### Fortran：正式计算

```bash
cd fortran
fpm run --list                 # 查看可运行的程序
fpm run --target qwz_band      # 示例：双层 QWZ 能带
fpm test                      # 运行测试
```

具体任务在 `fortran/app/` 中选择；物理模型和公共算法放在 `fortran/src/`。

### Python：BM benchmark 示例

```bash
source python/.venv/bin/activate
./python/scripts/tbg/run_tbg_bands.py
```

该命令依次计算 Python 能带、保存数据、绘图，不会自动运行 Fortran 或完成两者的比较。
默认数据为 `results/tbg_bm_bands.npz`，图片为 `figure/tbg/bm_bands.png`。

已有数据、只想重新绘图时，激活同一环境后运行：

```bash
./plot/tbg/tbg_bands.py
```

参数设置、分步运行和自定义输出路径见 [Python 使用说明](python/README.md)。

## 交叉验证约定

- 比较前对齐模型参数、基底与层/谷约定、单位、动量采样、截断和能量零点；不要直接使用两套程序各自的默认值。
- 对哈密顿量、本征值或目标物理量做数值误差比较，不仅凭图形判断；逐项比较矩阵前还需统一基底顺序与相位约定。
- 分别保留 Fortran 与 Python 的结果及参数记录。两者一致后，仍需检查截断、采样等数值收敛性，再用于正式分析。

许可证：MIT。
