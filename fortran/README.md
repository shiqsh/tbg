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