# fortran — TBG Fortran 数值实现

TBG 连续模型的 Fortran 数值实现，由 `fpm` 管理。

```bash
fpm build   # 编译项目
fpm test    # 运行环境自检（stdlib / OpenMP / LAPACK）
fpm run     # 运行默认 app
```

- `src/` — 可复用物理模块（规划：`tbg_constants` / `tbg_geometry` / `tbg_basis` / `tbg_hamiltonian` / `tbg_kpath` / `tbg_analytic`）
- `app/` — 可执行程序（`tbg/` 无应变、`heterostrain/` 含应变、`orbital-magnetization/` 轨道磁化 toy model）
- `test/` — 测试与环境自检

当前为早期搭建阶段，`src/` 与 `app/` 中多为占位代码。详见根目录 [README](../README.md)。
