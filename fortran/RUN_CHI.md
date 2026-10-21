# TBG chi：安全启动、逐点保存与续算

## 2026-10-20 的恢复

旧进程卡在第 18 点的 CHI 计算。已从其内存提取并分别读取两次，确认前 17 点的
`Mx_OM, My_OM, chi_fix_x/y, chi_fs_x/y, chi_tot_x/y` 逐字节一致；
每点的 `chi_tot = chi_fix + chi_fs` 也通过验证。

备份目录：`../data/tbg/fortran/susceptibility/recovery-20261020/`。

- `recovered-points-01-17.dat`：10 列可读结果，文件头列出列名；不是原来的最终 14 列格式。
- `recovered.checkpoint`：用于续算的精确参数、完整扫描网格和 17 行已完成结果。
- `snapshot.json`、`snapshot-verify.json`：两次读取的原始数组字节及内存快照。
- `validation.json`：校验结果与文件 SHA-256。
- `original-tbg_strain_chi`、`original-main.f90`：恢复时备份的可执行文件和主程序源文件。

第 18 点未完成，不能保留为一个已完成点。`Mx_from_chi/My_from_chi` 当时尚未计算，
没有将它们的未初始化内存当作结果；续算结束后仍由原来的 `M_FROM_CHI_CUTOFF` 重构。
内存快照是数值结果备份，不是可以还原全部进程状态的 core dump。

## 日常启动

从任意目录运行，输出目录必须存在，输出文件名必须是新的：

```sh
zsh /Users/shiqsh/Workspace/tbg/scripts/run_tbg_chi_safe.sh \
  --output /Users/shiqsh/Workspace/tbg/data/tbg/fortran/susceptibility/my-new-run.dat
```

脚本只为这次运行设置线程环境，不改系统设置或 `.zshrc`：

- 默认 10 个外层 OpenMP 线程，可以用 `TBG_THREADS=4 zsh ...` 调整。
- `OMP_MAX_ACTIVE_LEVELS=1`：不再嵌套创建 OpenMP 线程组。
- `VECLIB_MAXIMUM_THREADS=1`：限制 Apple vecLib 内部线程，避免与外层并行叠加。
- `caffeinate -i`：计算期间防止空闲导致的系统睡眠；屏幕仍可关闭。
- `/usr/bin/time -l`：结束时报告时间和内存。

不要在计算时主动选择系统“睡眠”或关机。防空闲睡眠并不是防止一切停机。
原问题高度疑似 OpenMP/Accelerate 同步卡住；短测试通过不等于证明永远不会再发生。

## 保存与续算

每个扫描点的 CHI 和 OM 都完成后，立即向 `输出文件.dat.checkpoint` 写入全部八个
结果并 `flush`，然后才打印该点完成。`flush` 让数据离开 Fortran 输出缓冲区，
不等同于断电情况下的磁盘持久化保证。

最终 `.dat` 仍使用原来的 14 列格式；在全扫描及参考 OM 重构完成前，它可能只有文件头。
判断进度应查看日志和 `.checkpoint`，不要仅看最终 `.dat` 大小。

续算时，输入原 checkpoint，输出使用另一个新文件名：

```sh
zsh /Users/shiqsh/Workspace/tbg/scripts/run_tbg_chi_safe.sh \
  --resume /absolute/path/previous.dat.checkpoint \
  --output /absolute/path/continued.dat
```

续算会验证维度、索引范围、theta/d/mu/eta/prefactor/delta 和两条参数网格，并复制已有点到
新的 checkpoint，再计算缺失的点。**必须保持哈密顿量、物理常数及核函数代码一致**；
checkpoint 不包含源代码指纹，无法自动检测模型公式的修改。

参数不匹配、数据行不完整或不连续时会报错，不会静默跳过。原文件从不覆盖。
若最后一行因外部终止而不完整，先保留原件，再检查副本，不要盲目删行。

诊断选项：

- `--resume FILE --check-resume`：只读取验证，既不计算也不创建文件。
- `--stop-after K`：最多再完成 K 个新点后正常退出并保留 checkpoint。

## 本次续算

最初用恢复备份 `recovered.checkpoint` 从 18/21 开始，已计算并保存到 20/21。
随后切换到显式 `KeepAlive=false`、`LaunchOnlyOnce=true` 的后台配置，从 21/21 继续，
避免 `launchctl submit` 的默认保活重启行为。前一阶段的日志和 checkpoint 全部保留。
新输出目录是 `../data/tbg/fortran/susceptibility/runs/20261020-resumed-18to21/`：

- `progress.log`：**当前**标准输出和逐阶段进度。
- `diagnostics.log`：**当前**诊断、时间统计和错误输出。
- `resumed-once.dat.checkpoint`：当前每完成一点就追加并 flush 的数据。
- `resumed-once.dat`：全扫描完成后的最终 14 列结果。
- `saved-points-01-20.checkpoint`：切换时的 20 点独立备份。
- `run-continued.log`、`resumed.dat.checkpoint`：前一阶段记录，保存到第 20 点。
- `resumed.dat`：前一阶段的最终输出占位文件，只有文件头，不是最终结果。

这是通过 macOS `launchctl` 启动的一次性后台计算，标签为
`local.tbg.chi.resume-once-20261020`；不是定时任务，不会自动重复计算。
原终端不再显示新程序的进度。可另外开终端查看：

```sh
tail -f /Users/shiqsh/Workspace/tbg/data/tbg/fortran/susceptibility/runs/20261020-resumed-18to21/progress.log
```

这里按 Ctrl+C 只退出 `tail`，不会停止计算。不要重复运行同一续算命令。
任务完成后，可用 `launchctl remove local.tbg.chi.resume-once-20261020` 清理已退出的任务登记（如果还在）；
**运行中执行 remove 会停止它**，请先确认最终结果已保存。
