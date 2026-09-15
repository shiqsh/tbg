#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt 

# =====
# 输出路径
# =====
ROOT = Path(__file__).resolve().parents[2]

FIGURE_DIR = (ROOT 
              / "figure" 
              / "qwz" 
              / "fortran"
              / "om")
FIGURE_DIR.mkdir(parents=True, exist_ok=True) # 如果目录不存在，则创建目录

# =====
# 数据文件
# =====

name = f'qwz_om_m_1.000_tperp_0.200_mu_-1.000'

data_file = (ROOT 
             / "data" 
             / "qwz" 
             / "fortran"
             / "om" 
             /f"{name}.dat")

# =====
# 参数
# =====
param = {}

with open(data_file, "r") as f:
    for line in f:

        if not line.startswith("#"):
            break

        if "=" in line:
            key, value = line[1:].split("=", 1)
            param[key.strip()] = float(value) # 

# 假设当前读到的行是：line = "# m = 1.0\n"
# line[1:]  → 去掉第一个字符 '#'，得到 " m = 1.0\n"
# .split("=", 1)  → 以 '=' 为分隔符，最多分割一次，
#                   返回列表 [" m ", " 1.0\n"]
# key, value = ... → 将列表解包：
#                    key = " m "
#                    value = " 1.0\n"
# key.strip()  → 去除 key 两端的空白，得到 "m"
# float(value) → 将 value 转换为浮点数（自动忽略空白和换行），得到 1.0


m      = param["m"]
t_perp = param["t_perp"]
mu     = param["mu"]
d      = param["d"]
delta  = param["delta"]
nk     = int(param["nk"])


# =====
# 读取数据
# =====
data = np.loadtxt(data_file)

lambda_list = data[:, 0]
Mx = data[:, 1]
My = data[:, 2]

# =====
# 绘图
# =====
fig, ax = plt.subplots(figsize=(8,8))

ax.plot(lambda_list, Mx, marker="o", label=r"$M_x$")
ax.plot(lambda_list, My, marker='o', label=r"$M_y$")

ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$M$")

ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True) # 也可以单独是x, y

ax.set_title(
    rf"Bilayer QWZ Orbital Magnetization "
    rf"($m={m:.1f}$, $t_\perp={t_perp:.1f}$, "
    rf"$\mu={mu:.1f}$, $\delta={delta:.2e}$, $N_k={nk}$)"
)

ax.grid(True)
ax.legend()

fig.tight_layout()

# -----
# 输出
# -----
figure_name = (
    f"qwz_om_fortran"
    f"_m_{m:.3f}"
    f"_tperp_{t_perp:.3f}"
    f"_mu_{mu:.3f}"
)

fig.savefig(
    FIGURE_DIR / f"{figure_name}.pdf",
    dpi=300,
    bbox_inches="tight"
    )

plt.show()