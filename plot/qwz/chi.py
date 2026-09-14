#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==========
# parameters
# ==========
m = 1.0
t_perp = 0.2
mu = -1.0
eta = 0.03
nk = 251

# ==========
# parameter strings
# ==========
s_m = f"{m:.3f}"
s_tperp = f"{t_perp:.3f}"
s_mu = f"{mu:.3f}"
s_eta = f"{eta:.3f}"
s_nk = f"{nk:d}"

# ==========
# project root
# ==========
ROOT = Path(__file__).resolve().parents[2]

# ==========
# Data paths
# ==========
CHI_DATA_DIR = (
    ROOT
    / "data"
    / "qwz"
    / "fortran"
    / "susceptibility"
)

OM_DATA_DIR = (
    ROOT
    / "data"
    / "qwz"
    / "fortran"
    / "om"
)

# ==========
# Figure root
# ==========
FIGURE_DIR = (
    ROOT
    / "figure"
    / "qwz"
    / "fortran"
    / "susceptibility"
)

# ==========
# Parameter-dependent figure directory
# ==========
run_name = (
    f"m_{s_m}"
    f"_tperp_{s_tperp}"
    f"_mu_{s_mu}"
    f"_eta_{s_eta}"
    f"_nk_{s_nk}"
)

RUN_FIGURE_DIR = FIGURE_DIR / run_name

RUN_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# ==========
# Input files
# ==========
chi_filename = (
    f"qwz_chi"
    f"_m_{s_m}"
    f"_tperp_{s_tperp}"
    f"_mu_{s_mu}"
    f"_eta_{s_eta}"
    f"_nk_{s_nk}"
    f".dat"
)

om_filename = (
    f"qwz_om"
    f"_m_{s_m}"
    f"_tperp_{s_tperp}"
    f"_mu_{s_mu}"
    f".dat"
)

# ==========
# Full input path
# ==========
chi_path = CHI_DATA_DIR / chi_filename
om_path = OM_DATA_DIR / om_filename

# ==========
# Check
# ==========
print("ROOT =", ROOT)
print("chi_path =", chi_path)
print("om_path =", om_path)
print("figure path =", RUN_FIGURE_DIR)

# ==========
# 读取 susceptibility 数据
# ==========
chi_data = np.loadtxt(chi_path)

# ----- check ----- 
if chi_data.ndim != 2:
    raise ValueError(
        "Susceptibility data must be a two-dimensional array"
    )

if chi_data.shape[1] != 9:
    raise ValueError(
        f"Expected 9 columns in susceptibility file, "
        f"but got {chi_data.shape[1]}"
    )

# ----- input data ----- 
lambda_list = chi_data[:, 0]
Mx_from_chi = chi_data[:, 1]
My_from_chi = chi_data[:, 2]
chi_fix_x_list = chi_data[:, 3]
chi_fix_y_list = chi_data[:, 4]
chi_fs_x_list = chi_data[:, 5]
chi_fs_y_list = chi_data[:, 6]
chi_tot_x_list = chi_data[:, 7]
chi_tot_y_list = chi_data[:, 8]

# ==========
# 读取OM数据
# ==========
om_data = np.loadtxt(om_path)

if om_data.ndim != 2:
    raise ValueError(
        "OM data must be a two-dimensional array"
    )

if om_data.shape[1] != 3:
    raise ValueError(
        f"Expected 3 columns in OM file, "
        f"but got {om_data.shape[1]}"
    )

lambda_om = om_data[:, 0]
Mx_list = om_data[:, 1]
My_list = om_data[:, 2]

# ==========
# 检查lambda 网格
# ==========
# 检查两个 lambda 数组的形状是否一致（例如长度、维度）
# lambda_list 和 lambda_om 分别来自磁化率计算和轨道磁化计算
if lambda_list.shape != lambda_om.shape:
    # 形状不同则无法逐点比较，抛出异常
    raise ValueError(
        "Susceptibility and OM lambda arrays "
        "have different dimensions"
    )

# 形状相同的情况下，再检查两个数组的数值是否在容差范围内完全一致
# np.allclose 会逐元素比较，默认相对容差 1e-05，绝对容差 1e-08
if not np.allclose(lambda_list, lambda_om):
    # 数值不一致说明两次计算使用了不同的 lambda 网格，抛出异常
    raise ValueError(
        "Susceptibility and OM calculations "
        "use different lambda meshes"
    )

# ==========
# 公共标题参数
# ==========
parameter_title = (
    rf"$m={m}$, "
    rf"$t_\perp={t_perp}$, "
    rf"$\mu={mu}$, "
    rf"$\eta={eta}$, "
    rf"$N_k={nk}$"
)


# ====================
# 绘图
# ====================

# ----------
# chi_x
# ----------
fig, ax = plt.subplots(
    figsize=(8, 8)
)

ax.plot(
    lambda_list,
    chi_fix_x_list,
    marker="o",
    label=r"$\chi^{\mathrm{FIX}}_{\lambda,x}$"
)

ax.plot(
    lambda_list,
    chi_fs_x_list,
    marker="o",
    label=r"$\chi^{\mathrm{FS}}_{\lambda,x}$"
)

ax.plot(
    lambda_list,
    chi_tot_x_list,
    marker="o",
    label=r"$\chi^{\mathrm{TOT}}_{\lambda,x}$"
)

ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0),
    useMathText=True
)

ax.set_xlabel(
    r"$\lambda$"
)

ax.set_ylabel(
    r"$\chi_{\lambda,x}$"
)

ax.set_title(
    "Bilayer QWZ Susceptibility\n"
    + parameter_title
)

ax.grid(True)

ax.legend()

fig.tight_layout()

fig.savefig(
    RUN_FIGURE_DIR / "chi_x.pdf",
    dpi=300,
    bbox_inches="tight"
)

# ----------
# chi_y
# ----------
fig, ax = plt.subplots(
    figsize=(8, 8)
)

ax.plot(
    lambda_list,
    chi_fix_y_list,
    marker="o",
    label=r"$\chi^{\mathrm{FIX}}_{\lambda,y}$"
)

ax.plot(
    lambda_list,
    chi_fs_y_list,
    marker="o",
    label=r"$\chi^{\mathrm{FS}}_{\lambda,y}$"
)

ax.plot(
    lambda_list,
    chi_tot_y_list,
    marker="o",
    label=r"$\chi^{\mathrm{TOT}}_{\lambda,y}$"
)

ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0),
    useMathText=True
)

ax.set_xlabel(
    r"$\lambda$"
)
ax.set_ylabel(
    r"$\chi_{\lambda,y}$"
)

ax.set_title(
    "Bilayer QWZ Susceptibility\n"
    + parameter_title
)

ax.grid(True)

ax.legend()

fig.tight_layout()

fig.savefig(
    RUN_FIGURE_DIR / "chi_y.pdf",
    dpi=300,
    bbox_inches="tight"
)

# ----------
# Benchmark Mx
# direct Mx vs. Mx(lambda_0) + integral chi_x dlambda
# ----------
fig, ax = plt.subplots(
    figsize=(8, 8)
)

ax.plot(
    lambda_om,
    Mx_list,
    marker="o",
    label=r"$M_x^{\mathrm{direct}}$"
)

ax.plot(
    lambda_list,
    Mx_from_chi,
    marker="o",
    label=(
        r"$M_x(\lambda_0)"
        r"+\int_{\lambda_0}^{\lambda}"
        r"\chi_{\lambda',x}^{\mathrm{TOT}}"
        r"\,d\lambda'$"
    )
)

ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0),
    useMathText=True
)

ax.set_xlabel(
    r"$\lambda$"
)

ax.set_ylabel(
    r"$M_x$"
)

ax.set_title(
    r"Bilayer QWZ: $M_x$ vs. integrated susceptibility"
    "\n"
    + parameter_title
)

ax.grid(True)

ax.legend()

fig.tight_layout()

fig.savefig(
    RUN_FIGURE_DIR / "benchmark_Mx.pdf",
    dpi=300,
    bbox_inches="tight"
)

# ----------
# Benchmark My
# direct My vs. My(lambda_0) + integral chi_y dlambda
# ----------
fig, ax = plt.subplots(
    figsize=(8, 8)
)

ax.plot(
    lambda_om,
    My_list,
    marker="o",
    label=r"$M_y^{\mathrm{direct}}$"
)

ax.plot(
    lambda_list,
    My_from_chi,
    marker="o",
    label=(
        r"$M_y(\lambda_0)"
        r"+\int_{\lambda_0}^{\lambda}"
        r"\chi_{\lambda',y}^{\mathrm{TOT}}"
        r"\,d\lambda'$"
    )
)

ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0),
    useMathText=True
)

ax.set_xlabel(
    r"$\lambda$"
)

ax.set_ylabel(
    r"$M_y$"
)

ax.set_title(
    r"Bilayer QWZ: $M_y$ vs. integrated susceptibility"
    "\n"
    + parameter_title
)

ax.grid(True)

ax.legend()

fig.tight_layout()

fig.savefig(
    RUN_FIGURE_DIR / "benchmark_My.pdf",
    dpi=300,
    bbox_inches="tight"
)


# ====================
# error check
# ====================
# ----------
# Mx error
# ----------
M_error_x = Mx_list - Mx_from_chi

max_abs_error_x = np.max(
    np.abs(M_error_x)
)

rms_error_x = np.sqrt(
    np.mean(M_error_x**2)
)


print()

print(
    f"max |Mx - Mx_from_chi| = "
    f"{max_abs_error_x:.8e}"
)

print(
    f"RMS error x = "
    f"{rms_error_x:.8e}"
)

# ----- Mx的无量纲误差 ----- 
M_scale_x = (
    np.max(Mx_list)
    - np.min(Mx_list)
)


if M_scale_x > 0.0:

    normalized_error_x = (
        max_abs_error_x
        / M_scale_x
    )
else:

    normalized_error_x = np.nan

print(
    f"normalized max error x = "
    f"{normalized_error_x:.8e}"
)


# ----------
# My error
# ----------
M_error_y = My_list - My_from_chi

max_abs_error_y = np.max(
    np.abs(M_error_y)
)

rms_error_y = np.sqrt(
    np.mean(M_error_y**2)
)


print()

print(
    f"max |My - My_from_chi| = "
    f"{max_abs_error_y:.8e}"
)

print(
    f"RMS error y = "
    f"{rms_error_y:.8e}"
)

# ----- My 的无量纲误差 ----- 
M_scale_y = (
    np.max(My_list)
    - np.min(My_list)
)


if M_scale_y > 0.0:

    normalized_error_y = (
        max_abs_error_y
        / M_scale_y
    )

else:

    normalized_error_y = np.nan

print(
    f"normalized max error y = "
    f"{normalized_error_y:.8e}"
)

# ====================
# show all figures 
# ====================
plt.show()






