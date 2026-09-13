#!/usr/bin/env python3

from pathlib import Path

import numpy as np 
import matplotlib.pyplot as plt

from qwz_om import OM, H_BILAYER_QWZ, sigma_0, sigma_1, sigma_2, sigma_3

ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = (
    ROOT
    / "data" 
    / "orbital-magnetization"
    / "qwz"
    / "python"
    / "susceptibility"
)

FIGURE_DIR = (
    ROOT 
    / "figure" 
    / "orbital-magnetization"
    / "qwz"
    / "python"
    / "susceptibility"
)

# make sure your file exist
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# ==========
# z-position operator
# ==========
def Z_QWZ(d):
    
    z = (d/2.0) * np.kron(sigma_3, sigma_0)

    return z

# ==========
# Derivatives of bilayer QWZ model Hamiltonian
# dH_dkx, dH_dky, dH_d_lambda
# =========
def H_DERIVATIVES(k, lambda_):

    kx, ky = k
    kxt = kx + lambda_
    kxb = kx - lambda_

    zero = np.zeros((2,2), dtype=complex)

    dH_dkx = np.block([
        [np.cos(kxt)*sigma_1-np.sin(kxt)*sigma_3, zero],
        [zero, np.cos(kxb)*sigma_1-np.sin(kxb)*sigma_3]
    ])

    dH_dky = np.block([
        [np.cos(ky)*sigma_2 - np.sin(ky)*sigma_3, zero],
        [zero, np.cos(ky)*sigma_2 - np.sin(ky)*sigma_3]
    ])

    dH_dlambda = np.block([
        [np.cos(kxt)*sigma_1 - np.sin(kxt)*sigma_3, zero],
        [zero, -np.cos(kxb)*sigma_1 + np.sin(kxb)*sigma_3]
    ])

    return dH_dkx, dH_dky, dH_dlambda

# ==========
# Gaussian representation of Dirac delta function
# delta_eta(E - mu) = exp[-((E -mu)/eta)^2]/(sqrt(pi)*eta)
# ==========
def DELTA_GAUSSIAN(E, mu, eta):

    if eta <= 0.0:
        raise ValueError("eta must be positive")

    delta = (
        np.exp(-((E - mu)/eta)**2)/(np.sqrt(np.pi)*eta)
    )

    return delta

# ==========
# FS contribution
# chi_FS(lambda, alpha)
# ==========
def CHI_FS(klist, dk, z, eta, lambda_, m, t_perp, mu, prefactor):

    sum_Gx = 0.0
    sum_Gy = 0.0

    # -----
    # BZ integration
    # -----
    for kx in klist:
        for ky in klist:

            k = [kx, ky]

            # -----Hamiltonian-----
            H = H_BILAYER_QWZ(k, lambda_, m, t_perp)
            # -----eigenvalue/vectors-----
            E, V = np.linalg.eigh(H)
            Vdag = V.conj().T   # 共轭转置，dagger, 这里是N times N 的矩阵做共轭转置

            # -----
            # dH_dkx, dH_dky, dH_dlambda
            # -----
            dH_dkx, dH_dky, dH_dlambda = H_DERIVATIVES(k, lambda_)

            # -----
            # 将算符转换为表示矩阵
            # -----
            Hx_e = Vdag @ dH_dkx @ V
            Hy_e = Vdag @ dH_dky @ V
            Hlambda_e = Vdag @ dH_dlambda @ V
            z_e = Vdag @ z @ V

            Gx = 0.0
            Gy = 0.0

            # -----
            # sum over m
            # -----
            for mm in range(len(E)):

                delta_E = DELTA_GAUSSIAN(E[mm], mu, eta)

                # ----- check 求导矩阵元的对角项是否为实数 -----
                Hx_mm = Hx_e[mm, mm]
                Hy_mm = Hy_e[mm, mm]
                Hlambda_mm = Hlambda_e[mm, mm]

                tol = 1.0e-12

                if abs(np.imag(Hx_mm)) > tol:
                    raise ValueError("<m|Hx|m> is not real")
                if abs(np.imag(Hy_mm)) > tol:
                    raise ValueError("<m|Hy|m> is not real")
                if abs(np.imag(Hlambda_mm)) > tol:
                    raise ValueError("<m|Hlambda|m> is not real")

                # ----- sum over n != mm -----
                A_x = 0.0
                A_y = 0.0
                A_lambda = 0.0

                for n in range(len(E)):

                    if n==mm:
                        continue

                    # ---- sum_{n!=m} Re[braket{m|Hlambda|n}braket{n|z|m}] -----
                    # ----- lambda -----
                    A_lambda += np.real(
                        Hlambda_e[mm, n] * z_e[n, mm]
                    )
                    # ----- if beta = x -----
                    A_x += np.real(
                        Hx_e[mm, n] * z_e[n, mm]
                    )
                    # ----- if beta = y -----
                    A_y += np.real(
                        Hy_e[mm, n] * z_e[n, mm]
                    )

                # ----- 回到对mm的求和 -----
                Gx += delta_E * np.real(
                    Hx_mm * A_lambda - Hlambda_mm * A_x
                )

                Gy += delta_E * np.real(
                    Hy_mm * A_lambda - Hlambda_mm * A_y
                )

            sum_Gx += Gx
            sum_Gy += Gy

    # -----
    # Integration Weight
    # -----
    weight = dk**2/(2.0*np.pi)**2

    # -----
    # chi_FS = + prefactor epsilon_{alpha beta} G_beta
    # -----
    chi_fs_x = (
        prefactor * weight * sum_Gy
    )

    chi_fs_y = (
        -prefactor * weight * sum_Gx
    )

    return chi_fs_x, chi_fs_y


# ==========
# FIX contribution 
# CHI_FIXED(lambda, alpha)
# Attention: directly use braket{m|{Hlambda, z}|n}
# ==========
def CHI_FIXED(klist, dk, z, lambda_, m, t_perp, mu, prefactor):

    sum_Fx = 0.0
    sum_Fy = 0.0

    for kx in klist:
        for ky in klist:

            k = [kx, ky]

            H = H_BILAYER_QWZ(k, lambda_, m, t_perp)

            E, V = np.linalg.eigh(H)
            Vdag = V.conj().T

            # -----
            # 构造微分表示矩阵
            # ------
            dH_dx, dH_dy, dH_dlambda = H_DERIVATIVES(k, lambda_)

            Zlambda = z @ dH_dlambda + dH_dlambda @ z
            Zx = z @ dH_dx + dH_dx @ z
            Zy = z @ dH_dy + dH_dy @ z

            Hx_e = Vdag @ dH_dx @ V
            Hy_e = Vdag @ dH_dy @ V
            Hlambda_e = Vdag @ dH_dlambda @ V

            Zlambda_e = Vdag @ Zlambda @ V
            Zx_e = Vdag @ Zx @ V
            Zy_e = Vdag @ Zy @ V

            # -----
            # T = 0
            # f_n(1-f_m) != 0 only when
            # n : occupied, m : unoccupied
            # -----
            occ = np.where(E < mu)[0]
            emp = np.where(E >= mu)[0]

            Fx = 0.0
            Fy = 0.0

            # ----------
            # sum over n, m
            # ----------
            for n in occ:
                for mm in emp:

                    # ----- 分母 -----
                    denominator = E[n] - E[mm]

                    # ----- beta = x -----
                    term_x = (
                        Hx_e[n, mm] * Zlambda_e[mm, n] - Hlambda_e[n, mm] * Zx_e[mm, n]
                    )

                    # ----- beta = y -----
                    term_y = (
                        Hy_e[n, mm] * Zlambda_e[mm, n] - Hlambda_e[n, mm] * Zy_e[mm, n]
                    )

                    Fx += np.real(
                        term_x / denominator
                    )

                    Fy += np.real(
                        term_y / denominator
                    )

            sum_Fx += Fx
            sum_Fy += Fy

    # -----
    # 积分权重
    # -----
    weight = dk**2/(2.0*np.pi)**2

    # ----------
    # chi_FIX = - prefactor epsilon_{alpha beta} F_beta
    # ----------
    chi_fix_x = (
        -prefactor * weight * sum_Fy
    )

    chi_fix_y = (
        + prefactor * weight * sum_Fx
    )

    return chi_fix_x, chi_fix_y


# ==========
# 主程序
# ==========
if __name__=="__main__":

    # ----------
    # Common parameters
    # ----------
    common_param = {
        "m": 1.0,
        "t_perp": 0.2,
        "mu": 1.0,
        "d": 1.0,
        "nk": 251,
        "prefactor": 1.0
    }

    eta = 0.0258

    klist = np.linspace(-np.pi, np.pi, common_param["nk"], endpoint=False)
    dk = 2.0*np.pi/common_param["nk"]
    z = Z_QWZ(common_param["d"])

    # ----------
    # Parameters actually passed to
    # both CHI_FIXED and CHI_FS
    # ----------
    chi_param = {
        "m": common_param["m"],
        "t_perp": common_param["t_perp"],
        "mu": common_param["mu"],
        "prefactor": common_param["prefactor"]
    }

    # ----- lambda scan -----
    lambda_list = np.linspace(0.0, 0.2, 21, endpoint=True)

    # -----
    # result arrays
    # -----
    chi_fix_x_list = np.empty_like(lambda_list)
    chi_fix_y_list = np.empty_like(lambda_list)
    chi_fs_x_list = np.empty_like(lambda_list)
    chi_fs_y_list = np.empty_like(lambda_list)
    chi_tot_x_list = np.empty_like(lambda_list)
    chi_tot_y_list = np.empty_like(lambda_list)

    # ----- Orbital magnetization ----- 
    Mx_list = np.empty_like(lambda_list)
    My_list = np.empty_like(lambda_list)

    # ----------
    # 计算Susceptibility
    # ----------
    for i, lambda_ in enumerate(lambda_list):
        # ----- FIX contribution -----
        chi_fix_x, chi_fix_y = CHI_FIXED(
            klist, dk, z, lambda_, **chi_param
        )
        # ----- FS contribution -----
        chi_fs_x, chi_fs_y = CHI_FS(
            klist, dk, z, eta, lambda_, **chi_param
        )

        # ----- total -----
        chi_tot_x = chi_fix_x + chi_fs_x
        chi_tot_y = chi_fix_y + chi_fs_y

        # ----- store result -----
        chi_fix_x_list[i] = chi_fix_x
        chi_fix_y_list[i] = chi_fix_y
        chi_fs_x_list[i] = chi_fs_x
        chi_fs_y_list[i] = chi_fs_y
        chi_tot_x_list[i] = chi_tot_x
        chi_tot_y_list[i] = chi_tot_y

        # ----- 计算轨道磁化 ----- 
        Mx, My, is_insulator, nocc_min, nocc_max, fermi_distance = OM(lambda_, **common_param)

        Mx_list[i] = Mx
        My_list[i] = My

        # ----------
        # print
        # ----------
        print(
            f"lambda = {lambda_:7.3f}   "
            f"FIX = ({chi_fix_x: .8e}, {chi_fix_y: .8e})   "
            f"FS = ({chi_fs_x: .8e}, {chi_fs_y: .8e})   "
            f"TOT = ({chi_tot_x: .8e}, {chi_tot_y: .8e})"
        )

    # ==========
    # Reconstruct M from susceptibility
    # ==========
    My_from_chi = np.empty_like((2, len(lambda_list)))

    M_from_chi = np.empty((2, len(lambda_list)))

    # 初值
    M_from_chi[0, 0] = Mx_list[0] # Mx
    M_from_chi[1, 0] = My_list[0] # My

    for i in range(1, len(lambda_list)):

        delta_lambda = lambda_list[i] - lambda_list[i-1]

        chi_mid_x = 0.5 * (chi_tot_x_list[i] + chi_tot_x_list[i-1])
        chi_mid_y = 0.5 * (chi_tot_y_list[i] + chi_tot_y_list[i-1])

        M_from_chi[0, i] = (M_from_chi[0, i-1] + chi_mid_x * delta_lambda)
        M_from_chi[1, i] = (M_from_chi[1, i-1] + chi_mid_y * delta_lambda)


    # ----------
    # 误差分析，针对nk和eta
    # ----------
    M_error_y = My_list - M_from_chi[1, :]
    max_abs_error_y = np.max(np.abs(M_error_y))

    rms_error_y = np.sqrt(np.mean(M_error_y**2))

    print(
        f"max |My - My_from_chi| = {max_abs_error_y:.8e}"
    )
    print(
        f"RMS error = {rms_error_y:.8e}"
    )
    # ----- Dimensionless quantity -----    
    M_scale = np.max(My_list) - np.min(My_list)

    normalized_error_y = max_abs_error_y / M_scale

    print(
        f"normalized max error = {normalized_error_y:.8e}"
    )


    # ==========
    # Output data，输出数据也是主程序的一部分
    # ==========

    name  = (
        f"qwz_chi"
        f"_m_{common_param['m']:.3f}"
        f"_tperp_{common_param['t_perp']:.3f}"
        f"_mu_{common_param['mu']:.3f}"
        f"_eta_{eta:.3f}"
        f"_nk_{common_param['nk']}"
    )

    RUN_FIGURE_DIR = FIGURE_DIR / name
    RUN_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    data_file = DATA_DIR / f"{name}.dat"

    with open(data_file, "w") as of: # output file

        # ----- 输出文件中的参数注释 -----
        of.write(f"# m = {common_param['m']:.16e}\n")
        of.write(f"# t_perp = {common_param['t_perp']:.16e}\n")
        of.write(f"# mu = {common_param['mu']:.16e}\n")
        of.write(f"# d = {common_param['d']:.16e}\n")
        of.write(f"# nk = {common_param['nk']}\n")
        of.write(f"# eta = {eta:.16e}\n")
        of.write(f"# prefactor = {common_param['prefactor']:.16e}\n")
        
        # ----- column information(表头)-----
        of.write(
            "# columns: lambda "
            "Mx My " 
            "Mx_from_chi My_from_chi " # 加入空格 
            "chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y\n"
            )

        for i in range(len(lambda_list)):

            of.write(
                f"{lambda_list[i]:32.16e}"
                f"{Mx_list[i]:32.16e}"
                f"{My_list[i]:32.16e}"
                f"{M_from_chi[0, i]:32.16e}"
                f"{M_from_chi[1, i]:32.16e}"
                f"{chi_fix_x_list[i]:32.16e}"
                f"{chi_fix_y_list[i]:32.16e}"
                f"{chi_fs_x_list[i]:32.16e}"
                f"{chi_fs_y_list[i]:32.16e}"
                f"{chi_tot_x_list[i]:32.16e}"
                f"{chi_tot_y_list[i]:32.16e}\n"
            )


    # ==========
    # 绘图
    # ==========
    # ----------
    # plot chi_x
    # ----------
    fig, ax = plt.subplots(figsize = (8,8))

    ax.plot(lambda_list, chi_fix_x_list, marker="o", label = r"$\chi^{\mathrm{FIX}}_{\lambda, x}$")
    ax.plot(lambda_list, chi_fs_x_list, marker="o", label = r"$\chi^{\mathrm{FS}}_{\lambda, x}$")
    ax.plot(lambda_list, chi_tot_x_list, marker="o", label = r"$\chi^{\mathrm{TOT}}_{\lambda, x}$")

    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\chi_{\lambda,x}$")

    ax.set_title(
        rf"Bilayer QWZ Susceptibility "
        rf"($m={common_param['m']}$, "
        rf"$t_\perp={common_param['t_perp']}$, "
        rf"$\mu={common_param['mu']}$, "
        rf"$\eta={eta}$, "
        rf"$N_k={common_param['nk']}$)"
        )

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
    RUN_FIGURE_DIR / f"{name}_x.pdf",
        dpi=300,
        bbox_inches="tight"
        )


    # ----------
    # plot chi_y
    # ----------
    fig, ax = plt.subplots(figsize=(8, 8))

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

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\chi_{\lambda,y}$")

    ax.set_title(
        rf"Bilayer QWZ Susceptibility "
        rf"($m={common_param['m']}$, "
        rf"$t_\perp={common_param['t_perp']}$, "
        rf"$\mu={common_param['mu']}$, "
        rf"$\eta={eta}$, "
        rf"$N_k={common_param['nk']}$)"
        )

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RUN_FIGURE_DIR / f"{name}_y.pdf",
        dpi=300,
        bbox_inches="tight"
        )

    # -----
    # Bechmark
    # M vs sum over chi_{i} delta{lambda_i}
    # -----
    fig, ax = plt.subplots(figsize=(8,8))

    ax.plot(lambda_list, My_list, marker="o", label=r"$M_y $")
    ax.plot(lambda_list, M_from_chi[1,:], marker="o", label=r"$M_y(0)+\int_0^\lambda\chi_{\lambda',y}^{\mathrm{TOT}}d\lambda'$")

    ax.ticklabel_format(axis="y", style="sci",scilimits=(0,0), useMathText=True)

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$M_y$")

    ax.set_title(
        rf"Bilayer QWZ: $M_y$ vs. integrated susceptibility "
        rf"($m={common_param['m']}$, "
        rf"$t_\perp={common_param['t_perp']}$, "
        rf"$\mu={common_param['mu']}$, "
        rf"$\eta={eta}$, "
        rf"$N_k={common_param['nk']}$)"
    )

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RUN_FIGURE_DIR / f"{name}_benchmark_my.pdf",
        dpi = 300,
        bbox_inches="tight"
        )

    # ----- show all figures ----- 
    plt.show()

                


    

    

        
        







