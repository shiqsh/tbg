#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ==========
# Pauli matrices
# ==========
sigma_0 = np.eye(2,dtype=complex)
sigma_1 = np.array([[0,1],[1,0]],dtype=complex)
sigma_2 = np.array([[0,-1j],[1j,0]],dtype=complex)
sigma_3 = np.array([[1,0],[0,-1]],dtype=complex)

# ==========
# QWZ Hamiltonian
# ==========
def H_QWZ(k, m):
    kx, ky = k
    H = ( np.sin(kx)*sigma_1 + np.sin(ky)*sigma_2 + (m + np.cos(kx) + np.cos(ky))*sigma_3 )
    return H

# ==========
# Bilayer QWZ Hamiltonian
# ==========
def H_BILAYER_QWZ(k, lambda_, m, t_perp):
    kx, ky = k
    H_top = H_QWZ((kx + lambda_, ky), m)
    H_bottom = H_QWZ((kx - lambda_, ky), m)


    H = np.block([
        [H_top, t_perp*sigma_0],
        [t_perp*sigma_0, H_bottom]
    ])
    return H


# ==========
# Z operator 
# Z = {H- mu, z}, z = d/2 tau_z otimes sigma_0
# ==========
def BUILD_Z(H, mu, d):

    I = np.eye(4, dtype=complex)
    tau_z = np.array([[1,0],[0,-1]],dtype=complex)

    z = (d/2.0) * np.kron(tau_z, sigma_0)

    H_mu = H - mu * I

    Z = H_mu @ z + z @ H_mu

    return Z


# ==========
# Occupied projector
# ========== 
def PROJECTOR(H, mu):
    E, V = np.linalg.eigh(H)

    occ = E < mu
    V_occ = V[:, occ]
    P = V_occ @ V_occ.conj().T

    nocc = np.sum(occ)
    
    return P,E, nocc

# ==========
# Inplane Orbital Magnetization
# ==========
def OM(lambda_, m=1.0, t_perp=0.2, mu=0.0, d=1.0, nk=101, prefactor=1.0): # 这里的prefactor 实际上是e/hbar c

    klist = np.linspace(-np.pi, np.pi, nk, endpoint=False) # Endpoint  = False 避免了取到两个相同的点，保证了k点的唯一性， 这里的取值为 -pi, -pi + 2*pi/nk, ..., pi + 2*(nk)*pi/nk

    dk = 2.0 * np.pi / nk

    P_grid = np.empty((nk, nk, 4, 4), dtype=complex)
    H_grid = np.empty((nk, nk, 4, 4), dtype=complex)

    E_grid = np.empty((nk, nk, 4), dtype=float)
    nocc_grid = np.empty((nk, nk), dtype=int)


    # ----------
    # Construct the projectors over the whole BZ 
    # ----------
    for i, kx in enumerate(klist):
        for j, ky in enumerate(klist):

            H = H_BILAYER_QWZ([kx, ky], lambda_, m, t_perp)

            P, E, nocc = PROJECTOR(H, mu)

            H_grid[i, j] = H
            P_grid[i, j] = P
            E_grid[i, j] = E
            nocc_grid[i, j] = nocc

    # ----------
    # Diagnose Insulator or Metal (just for diagnosis)
    # ----------
    nocc_min = np.min(nocc_grid)
    nocc_max = np.max(nocc_grid)

    is_insulator = (nocc_min == nocc_max)

    fermi_distance = np.min(np.abs(E_grid - mu))

    #---------
    # Finite difference of P
    #---------
    dP_dx = (np.roll(P_grid, -1, axis=0) - np.roll(P_grid, 1, axis=0)) / (2.0 * dk)
    dP_dy = (np.roll(P_grid, -1, axis=1) - np.roll(P_grid, 1, axis=1)) / (2.0 * dk)

    #-----
    # Empty projector
    #-----
    # I = np.eye(4, dtype=complex)
    # Q_grid = I[None, None, :, :] - P_grid

    # ---------
    # Integration over the BZ
    # ---------
    Mx = 0.0
    My = 0.0

    I = np.eye(4, dtype=complex)

    for i in range(nk):
        for j in range(nk):

            P = P_grid[i, j]
            Q = I - P
            Z = BUILD_Z(H_grid[i, j], mu, d)

            # -----
            # Fx = Re Tr [P(dP/dx) Q Z]
            # Fy = Re Tr [P(dP/dy) Q Z]
            # -----
            Fx = np.real(np.trace(P @ dP_dx[i, j] @ Q @ Z))
            Fy = np.real(np.trace(P @ dP_dy[i, j] @ Q @ Z))

            # -----
            # Mx = - Fx, My = Fy
            # += 就是累加，实际上也就是积分了
            # -----
            Mx += -prefactor * Fy 
            My += prefactor * Fx

            # -----------
            # Integration weight
            # -----------
            weight = dk**2 / (2.0 * np.pi)**2
            Mx *= weight
            My *= weight

            return (Mx, My, is_insulator, nocc_min, nocc_max, fermi_distance)


# ==========
# Main function
# ==========
if __name__ == "__main__":
    m = 1.0
    t_perp = 0.2
    mu = 0.0
    d = 1.0
    nk = 101
    prefactor = 1.0

    lambda_list = np.linspace(0.0, 0.2, 20, endpoint=False)

    Mx_list = np.empty_like(lambda_list)
    My_list = np.empty(len(lambda_list))

    for i, lambda_ in enumerate(lambda_list):

        (Mx, My, is_insulator, nocc_min, nocc_max, fermi_distance) = OM(lambda_, m=m, t_perp=t_perp, mu=mu, d=d, nk=nk, prefactor=prefactor)

        Mx_list[i] = Mx
        My_list[i] = My

        if is_insulator:
            system_type = "Insulator"
        else:
            system_type = "Metal"

        print(
            f"lambda = {lambda_:7.3f}   "
            f"Mx = {Mx: .16e}   "
            f"My = {My: .16e}   "
            f"{system_type:9s}   "
            f"Nocc = [{nocc_min}, {nocc_max}]   "
            f"min|E-mu| = {fermi_distance:.6e}"
        )

    # ---------
    # Plot Mx and My as a function of lambda
    # ---------
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(lambda_list, Mx_list, marker='o', label=r"$M_x$")
    ax.plot(lambda_list, My_list, marker='o', label=r"$M_y$")

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$M$")
    ax.set_title(rf"Bilayer QWZ Orbital Magnetization ($m={m}$, $t_\perp={t_perp}$, $\mu={mu}$, $d={d}$)" ) # 这里rf， r: raw string, f: format string

    ax.grid(True)
    ax.legend()

    fig.tight_layout() # 请注意这里是fig，而不是plt
    plt.show() # plt 仍然是 pyplot 模块









    
 
    

