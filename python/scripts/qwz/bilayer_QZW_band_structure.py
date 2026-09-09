#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


# =========================
# parameters
# =========================

m = 1.0
lam = 0.025
tp = 0.2


# =========================
# Pauli matrices
# =========================

sigma0 = np.eye(2,dtype=complex)
sigma1 = np.array([[0,1],[1,0]],dtype=complex)
sigma2 = np.array([[0,-1j],[1j,0]],dtype=complex)
sigma3 = np.array([[1,0],[0,-1]],dtype=complex)


# =========================
# QWZ Hamiltonian
# =========================

def qwz(kx,ky):
    dx = np.sin(kx)
    dy = np.sin(ky)
    dz = m + np.cos(kx)+np.cos(ky)

    return dx*sigma1 + dy*sigma2 + dz*sigma3


# =========================
# bilayer QWZ
# =========================

def bilayer_qwz(kx,ky):
    ht = qwz(kx+lam,ky)
    hb = qwz(kx-lam,ky)

    H = np.zeros((4,4),dtype=complex)


    H[:2,:2]=ht
    H[2:,2:]=hb

    H[:2,2:]=tp*sigma0
    H[2:,:2]=tp*sigma0

    return H


# =========================
# k path
# Gamma -> X -> M -> Gamma
# =========================

# ----------
# high symmetric point
# ----------
Gamma = np.array([0.0, 0.0])
X = np.array([np.pi, 0.0])
M = np.array([np.pi, np.pi])

N = 100

def segment(start, end, n, endpoint=False):
    t = np.linspace(0.0, 1.0, n, endpoint=endpoint)
    return start[None, :] + t[:, None] * (end - start)[None, :]

k_path = np.vstack([
    segment(Gamma, X, N, endpoint=False),
    segment(X, M, N, endpoint=False),
    segment(M, Gamma, N, endpoint=True),
])

bands = []

for kx, ky in k_path:
    H = bilayer_qwz(kx, ky)
    E = np.linalg.eigvalsh(H)
    bands.append(E)

bands = np.array(bands)

# np.diff 相邻元素后面一个减去前面一个
delta_k = np.diff(k_path, axis=0) # axis=0 ：沿着“行的方向”往下走，也就是ky - kx
# np.linalg.norm 对每一行中的两个元素 \(x,y\) 求 norm
dk = np.linalg.norm(delta_k, axis=1) # axis = 1: ：沿着“列的方向”横着走, 实际上求的是｜delta ki｜
# np.cumsum : cumlative sum 累积求和：每一步都保存了下来， 
k_distance = np.concatenate(([0.0], np.cumsum(dk)))

fig, ax = plt.subplots(figsize=(6, 6))

for n in range(4):
    ax.plot(k_distance, bands[:, n])

tick_indices = [0, N, 2*N, 3*N - 1]

ax.set_xticks(
    k_distance[tick_indices],
    [r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$"],
)

ax.set_ylabel("Energy")
ax.set_title("Bilayer QWZ dispersion")
ax.grid(True)

fig.tight_layout() # 进行排版
plt.show()