"""Bistritzer–MacDonald 模型：供其他计算程序导入的核心模块。

本模块描述单谷、单自旋的连续模型，保持原 tbg.f90 的动量和层间耦合约定。
只包含模型、动量路径和本征值计算；实际运行任务位于 python/scripts/。

单位约定
--------
长度：nm；动量：nm^-1；哈密顿量和本征值：meV；输入转角：度。
原 Fortran 的 vF = 5.944e-10 eV*m 实际表示 hbar*v_F，
换算后为 594.4 meV*nm。因此 hbar*v_F*k 已经是能量，不能再乘一次 hbar。

基底约定
--------
|layer, n1, n2, sublattice>，先下层，后上层。
每层先遍历 n1，再遍历 n2；每个平面波内依次排列 A、B 子晶格。
设 c=cutoff，N_G=(2c+1)^2，则每层有 2*N_G 个态，总维数为 4*N_G。

与 Fortran 的主要对应关系
-------------------------
BMModel.hamiltonian       <-> HTBG（本版本取 te=theta）
BMModel._build_interlayer <-> T1/T2/T3、ELEMENTS、ASMBL 的层间部分
make_k_path               <-> KPATH
BMModel.layer_z_diagonal  <-> Z_TBG
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BMParameters:
    """集中存放一组模型参数。

    dataclass 自动生成构造函数；frozen=True 防止直接修改参数字段。
    原因是模型会缓存倒格矢和层间耦合：改变参数时应创建新的 BMModel，
    让这些缓存一起重算，而不是只改一个 w0 或 theta 字段。

    默认 1.09 度是保留原 hbar*v_F 时的近魔角示例，不代表普适魔角。
    w0=w1=110 meV 时，中心两带仍存在有限色散。
    """

    theta_deg: float = 1.09                  # 两层相对转角，单位：度
    w0_mev: float = 110.0                   # AA/BB，即同子晶格层间耦合
    w1_mev: float = 110.0                   # AB/BA，即不同子晶格层间耦合
    cutoff: int = 4                        # -cutoff <= n1,n2 <= cutoff
    lattice_nm: float = np.sqrt(3) * 0.142  # 晶格常数；0.142 nm 是碳碳键长
    hbar_vf_mev_nm: float = 594.4            # hbar*v_F，单位：meV*nm

    def __post_init__(self):
        """构造参数对象后，检查会使几何或矩阵维数失效的输入。"""
        scalars = (self.theta_deg, self.w0_mev, self.w1_mev,
                   self.lattice_nm, self.hbar_vf_mev_nm)
        if not np.all(np.isfinite(scalars)):
            raise ValueError('模型参数必须是有限实数')
        if self.theta_deg <= 0 or self.lattice_nm <= 0 or self.hbar_vf_mev_nm <= 0:
            raise ValueError('转角、晶格常数和 hbar*v_F 必须为正数')
        if not isinstance(self.cutoff, int) or self.cutoff < 0:
            raise ValueError('cutoff 必须为非负整数')


class BMModel:
    """给定参数后，提供任意动量处的 H(k)。

    初始化时预先生成平面波基底、动量偏移和层间耦合。
    hamiltonian(k) 只需更新各平面波上的 2x2 Dirac 块。

    常用属性
    --------
    dimension : 总矩阵维数，cutoff=4 时为 324。
    pairs     : 形状 (N_G, 2)，每行为整数对 (n1,n2)。
    G         : 形状 (N_G, 2)，每行为 n1*b1+n2*b2，单位 nm^-1。
    b1, b2    : 两个莫尔倒格矢，单位 nm^-1。
    """

    def __init__(self, parameters=BMParameters()):
        self.parameters = parameters
        self.theta = np.deg2rad(parameters.theta_deg)

        # 两层 Dirac 点的间距 k_theta = 2*K_D*sin(theta/2)，
        # 其中 K_D=4*pi/(3*a)。q1、b1、b2 的方向与原 Fortran 完全一致。
        self.k_theta = 8 * np.pi * np.sin(self.theta / 2) / (3 * parameters.lattice_nm)
        self.q1 = self.k_theta * np.array([0.0, -1.0])
        self.b1 = self.k_theta * np.array([np.sqrt(3) / 2, -1.5])
        self.b2 = self.k_theta * np.array([np.sqrt(3) / 2, 1.5])

        # meshgrid(indexing='ij') 保持 n1 在外层、n2 在内层的遍历顺序。
        # stack 后形状为 (2c+1,2c+1,2)，再整理成逐行的 (n1,n2) 列表。
        n = np.arange(-parameters.cutoff, parameters.cutoff + 1)
        n1_grid, n2_grid = np.meshgrid(n, n, indexing='ij')
        self.pairs = np.stack([n1_grid, n2_grid], axis=-1).reshape(-1, 2)
        self.G = self.pairs @ np.array([self.b1, self.b2])
        self.n_g = len(self.G)
        self.dimension = 4 * self.n_g

        # 下层 p_b=k-q1-G，上层 p_t=k-G。
        # offsets 的前 N_G 行属于下层，后 N_G 行属于上层，形状 (2*N_G,2)。
        self._offsets = np.vstack([-self.q1 - self.G, -self.G])

        # 两层 Dirac 矩阵分别旋转 +theta/2、-theta/2。
        # 相位 e^(-i*phi) 乘在 AB 元素上；BA 元素必须取其复共轭。
        angles = np.repeat([self.theta / 2, -self.theta / 2], self.n_g)
        self._phases = np.exp(-1j * angles)

        # Python 从 0 开始：A 子晶格索引为 0,2,4,...，B 索引为 A+1。
        self._a_indices = 2 * np.arange(2 * self.n_g)
        self._h_interlayer = self._build_interlayer()

    def _build_interlayer(self):
        """构造与 k 无关的层间部分 [[0,U],[U^dagger,0]]，单位 meV。

        U 的行对应下层，列对应上层，每个 (n1,n2) 占据两个子晶格态。
        字典 index 将整数对映射到单层平面波编号 g：
            g = (n1+c)*(2c+1) + (n2+c)。
        完整基底索引为 I = 2*(layer*N_G+g)+s，
        其中 layer=0/1 表示下/上层，s=0/1 表示 A/B。
        """
        p = self.parameters
        index = {tuple(pair): i for i, pair in enumerate(self.pairs)}
        U = np.zeros((2 * self.n_g, 2 * self.n_g), dtype=complex)

        # 每个连接由 (上层 n1 偏移，上层 n2 偏移，子晶格耦合相位) 指定。
        # T1: (n1,n2) -> (n1,n2)
        # T2: (n1,n2) -> (n1,n2-1)
        # T3: (n1,n2) -> (n1+1,n2)
        links = [(0, 0, 0), (0, -1, 2 * np.pi / 3), (1, 0, -2 * np.pi / 3)]
        for dn1, dn2, phase in links:
            # T(phi)=w0*I+w1*(cos(phi)*sigma_x+sin(phi)*sigma_y)。
            # 用复指数表示非对角元，可以直接看出 AB 与 BA 的共轭关系。
            z = p.w1_mev * np.exp(-1j * phase)
            T = np.array([[p.w0_mev, z], [z.conjugate(), p.w0_mev]])
            for i, (n1, n2) in enumerate(self.pairs):
                j = index.get((n1 + dn1, n2 + dn2))
                if j is not None:
                    # i 对应的 A/B 行是 2i,2i+1；j 对应的列同理。
                    # 找不到 j 表示连接超出平面波截断，直接略去，不做周期回绕。
                    U[2*i:2*i+2, 2*j:2*j+2] = T

        zero = np.zeros_like(U)
        # 另一层间块必须是 U^dagger，而不是 U.T；这里显式保证厄米性。
        return np.block([[zero, U], [U.conj().T, zero]])

    def hamiltonian(self, k_nm):
        """构造一个动量点上的完整哈密顿量。

        Parameters
        ----------
        k_nm : array-like, shape (2,)
            [kx,ky]，单位 nm^-1；不是沿高对称路径的累积距离。

        Returns
        -------
        H : complex ndarray, shape (dimension,dimension)
            单位 meV。基底顺序固定，H=H^dagger。
            每次返回独立数组，可以在调用方添加势能等项。

        单个 Dirac 块为 h(p,phi)=[[0,z],[z*,0]]，其中
            z = hbar*v_F*(px-i*py)*exp(-i*phi)。
        它等价于原 HSLG 的极坐标写法，且在 p=0 时无需处理 atan2。
        """
        k = np.asarray(k_nm, dtype=float)
        if k.shape != (2,) or not np.all(np.isfinite(k)):
            raise ValueError('k 必须是形状为 (2,) 的有限二维动量')

        # NumPy 广播：一个 (2,) 的 k 与所有 (2*N_G,2) 动量偏移相加。
        momenta = k + self._offsets
        z = (self.parameters.hbar_vf_mev_nm
             * (momenta[:, 0] - 1j * momenta[:, 1]) * self._phases)

        H = self._h_interlayer.copy()  # 避免本次更新改变缓存，污染下一个 k 点。
        a = self._a_indices
        H[a, a + 1] = z              # 同一层、同一 G 内的 A -> B 元素
        H[a + 1, a] = z.conj()       # 对应的 B -> A 元素
        return H

    def layer_z_diagonal(self, separation_nm=0.335):
        """返回 z 算符的对角元素，shape (dimension,)，单位 nm。

        下层取 -d/2，上层取 +d/2；每层有 2*N_G 个子晶格态。
        只保存对角线可以节省内存。需要完整矩阵时使用 np.diag(z)。
        separation_nm 是层间距，与面内碳碳键长 0.142 nm 不同。
        """
        return np.repeat([-separation_nm / 2, separation_nm / 2], 2 * self.n_g)


def make_k_path(model, points_per_segment=60):
    """生成原 Fortran 中的 K_b -> Gamma -> M -> K_t 路径。

    points_per_segment 指每段的区间数，总动量点数为 3*n+1。
    返回 k_points (3*n+1,2)、distance (3*n+1,)、ticks (4,)，单位均为 nm^-1。
    distance 是沿路径的累积弧长，ticks 用作绘图中四个高对称点的位置。
    """
    if not isinstance(points_per_segment, int) or points_per_segment < 1:
        raise ValueError('points_per_segment 必须为正整数')
    vertices = model.k_theta * np.array([
        [0.0, -1.0],                         # K_b
        [np.sqrt(3) / 2, -0.5],              # Gamma：此动量约定下不在原点
        [0.0, -0.5],                        # M
        [0.0, 0.0],                         # K_t
    ])

    # 每段包含起点、不包含终点，避免相邻段重复收录同一高对称点。
    # 最后再补上整条路径的终点 K_t。
    segments = [np.linspace(start, end, points_per_segment, endpoint=False)
                for start, end in zip(vertices[:-1], vertices[1:])]
    k_points = np.vstack([*segments, vertices[-1:]])
    step_lengths = np.linalg.norm(np.diff(k_points, axis=0), axis=1)
    distance = np.r_[0.0, np.cumsum(step_lengths)]
    ticks = distance[np.arange(4) * points_per_segment]
    return k_points, distance, ticks


def solve_bands(model, k_points):
    """计算一组动量上的所有能量，返回 shape (N_k,dimension)，单位 meV。

    使用厄米矩阵专用的 eigvalsh，结果按能量从低到高排列。
    这里只计算能量、不计算本征矢；需要本征矢时可调用：
        energies, vectors = np.linalg.eigh(model.hamiltonian(k))
    vectors 的每一列是一条能带的本征矢。

    本函数保留原始能量，不平移、不缩放；能量零点由绘图脚本统一选择。
    """
    points = np.asarray(k_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) == 0:
        raise ValueError('k_points 必须是非空的 (N_k,2) 数组')
    return np.array([np.linalg.eigvalsh(model.hamiltonian(k)) for k in points])
