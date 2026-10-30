# TBG 面内轨道磁化：原始巨势独立核验

日期：2026-10-30。所有磁矩单位为 **mu_B / moire cell，单自旋、单谷**。

## 判断

当前 `om_cutoff` 的结果与有限 BM 矩阵的**未正则化全占据带巨势导数不相等**。
差异的主要来源是有限平面波截断留下的 BZ 边界项，不是电荷符号、Bohr 磁子换算或缺少简并因子。
因此不能把原始巨势的数值直接视作物理标准答案，也不能据此宣称现有 `ncut=100` 已得到收敛的绝对轨道磁化。

原有 ket 和 projector 实现彼此一致。完整带空间中的公式，与热力学公式在 BZ 周期性及相应积分成立时等价；本次推导和独立边界积分确认了有限截断下二者之间的关系。现有 `ncut=100` 还有可见的带窗口依赖，应在确定高能正则化/参考态后进一步检查 `tr`、`ncut`、`nk`。

## 模型、程序和参数

- 保留原有 `tbg_strain_om.f90`、`om_cutoff.f90` 和物理参数不变；只新增独立磁场模块、对照程序和回归测试。
- 主对照固定 theta=1.4°、epsilon=0.1%、phi=0°、mu=0.5 meV、d=3.35e-10 m。
- 与现有生产数据第一行一致的网格为 nk=40；tr=4 对应 N=324，完整带窗口 ncut=161。
- 补充计算使用 nk=8、12；B0=0.1、0.01、0.001 T；tr=3、4、5；phi=45°及零应变。
- 工作区计算期间有其他修改，包含 tr 从 4 变成 5。本次使用独立源码快照，并在每份结果首行打印实际 tr、N、nk；没有改写工作区的 constants.f90。
- `ncut` 是能带求和窗口，`tr` 是哈密顿量平面波基底截断，两者不可混用。
- 所有数值使用工程自身的 hbar、m_e、eV，不另换物理常数。

## 从巨势得到的独立公式

取电子电荷 -e，e>0，A=(z By,-z Bx,0)，则 curl A=B。
当前 BM 模型对 k 线性，层间隧穿不依赖 k，且 [z,Hx]=[z,Hy]=0，因此

```math
H(B)=H(0)+\frac{e}{\hbar}z(B_y H_x-B_x H_y),\qquad H_i=\partial_{k_i}H.
```

用户提供的磁场耦合符号正确。新增模块采用反对易形式保持厄米性；在本模型中它与逐行乘 z 完全相同。

```math
\Omega(B)=\left\langle\sum_n\min(E_n(B)-\mu,0)\right\rangle_k,
\qquad M_i=-\partial_{B_i}\Omega/\mu_B.
```

独立的 Hellmann–Feynman 校验为

```math
\mathbf M_\Omega=C\left\langle\sum_{n\in\mathrm{occ}}
  \big(\mathrm{Re}\langle n|zH_y|n\rangle,
      -\mathrm{Re}\langle n|zH_x|n\rangle\big)\right\rangle_k,
\qquad C=2m_e/\hbar^2.
```

这里 k 平均为 `sum / nk**2`，对应每个 moire 原胞；无需再乘或除 BZ 面积。

## 为什么与原有公式不同

定义完整带空间中的量

```math
F_i=\mathrm{Re}\sum_{n\in\mathrm{occ},m\in\mathrm{empty}}
\frac{E_n+E_m-2\mu}{E_n-E_m}(H_i)_{nm}z_{mn},
\quad D_i=\sum_{n\in\mathrm{occ}}\mathrm{Re}\langle n|zH_i|n\rangle,
\quad G=\sum_{n\in\mathrm{occ}}(E_n-\mu)z_{nn}.
```

对 G 微分，利用本征态微扰公式，先把占据—占据矩阵元成对求和，可得

```math
\partial_{k_i}G=D_i+F_i.
```

占据变化的分布项含有 `(E_n-mu) delta(E_n-mu)`，因此不额外贡献该关系。需要完整本征态空间；简并处应理解为投影算符/适当极限。

代码返回 `M_formula=C*(-<Fy>,<Fx>)`，所以

```math
\mathbf M_{\rm formula}^{\rm full}
=\mathbf M_\Omega-\mathbf M_{\rm boundary},\qquad
\mathbf M_{\rm boundary}=C\langle(\partial_{k_y}G,-\partial_{k_x}G)\rangle_k.
```

若 G 在 BZ 上周期，则边界项为零。固定的有限 G 向量矩形截断不能在 `k -> k+b_i` 下实现完整基底重标记，高能态使 G 在相对边界处不匹配；故当前边界项很大。程序通过四条 BZ 边上的 G 独立积分测量它，**不是把两种磁化相减定义为边界项**。有限 nk 下的差额还包含二维和一维求积误差。

`ncut=100` 只保留占据—空态配对的部分项，P+Q 不等于完整单位阵。即使在 Omega 中仅保留相同编号的能带，其导数也不是这个双重截断公式：被保留本征态的导数仍涉及窗口外的态。因此不应通过简单截掉 Omega 的深能带，强行使两种算法“口径一致”。

## 已观察到的核验结果

主对照（tr=4、nk=40、B0=0.01 T、phi=0°）：

| 方法 | Mx | My |
|---|---:|---:|
| 原程序 ket，ncut=100 | -0.02724630203 | -0.29260155831 |
| 全带公式，ncut=161 | -0.03994833213 | -0.29158017206 |
| 原始巨势中心差分 | 22.54764153748 | -0.49183896519 |
| 零场 Hellmann–Feynman | 22.54764138726 | -0.49183898956 |
| 独立 BZ 边界积分 | 22.58756991561 | -0.20025648634 |
| 巨势差分减去边界项 | -0.03992837814 | -0.29158247885 |

原始巨势差分与 HF 最大相差 1.51e-7；减去独立边界项后，与全带公式最大相差 2.00e-5。
后者由 nk=8 的约 7.28e-4 降至 nk=40 的约 1.98e-5（以 HF 核验），支持边界关系和求积误差解释。
`ncut=100 -> 161` 的 Mx 改变约 0.01270；这是对带窗口的敏感性，**不能直接称为相对于物理真值的误差**。

- nk=8、tr=4、phi=0°：原有 ket 为 `(-0.0427371667642,-0.287873607850)`；原有 projector（delta=1e-5 |b1|）为 `(-0.0427371663670,-0.287873607772)`，最大分量差约 4e-10。
- nk=12、tr=4、phi=0°：B0=0.1 T 时，Omega 差分与 HF 的最大分量差约 8.7e-7；B0=0.001 T 时约 4.9e-7。极小 B 受大数相减影响，不是越小越好。
- nk=12、tr=3、4、5 时，原始 HF 的 Mx 分别约 17.5077、22.5483、27.5830；对应边界项约 17.5584、22.5876、27.6170。增大 tr 并未消除原始巨势中的大背景。
- 零应变、tr=4、nk=12 时，原始 Omega 差分 Mx 约 22.5671，边界项约 22.5817。这个明显的大残余支持其截断背景解释。有限矩形基底及有限网格也可能在公式中留下较小的对称性残余。
- phi=45°、tr=4、nk=12 时，Omega 差分约 `(23.80840987,0.39344765)`，与 HF 在 3e-8 内吻合；减去独立边界项后约 `(0.28858915,0.00852508)`，全带公式为 `(0.28819237,0.00850916)`，残差与求积误差一致。

完整机器精度结果、各次计算的四个 Omega 值和运行时间见本目录 `.txt`。`summary.csv`、`comparison.png` 由绘图脚本从成功完成的日志生成。

## 物理解释的范围

本核验确认磁场耦合、单位、公式的符号，以及 ket/projector 实现一致性；也定位了“全带原始巨势 vs 公式”的大差异。它**不构成绝对物理磁化已经收敛的证明**。

若要报告绝对单谷磁化，应先明确对连续模型 Dirac 海的高能正则化，或者使用具有周期完整 BZ 的晶格模型，并在同一规定下比较磁场法与公式。若研究同一个哈密顿量的掺杂变化，可考虑比较 `M(mu)-M(mu_ref)`，以消去相同的深能带背景，但仍须检查参考态、BZ 边界和网格误差；不同应变之间相减不能自动保证相同背景被消除。

只求两个平带的磁矩是另一项定义清楚的量，本程序输出 `HF_active_only` 供区分。它既不是原始全海巨势结果，也不是 `ncut=0` 的占据—空态双重截断公式。

当前模型是单谷、单自旋。总实验磁化还要结合不同谷和自旋的占据；时间反演对称、等占据双谷的磁化相互抵消，不能直接把单谷结果乘 4。磁场耦合及谷分辨磁化的背景可参见 [Antebi, Stern & Berg, arXiv:2112.14785，§II、§IV.1](https://arxiv.org/html/2112.14785v1)。该文模型参数与本工程不同；此处所有数值来自本地代码，不是文献数值。

## 复现

在项目根目录执行，快照编译不会修改生产 constants：

```sh
python3 scripts/build_om_magnetic_audit.py --tr 4
OMP_NUM_THREADS=2 OMP_STACKSIZE=128M VECLIB_MAXIMUM_THREADS=1 \
  ./fortran/build/om_magnetic_audit_tr4_accelerate/compare_om_magnetic \
  40 0.01 0 0.1 0.5 1
```

六个位置参数分别为 `nk B0_T phi_deg epsilon_percent mu_meV verify`。B0=0 只跳过四次有限磁场对角化；verify=1 额外调用原有 ket 例程核验，verify=2 再调用原有 projector（适合小网格）。theta 固定 1.4°。程序对 mu 超出两个 active bands 的占据范围会停止，避免沿用原例程的隐含占据假设。

也可从 fortran 目录使用标准 fpm 入口，注意此时采用工作区当前 tr：

```sh
OMP_NUM_THREADS=2 OMP_STACKSIZE=128M VECLIB_MAXIMUM_THREADS=1 \
  fpm run --target compare_om_magnetic --profile release -- 12 0.01 0 0.1 0.5 1
fpm test --target check_om_magnetic
```

本次实际编译使用 gfortran 16.1.0 和本地 stdlib 依赖缓存，主要使用 Accelerate。源码快照及 SHA-256 清单位于各自 `fortran/build/om_magnetic_audit_tr*_accelerate` 目录。测试也以独立快照方式编译，启用 `-fcheck=all`；测试日志保存在 `test_result.txt`。
