program TBG_STRAIN_OM_MU
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF_MU
    use stdlib_math, only : linspace
    implicit none

    integer, parameter :: N = 2*NTBG

    integer :: nk, ncut
    integer :: imu
    integer :: fileunit

    real(8) :: theta, strain_epsilon, strain_phi
    real(8) :: d
    real(8) :: delta, delta_h
    real(8) :: prefactor

    real(8), allocatable :: mu_list(:)
    real(8), allocatable :: Mx_list(:), My_list(:)

    real(8) :: bs(2,2)
    real(8) :: ham_param(3)
    real(8) :: zdiag(N)

    character(len=256) :: filename
    character(len=24) :: s_theta, s_epsilon, s_phi
    character(len=24) :: s_mu_min, s_mu_max
    character(len=16) :: s_nk, s_ncut

    ! -----
    ! 固定参数
    ! -----
    d = 3.35d-10
    theta = (1.40d0)*PI/180.d0

    strain_epsilon = 1.d-3
    strain_phi = (90.d0)*PI/180.d0

    prefactor = 2.d0*m_e/hbar**2

    ! -----
    ! 可变参数
    ! -----
    delta = 1.d-3 ! 中心差分步长
    ncut = 100 ! 对称截断数
    nk = 120 ! BZ mesh

    ! ----------
    ! 化学势
    ! ----------
    mu_list = linspace(-20.d0, 20.d0, 81)*(1.d-3)*eV

    ! ----- allocate the oribtal moment arrays ----- 
    allocate(Mx_list(size(mu_list)), My_list(size(mu_list)))

    ! =====
    ! 构造哈密顿量
    ! =====
    ! ----- Hamiltonian parameters ----- 
    ham_param = [theta, strain_epsilon, strain_phi]
    ! ----- moire reciprocal lattice vectors ----- 
    call GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs)
    ! ----- 中心差分步长 ----- 
    delta_h = delta*(norm2(bs(:,1)))

    ! ----------
    ! 层间位置算符
    ! ----------
    call Z_TBG(d, zdiag)

    ! ----------
    ! 字符串参数
    ! ----------
    write(s_theta, '(F8.3)') theta*180.d0/PI
    write(s_epsilon, '(F8.3)') strain_epsilon*100.d0
    write(s_phi, '(F8.3)') strain_phi*180.d0/PI

    write(s_mu_min, '(F8.3)') mu_list(1)/(1.d-3*eV)
    write(s_mu_max, '(F8.3)') mu_list(size(mu_list))/(1.d-3*eV)

    write(s_nk, '(I0)') nk
    write(s_ncut, '(I0)') ncut

    ! ==========
    ! 输出文件
    ! ==========
    ! ----- 文件名 ----- 
    filename = '../data/tbg/fortran/om/mu/tbg_om_mu_theta_'//trim(adjustl(s_theta))// &
               '_epsilon_'//trim(adjustl(s_epsilon))// &
               '_phi_'//trim(adjustl(s_phi))// &
               '_mu_'//trim(adjustl(s_mu_min))//'to'//trim(adjustl(s_mu_max))// &
               '_nk_'//trim(adjustl(s_nk))// &
               '_ncut_'//trim(adjustl(s_ncut))//'.dat'

    ! ----- open file ----- 
    open(newunit=fileunit, file = trim(filename), status = 'replace')

    ! ----- Output header -----*
    write(fileunit,'(A)') '# TBG in-plane orbital magnetization vs chemical potential'

    write(fileunit,'(A,F12.6)') '# theta_deg = ', theta*180.d0/PI
    write(fileunit,'(A,F12.6)') '# epsilon_percent = ', strain_epsilon*100.d0
    write(fileunit,'(A,F12.6)') '# phi_deg = ', strain_phi*180.d0/PI

    write(fileunit,'(A,I0)') '# nk = ', nk
    write(fileunit,'(A,I0)') '# ncut = ', ncut
    write(fileunit,'(A,ES16.8)') '# delta = ', delta

    write(fileunit,'(A,2F12.6)') '# mu_range_meV = ', &
        mu_list(1)/(1.d-3*eV), mu_list(size(mu_list))/(1.d-3*eV)
    write(fileunit,'(A,I0)') '# n_mu = ', size(mu_list)

    write(fileunit,'(A)') '# mu(meV) Mx(mu_B/cell) My(mu_B/cell)'

    ! ----- 计算轨道磁化 ----- 
    call CALCULATE_OM_CUTOFF_MU(N, nk, ncut, bs(:,1), bs(:,2), delta_h, mu_list, prefactor, zdiag, ham_param, HAM_STRAIN,&
                                Mx_list, My_list)

    ! ----- 输出 ----- 
    do imu = 1, size(mu_list)

        write(fileunit, '(F10.4, 2ES32.16E3)') &
            mu_list(imu)/(1.d-3*eV), Mx_list(imu), My_list(imu)

    enddo

    close(fileunit)
    write(*,'(A,A)') 'Results saved to: ', trim(filename)

end program TBG_STRAIN_OM_MU