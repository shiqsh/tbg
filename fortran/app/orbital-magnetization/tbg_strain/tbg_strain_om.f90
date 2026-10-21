program TBG_STRAIN_OM
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF
    use stdlib_math, only : linspace
    implicit none

    integer, parameter :: N = 2*NTBG

    integer :: nk, ncut
    integer :: iepsilon, iphi
    integer :: ie_start, ie_end, ip_start, ip_end
    integer :: fileunit

    real(8) :: theta, strain_epsilon, strain_phi
    real(8) :: mu, d
    real(8) :: delta, delta_h
    real(8) :: prefactor

    real(8), allocatable :: epsilon_list(:)
    real(8), allocatable :: phi_list(:)
    
    real(8) :: bs(2,2)
    real(8) :: ham_param(3)
    real(8) :: zdiag(N)

    real(8) :: Mx, My

    character(len=256) :: filename
    character(len = 16) ::  s_epsilon, s_phi
    character(len=16) :: s_nk, s_ncut, s_mu, s_theta

    ! -----
    ! 固定参数
    ! -----
    d = 3.35d-10
    theta = (1.40d0)*PI/180.d0
    prefactor = 2.d0*m_e/hbar**2

    ! -----
    ! 可变参数
    ! -----
    mu = 0.d0 
    delta = 1.d-3 ! 中心差分步长
    ncut = 80 ! 对称截断数
    nk = 80 ! BZ mesh

    ! -----
    ! strian paramter
    ! -----
    epsilon_list = linspace(0.d0, 0.d0, 1)
    phi_list = linspace(0.d0, PI, 181)

    ! =====
    ! 选择扫描范围
    ! start = end : 固定参数
    ! start < end : 扫描参数
    ! =====
    ! ----------
    ! 固定phi = 90 deg, 扫描epsilon
    ! ----------
    ie_start = 1
    ie_end = size(epsilon_list)

    ip_start = 31
    ip_end = 31
    ! ----- 检查索引 ----- 
    if (ie_start < 1 .or. ie_end > size(epsilon_list) &
        .or. ie_start > ie_end) then
        error stop "Invalid epsilon index range"
    endif

    if (ip_start < 1 .or. ip_end > size(phi_list) &
        .or. ip_start > ip_end) then
        error stop "Invalid phi index range"
    endif

    ! ----------
    ! 层间位置算符
    ! ----------
    call Z_TBG(d, zdiag)

    ! ----------
    ! 字符串参数
    ! ----------
    write(s_theta, '(F8.3)') theta*180.d0/PI
    write(s_mu, '(F8.3)') mu/eV*1000.d0 ! meV输出
    write(s_nk, '(I0)') nk
    write(s_ncut, '(I0)') ncut

    ! ----- strain epsilon ----- 
    if (ie_start == ie_end) then
        write(s_epsilon, '(F8.3)') epsilon_list(ie_start)*100.d0
    else
        write(s_epsilon, '(A,I0,A,I0)') 'scan_', ie_start, 'to', ie_end
    endif

    ! ----- strain phi ----- 
    if (ip_start == ip_end) then
        write(s_phi, '(F8.3)') phi_list(ip_start)*180.d0/PI
    else
        write(s_phi, '(A,I0,A,I0)') 'scan_', ip_start, 'to', ip_end
    endif

    ! ==========
    ! 输出文件
    ! ==========
    ! ----- file name ----- 
    filename = '../data/tbg/fortran/om/tbg_om_theta'//trim(adjustl(s_theta))// &
               '_epsilon_'//trim(adjustl(s_epsilon))// &
               '_phi_'//trim(adjustl(s_phi))// &
               '_mu_'//trim(adjustl(s_mu))// &
               '_nk_'//trim(adjustl(s_nk))// &
               '_ncut_'//trim(adjustl(s_ncut))//'.dat'

    ! ----- open file ----- 
    open(newunit=fileunit, file = trim(filename), status = 'replace')

    ! ----- output the header ----- 
    write(fileunit,'(A)') '# TBG under heterostrain in-plane orbital magnetization'
    write(fileunit,'(A,F12.6)') '# theta_deg = ', theta*180.d0/PI
    write(fileunit,'(A,F12.6)') '# mu_meV = ', mu/eV*1000.d0
    write(fileunit,'(A,I0)') '# nk = ', nk
    write(fileunit,'(A,I0)') '# ncut = ', ncut
    write(fileunit,'(A,ES16.8)') '# delta = ', delta

    ! ----- fixed or scan epsilon ----- 
    if (ie_start == ie_end) then
        write(fileunit,'(A,F12.6)') &
            '# epsilon_fixed_percent = ', epsilon_list(ie_start)*100.d0
    else
        write(fileunit,'(A,2F12.6)') &
            '# epsilon_scan_percent = ', &
            epsilon_list(ie_start)*100.d0, epsilon_list(ie_end)*100.d0
    endif
    ! ----- fixed or scan phi ----- 
    if (ip_start == ip_end) then
        write(fileunit,'(A,F12.6)') &
            '# phi_fixed_deg = ', phi_list(ip_start)*180.d0/PI
    else
        write(fileunit,'(A,2F12.6)') &
            '# phi_scan_deg = ', &
            phi_list(ip_start)*180.d0/PI, phi_list(ip_end)*180.d0/PI
    endif

    ! ----- Data columns ----- 
    write(fileunit,'(A)') '# epsilon(%) phi(deg) Mx(mu_B/cell) My(mu_B/cell)'

    ! ----- output data ----- 
    do iepsilon = ie_start, ie_end
        do iphi = ip_start, ip_end 
            ! ----- current strian parameter ----- 
            strain_epsilon = epsilon_list(iepsilon)
            strain_phi = phi_list(iphi)
            ! ----- Hamiltonian parameters ----- 
            ham_param = [theta, strain_epsilon, strain_phi]
            ! ----- moire reciprocal lattice vectors ----- 
            call GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs)
            ! ----- 中心差分步长 ----- 
            delta_h = delta*(norm2(bs(:,1)))

            call CALCULATE_OM_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), delta_h, mu, prefactor, zdiag, ham_param, HAM_STRAIN, Mx, My)

            write(fileunit, '(4ES32.16E3)') strain_epsilon*100.d0, strain_phi*180.d0/PI, Mx, My

            ! ----- screen output ----- 
            ! write(*,'(A,F8.4,A,F8.2,A,2ES16.6)') 'epsilon(%) = ', strain_epsilon*100.d0, ' phi(deg) = ', strain_phi*180.d0/PI, ' Mx, My = ', Mx, My
        
        enddo
    enddo

    close(fileunit)

    write(*,'(A,A)') 'Results saved to: ', trim(filename)

end program TBG_STRAIN_OM