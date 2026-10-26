program TBG_STRAIN_OM
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN, HAM_DERIVATIVES_STRAIN_PHI
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF, CALCULATE_OM_CUTOFF_KET
    use stdlib_math, only : linspace
    ! 新增：显式引用标准输出单元；flush 后，屏幕与 tee 保存的日志可以及时显示进度。
    use, intrinsic :: iso_fortran_env, only : output_unit
    ! 新增：墙钟计时包含实际等待时间；线程数用于在日志中核对本次运行配置。
    use omp_lib, only : omp_get_wtime, omp_get_max_threads
    implicit none

    integer, parameter :: N = 2*NTBG

    integer :: nk, ncut
    integer :: iepsilon, iphi
    integer :: ie_start, ie_end, ip_start, ip_end
    integer :: fileunit

    ! 新增：记录扫描点序号、总点数，以及每点和整次计算的墙钟耗时。
    integer :: ipoint, npoints
    real(8) :: run_started, point_started, point_elapsed

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

    ! 新增：从参数初始化开始计时，最终与每个扫描点的耗时一起打印。
    run_started = omp_get_wtime()

    ! -----
    ! 固定参数
    ! -----
    d = 3.35d-10
    theta = (1.40d0)*PI/180.d0
    prefactor = 2.d0*m_e/hbar**2

    ! -----
    ! 可变参数
    ! -----
    mu = (0.5d0)*1.d-3*eV
    delta = 1.d-3 ! 中心差分步长
    ncut = 100 ! 对称截断数
    nk = 40 ! BZ mesh,\

    ! -----
    ! strian paramter
    ! -----
    epsilon_list = linspace(0.d0, 2.d-3, 21)
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

    ip_start = 21
    ip_end = 21
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
    ! 新增：与 CHI 程序一样使用绝对目录，方便从日志直接定位数值结果文件。
    filename = '/Users/shiqsh/Workspace/tbg/data/tbg/fortran/om/tbg_om_theta'//trim(adjustl(s_theta))// &
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

    ! 新增：立即刷新文件头，让计算尚未完成时也能检查本次输出文件和参数。
    flush(fileunit)

    ! 新增：统计全部扫描点；原有的单参数或双参数扫描顺序均保持不变。
    npoints = (ie_end-ie_start+1)*(ip_end-ip_start+1)
    ipoint = 0

    ! 新增：日志由 tee 同时显示和保存，本程序不再重复打开同一个 progress.log。
    ! 运行示例：frt tbg_strain_om 2>&1 | tee /Users/shiqsh/Workspace/tbg/data/tbg/fortran/om/progress.log
    ! 新增：启动时输出算法、数值参数、物理参数和路径，便于核对是否运行了预期配置。
    write(*,'(A)') 'TBG strain OM: ket calculation starting'
    write(*,'(A,I0,A,I0,A,I0)') 'N = ', N, ', nk = ', nk, ', ncut = ', ncut
    write(*,'(A,F10.4,A,F10.4)') 'theta(deg) = ', theta*180.d0/PI, ', mu(meV) = ', mu/eV*1000.d0
    write(*,'(A,ES12.4)') 'delta = ', delta
    write(*,'(A,I0)') 'OpenMP max threads = ', omp_get_max_threads()
    write(*,'(A,I0)') 'Total scan points = ', npoints
    write(*,'(A,A)') 'Data file: ', trim(filename)
    ! 新增：说明一行结果对应完整的 BZ 积分，防止将首行尚未产生误判为写入失败。
    write(*,'(A,I0,A)') 'Each scan point integrates ', nk*nk, ' k points before saving one data row.'
    ! 新增：即使标准输出经管道交给 tee，也立即显示上面的启动信息。
    flush(output_unit)

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

            ! 新增：与 CHI 相同，在耗时计算开始前输出点序号、当前应变参数并立即刷新。
            ipoint = ipoint + 1
            point_started = omp_get_wtime()
            write(*,'(A,I0,A,I0,A,F10.5,A,F10.3,A)') 'point ', ipoint, '/', npoints, &
                ' epsilon(%) = ', strain_epsilon*100.d0, ' phi(deg) = ', strain_phi*180.d0/PI, ': OM starting'
            flush(output_unit)

            ! call CALCULATE_OM_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), delta_h, mu, prefactor, zdiag, ham_param, HAM_STRAIN, Mx, My)
            call CALCULATE_OM_CUTOFF_KET(N, nk, ncut, bs(:,1), bs(:,2), mu, prefactor,zdiag, ham_param, &
                                        HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, Mx, My)

            write(fileunit, '(4ES32.16E3)') strain_epsilon*100.d0, strain_phi*180.d0/PI, Mx, My

            ! 新增：每个扫描点完成后立即刷新数值行，再报告 saved，保留已完成的结果。
            flush(fileunit)
            ! 新增：计时覆盖本点的 OM 计算和结果写入，便于检查实际运行速度。
            point_elapsed = omp_get_wtime()-point_started

            ! ----- screen output ----- 
            ! 新增：逐点输出与数据文件对应的 Mx、My，便于在屏幕及日志中交叉核对。
            write(*,'(A,I0,A,I0,A,F10.5,A,F10.3,2(A,ES32.16))') &
                'point ', ipoint, '/', npoints, &
                ' epsilon(%) = ', strain_epsilon*100.d0, ' phi(deg) = ', strain_phi*180.d0/PI, &
                ' Mx = ', Mx, ' My = ', My
            ! 新增：只有数值文件刷新成功后才打印保存完成信息，并给出每点耗时。
            write(*,'(A,F12.2,A)') '  saved; point wall time = ', point_elapsed, ' s'
            flush(output_unit)
        
        enddo
    enddo

    close(fileunit)

    write(*,'(A,A)') 'Results saved to: ', trim(filename)
    ! 新增：正常结束时报告全部扫描完成及总耗时，便于从日志确认本次任务完成。
    write(*,'(A,I0,A,I0,A,F12.2,A)') 'Completed ', ipoint, '/', npoints, &
        ' points; total wall time = ', omp_get_wtime()-run_started, ' s'
    flush(output_unit)

end program TBG_STRAIN_OM
