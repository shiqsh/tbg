program TBG_STRAIN_CHI_NCUT
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN, &
                        HAM_DERIVATIVES_STRAIN_EPSILON, &
                        HAM_DERIVATIVES_STRAIN_PHI
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF_KET
    use SUSCEPTIBILITY_CUTOFF, only : CALCULATE_CHI_CUTOFF
    use stdlib_math, only : linspace
    use, intrinsic :: iso_fortran_env, only : output_unit !明确指定输出到标准输出
    use omp_lib, only : omp_set_max_active_levels, omp_get_wtime !计算程序实际运行时间的函数

    ! ==========
    ! 参数声明
    ! ==========
    integer, parameter :: N = 2*NTBG
    ! ----- ncut 数组 ----- 
    integer, allocatable :: ncut_list(:)
    ! ----- 确定应力参数 -----
    character(len=*), parameter :: response_parameter = 'strain_epsilon'

    integer :: nk, ncut
    integer :: icut, npoints
    integer :: fileunit, ios

    ! ----- 几何参数 ----- 
    real(8) :: theta, strain_epsilon, strain_phi
    real(8) :: mu, eta, d, prefactor

    ! ----- BZ和哈密顿量参数 -----
    real(8) :: bs(2, 2)
    real(8) :: ham_param(3)
    real(8) :: zdiag(N)

    ! ----- OM 和 极化率 ----- 
    real(8) :: Mx, My
    real(8) :: chi_fix_x, chi_fix_y
    real(8) :: chi_fs_x, chi_fs_y
    real(8) :: chi_tot_x, chi_tot_y

    ! ----- 计时变量 ----- 
    real(8) :: run_started, point_started, chi_elapsed

    ! ----- 输出文件名相关字符串 -----
    character(len=512) :: filename
    character(len=24) :: s_theta, s_epsilon, s_phi
    character(len=24) :: s_mu, s_eta, s_nk
    character(len=8) :: date_tag
    character(len=10) :: time_tag

    ! ==========
    ! 参数赋值
    ! ==========
    d              = 3.35d-10
    theta          = 1.4d0*PI/180.d0
    strain_epsilon = 1.d-3             
    strain_phi     = 20.d0*PI/180.d0

    ! ----- 化学势，高斯展宽与OM系数 -----
    mu             = 0.5d-3*eV        ! 0.5 meV
    eta            = 0.25d-3*eV
    prefactor = 2.d0*m_e/hbar**2

    ! ----- BZ mesh and cut ----- 
    nk = 40
    ncut_list = nint(linspace(20.d0, 160.d0, 8))

    ! ==========
    ! 输入参数与截断检查
    ! ==========
    if (N < 2 .or. mod(N,2) /= 0) error stop 'TBG_STRAIN_CHI_NCUT: invalid N'
    if (nk <= 0)  error stop 'TBG_STRAIN_CHI_NCUT: invalid nk'
    if (eta <= 0.d0) error stop 'TBG_STRAIN_CHI_NCUT: invalid eta'

    ! ----- ncut 检查 ----- 
    if (size(ncut_list) == 0) error stop 'TBG_STRAIN_CHI_NCUT: empty ncut_list' 
    if (any(ncut_list < 0) .or. any(ncut_list > N/2 - 1)) then
        error stop 'TBG_STRAIN_CHI_NCUT :  invalid ncut'
    end if

    ! ----- 极化率响应参数检查 -----
    if (response_parameter /= 'strain_epsilon' .and. &
        response_parameter /= 'strain_phi') then
        error stop 'TBG_STRAIN_CHI_NCUT: invalid response parameter'
    endif
    
    ! ==========
    ! 构造哈密顿量和几何
    ! ==========
    ! ----- 禁止 OpenMP 嵌套并行 ----- 
    ! CHI 和 ket OM 的子程序各自在BZ k-mesh 上并行
    call omp_set_max_active_levels(1)

    ! ----- 哈密顿量参数 ----- 
    ham_param = [theta, strain_epsilon, strain_phi]

    ! ----- BZ geometry ----- 
    call GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs )

    ! ----- 层间位置算符 -----
    call Z_TBG(d, zdiag)

    ! ==========
    ! 构造输出文件路径
    ! ==========
    ! ----- 将数值参数转换为字符串 -----
    write(s_theta,   '(F8.3)') theta*180.d0/PI
    write(s_epsilon, '(F8.4)') strain_epsilon*100.d0
    write(s_phi,     '(F8.3)') strain_phi*180.d0/PI
    write(s_mu,      '(F8.3)') mu/eV*1000.d0
    write(s_eta,     '(F8.3)') eta/eV*1000.d0
    write(s_nk,      '(I0)') nk

    ! ----- 添加日期和时间标签，避免覆盖已有文件 -----  
    call date_and_time(date=date_tag, time=time_tag)

    ! ----- 构造绝对输出路径 -----
    filename = '/Users/shiqsh/Workspace/tbg/data/tbg/fortran/susceptibility/ncut/'// &
               'tbg_chi_ncut_'//trim(response_parameter)// &
               '_theta'//trim(adjustl(s_theta))// &
               '_epsilon_'//trim(adjustl(s_epsilon))// &
               '_phi_'//trim(adjustl(s_phi))// &
               '_mu_'//trim(adjustl(s_mu))// &
               '_eta_'//trim(adjustl(s_eta))// &
               '_nk_'//trim(adjustl(s_nk))//'.dat'

    ! ----- 新建文件，禁止覆盖 -----
    open(newunit=fileunit, file=trim(filename), &
        status='new', action='write', iostat=ios)

    if (ios /= 0) then
        error stop 'TBG_STRAIN_CHI_NCUT: cannot create output file'
    endif

    ! ==========
    ! 输出文件头
    ! ==========
    ! ----- 物理参数 -----
    write(fileunit,'(A)') '# TBG fixed-strain scan over ncut'
    write(fileunit,'(A,F12.6)') '# theta_deg = ', theta*180.d0/PI
    write(fileunit,'(A,F12.6)') '# epsilon_percent = ', strain_epsilon*100.d0
    write(fileunit,'(A,F12.6)') '# phi_deg = ', strain_phi*180.d0/PI
    write(fileunit,'(A,F12.6)') '# mu_meV = ', mu/eV*1000.d0
    write(fileunit,'(A,F12.6)') '# eta_meV = ', eta/eV*1000.d0
    ! ----- 数值参数 -----
    write(fileunit,'(A,I0)') '# nk = ', nk
    write(fileunit,'(A,*(I0,1X))') '# ncut_list = ', ncut_list
    ! ----- 极化率所对应的求导参数 -----
    write(fileunit,'(A,A)') '# derivative_parameter = ', response_parameter
    ! ----- 物理量单位 -----
    write(fileunit,'(A)') '# OM_unit = mu_B/cell'

    if (response_parameter == 'strain_epsilon') then
        write(fileunit,'(A)') '# chi_unit = mu_B/cell per unit strain_epsilon'
    else
        write(fileunit,'(A)') '# chi_unit = mu_B/cell per radian'
    endif

    ! ----- 输出数据列 -----
    write(fileunit,'(A)') &
        '# columns: ncut Mx_OM My_OM chi_fix_x chi_fix_y '// &
        'chi_fs_x chi_fs_y chi_tot_x chi_tot_y'

    ! ----- 建立缓冲区 ----- 
    flush(fileunit)

    ! ==========
    ! 进行扫描
    ! ==========
    ! ----- 总扫描点数 ----- 
    npoints = size(ncut_list)

    ! ----- 开始计时 ----- 
    run_started = omp_get_wtime()

    ! ----- 输出初始信息 -----
    write(*,'(A,A)') 'Output file: ', trim(filename)
    write(*,'(A,A)') 'Derivative parameter: ', response_parameter
    write(*,'(A,I0,A,I0)') 'nk = ', nk, ', number of ncut points = ', npoints

    flush(output_unit)

    ! ----------
    ! 开始循环
    ! ----------
    do icut = 1, size(ncut_list)
        ! ----- 当前 cutoff -----
        ncut = ncut_list(icut)
        ! ----- 当前扫描点计时 -----
        point_started = omp_get_wtime()
        ! ----- 输出计算开始信息 -----
        write(*,'(A,I0,A,I0,A,I0,A)') &
            'point ', icut, '/', npoints, ', ncut = ', ncut, ': CHI starting'
        flush(output_unit)

        ! -----
        ! 计算极化率
        ! -----
        select case (response_parameter)

        case('strain_epsilon')
            ! ----- H_lambda = dH/d(epsilon) -----
            call CALCULATE_CHI_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), &
                                     mu, eta, prefactor, zdiag, ham_param, &
                                     HAM_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, &
                                     chi_fix_x, chi_fix_y, &
                                     chi_fs_x, chi_fs_y, &
                                     chi_tot_x, chi_tot_y)
        case ('strain_phi')
            ! ----- H_lambda = dH/d(phi) -----
            call CALCULATE_CHI_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), &
                                        mu, eta, prefactor, zdiag, ham_param, &
                                        HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, &
                                        chi_fix_x, chi_fix_y, &
                                        chi_fs_x, chi_fs_y, &
                                        chi_tot_x, chi_tot_y)
        end select

        ! ----- CHI 计算耗时 -----
        chi_elapsed = omp_get_wtime()-point_started

        write(*,'(A,I0,A,F12.2,A)') &
            'ncut = ', ncut, ': CHI done after ', chi_elapsed, &
            ' s; ket OM starting'
        flush(output_unit)

        ! ----------
        ! 直接计算轨道磁化
        ! ----------
        ! ----- ket OM 只使用 Hx = dH/dkx、Hy = dH/dky -----
        ! ----- H_lambda 不参与直接 OM 的计算 -----
        call CALCULATE_OM_CUTOFF_KET(N, nk, ncut, bs(:,1), bs(:,2), &
                                     mu, prefactor, zdiag, ham_param, &
                                     HAM_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, &
                                     Mx, My)

        ! ----------
        ! 保存结果
        ! ----------
        write(fileunit,'(I5,1X,8(ES26.16E3,1X))',iostat=ios) &
        ncut, Mx, My, &
        chi_fix_x, chi_fix_y, &
        chi_fs_x, chi_fs_y, &
        chi_tot_x, chi_tot_y

        if (ios /= 0) then
            error stop 'TBG_STRAIN_CHI_NCUT: output write failed'
        endif

        ! ----- 每完成一个 ncut，立即刷新数据文件 -----
        flush(fileunit)

        ! ==========
        ! 输出屏幕日志
        ! ==========
        ! ----- 直接轨道磁化 -----
        write(*,'(A,I0,A,2(ES16.6,1X))') &
            'ncut = ', ncut, ', OM (Mx, My) = ', Mx, My
        ! ----- FIX contribution -----
        write(*,'(A,2(ES16.6,1X))') &
            '           FIX (x, y) = ', chi_fix_x, chi_fix_y
        ! ----- FS contribution -----
        write(*,'(A,2(ES16.6,1X))') &
            '           FS  (x, y) = ', chi_fs_x, chi_fs_y
        ! ----- TOTAL susceptibility -----
        write(*,'(A,2(ES16.6,1X))') &
            '           TOT (x, y) = ', chi_tot_x, chi_tot_y
        ! ----- 当前扫描点总耗时 -----
        write(*,'(A,F12.2,A)') &
            '  saved; point wall time = ', &
            omp_get_wtime()-point_started, ' s'

        flush(output_unit)

    enddo

    ! ==========
    ! 计算结束
    ! ==========
    ! ----- 关闭输出文件 ----- 
    close(fileunit)
    ! ----- 打印输出路径 -----
    write(*,'(A,A)') 'All cutoff points saved to: ', trim(filename)
    ! ----- 打印总耗时 ----- 
    write(*,'(A,F12.2,A)') &
        'Total wall time = ', omp_get_wtime()-run_started, ' s'

    flush(output_unit)

end program TBG_STRAIN_CHI_NCUT