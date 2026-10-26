program TBG_STRAIN_CHI
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON,HAM_DERIVATIVES_STRAIN_PHI
    use OM_CUTOFF, only :  CALCULATE_OM_CUTOFF_KET
    use SUSCEPTIBILITY_CUTOFF, only : CALCULATE_CHI_CUTOFF, M_FROM_CHI_CUTOFF, M_FROM_CHI_CUTOFF_KET
    use stdlib_math, only : linspace
    use CHI_CHECKPOINT, only : LOAD_CHI_CHECKPOINT, CREATE_CHI_CHECKPOINT, APPEND_CHI_POINT
    use, intrinsic :: iso_fortran_env, only : output_unit
    use omp_lib, only : omp_set_max_active_levels, omp_get_wtime
    implicit none

    integer, parameter :: N = 2*NTBG

    integer :: nk, ncut
    integer :: iepsilon, iphi
    integer :: ilambda, nlambda
    integer :: ie_start, ie_end
    integer :: ip_start, ip_end
    integer :: fileunit

    integer :: checkpoint_unit, completed, iarg, ios, stop_after
    real(8), allocatable :: saved_results(:,:)
    real(8) :: point_started
    character(len=1024) :: resume_path, output_path, arg
    logical :: check_resume

    real(8) :: theta
    real(8) :: strain_epsilon, strain_phi
    real(8) :: mu, d
    real(8) :: eta, prefactor
    real(8) :: delta, delta_h, delta_h_ref
    real(8) :: bs(2,2), bs_ref(2,2)
    real(8) :: ham_param(3), ham_param_ref(3)
    real(8) :: zdiag(N)

    real(8), allocatable :: epsilon_list(:), phi_list(:) 
    real(8), allocatable :: chi_fix_x(:), chi_fix_y(:)
    real(8), allocatable :: chi_fs_x(:), chi_fs_y(:)
    real(8), allocatable :: chi_tot_x(:), chi_tot_y(:)
    real(8), allocatable :: Mx_OM(:), My_OM(:) ! OM from cutoff
    real(8), allocatable :: Mx_from_chi(:), My_from_chi(:)

    character(len=24) :: scan_parameter ! 'strain_epsilon' or 'strain_phi'
    character(len = 64) :: chi_unit
    character(len = 512) :: filename    
    character(len = 24) :: s_theta, s_epsilon, s_phi    
    character(len = 24) :: s_mu, s_eta, s_nk, s_ncut

    ! =====
    ! 参数赋值
    ! =====
    ! ----------
    ! 固定参数
    ! ----------
    d = 3.35d-10
    theta = (1.4d0)*PI/180.d0
    prefactor = 2.d0*m_e/hbar**2

    ! ----------
    ! 可变参数
    ! ----------
    ! ----- 化学势与高斯展宽 ----- 
    mu = (0.5d0)*1.d-3*eV
    eta = (0.25d0)*1.d-3*eV  

    ! ----- 截断与mesh ----- 
    delta = 1.d-4 ! 中心差分步长
    ncut = 120
    nk = 40! 这里需要用奇数 for mu = 0, but odd nk does not satisfy the C3 symmetry, 

    ! BLAS thread limiting and idle-sleep protection are supplied by run_tbg_chi_safe.sh.
    ! Keep nested OpenMP regions serial; the outer k mesh remains parallel.
    call omp_set_max_active_levels(1)
    resume_path = ''
    output_path = ''
    stop_after = huge(0)
    check_resume = .false.
    iarg = 1
    do while (iarg <= command_argument_count())
        call get_command_argument(iarg, arg)
        select case (trim(arg))
        case ('--check-resume')
            check_resume = .true.
        case ('--resume', '--output', '--stop-after')
            if (iarg == command_argument_count()) error stop 'Missing command-line option value'
            iarg = iarg + 1
            select case (trim(arg))
            case ('--resume')
                call get_command_argument(iarg, resume_path)
            case ('--output')
                call get_command_argument(iarg, output_path)
                if (len_trim(output_path) > len(filename)) error stop 'Output path too long'
            case ('--stop-after')
                call get_command_argument(iarg, arg)
                read(arg,*,iostat=ios) stop_after
                if (ios /= 0) error stop 'Invalid --stop-after value'
                if (stop_after < 0) error stop '--stop-after must be nonnegative'
            end select
        case default
            error stop 'Options: --resume FILE --output FILE --check-resume --stop-after NEW_POINTS'
        end select
        iarg = iarg + 1
    enddo

    ! ----------
    ! 维度检查
    ! ----------
    if (N < 2 .or. mod(N,2) /= 0) error stop 'TBG_STRAIN_CHI: invalid N'
    if (ncut < 0 .or. ncut > N/2-1) error stop 'TBG_STRAIN_CHI: invalid ncut'
    if (nk < 1) error stop 'TBG_STRAIN_CHI: invalid nk'
    if (eta <= 0.d0 .or. delta <= 0.d0) error stop 'TBG_STRAIN_CHI: invalid eta or delta'

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

    ! ----------
    ! 固定 epsilon 0.1%, 扫描phi
    ! ----------
    !   ie_start = 11
    !   ie_end = 11
    !   ip_start = 1
    !   ip_end = size(phi_list)

    ! ----- 检查索引 ----- 
    if (ie_start < 1 .or. ie_end > size(epsilon_list) .or. ie_start > ie_end) then
        error stop 'TBG_STRAIN_CHI: invalid epsilon index range'
    endif
    if (ip_start < 1 .or. ip_end > size(phi_list) .or. ip_start > ip_end) then
        error stop 'TBG_STRAIN_CHI: invalid phi index range'
    endif

    ! ==========
    ! 应力参数扫描检查
    ! ==========
    ! -----
    ! 扫描检查
    ! 有且只能有一个参数进行扫描，两个参数扫描即报错
    ! -----
    if ((ie_start < ie_end) .eqv. (ip_start < ip_end)) then
        error stop 'TBG_STRAIN_CHI: scan exactly one of strain_epsilon and strain_phi'
    endif
    ! -----
    ! 确定扫描参数到底是strain_epsilon 还是 strain_phi
    ! 通过比较索引的大小即可进行比较
    ! -----
    if (ie_start < ie_end) then
        scan_parameter = 'strain_epsilon'
        nlambda = ie_end-ie_start+1
        chi_unit = 'mu_B/cell per unit strain_epsilon'
    else
        scan_parameter = 'strain_phi'
        nlambda = ip_end-ip_start+1
        chi_unit = 'mu_B/cell per radian'
    endif

    ! =====
    ! 分配数组
    ! =====
    allocate(chi_fix_x(nlambda), chi_fix_y(nlambda))
    allocate(chi_fs_x(nlambda), chi_fs_y(nlambda))
    allocate(chi_tot_x(nlambda), chi_tot_y(nlambda))
    allocate(Mx_OM(nlambda), My_OM(nlambda))
    allocate(Mx_from_chi(nlambda), My_from_chi(nlambda))
    allocate(saved_results(8,nlambda))
    saved_results = 0.d0
    completed = 0
    if (len_trim(resume_path) > 0) then
        call LOAD_CHI_CHECKPOINT(trim(resume_path), scan_parameter, &
            [N,nk,ncut,nlambda,ie_start,ie_end,ip_start,ip_end,size(epsilon_list),size(phi_list)], &
            [theta,d,mu,eta,prefactor,delta], epsilon_list, phi_list, saved_results, completed)
        Mx_OM(:completed) = saved_results(1,:completed)
        My_OM(:completed) = saved_results(2,:completed)
        chi_fix_x(:completed) = saved_results(3,:completed)
        chi_fix_y(:completed) = saved_results(4,:completed)
        chi_fs_x(:completed) = saved_results(5,:completed)
        chi_fs_y(:completed) = saved_results(6,:completed)
        chi_tot_x(:completed) = saved_results(7,:completed)
        chi_tot_y(:completed) = saved_results(8,:completed)
        write(*,'(A,I0,A,I0,A)') 'Recovered ', completed, '/', nlambda, ' completed points.'
    endif
    if (check_resume) then
        if (len_trim(resume_path) == 0) error stop '--check-resume requires --resume FILE'
        write(*,'(A)') 'Checkpoint parameters and results validated; no calculation or files written.'
        stop
    endif

    ! Initialize the reference even when the scan resumes after point 1.
    ham_param_ref = [theta, epsilon_list(ie_start), phi_list(ip_start)]
    call GEOMETRY_STRAIN(theta, epsilon_list(ie_start), phi_list(ip_start), bs_ref)
    delta_h_ref = delta*norm2(bs_ref(:,1))

    ! ----------
    ! 层间位置算符
    ! ----------
    call Z_TBG(d, zdiag)

    ! ==========
    ! 输出文件名
    ! ==========
    ! ----------
    ! 字符串参数
    ! ----------
    write(s_theta, '(F8.3)') theta*180.d0/PI
    write(s_mu, '(F8.3)') mu/eV*1000.d0 ! meV输出
    write(s_eta, '(F8.3)') eta/eV*1000.d0 !meV输出
    write(s_nk, '(I0)') nk
    write(s_ncut, '(I0)') ncut

    ! ----- strain epsilon label -----
    if (ie_start == ie_end) then
        write(s_epsilon, '(F8.4)') epsilon_list(ie_start)*100.d0
    else
        write(s_epsilon, '(A,I0,A,I0)') 'scan', ie_start, 'to', ie_end
    endif

    ! ----- strain phi label ----- 
    if (ip_start == ip_end) then
        write(s_phi, '(F8.3)') phi_list(ip_start)*180.d0/PI
    else
        write(s_phi, '(A,I0,A,I0)') 'scan_', ip_start, 'to', ip_end
    endif

    ! ----------
    !  构造输出路径 
    ! 请使用绝对路径，这样方便建立logo文件
    ! ----------
    filename = '/Users/shiqsh/Workspace/tbg/data/tbg/fortran/susceptibility/'// &
               'tbg_chi_'//trim(scan_parameter)// &
               '_theta'//trim(adjustl(s_theta))// &
               '_epsilon_'//trim(adjustl(s_epsilon))// &
               '_phi_'//trim(adjustl(s_phi))// &
               '_mu_'//trim(adjustl(s_mu))// &
               '_eta_'//trim(adjustl(s_eta))// &
               '_nk_'//trim(adjustl(s_nk))// &
               '_ncut_'//trim(adjustl(s_ncut))//'.dat'

    if (len_trim(output_path) > 0) filename = trim(output_path)
    ! Never overwrite the old output or the recovered checkpoint.
    open(newunit=fileunit, file=trim(filename), status='new', action='write', iostat=ios)
    if (ios /= 0) error stop 'Cannot create output; use --output with a NEW filename'
    call CREATE_CHI_CHECKPOINT(trim(filename)//'.checkpoint', scan_parameter, &
        [N,nk,ncut,nlambda,ie_start,ie_end,ip_start,ip_end,size(epsilon_list),size(phi_list)], &
        [theta,d,mu,eta,prefactor,delta], epsilon_list, phi_list, checkpoint_unit)
    do ilambda = 1, completed
        call APPEND_CHI_POINT(checkpoint_unit, ilambda, saved_results(:,ilambda))
    enddo
    write(*,'(A,A)') 'Per-point checkpoint: ', trim(filename)//'.checkpoint'

    ! ----- output the header ----- 
    ! ----- Physical parameters -----
    write(fileunit,'(A,F10.4)') '# theta_deg = ', theta*180.d0/PI
    write(fileunit,'(A,F10.4)') '# mu_meV = ', mu/eV*1000.d0
    write(fileunit,'(A,F10.4)') '# eta_meV = ', eta/eV*1000.d0
    ! ----- Strain parameters -----
    write(fileunit,'(A,2F12.6)') '# epsilon_percent (start, end) = ', &
        epsilon_list(ie_start)*100.d0, epsilon_list(ie_end)*100.d0
    write(fileunit,'(A,2F12.6)') '# phi_deg (start, end) = ', &
        phi_list(ip_start)*180.d0/PI, phi_list(ip_end)*180.d0/PI
    write(fileunit,'(A,A)') '# derivative_parameter = ', trim(scan_parameter)
    ! ----- Numerical parameters -----
    write(fileunit,'(A,I0,A,I0)') '# nk = ', nk, ', ncut = ', ncut
    write(fileunit,'(A,ES12.4)') '# OM_delta = ', delta
    ! ----- Units -----
    write(fileunit,'(A)') '# OM_unit = mu_B/cell'
    write(fileunit,'(A,A)') '# chi_unit = ', trim(chi_unit)
    ! ----- Data columns -----
    write(fileunit,'(A)') &
        '# columns: epsilon_percent phi_deg Mx_OM My_OM Mx_from_chi My_from_chi '// &
        'chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y delta_Mx delta_My'
    flush(fileunit)
    

    ! =====
    ! 计算极化率对应的轨道磁矩与pristine轨道磁矩
    ! =====
    do ilambda = completed+1, nlambda
        if (ilambda-completed > stop_after) then
            close(checkpoint_unit)
            close(fileunit)
            write(*,'(A)') 'Stopped after requested new points; resume from the checkpoint with a new output name.'
            stop
        endif
        point_started = omp_get_wtime()
        write(*,'(A,I0,A,I0,A)') 'point ', ilambda, '/', nlambda, ': CHI starting'
        flush(output_unit)
        ! ----- 确定扫描索引 ----- 
        if (scan_parameter == 'strain_epsilon') then
            iepsilon = ie_start+ilambda-1
            iphi = ip_start
        else
            iepsilon = ie_start
            iphi = ip_start+ilambda-1
        endif
        ! ----- 当前哈密顿量与几何 ----- 
        strain_epsilon = epsilon_list(iepsilon)
        strain_phi = phi_list(iphi)
        ham_param = [theta, strain_epsilon, strain_phi]
        call GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs)
        delta_h = delta*norm2(bs(:,1))

        ! ----- 计算轨道磁化 ----- 
        if (scan_parameter == 'strain_epsilon') then
            call CALCULATE_CHI_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), mu, eta, prefactor, &
                                      zdiag, ham_param, HAM_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, &
                                      chi_fix_x(ilambda), chi_fix_y(ilambda), &
                                      chi_fs_x(ilambda), chi_fs_y(ilambda), &
                                      chi_tot_x(ilambda), chi_tot_y(ilambda))
        else
            call CALCULATE_CHI_CUTOFF(N, nk, ncut, bs(:,1), bs(:,2), mu, eta, prefactor, &
                                      zdiag, ham_param, HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, &
                                      chi_fix_x(ilambda), chi_fix_y(ilambda), &
                                      chi_fs_x(ilambda), chi_fs_y(ilambda), &
                                      chi_tot_x(ilambda), chi_tot_y(ilambda))
        endif

        ! ----- 直接计算轨道磁矩 ----- 
        write(*,'(A,I0,A,F12.2,A)') 'point ', ilambda, ': CHI done after ', &
            omp_get_wtime()-point_started, ' s; OM starting'
        flush(output_unit)
        call CALCULATE_OM_CUTOFF_KET(N, nk, ncut, bs(:,1), bs(:,2), mu, prefactor, zdiag, ham_param, &
                                    HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, Mx_OM(ilambda), My_OM(ilambda))

        ! Flush all eight computed quantities before printing "point completed".
        call APPEND_CHI_POINT(checkpoint_unit, ilambda, &
            [Mx_OM(ilambda),My_OM(ilambda),chi_fix_x(ilambda),chi_fix_y(ilambda), &
             chi_fs_x(ilambda),chi_fs_y(ilambda),chi_tot_x(ilambda),chi_tot_y(ilambda)])

        ! -----
        ! 在屏幕上进行打印，确认算到哪个参数了  
        ! -----
        write(*,'(A,I0,A,I0,A,F10.5,A,F10.3,2(A, ES32.16))') &
            'point ', ilambda, '/', nlambda, &
            ' epsilon(%) = ', strain_epsilon*100.d0, &
            ' phi(deg) = ', strain_phi*180.d0/PI, &
            ' Mx = ', Mx_OM(ilambda), &
            ' My = ', My_OM(ilambda)
        write(*,'(A,F12.2,A)') '  saved; point wall time = ', omp_get_wtime()-point_started, ' s'
        flush(output_unit)

    enddo
    close(checkpoint_unit)
    write(*,'(A)') 'All scan points saved. Calculating reference OM and reconstructing M from chi...'
    flush(output_unit)

    ! ==========
    ! 计算重构的轨道磁矩
    ! ==========
    if (scan_parameter == 'strain_epsilon') then
        call M_FROM_CHI_CUTOFF_KET(N, nk, ncut, bs_ref(:,1), bs_ref(:,2), mu, prefactor, & 
                               zdiag, ham_param_ref, HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, &
                               epsilon_list(ie_start:ie_end), chi_tot_x, chi_tot_y, &
                               Mx_from_chi, My_from_chi)
    else
        call M_FROM_CHI_CUTOFF_KET(N, nk, ncut, bs_ref(:,1), bs_ref(:,2), mu, prefactor, &
                               zdiag, ham_param_ref, HAM_STRAIN,  HAM_DERIVATIVES_STRAIN_PHI, &
                               phi_list(ip_start:ip_end), chi_tot_x, chi_tot_y, &
                               Mx_from_chi, My_from_chi)
    endif
    ! if (scan_parameter == 'strain_epsilon') then
    !     call M_FROM_CHI_CUTOFF(N, nk, ncut, bs_ref(:,1), bs_ref(:,2), delta_h_ref, mu, prefactor, & 
    !                            zdiag, ham_param_ref, HAM_STRAIN, &
    !                            epsilon_list(ie_start:ie_end), chi_tot_x, chi_tot_y, &
    !                            Mx_from_chi, My_from_chi)
    ! else
    !     call M_FROM_CHI_CUTOFF(N, nk, ncut, bs_ref(:,1), bs_ref(:,2), delta_h_ref, mu, prefactor, &
    !                            zdiag, ham_param_ref, HAM_STRAIN, &
    !                            phi_list(ip_start:ip_end), chi_tot_x, chi_tot_y, &
    !                            Mx_from_chi, My_from_chi)
    ! endif


    ! ==========
    ! 输出结果
    ! ==========
    do ilambda = 1, nlambda

        ! ----- 确定参数与索引 ----- 
        if (scan_parameter == 'strain_epsilon') then
            iepsilon = ie_start+ilambda-1
            iphi = ip_start
        else
            iepsilon = ie_start
            iphi = ip_start+ilambda-1
        endif

        write(fileunit,'(14(ES32.16E3,1X))') &
            epsilon_list(iepsilon)*100.d0, &
            phi_list(iphi)*180.d0/PI, &
            Mx_OM(ilambda), My_OM(ilambda), &
            Mx_from_chi(ilambda), My_from_chi(ilambda), &
            chi_fix_x(ilambda), chi_fix_y(ilambda), &
            chi_fs_x(ilambda), chi_fs_y(ilambda), &
            chi_tot_x(ilambda), chi_tot_y(ilambda), &
            Mx_from_chi(ilambda)-Mx_OM(ilambda), &
            My_from_chi(ilambda)-My_OM(ilambda)

    enddo

    close(fileunit)
    write(*,'(A,A)') 'Results saved to: ', trim(filename)

end program TBG_STRAIN_CHI

        










