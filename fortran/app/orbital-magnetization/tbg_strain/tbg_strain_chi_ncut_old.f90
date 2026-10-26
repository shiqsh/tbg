program TBG_STRAIN_CHI_NCUT
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN, &
                      HAM_DERIVATIVES_STRAIN_EPSILON, HAM_DERIVATIVES_STRAIN_PHI
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF_KET
    use SUSCEPTIBILITY_CUTOFF, only : CALCULATE_CHI_CUTOFF, M_FROM_CHI_CUTOFF
    use stdlib_math, only : linspace
    use, intrinsic :: iso_fortran_env, only : output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use omp_lib, only : omp_set_max_active_levels, omp_get_wtime
    implicit none

    integer, parameter :: N = 2*NTBG
    integer :: nk, ncut_first, ncut_last, ncut_step, ncut, ncuts, icut
    integer :: npath, ilambda, iarg, ios, fileunit, pathunit
    real(8) :: theta, strain_phi, epsilon_fixed, epsilon_ref
    real(8) :: mu, eta, d, delta, prefactor, delta_h_ref
    real(8) :: bs(2,2), bs_fixed(2,2), bs_ref(2,2)
    real(8) :: ham_param(3), ham_param_fixed(3), ham_param_ref(3), zdiag(N)
    real(8) :: Mx_OM, My_OM, point_started, total_started
    real(8), allocatable :: epsilon_path(:)
    real(8), allocatable :: chi_fix(:,:), chi_fs(:,:), chi_tot(:,:), M_from_chi(:,:)
    character(len=32) :: s_theta, s_epsilon, s_phi, s_mu, s_eta, s_nk
    character(len=32) :: s_first, s_last, s_step, s_npath
    character(len=1024) :: filename, path_filename, output_path, arg, option, io_message
    logical :: check_only, exists

    ! ================================================================
    ! 参数区：固定物理点，只扫描 ncut。
    ! epsilon 是无量纲应变：5.d-4 = 0.05%，不是 0.05。
    ! ================================================================
    theta = 1.4d0*PI/180.d0
    strain_phi = 20.d0*PI/180.d0
    epsilon_fixed = 5.d-4
    mu = 0.5d-3*eV
    eta = 0.25d-3*eV
    d = 3.35d-10
    prefactor = 2.d0*m_e/hbar**2

    nk = 40
    ncut_first = 20
    ncut_last = 160
    ncut_step = 20

    ! 重构 OM 仍沿 epsilon 积分，不能对 ncut 积分 chi。
    ! 0 -> 0.05%，6 个节点：0, 0.01%, ..., 0.05%。
    ! 该间距与原程序 linspace(0, 0.002, 21) 的前 6 点相同。
    epsilon_ref = 0.d0
    npath = 6
    delta = 1.d-3              ! 原 projector 参考 OM 的相对差分步长

    ! ================================================================
    ! 可选命令行参数；不指定时使用上面的参数区。
    ! --check-only 仅显示参数，不计算、不创建结果文件。
    ! 其余选项也可用于小规模验证，不必修改源文件。
    ! ================================================================
    check_only = .false.
    output_path = ''
    iarg = 1
    do while (iarg <= command_argument_count())
        call get_command_argument(iarg, option)
        select case (trim(option))
        case ('--check-only')
            check_only = .true.
        case ('--output', '--nk', '--path-points', '--ncut-last')
            if (iarg == command_argument_count()) error stop 'Missing command-line option value'
            iarg = iarg+1
            call get_command_argument(iarg, arg)
            ios = 0
            select case (trim(option))
            case ('--output')
                output_path = trim(arg)
                if (len_trim(output_path) == 0) error stop 'Empty output path'
            case ('--nk')
                read(arg,*,iostat=ios) nk
            case ('--path-points')
                read(arg,*,iostat=ios) npath
            case ('--ncut-last')
                read(arg,*,iostat=ios) ncut_last
            end select
            if (ios /= 0) error stop 'Invalid integer option value'
        case default
            error stop 'Options: --check-only --output FILE --nk N --path-points N --ncut-last N'
        end select
        iarg = iarg+1
    enddo

    if (N < 2 .or. mod(N,2) /= 0) error stop 'Invalid Hamiltonian dimension'
    if (nk < 1 .or. npath < 1) error stop 'nk and path-points must be positive'
    if (eta <= 0.d0 .or. delta <= 0.d0) error stop 'eta and delta must be positive'
    if (ncut_step <= 0) error stop 'ncut_step must be positive'
    if (ncut_first < 0 .or. ncut_last < ncut_first .or. ncut_last > N/2-1) &
        error stop 'Invalid ncut range: require 0 <= first <= last <= N/2-1'
    if (mod(ncut_last-ncut_first,ncut_step) /= 0) error stop 'ncut_last must lie on the ncut grid'
    if (epsilon_fixed == epsilon_ref) then
        npath = 1
    else if (npath < 2) then
        error stop 'A nonzero epsilon integration interval needs at least two path points'
    endif
    ncuts = (ncut_last-ncut_first)/ncut_step+1
    call omp_set_max_active_levels(1)

    allocate(epsilon_path(npath))
    allocate(chi_fix(2,npath), chi_fs(2,npath), chi_tot(2,npath), M_from_chi(2,npath))
    if (npath == 1) then
        epsilon_path(1) = epsilon_fixed
    else
        epsilon_path = linspace(epsilon_ref,epsilon_fixed,npath)
        epsilon_path(1) = epsilon_ref
        epsilon_path(npath) = epsilon_fixed
    endif
    ham_param_fixed = [theta,epsilon_fixed,strain_phi]
    ham_param_ref = [theta,epsilon_ref,strain_phi]
    call GEOMETRY_STRAIN(theta,epsilon_fixed,strain_phi,bs_fixed)
    call GEOMETRY_STRAIN(theta,epsilon_ref,strain_phi,bs_ref)
    delta_h_ref = delta*norm2(bs_ref(:,1))
    call Z_TBG(d,zdiag)

    write(*,'(A,F8.3,A,F10.5)') 'Fixed phi (deg) = ',strain_phi*180.d0/PI, &
                               ', epsilon (%) = ',epsilon_fixed*100.d0
    write(*,'(A,I0,A,I0,A,I0)') 'N = ',N,', nk = ',nk,', maximum allowed ncut = ',N/2-1
    write(*,'(A,I0,A,I0,A,I0,A,I0)') 'ncut: ',ncut_first,' to ',ncut_last,' step ',ncut_step, &
                                    '; number of cutoffs = ',ncuts
    write(*,'(A,I0,A)') 'Each cutoff uses ',npath,' epsilon path points; chi is dM/d(epsilon).'
    write(*,'(A)') 'Direct OM: KET. Reconstruction reference: original PROJECTOR finite difference.'
    flush(output_unit)
    if (check_only) then
        write(*,'(A)') 'Parameters checked; no numerical integration or output files created.'
        stop
    endif

    ! ================================================================
    ! 每个 ncut 一行最终数据；另存积分路径上的 chi，避免只存在内存中。
    ! status=new：绝不覆盖已有结果。需要重跑时用 --output 指定新文件。
    ! ================================================================
    write(s_theta,'(F9.3)') theta*180.d0/PI
    write(s_epsilon,'(F9.4)') epsilon_fixed*100.d0
    write(s_phi,'(F9.3)') strain_phi*180.d0/PI
    write(s_mu,'(F9.3)') mu/eV*1000.d0
    write(s_eta,'(F9.3)') eta/eV*1000.d0
    write(s_nk,'(I0)') nk
    write(s_first,'(I0)') ncut_first
    write(s_last,'(I0)') ncut_last
    write(s_step,'(I0)') ncut_step
    write(s_npath,'(I0)') npath
    filename = '/Users/shiqsh/Workspace/tbg/data/tbg/fortran/susceptibility/'// &
        'tbg_chi_ncut_theta'//trim(adjustl(s_theta))//'_epsilon_'//trim(adjustl(s_epsilon))// &
        '_phi_'//trim(adjustl(s_phi))//'_mu_'//trim(adjustl(s_mu))//'_eta_'//trim(adjustl(s_eta))// &
        '_nk_'//trim(s_nk)//'_ncut_'//trim(s_first)//'to'//trim(s_last)//'step'//trim(s_step)// &
        '_npath_'//trim(s_npath)//'.dat'
    if (len_trim(output_path) > 0) filename = trim(output_path)
    if (len_trim(filename)+len('.path.dat') > len(path_filename)) error stop 'Output path too long'
    path_filename = trim(filename)//'.path.dat'
    inquire(file=trim(filename),exist=exists)
    if (exists) error stop 'Output already exists; choose a NEW --output filename'
    inquire(file=trim(path_filename),exist=exists)
    if (exists) error stop 'Path output already exists; choose a NEW --output filename'
    open(newunit=fileunit,file=trim(filename),status='new',action='write',iostat=ios,iomsg=io_message)
    call check_io(ios,io_message)
    open(newunit=pathunit,file=trim(path_filename),status='new',action='write',iostat=ios,iomsg=io_message)
    call check_io(ios,io_message)
    call write_header(fileunit)
    call write_header(pathunit)
    write(fileunit,'(A)') '# columns: ncut epsilon_percent phi_deg Mx_OM My_OM Mx_from_chi My_from_chi '// &
        'chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y delta_Mx delta_My'
    write(pathunit,'(A)') '# columns: ncut epsilon_percent phi_deg '// &
        'chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y'
    flush(fileunit)
    flush(pathunit)
    write(*,'(A,A)') 'Summary output: ',trim(filename)
    write(*,'(A,A)') 'Integration-path output: ',trim(path_filename)
    total_started = omp_get_wtime()

    ! ================================================================
    ! 外层：扫描 ncut。内层：仅为重构 OM 计算辅助 epsilon 积分路径。
    ! ================================================================
    do ncut = ncut_first,ncut_last,ncut_step
        icut = (ncut-ncut_first)/ncut_step+1
        point_started = omp_get_wtime()
        write(*,'(A,I0,A,I0,A,I0)') 'Cutoff ',icut,'/',ncuts,': ncut = ',ncut
        flush(output_unit)

        do ilambda = 1,npath
            ham_param = [theta,epsilon_path(ilambda),strain_phi]
            ! epsilon 变化时，必须更新应变后的倒格矢，不能一直使用终点 bs_fixed。
            call GEOMETRY_STRAIN(theta,epsilon_path(ilambda),strain_phi,bs)
            call CALCULATE_CHI_CUTOFF(N,nk,ncut,bs(:,1),bs(:,2),mu,eta,prefactor, &
                zdiag,ham_param,HAM_STRAIN,HAM_DERIVATIVES_STRAIN_EPSILON, &
                chi_fix(1,ilambda),chi_fix(2,ilambda),chi_fs(1,ilambda),chi_fs(2,ilambda), &
                chi_tot(1,ilambda),chi_tot(2,ilambda))
            if (.not.all(ieee_is_finite([chi_fix(:,ilambda),chi_fs(:,ilambda),chi_tot(:,ilambda)]))) &
                error stop 'Nonfinite chi; this path point was not saved'
            write(pathunit,'(I6,1X,8(ES26.17E3,1X))',iostat=ios,iomsg=io_message) &
                ncut,epsilon_path(ilambda)*100.d0,strain_phi*180.d0/PI, &
                chi_fix(:,ilambda),chi_fs(:,ilambda),chi_tot(:,ilambda)
            call check_io(ios,io_message)
            flush(pathunit,iostat=ios,iomsg=io_message)
            call check_io(ios,io_message)
            write(*,'(A,I0,A,I0,A,F10.5,A)') '  chi path ',ilambda,'/',npath, &
                ', epsilon (%) = ',epsilon_path(ilambda)*100.d0,' saved'
            flush(output_unit)
        enddo

        ! 仅在固定物理点计算 direct OM，不在每个路径点重复计算它。
        call CALCULATE_OM_CUTOFF_KET(N,nk,ncut,bs_fixed(:,1),bs_fixed(:,2),mu,prefactor, &
            zdiag,ham_param_fixed,HAM_STRAIN,HAM_DERIVATIVES_STRAIN_PHI,Mx_OM,My_OM)

        ! 每个 ncut 都重新计算同一参考物理点的 M_ref；不能跨 ncut 共用 M_ref。
        ! 保持附件的原组合：reference OM 用 projector，direct OM 用 ket。
        call M_FROM_CHI_CUTOFF(N,nk,ncut,bs_ref(:,1),bs_ref(:,2),delta_h_ref,mu,prefactor, &
            zdiag,ham_param_ref,HAM_STRAIN,epsilon_path,chi_tot(1,:),chi_tot(2,:), &
            M_from_chi(1,:),M_from_chi(2,:))
        if (.not.all(ieee_is_finite([Mx_OM,My_OM,M_from_chi(:,npath)]))) &
            error stop 'Nonfinite OM; this cutoff was not saved'

        ! 只输出固定 epsilon/phi 的结果：每个 ncut 一行，共 8 行。
        write(fileunit,'(I6,1X,14(ES26.17E3,1X))',iostat=ios,iomsg=io_message) &
            ncut,epsilon_fixed*100.d0,strain_phi*180.d0/PI,Mx_OM,My_OM,M_from_chi(:,npath), &
            chi_fix(:,npath),chi_fs(:,npath),chi_tot(:,npath), &
            M_from_chi(1,npath)-Mx_OM,M_from_chi(2,npath)-My_OM
        call check_io(ios,io_message)
        flush(fileunit,iostat=ios,iomsg=io_message)
        call check_io(ios,io_message)
        write(*,'(A,I0,2(A,ES18.10),A,F10.2,A)') '  ncut=',ncut,' Mx=',Mx_OM,' My=',My_OM, &
            '; saved after ',omp_get_wtime()-point_started,' s'
        flush(output_unit)
    enddo

    close(pathunit,iostat=ios,iomsg=io_message)
    call check_io(ios,io_message)
    close(fileunit,iostat=ios,iomsg=io_message)
    call check_io(ios,io_message)
    write(*,'(A,F12.2,A)') 'All cutoffs completed in ',omp_get_wtime()-total_started,' s'
    write(*,'(A,A)') 'Results saved to: ',trim(filename)

contains

    subroutine check_io(status,message)
        integer, intent(in) :: status
        character(*), intent(in) :: message
        if (status /= 0) then
            write(*,'(A)') trim(message)
            error stop 'NCUT scan: file I/O failed'
        endif
    end subroutine check_io

    subroutine write_header(unit)
        integer, intent(in) :: unit
        write(unit,'(A)') '# TBG fixed-point ncut convergence scan'
        write(unit,'(A,F12.6)') '# theta_deg = ',theta*180.d0/PI
        write(unit,'(A,F12.6)') '# fixed_phi_deg = ',strain_phi*180.d0/PI
        write(unit,'(A,ES26.17E3)') '# fixed_epsilon = ',epsilon_fixed
        write(unit,'(A,ES26.17E3)') '# reference_epsilon = ',epsilon_ref
        write(unit,'(A,ES26.17E3)') '# interlayer_d_m = ',d
        write(unit,'(A,ES26.17E3)') '# prefactor = ',prefactor
        write(unit,'(A,F12.6,A,F12.6)') '# mu_meV = ',mu/eV*1000.d0,', eta_meV = ',eta/eV*1000.d0
        write(unit,'(A,I0,A,I0,A,I0)') '# N = ',N,', nk = ',nk,', epsilon_path_points = ',npath
        write(unit,'(A,3(I0,1X))') '# ncut_first last step = ',ncut_first,ncut_last,ncut_step
        write(unit,'(A,ES26.17E3)') '# reference_OM_delta = ',delta
        write(unit,'(A)') '# direct_OM = KET; reference_OM = PROJECTOR; integration = epsilon trapezoid'
        write(unit,'(A)') '# OM_unit = mu_B/cell; chi_unit = mu_B/cell per unit strain_epsilon'
        write(unit,'(A)') '# delta_M = M_from_chi - M_OM; epsilon_percent columns use epsilon*100'
    end subroutine write_header

end program TBG_STRAIN_CHI_NCUT
