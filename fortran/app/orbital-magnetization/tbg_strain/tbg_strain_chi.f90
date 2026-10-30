program TBG_STRAIN_CHI
    use CONSTANTS, only : PI, m_e, hbar, NTBG, eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN, GEOMETRY_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, HAM_DERIVATIVES_STRAIN_PHI
    use OM_CUTOFF, only :  CALCULATE_OM_CUTOFF_KET
    use SUSCEPTIBILITY_CUTOFF, only : CALCULATE_CHI_CUTOFF, M_FROM_CHI_CUTOFF, M_FROM_CHI_CUTOFF_KET
    use CHI_CHECKPOINT, only : LOAD_CHI_CHECKPOINT, CREATE_CHI_CHECKPOINT, APPEND_CHI_POINT
    use stdlib_math, only : linspace
    use, intrinsic :: iso_fortran_env, only : output_unit
    use omp_lib, only : omp_set_max_active_levels, omp_get_wtime
    implicit none

    integer, parameter :: N = 2*NTBG

    integer :: nk, ncut
    integer :: iepsilon, iphi, ilambda, nlambda
    integer :: ie_start, ie_end, ip_start, ip_end
    integer :: fileunit, checkpoint_unit
    integer :: completed, iarg, ios, stop_after

    real(8), allocatable :: saved_results(:,:)
    real(8) :: point_started

    character(len=1024) :: resume_path, output_path, arg, filename
    character(len = 1024) :: io_message
    character(len = 24) :: scan_parameter, checkpoint_mode
    character(len = 64) :: chi_unit
    character(len = 24) :: s_theta, s_epsilon, s_phi
    character(len = 24) :: s_mu, s_eta, s_nk, s_ncut
    character(len = 8) :: date_stamp
    character(len = 10) :: time_stamp

    logical :: check_resume, file_exists, checkpoint_exists
    
    ! ==========
    ! 哈密顿量参数与物理量
    ! ==========
    real(8) :: theta, strain_epsilon, strain_phi
    real(8) :: mu ,d ,eta, prefactor, delta
    real(8) :: bs(2,2), bs_ref(2, 2)
    real(8) :: ham_param(3), ham_param_ref(3)
    real(8) :: zdiag(N)

    real(8) :: chi0_fix_x, chi0_fix_y
    real(8) :: chi0_fs_x, chi0_fs_y
    real(8) :: chi0_tot_x, chi0_tot_y


    ! =====
    ! 输出物理量
    ! =====
    ! ----- 重构残差 -----
    real(8) :: residual_x, residual_y
    ! ----- 扫描数组 -----
    real(8), allocatable :: epsilon_list(:), phi_list(:)
    ! ----- susceptibility 数组 -----
    real(8), allocatable :: chi_fix_x(:), chi_fix_y(:)
    real(8), allocatable :: chi_fs_x(:), chi_fs_y(:)
    real(8), allocatable :: chi_tot_x(:), chi_tot_y(:)
    ! ----- OM 数组 -----
    real(8), allocatable :: Mx_OM(:), My_OM(:)
    real(8), allocatable :: Mx_from_chi(:), My_from_chi(:)


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

    ! ----- 禁止并行嵌套; 外层k 网格仍可并行 ----- 
    call omp_set_max_active_levels(1)

    ! ==========
    ! 命令行参数
    ! 支持以下四个选项 
    ! --resume FILE
    !     从指定的 checkpoint 文件恢复已完成的计算。
    ! --output FILE
    !     指定本次计算的数据输出文件路径。
    ! --check-resume
    !     只检查 checkpoint 是否有效，不进行新的计算。
    ! --stop-after NEW_POINTS
    !     限制本次最多新计算多少个参数点。
    !
    ! 如果运行程序时不输入任何附加参数，则使用下方的
    ! 默认设置，从头计算全部参数点
    ! ==========
    ! ----------
    ! .1 初始化命令行选项的默认值
    ! ----------
    ! 默认不读取已有的 checkpoint 文件。
    ! 空字符串表示用户没有指定恢复路径。
    resume_path = ''
    ! 默认不使用外部指定的输出文件名。
    ! 后续程序将按照 theta、epsilon、phi、mu 等参数
    ! 自动构造输出文件名。
    output_path = ''
    ! huge(0) 返回默认整数类型能够表示的最大正整数。
    ! 将 stop_after 初始化为最大整数，相当于默认
    ! 不限制本次新增计算的参数点数量。
    stop_after = huge(0)
    ! 默认不进入 checkpoint 检查模式。
    ! 因此程序将正常执行后续物理计算。
    check_resume = .false.

    ! ----------
    ! .2 逐个读取终端传入的命令行参数
    ! ----------
    ! iarg 表示当前正在处理第几个命令行参数。
    ! Fortran 的命令行参数从第 1 个开始编号。
    iarg = 1

    ! command_argument_count() 是 Fortran 内置函数，
    ! 返回用户在终端传入的命令行参数数量。
    !
    ! 例如：
    !     ./tbg_strain_chi --stop-after 3
    !
    ! 一共有两个参数：
    !     第 1 个参数：--stop-after
    !     第 2 个参数：3
    !
    ! 因此 command_argument_count() 返回 2。
    !
    ! 当 iarg 超过参数数量时，结束读取循环。
    do while (iarg <= command_argument_count())

        ! get_command_argument 是 Fortran 内置子程序。
        ! 将第 iarg 个命令行参数读取到字符串 arg 中。
        !
        ! 例如：
        !     iarg = 1
        !     arg = '--stop-after'
        call get_command_argument(iarg, arg)
        ! ----------
        ! .3 根据读取到的参数名称，选择对应的操作
        ! ----------
        ! trim(arg) 删除字符串末尾的空格。
        ! select case 根据 arg 的字符串内容，
        ! 执行对应的 case 分支。
        select case (trim(arg))
            ! -----
            ! 选项A： --check-resume
            ! -----
            case ('--check-resume')
                ! 如果用户输入 --check-resume，
                ! 就将检查模式设置为 .true.。
                ! 主程序后面会根据 check_resume 的值，
                ! 决定是否仅验证 checkpoint 然后停止。
                ! 这个选项不需要额外的参数值。
                check_resume = .true.
            ! -----
            ! 选项B ： 需要额外参数值的三个选项
            ! -----
            case ('--resume', '--output', '--stop-after')
                ! 这三个选项有一个共同特点：
                ! 选项名称后面必须再提供一个值。
                ! 例如：
                !     --resume old.dat.checkpoint
                !     --output result.dat
                !     --stop-after 3
                ! 如果当前 iarg 已经等于参数总数，
                ! 说明当前选项是最后一个参数，
                ! 它后面没有提供所需要的值。
                ! 例如只输入 --stop-after，
                ! 却没有输入需要计算的点数。
                if (iarg == command_argument_count()) then
                    error stop 'TBG_STRAIN_CHI: missing command-line option value'
                endif

                ! 将参数索引增加 1。
                ! 例如：
                !     第 1 个参数是 --stop-after
                !     第 2 个参数是 3
                ! 现在将 iarg 从 1 增加到 2，
                ! 以便读取后面的参数值。
                iarg = iarg + 1
                ! 根据刚才读取的选项名称，
                ! 判断应该如何处理紧随其后的参数值。
                ! 注意：此时 arg 仍保存选项名称，
                ! 例如 '--stop-after'。
                select case (trim(arg))
                    ! -----
                    ! B1. 指定 checkpoint 恢复路径
                    ! -----
                    case ('--resume')
                        ! 读取第 iarg 个参数，
                        ! 并将其作为 checkpoint 文件路径。
                        ! 例如：
                        !     --resume old.dat.checkpoint
                        ! 得到：
                        !     resume_path = 'old.dat.checkpoint'
                        call get_command_argument(iarg, resume_path)
                    ! -----
                    ! B2. 指定输出文件路径
                    ! -----
                    case ('--output')
                        ! 读取用户指定的输出文件名。
                        ! 例如：
                        !     --output result.dat
                        ! 得到：
                        !     output_path = 'result.dat'
                        ! 主程序后面将使用这个路径，
                        ! 替代自动生成的默认输出文件名。
                        call get_command_argument(iarg, output_path)
                    ! -----
                    ! B3. 限制本次新计算的参数点数
                    ! -----
                    case ('--stop-after')
                        ! 读取参数值到字符串 arg。
                        ! 例如：
                        !     --stop-after 3
                        ! 此时 arg = '3'。
                        call get_command_argument(iarg, arg)       
                        ! 命令行参数读取后属于字符串类型，
                        ! 而 stop_after 是 integer 类型。
                        ! 因此需要通过内部 read，
                        ! 将字符串转换为整数。
                        ! iostat=ios 用于检查读取是否成功：
                        !     ios = 0：成功
                        !     ios /= 0：转换失败
                        read(arg,*,iostat=ios) stop_after
                        ! 如果用户输入的不是有效整数，
                        ! 例如 --stop-after abc，
                        ! 就终止程序并报告错误。
                        if (ios /= 0) then
                            error stop 'TBG_STRAIN_CHI: invalid --stop-after'
                        endif
                        ! 不允许设置负数计算点。
                        !
                        ! 例如 --stop-after -3 不合法。
                        if (stop_after < 0) then
                            error stop 'TBG_STRAIN_CHI: negative --stop-after'
                        endif
            
                    end select
            ! ----------
            ! 选项 C：无法识别的命令行参数
            ! ----------
            case default
                ! 如果用户传入了不属于上述四种选项的参数，
                ! 说明命令行输入不符合程序要求。
                ! 例如： --unknown
                ! 此时终止程序，同时显示支持的选项
                error stop &
            'Options: --resume FILE --output FILE --check-resume --stop-after NEW_POINTS'

        end select
        ! ----------
        ! .4 移动到下一个命令行参数
        ! ----------
        ! 这里的 iarg = iarg + 1
        ! 与前面 case 分支中的 iarg = iarg + 1
        ! 具有不同作用。
        ! 前面的增加：用于跳过选项名称，读取选项对应的值。
        ! 这里的增加：用于进入下一轮循环，处理下一个选项。
        ! 例如：--stop-after 3 --output result.dat
        ! 处理顺序：
        !     iarg=1  读取 --stop-after
        !     iarg=2  读取 3
        !     iarg=3  读取 --output
        !     iarg=4  读取 result.dat
        !     iarg=5  循环结束
        iarg = iarg + 1
    enddo

    ! ==========
    ! 5. 数值参数检查与扫描列表
    ! ==========
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
        checkpoint_mode = 'strain_epsilon_chi0'
        nlambda = ie_end-ie_start+1
        chi_unit = 'mu_B/cell per unit strain_epsilon'
        ! Mx_OM(1) 必须对应 epsilon=0。
        if (abs(epsilon_list(ie_start)) > 1.d-14) then
            error stop 'TBG_STRAIN_CHI: epsilon scan must start from zero'
        endif
    else
        scan_parameter = 'strain_phi'
        checkpoint_mode = 'strain_phi'
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
    completed = 0 ! 表示已经完成计算的扫描点数

    ! ==========
    ! 读取旧checkpoint, 实现断点续算
    ! ==========
    ! ----------
    ! .1 判断用户是否指定了恢复路径
    ! ----------
    ! len_trim 返回去除尾部空格后的字符串长度。
    ! 如果 resume_path 不是空字符串，
    ! 说明用户提供了 --resume FILE。
    ! 例如：resume_path = 'previous.dat.checkpoint'
    ! 则 len_trim(resume_path) > 0。
    if (len_trim(resume_path) > 0) then
        ! ----------
        ! .2 加载 checkpoint
        ! ----------
        ! LOAD_CHI_CHECKPOINT 读取已有的 checkpoint，
        ! 同时检查它是否与当前物理参数、数值参数匹配。
        ! 参数分成几个部分：
        ! 1. trim(resume_path) 要读取的 checkpoint 文件。
        ! 2. scan_parameter 扫描类型：epsilon 或 phi。
        ! 3. 整数参数数组 N、nk、ncut、扫描范围和数组长度。
        ! 4. 实数参数数组 theta、d、mu、eta、prefactor、delta。
        ! 5. epsilon_list、phi_list 当前使用的完整扫描网格。
        ! 6. saved_results 输出：恢复的物理量。
        ! 7. completed 输出：已经完成的参数点数量。
        call LOAD_CHI_CHECKPOINT(trim(resume_path), trim(checkpoint_mode), &
        [N,nk,ncut,nlambda,ie_start,ie_end,ip_start,ip_end, &
         size(epsilon_list),size(phi_list)], &
        [theta,d,mu,eta,prefactor,delta], &
        epsilon_list, phi_list, saved_results, completed)
        ! ----------
        ! .3 将 checkpoint 的数据还原到物理量数组
        ! ----------
        Mx_OM(:completed) = saved_results(1,:completed)
        My_OM(:completed) = saved_results(2,:completed)
        chi_fix_x(:completed) = saved_results(3,:completed)
        chi_fix_y(:completed) = saved_results(4,:completed)
        chi_fs_x(:completed) = saved_results(5,:completed)
        chi_fs_y(:completed) = saved_results(6,:completed)
        chi_tot_x(:completed) = saved_results(7,:completed)
        chi_tot_y(:completed) = saved_results(8,:completed)
        ! ----------
        ! 打印恢复进度
        ! ----------
        write(*,'(A,I0,A,I0,A)') 'Recovered ', completed, '/', nlambda, ' completed points.' ! like recovered 15/21 completed points
    endif

    ! ==========
    ! 仅检查checkpoint 的模式
    ! ==========
    ! 如果用户输入 --check-resume, check_resume 就被设为 .true.。
    if (check_resume) then
        ! check-resume 必须与 resume 配合使用。i.e.如果没有提供 checkpoint 路径，就报错。
        if (len_trim(resume_path) == 0) error stop '--check-resume requires --resume FILE'
            ! 说明前面的数据读取和参数检查已经通过
            write(*,'(A)') 'Checkpoint parameters and results validated; no calculation or files written.'
        stop
    endif


    ! ==========
    ! 构造哈密顿量
    ! ==========
    call Z_TBG(d, zdiag)

    if (scan_parameter == "strain_epsilon") then
        ! ----- set epsilon = 0.d0 ----- 
        ham_param_ref = [theta, 0.d0, phi_list(ip_start)]
        call GEOMETRY_STRAIN(theta, 0.d0, phi_list(ip_start), bs_ref)
    else
        ham_param_ref = [theta, epsilon_list(ie_start), phi_list(ip_start)]
        call GEOMETRY_STRAIN(theta, epsilon_list(ie_start), phi_list(ip_start), bs_ref)
    end if


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
    ! 输出时间
    ! ----------
    call date_and_time(date=date_stamp, time=time_stamp)

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

    ! ----- 命令行指定的路径优先 -----
    if (len_trim(output_path) > 0) then
        filename = trim(output_path)
    endif
    ! ----- 防止文件被覆盖 ----- 
    inquire(file=trim(filename), exist=file_exists)
    inquire(file=trim(filename)//'.checkpoint', exist=checkpoint_exists)

    if (file_exists .or. checkpoint_exists) then
        write(*,'(A)') 'Output or checkpoint already exists: '
        write(*,'(A)') trim(filename)
        error stop 'TBG_STRAIN_CHI: use a NEW --output filename'
    endif

    ! ----------
    ! 建立输出文件
    ! ----------
    io_message = ''
    open(newunit=fileunit, file=trim(filename), status='new', action='write', &
        iostat=ios, iomsg=io_message)

    if (ios /= 0) then
        write(*,'(A)') 'Cannot create output file: '//trim(filename)
        write(*,'(A,I0)') 'IOSTAT = ', ios
        write(*,'(A)') 'IOMSG = '//trim(io_message)
        error stop 'TBG_STRAIN_CHI: cannot create output file'
    endif

    ! ----------
    ! 建立新checkpoint
    ! 不覆盖已有的checkpoint
    ! ----------
    call CREATE_CHI_CHECKPOINT( &
        trim(filename)//'.checkpoint', trim(checkpoint_mode), &
        [N,nk,ncut,nlambda,ie_start,ie_end,ip_start,ip_end, &
         size(epsilon_list),size(phi_list)], &
        [theta,d,mu,eta,prefactor,delta], &
        epsilon_list, phi_list, checkpoint_unit)
    ! ----- 将已恢复的数据写入新的checkpoint ----- 
    do ilambda = 1, completed
        call APPEND_CHI_POINT( &
            checkpoint_unit, ilambda, saved_results(:,ilambda))
    enddo

    write(*,'(A,A)')'Per-point checkpoint: ', trim(filename)//'.checkpoint'


    ! ==========
    ! epsilon == 0 处的susceptibility
    ! ==========
    chi0_fix_x = 0.d0
    chi0_fix_y = 0.d0
    chi0_fs_x = 0.d0
    chi0_fs_y = 0.d0
    chi0_tot_x = 0.d0
    chi0_tot_y = 0.d0

    if (scan_parameter == 'strain_epsilon') then 
        if (completed > 0) then 
            ! ----------
            ! 续算：直接恢复已保存的 chi(0)
            ! ----------
            chi0_fix_x = chi_fix_x(1)
            chi0_fix_y = chi_fix_y(1)
            chi0_fs_x = chi_fs_x(1)
            chi0_fs_y = chi_fs_y(1)
            chi0_tot_x = chi_tot_x(1)
            chi0_tot_y = chi_tot_y(1)
        else
            ! ----------
            ! 首次计算：计算 epsilon=0 的 chi
            ! ----------
            write(*,'(A)') 'Calculating chi(epsilon=0) once...'
            flush(output_unit)
            call CALCULATE_CHI_CUTOFF( &
                N, nk, ncut, bs_ref(:,1), bs_ref(:,2), &
                mu, eta, prefactor, zdiag, ham_param_ref, &
                HAM_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, &
                chi0_fix_x, chi0_fix_y, &
                chi0_fs_x, chi0_fs_y, &
                chi0_tot_x, chi0_tot_y)
        endif
        write(*,'(A,2(1X,ES20.10))') 'chi_0_total(x,y) =', chi0_tot_x, chi0_tot_y
        flush(output_unit)
    endif


    ! ==========
    ! 输出表头
    ! ==========
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
    write(fileunit,'(A, I0, A,I0,A,I0)')'# N = ', N, '# nk = ', nk, ', ncut = ', ncut
    write(fileunit,'(A,ES12.4)') '# OM_delta = ', delta
    ! ----- Units -----
    write(fileunit,'(A)') '# OM_unit = mu_B/cell'
    write(fileunit,'(A,A)') '# chi_unit = ', trim(chi_unit)
    ! ----- 重构说明 ----- 
    if (scan_parameter == 'strain_epsilon') then
        write(fileunit,'(A)') '# chi_mode = chi(epsilon=0, fixed phi), repeated on every row'
        write(fileunit,'(A)') '# M_from_chi = chi(epsilon=0)*epsilon (NO numerical baseline)'
        write(fileunit,'(A)')  '# delta_M = M_from_chi - (M_OM(epsilon)-M_OM(0))'
        write(fileunit,'(A,2(1X,ES24.16E3))')  '# chi0_tot_xy =', chi0_tot_x, chi0_tot_y
    else
        write(fileunit,'(A)') '# chi_mode = chi_phi(epsilon fixed, phi current)'
        write(fileunit,'(A)') '# M_from_chi = M_OM(0) + chi(epsilon=0)*epsilon'
        write(fileunit,'(A)') '# delta_M = M_from_chi - M_OM(epsilon)'
    endif

    ! ---------- 数据列 ----------
    write(fileunit,'(A)') '# columns: epsilon_percent phi_deg Mx_OM My_OM '// &
        'Mx_from_chi My_from_chi '// &
        'chi_fix_x chi_fix_y chi_fs_x chi_fs_y '// &
        'chi_tot_x chi_tot_y delta_Mx delta_My'

    flush(fileunit)

    ! ==========
    ! 逐点进行计算
    ! ==========
    do ilambda = completed+1, nlambda
        ! ---------- 检查计算点数量限制 ----------
        if (ilambda-completed > stop_after) then
            close(checkpoint_unit)
            close(fileunit)
            write(*,'(A)') &
                'Stopped after requested new points; '// &
                'resume with a NEW --output.'
            stop
        endif

        point_started = omp_get_wtime()

        ! ----------
        ! 确定当前 epsilon、phi 索引
        ! ----------
        if (scan_parameter == 'strain_epsilon') then
            iepsilon = ie_start+ilambda-1
            iphi = ip_start
        else
            iepsilon = ie_start
            iphi = ip_start+ilambda-1
        endif
        ! ----------
        ! 当前应变参数
        ! ----------
        strain_epsilon = epsilon_list(iepsilon)
        strain_phi = phi_list(iphi)

        ham_param = [theta, strain_epsilon, strain_phi]
        ! ----------
        ! 当前 moire BZ
        ! ----------
        call GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs)

        ! ==========
        ! 计算极化率
        ! ==========
        if (scan_parameter == 'strain_epsilon') then
            ! ----------
            ! epsilon 扫描：
            ! 不再逐点重新计算 susceptibility。
            ! 全部使用 epsilon=0 的结果。
            ! ----------
            chi_fix_x(ilambda) = chi0_fix_x
            chi_fix_y(ilambda) = chi0_fix_y

            chi_fs_x(ilambda) = chi0_fs_x
            chi_fs_y(ilambda) = chi0_fs_y

            chi_tot_x(ilambda) = chi0_tot_x
            chi_tot_y(ilambda) = chi0_tot_y

        else
            ! ----------
            ! phi 扫描：
            ! 逐点计算 chi_phi
            ! ----------
            write(*,'(A,I0,A,I0,A)') &
                'point ', ilambda, '/', nlambda, ': CHI starting'

            flush(output_unit)

            call CALCULATE_CHI_CUTOFF( &
                N, nk, ncut, bs(:,1), bs(:,2), &
                mu, eta, prefactor, zdiag, ham_param, &
                HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, &
                chi_fix_x(ilambda), chi_fix_y(ilambda), &
                chi_fs_x(ilambda), chi_fs_y(ilambda), &
                chi_tot_x(ilambda), chi_tot_y(ilambda))

        endif

        ! ==========
        ! 直接计算OM
        ! ==========
        write(*,'(A,I0,A,I0,A)') 'point ', ilambda, '/', nlambda, ': OM starting'

        flush(output_unit)

        if (scan_parameter == 'strain_epsilon') then
            call CALCULATE_OM_CUTOFF_KET( &
                N, nk, ncut, bs(:,1), bs(:,2), &
                mu, prefactor, zdiag, ham_param, &
                HAM_STRAIN, HAM_DERIVATIVES_STRAIN_EPSILON, &
                Mx_OM(ilambda), My_OM(ilambda))
        else
            call CALCULATE_OM_CUTOFF_KET( &
                N, nk, ncut, bs(:,1), bs(:,2), &
                mu, prefactor, zdiag, ham_param, &
                HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, &
                Mx_OM(ilambda), My_OM(ilambda))
        endif

        ! ==========
        ! 单点 checkpoint
        ! ==========
        call APPEND_CHI_POINT(checkpoint_unit, ilambda, &
            [Mx_OM(ilambda), My_OM(ilambda), &
             chi_fix_x(ilambda), chi_fix_y(ilambda), &
             chi_fs_x(ilambda), chi_fs_y(ilambda), &
             chi_tot_x(ilambda), chi_tot_y(ilambda)])

        ! ==========
        ! 屏幕打印计算进度
        ! ==========
        write(*,'(A,I0,A,I0,A,F10.5,A,F10.3,2(A,ES24.14))') &
        'point ', ilambda, '/', nlambda, &
        ' epsilon(%) = ', strain_epsilon*100.d0, &
        ' phi(deg) = ', strain_phi*180.d0/PI, &
        ' Mx = ', Mx_OM(ilambda), &
        ' My = ', My_OM(ilambda)
 
        write(*,'(A,F12.2,A)') &
             '  checkpoint saved; wall time = ', &
             omp_get_wtime()-point_started, ' s'
 
        flush(output_unit)

    enddo

    close(checkpoint_unit)

    ! ==========
    ! 重构轨道磁化
    ! ==========
    if (scan_parameter == 'strain_epsilon') then
        Mx_from_chi = Mx_OM(1) + chi0_tot_x * (epsilon_list(ie_start:ie_end) - epsilon_list(ie_start))
        My_from_chi = My_OM(1) + chi0_tot_y * (epsilon_list(ie_start:ie_end) - epsilon_list(ie_start))
    else
        Mx_from_chi(1) = Mx_OM(1)
        My_from_chi(1) = My_OM(1)
        
        do ilambda = 2, nlambda
            Mx_from_chi(ilambda) = Mx_from_chi(ilambda-1) + 0.5d0 * (chi_tot_x(ilambda)+chi_tot_x(ilambda-1)) * &
                                    (phi_list(ip_start+ilambda-1) - phi_list(ip_start+ilambda-2))
            My_from_chi(ilambda) = My_from_chi(ilambda-1) + 0.5d0 * (chi_tot_y(ilambda)+chi_tot_y(ilambda-1)) * &
                                    (phi_list(ip_start+ilambda-1) - phi_list(ip_start+ilambda-2))
        enddo
    end if


    ! ==========
    ! 输出完整计算结果
    ! ==========
    do ilambda = 1, nlambda
        ! ----- 确定当前索引 ----- 
        if (scan_parameter == 'strain_epsilon') then
            iepsilon = ie_start+ilambda-1
            iphi = ip_start
        else
            iepsilon = ie_start
            iphi = ip_start+ilambda-1
        endif
        ! ----------
        ! 计算两种OM之间的残差
        ! ----------
        residual_x = Mx_from_chi(ilambda) - Mx_OM(ilambda)
        residual_y = My_from_chi(ilambda) - My_OM(ilambda)
        ! ----------
        ! 写入14列数据
        ! ----------
        write(fileunit,'(14(ES32.16E3,1X))') &
            epsilon_list(iepsilon)*100.d0, &
            phi_list(iphi)*180.d0/PI, &
            Mx_OM(ilambda), My_OM(ilambda), &
            Mx_from_chi(ilambda), My_from_chi(ilambda), &
            chi_fix_x(ilambda), chi_fix_y(ilambda), &
            chi_fs_x(ilambda), chi_fs_y(ilambda), &
            chi_tot_x(ilambda), chi_tot_y(ilambda), &
            residual_x, residual_y

    enddo

    ! ----------
    ! 关闭输出文件
    ! ----------
    close(fileunit) 
    ! ----------
    ! 在终端打印输出文件路径
    ! ----------
    write(*,'(A,A)') 'Results saved to: ', trim(filename)

end program TBG_STRAIN_CHI

        










