program dirac_energy
    use, intrinsic :: iso_fortran_env, only : error_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use CONSTANTS, only : PI, eV, tr, NTBG, w0, w1, vF
    use STRAIN, only : GEOMETRY_STRAIN, nu, alpha, gamma_vF
    use TBG, only : MOIRE_BZ_SCALE
    use DIRAC, only : dirac_candidate, FIND_DIRAC_POINTS
    implicit none

    real(8), parameter :: epsilon_max_percent = 0.2d0, phi_max_deg = 180.d0
    real(8) :: phi_fixed_deg = 0.d0, epsilon_fixed_percent = 0.1d0, theta_deg = 1.4d0
    real(8) :: theta, momentum_scale
    integer :: n_epsilon = 21, n_phi = 91
    character(:), allocatable :: mode, output_dir, arg

    mode = 'both'
    output_dir = '../data/tbg/fortran/dirac'
    if (command_argument_count() > 7) call fail('Too many arguments. Use --help for usage.')
    if (command_argument_count() >= 1) then
        arg = argument(1)
        if (arg == '--help' .or. arg == '-h') then
            call usage()
            stop
        end if
        mode = arg
    end if
    call read_real_argument(2, phi_fixed_deg)
    call read_real_argument(3, epsilon_fixed_percent)
    call read_real_argument(4, theta_deg)
    if (command_argument_count() >= 5) output_dir = argument(5)
    call read_integer_argument(6, n_epsilon)
    call read_integer_argument(7, n_phi)

    if (mode /= 'both' .and. mode /= 'epsilon' .and. mode /= 'phi') call fail('Mode must be both, epsilon, or phi.')
    if (phi_fixed_deg < 0.d0 .or. phi_fixed_deg > phi_max_deg) call fail('Fixed phi must be in [0, 180] degrees.')
    if (epsilon_fixed_percent < 0.d0 .or. epsilon_fixed_percent > epsilon_max_percent) &
        call fail('Fixed epsilon must be in [0, 0.2] percent.')
    if (theta_deg <= 0.d0 .or. theta_deg >= 180.d0) call fail('Theta must be in (0, 180) degrees.')
    if (n_epsilon < 2 .or. n_phi < 2) call fail('Each scan needs at least two points, including both endpoints.')
    if (len_trim(output_dir) == 0) call fail('Output directory cannot be empty.')
    call make_directory(output_dir)
    theta = theta_deg*PI/180.d0
    call MOIRE_BZ_SCALE(theta, momentum_scale)

    if (mode == 'both' .or. mode == 'epsilon') call scan(.true., n_epsilon)
    if (mode == 'both' .or. mode == 'phi') call scan(.false., n_phi)

contains

    subroutine scan(vary_epsilon, n_points)
        logical, intent(in) :: vary_epsilon
        integer, intent(in) :: n_points
        type(dirac_candidate), allocatable :: candidates(:)
        real(8) :: epsilon_percent, phi_deg, epsilon, phi, bs(2,2), q1(2), seeds(2,2), previous(2,2)
        logical :: previous_valid(2)
        integer :: unit, ios, i, branch, n_found, failures, stride
        character(32) :: count_tag
        character(512) :: io_message
        character(:), allocatable :: filename, scan_name

        write(count_tag, '(I0)') n_points
        if (vary_epsilon) then
            scan_name = 'epsilon'
            filename = output_dir//'/epsilon_0-0p2pct_phi_'//number_tag(phi_fixed_deg)// &
                'deg_theta_'//number_tag(theta_deg)//'deg_n'//trim(count_tag)//'.dat'
        else
            scan_name = 'phi'
            filename = output_dir//'/phi_0-180deg_epsilon_'//number_tag(epsilon_fixed_percent)// &
                'pct_theta_'//number_tag(theta_deg)//'deg_n'//trim(count_tag)//'.dat'
        end if
        open(newunit=unit, file=filename, status='replace', action='write', iostat=ios, iomsg=io_message)
        if (ios /= 0) call fail('Cannot open '//filename//': '//trim(io_message))
        call write_header(unit, scan_name, n_points)
        previous_valid = .false.
        previous = 0.d0
        failures = 0
        stride = max(1, (n_points - 1)/10)

        do i = 1, n_points
            epsilon_percent = epsilon_fixed_percent
            phi_deg = phi_fixed_deg
            if (vary_epsilon) then
                epsilon_percent = epsilon_max_percent*real(i - 1, 8)/real(n_points - 1, 8)
            else
                phi_deg = phi_max_deg*real(i - 1, 8)/real(n_points - 1, 8)
            end if
            epsilon = epsilon_percent/100.d0
            phi = phi_deg*PI/180.d0

            ! Only immediately preceding valid results become continuation seeds.
            ! An invalid branch restarts at its nominal cone for this geometry.
            if (any(previous_valid)) then
                call GEOMETRY_STRAIN(theta, epsilon, phi, bs, q1)
                seeds(:,1) = 0.d0
                seeds(:,2) = q1
                do branch = 1, 2
                    if (previous_valid(branch)) seeds(:,branch) = previous(:,branch)
                end do
                call FIND_DIRAC_POINTS(theta, epsilon, phi, candidates, n_found, seeds=seeds)
            else
                call FIND_DIRAC_POINTS(theta, epsilon, phi, candidates, n_found)
            end if

            do branch = 1, 2
                write(unit, '(2(ES24.16,1X),I1,1X,4(ES24.16,1X),2(I1,1X),I0)', &
                    iostat=ios, iomsg=io_message) &
                    epsilon_percent, phi_deg, branch, candidates(branch)%kx/momentum_scale, &
                    candidates(branch)%ky/momentum_scale, candidates(branch)%energy/eV*1.d3, &
                    candidates(branch)%gap/eV*1.d3, merge(1, 0, candidates(branch)%converged), &
                    merge(1, 0, candidates(branch)%valid), candidates(branch)%evaluations
                if (ios /= 0) call fail('Cannot write '//filename//': '//trim(io_message))
                previous_valid(branch) = candidates(branch)%valid
                if (previous_valid(branch)) previous(:,branch) = [candidates(branch)%kx, candidates(branch)%ky]
                if (.not. candidates(branch)%valid) failures = failures + 1
            end do
            if (i == 1 .or. i == n_points .or. mod(i - 1, stride) == 0) &
                write(*, '(A,1X,I0,A,I0,2(A,F9.4),A,I0)') trim(scan_name), i, '/', n_points, &
                    ' epsilon(%)=', epsilon_percent, ' phi(deg)=', phi_deg, ' valid cones=', n_found
        end do
        close(unit, iostat=ios, iomsg=io_message)
        if (ios /= 0) call fail('Cannot close '//filename//': '//trim(io_message))
        write(*, '(A)') 'Saved '//filename
        if (failures > 0) write(*, '(A,I0,A)') 'Warning: ', failures, ' invalid branch rows retained; inspect valid/gap columns.'
    end subroutine scan

    subroutine write_header(unit, scan_name, n_points)
        integer, intent(in) :: unit, n_points
        character(*), intent(in) :: scan_name
        write(unit, '(A)') '# Shifted Dirac-point energy; heterostrain BM model'
        write(unit, '(A,A,A,I0)') '# scan=', scan_name, ' n_points=', n_points
        write(unit, '(A,ES24.16)') '# theta_deg=', theta_deg
        write(unit, '(A,ES24.16)') '# fixed_phi_deg=', phi_fixed_deg
        write(unit, '(A,ES24.16)') '# fixed_epsilon_percent=', epsilon_fixed_percent
        write(unit, '(A,I0,A,I0)') '# tr=', tr, ' Hamiltonian_dimension=', 2*NTBG
        write(unit, '(A,ES24.16,A,ES24.16)') '# w0_meV=', w0/eV*1.d3, ' w1_meV=', w1/eV*1.d3
        write(unit, '(A,ES24.16)') '# hbar_vF_eV_m=', vF/eV
        write(unit, '(A,ES24.16,A,ES24.16,A,ES24.16)') '# nu=', nu, ' alpha_eV=', alpha/eV, &
            ' gamma_vF_eV=', gamma_vF/eV
        write(unit, '(A,ES24.16)') '# l_m^-1=MOIRE_BZ_SCALE(theta)=', momentum_scale
        write(unit, '(A)') '# epsilon_percent / 100 is dimensionless strain; theta and phi are converted to radians internally.'
        write(unit, '(A)') '# E_D=(E_conduction+E_valence)/2; gap=E_conduction-E_valence; both output in meV.'
        write(unit, '(A)') '# branch=1: cone continued from top/origin; branch=2: cone continued from bottom/q1.'
        write(unit, '(A)') '# Branch labels track trajectories, not permanent layer identities.'
        write(unit, '(A)') '# Default gap_tolerance_eV=1e-8; k_tolerance=1e-8 in units theta*|K_graphene|.'
        write(unit, '(A)') '# converged/valid are 1=true, 0=false; failed rows are retained. No division by epsilon.'
        write(unit, '(A)') '# epsilon_percent phi_deg branch kx/l ky/l E_D_meV gap_meV converged valid evaluations'
    end subroutine write_header

    function argument(index) result(value)
        integer, intent(in) :: index
        character(:), allocatable :: value
        integer :: length, stat
        call get_command_argument(index, length=length, status=stat)
        if (stat /= 0) call fail('Cannot read command-line argument.')
        allocate(character(length) :: value)
        call get_command_argument(index, value=value, status=stat)
        if (stat /= 0) call fail('Cannot read command-line argument.')
    end function argument

    subroutine read_real_argument(index, value)
        integer, intent(in) :: index
        real(8), intent(inout) :: value
        character(:), allocatable :: text
        integer :: stat
        if (command_argument_count() < index) return
        text = argument(index)
        read(text, *, iostat=stat) value
        if (stat /= 0) call fail('Invalid numeric argument: '//text)
        if (.not. ieee_is_finite(value)) call fail('Numeric arguments must be finite: '//text)
    end subroutine read_real_argument

    subroutine read_integer_argument(index, value)
        integer, intent(in) :: index
        integer, intent(inout) :: value
        character(:), allocatable :: text
        integer :: stat
        if (command_argument_count() < index) return
        text = argument(index)
        read(text, *, iostat=stat) value
        if (stat /= 0) call fail('Invalid integer argument: '//text)
    end subroutine read_integer_argument

    function number_tag(value) result(tag)
        real(8), intent(in) :: value
        character(:), allocatable :: tag
        character(32) :: buffer
        integer :: i
        write(buffer, '(F0.6)') value
        tag = trim(adjustl(buffer))
        do while (tag(len(tag):len(tag)) == '0')
            tag = tag(:len(tag)-1)
        end do
        if (tag(len(tag):len(tag)) == '.') tag = tag(:len(tag)-1)
        if (len(tag) == 0) tag = '0'
        if (tag(1:1) == '.') tag = '0'//tag
        do i = 1, len(tag)
            select case (tag(i:i))
            case ('.')
                tag(i:i) = 'p'
            case ('-')
                tag(i:i) = 'm'
            end select
        end do
    end function number_tag

    subroutine make_directory(path)
        character(*), intent(in) :: path
        character(:), allocatable :: quoted
        character(512) :: message
        integer :: i, cmdstat, exitstat
        ! POSIX single-quote escaping: an apostrophe becomes '"'"'.
        quoted = achar(39)
        do i = 1, len_trim(path)
            if (path(i:i) == achar(39)) then
                quoted = quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
            else
                quoted = quoted//path(i:i)
            end if
        end do
        quoted = quoted//achar(39)
        call execute_command_line('mkdir -p -- '//quoted, exitstat=exitstat, cmdstat=cmdstat, cmdmsg=message)
        if (cmdstat /= 0) call fail('Cannot create output directory: '//trim(message))
        if (exitstat /= 0) call fail('mkdir failed for '//path)
    end subroutine make_directory

    subroutine usage()
        write(*, '(A)') 'Usage: dirac_energy [both|epsilon|phi] [phi_fixed_deg] [epsilon_fixed_percent] '// &
            '[theta_deg] [output_dir] [n_epsilon] [n_phi]'
        write(*, '(A)') 'Defaults: both 0 0.1 1.4 ../data/tbg/fortran/dirac 21 91'
        write(*, '(A)') 'Inclusive scans: epsilon=0..0.2 percent; phi=0..180 degrees.'
    end subroutine usage

    subroutine fail(message)
        character(*), intent(in) :: message
        write(error_unit, '(A)') 'dirac_energy: '//message
        error stop 1
    end subroutine fail

end program dirac_energy
