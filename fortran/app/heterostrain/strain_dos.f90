program strain_dos
    use, intrinsic :: iso_fortran_env, only : error_unit, int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use CONSTANTS, only : PI, eV, tr, NTBG, w0, w1, vF
    use STRAIN, only : nu, alpha, gamma_vF
    use DOS, only : CAL_DOS_TRIANGLE
    implicit none

    real(8) :: theta_deg = 1.4d0, epsilon_percent = 0.1d0, phi_deg = 90.d0
    real(8) :: e_min_mev = -50.d0, e_max_mev = 50.d0, elapsed
    real(8), allocatable :: energy(:), density(:)
    integer :: nk = 100, n_intervals = 2000, unit, ios, i
    integer(int64) :: clock_start, clock_finish, clock_rate
    character(512) :: io_message
    character(:), allocatable :: output_dir, filename, arg

    output_dir = '../data/tbg/fortran/dos'
    if (command_argument_count() > 8) call fail('Too many arguments. Use --help for usage.')
    if (command_argument_count() >= 1) then
        arg = argument(1)
        if (arg == '--help' .or. arg == '-h') then
            call usage()
            stop
        end if
    end if
    call read_real_argument(1, theta_deg)
    call read_real_argument(2, epsilon_percent)
    call read_real_argument(3, phi_deg)
    call read_integer_argument(4, nk)
    call read_integer_argument(5, n_intervals)
    call read_real_argument(6, e_min_mev)
    call read_real_argument(7, e_max_mev)
    if (command_argument_count() >= 8) output_dir = argument(8)

    if (theta_deg <= 0.d0 .or. theta_deg >= 180.d0) call fail('Theta must be in (0, 180) degrees.')
    if (epsilon_percent < 0.d0 .or. epsilon_percent >= 100.d0) &
        call fail('Epsilon must be in [0, 100) percent.')
    if (phi_deg < 0.d0 .or. phi_deg > 180.d0) call fail('Phi must be in [0, 180] degrees.')
    if (nk < 1 .or. n_intervals < 1) call fail('nk and n_intervals must be positive integers.')
    if (e_max_mev <= e_min_mev) call fail('E_max_meV must be greater than E_min_meV.')
    if (len_trim(output_dir) == 0) call fail('Output directory cannot be empty.')
    call make_directory(output_dir)
    filename = output_dir//'/dos_theta'//number_tag(theta_deg)//'_epsilon'//number_tag(epsilon_percent)// &
        '_phi'//number_tag(phi_deg)//'_nk'//integer_tag(nk)//'_ne'//integer_tag(n_intervals)// &
        '_emin'//number_tag(e_min_mev)//'_emax'//number_tag(e_max_mev)//'.dat'

    allocate(energy(n_intervals + 1), density(n_intervals + 1))
    write(*, '(A,3(F9.4,A))') 'Calculating triangle DOS: theta=', theta_deg, ' deg, epsilon=', &
        epsilon_percent, ' %, phi=', phi_deg, ' deg'
    write(*, '(A,I0,A,I0,A,I0)') 'Grid cells: ', nk, ' x ', nk, '; energy intervals: ', n_intervals
    call system_clock(count=clock_start, count_rate=clock_rate)
    call CAL_DOS_TRIANGLE(theta_deg*PI/180.d0, epsilon_percent/100.d0, phi_deg*PI/180.d0, &
        e_min_mev*1.d-3*eV, e_max_mev*1.d-3*eV, n_intervals, energy, density, nkx=nk, nky=nk)
    call system_clock(count=clock_finish)
    if (.not. all(ieee_is_finite(energy)) .or. .not. all(ieee_is_finite(density))) &
        call fail('The DOS calculation returned a non-finite result.')
    elapsed = 0.d0
    if (clock_rate > 0_int64) elapsed = real(clock_finish - clock_start, 8)/real(clock_rate, 8)

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios, iomsg=io_message)
    if (ios /= 0) call fail('Cannot open '//filename//': '//trim(io_message))
    call write_header(unit)
    do i = 1, size(energy)
        write(unit, '(2(ES24.16,1X))', iostat=ios, iomsg=io_message) energy(i)/eV, density(i)*eV
        if (ios /= 0) call fail('Cannot write '//filename//': '//trim(io_message))
    end do
    close(unit, iostat=ios, iomsg=io_message)
    if (ios /= 0) call fail('Cannot close '//filename//': '//trim(io_message))
    write(*, '(A,F12.3,A)') 'Calculation wall time: ', elapsed, ' s'
    write(*, '(A)') 'Saved '//filename

contains

    subroutine write_header(unit)
        integer, intent(in) :: unit
        write(unit, '(A)') '# Heterostrain BM density of states'
        write(unit, '(A)') '# method=bin_averaged_analytic_triangle'
        write(unit, '(A,ES24.16)') '# theta_deg=', theta_deg
        write(unit, '(A,ES24.16)') '# epsilon_percent=', epsilon_percent
        write(unit, '(A,ES24.16)') '# phi_deg=', phi_deg
        write(unit, '(A,I0)') '# nkx=', nk
        write(unit, '(A,I0)') '# nky=', nk
        write(unit, '(A,I0)') '# n_intervals=', n_intervals
        write(unit, '(A,ES24.16)') '# E_min_meV=', e_min_mev
        write(unit, '(A,ES24.16)') '# E_max_meV=', e_max_mev
        write(unit, '(A,I0)') '# tr=', tr
        write(unit, '(A,I0)') '# Hamiltonian_dimension=', 2*NTBG
        write(unit, '(A,ES24.16)') '# w0_meV=', w0/eV*1.d3
        write(unit, '(A,ES24.16)') '# w1_meV=', w1/eV*1.d3
        write(unit, '(A,ES24.16)') '# hbar_vF_eV_m=', vF/eV
        write(unit, '(A,ES24.16)') '# nu=', nu
        write(unit, '(A,ES24.16)') '# alpha_eV=', alpha/eV
        write(unit, '(A,ES24.16)') '# gamma_vF_eV=', gamma_vF/eV
        write(unit, '(A,ES24.16)') '# calculation_wall_seconds=', elapsed
        write(unit, '(A)') '# sampling=reciprocal_parallelogram_[0,1]^2; two triangles per grid cell; all bands'
        write(unit, '(A)') '# normalization=states_per_eV_per_moire_cell_single_spin_single_valley'
        write(unit, '(A)') '# A full energy integral equals one state per band; no spin/valley factor of four.'
        write(unit, '(A)') '# Triangle density is already per energy; no extra inverse energy-grid spacing.'
        write(unit, '(A)') '# Analytic triangle CDF differences give bin averages; internal bin edges are node midpoints.'
        write(unit, '(A)') '# Endpoint bins have half width; trapezoidal integration equals the sum of bin state counts.'
        write(unit, '(A)') '# No Gaussian broadening, no smoothing, no shift of the energy origin.'
        write(unit, '(A)') '# energy_eV DOS_states_per_eV_per_cell_single_spin_single_valley'
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
        character(48) :: buffer
        integer :: i
        write(buffer, '(F0.6)') value
        tag = trim(adjustl(buffer))
        do while (tag(len(tag):len(tag)) == '0')
            tag = tag(:len(tag)-1)
        end do
        if (tag(len(tag):len(tag)) == '.') tag = tag(:len(tag)-1)
        if (len(tag) == 0) tag = '0'
        if (tag(1:1) == '.') tag = '0'//tag
        if (len(tag) >= 2) then
            if (tag(:2) == '-.') tag = '-0'//tag(2:)
        end if
        do i = 1, len(tag)
            select case (tag(i:i))
            case ('.')
                tag(i:i) = 'p'
            case ('-')
                tag(i:i) = 'm'
            end select
        end do
    end function number_tag

    function integer_tag(value) result(tag)
        integer, intent(in) :: value
        character(:), allocatable :: tag
        character(32) :: buffer
        write(buffer, '(I0)') value
        tag = trim(buffer)
    end function integer_tag

    subroutine make_directory(path)
        character(*), intent(in) :: path
        character(:), allocatable :: quoted
        character(512) :: message
        integer :: i, cmdstat, exitstat
        ! POSIX single-quote escaping: an apostrophe becomes '\"'\"'.
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
        write(*, '(A)') 'Usage: strain_dos [theta_deg] [epsilon_percent] [phi_deg] [nk] [n_intervals] '// &
            '[E_min_meV] [E_max_meV] [output_dir]'
        write(*, '(A)') 'Defaults: 1.4 0.1 90 100 2000 -50 50 ../data/tbg/fortran/dos'
        write(*, '(A)') 'nk counts grid cells along each reciprocal basis vector; nk+1 vertices per direction.'
        write(*, '(A)') 'n_intervals+1 inclusive energy points. Output is single-spin, single-valley states/eV/cell.'
    end subroutine usage

    subroutine fail(message)
        character(*), intent(in) :: message
        write(error_unit, '(A)') 'strain_dos: '//message
        error stop 1
    end subroutine fail

end program strain_dos
