program QWZ_OM
    use CONSTANTS, only : PI 
    use BILAYER_QWZ, only : HAM_QWZ, Z_BQWZ
    use OM, only : CALCULATE_OM
    use stdlib_math, only : linspace
    implicit none

    integer, parameter :: N = 4

    integer :: nk, fileunit, i
    real(8) :: m, t_perp
    real(8) :: mu, d, delta, delta_h
    real(8) :: prefactor
    real(8) :: b1(2), b2(2)
    real(8) :: ham_param(3)
    real(8) :: zdiag(N)
    real(8) :: lambda_list(21)
    character(len=256) :: filename
    character(len=16) :: s_m, s_tperp, s_mu

    real(8) :: Mx, My

    ! ==========
    ! Bilayer QWZ model
    ! ==========
    m = 1.d0
    t_perp =  0.2d0
    lambda_list = linspace(0.0, 0.2, 21)

    ! ----- chemical potential -----
    mu = 0.2d0
    d = 1.d0

    ! -----
    ! e/(hbar*c)
    ! -----
    prefactor = 1.d0

    ! -----
    ! Square-lattice reciprocal vectors
    ! k = u*b1 + v*b2. u,v in [-1/2,1/2)
    ! -----
    b1 = [2.d0*PI, 0.d0]
    b2 = [0.d0, 2.d0 * PI]

    ! -----
    ! Numerical parameter 
    ! delta : 差分步长，nk ： 沿着倒格矢基矢每一条边的分割数
    ! -----
    delta = 1.d-3
    nk = 101
    delta_h  = delta*(norm2(b1))

    ! ==========
    ! z-position operator
    ! ==========
    call Z_BQWZ(d, zdiag)


    ! ----- character string -----
    write(s_m, '(F8.3)')m
    write(s_tperp, '(F8.3)')t_perp
    write(s_mu, '(F8.3)')mu 

    ! ==========
    ! Output
    ! ==========
    ! -----file name-----
    filename = '../data/qwz/fortran/om/qwz_om_m_'//trim(adjustl(s_m))//'_tperp_'//trim(adjustl(s_tperp))//'_mu_'//trim(adjustl(s_mu))//'.dat'

    open(newunit=fileunit, file = trim(filename), status = 'replace')

    write(fileunit, '(A, ES32.16)') '# m = ', m
    write(fileunit, '(A, ES32.16)') '# t_perp = ', t_perp
    write(fileunit, '(A, ES32.16)') '# mu = ', mu
    write(fileunit, '(A, ES32.16)') '# d = ', d
    write(fileunit, '(A,ES32.16)') '# delta=', delta
    write(fileunit, '(A,I0)')      '# nk=', nk
    write(fileunit, '(A)') '# columns: Lambda Mx My'

    do i = 1, size(lambda_list)

        ham_param(1) = lambda_list(i)
        ! ----- Hamiltonian parameter -----
        ham_param = [lambda_list(i), m, t_perp]

        call CALCULATE_OM(4, nk, b1, b2, delta_h, mu, prefactor, zdiag, ham_param, HAM_QWZ, Mx, My)
        write(fileunit, '(3ES32.16)')lambda_list(i), Mx, My

    enddo

    close(fileunit)


end program QWZ_OM