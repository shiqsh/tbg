PROGRAM QWZ_CHI
    use CONSTANTS, only : PI
    use SUSCEPTIBILITY, only : CALCULATE_CHI, M_FROM_CHI
    use BILAYER_QWZ, only : HAM_QWZ, Z_BQWZ, HAM_DERIVATIVES_QWZ
    use stdlib_math, only : linspace
    implicit none

    integer, parameter :: N = 4
    integer, parameter :: nlambda = 21

    integer :: nk
    integer :: ilambda
    integer :: fileunit

    real(8) :: m, t_perp
    real(8) :: mu, d, eta, prefactor
    real(8) :: delta, delta_h
    real(8) :: b1(2), b2(2)
    real(8) :: ham_param(3), ham_param_ref(3)
    real(8) :: zdiag(N)
    real(8) :: lambda_list(nlambda)

    ! ----- 极化率 ----- 
    real(8) :: chi_fs_x(nlambda), chi_fs_y(nlambda)
    real(8) :: chi_fix_x(nlambda), chi_fix_y(nlambda)
    real(8) :: chi_tot_x(nlambda), chi_tot_y(nlambda)
    ! ----- M reconstruct from chi ----- 
    real(8) :: Mx_from_chi(nlambda)
    real(8) :: My_from_chi(nlambda)
    ! ----- output ----- 
    character(len=256) :: filename
    character(len=16) :: s_m, s_tperp, s_mu, s_eta, s_nk

    

    ! ==========
    ! Bilayer QWZ model
    ! ==========
    m = 1.d0
    t_perp =  0.2d0
    lambda_list = linspace(0.0, 0.2, nlambda)
    d = 1.d0

    ! ----- chemical potential -----
    mu = -1.0d0

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

    ! ----------
    ! central-difference step used by OM
    ! ----------
    delta = 1.d-5
    delta_h = delta*norm2(b1)

    ! =====
    ! 高斯函数参数展宽
    ! =====
    eta = 0.03

    ! =====
    ! BZ mesh 
    ! =====
    nk = 251

    ! =====
    ! z-position operator
    ! =====
    call Z_BQWZ(d, zdiag)


    ! ==========
    ! 计算极化率
    ! ==========
    do ilambda = 1, nlambda

        ham_param = [lambda_list(ilambda), m, t_perp]

        call CALCULATE_CHI(N, nk, b1, b2, mu, eta, prefactor, zdiag, ham_param, &
                            HAM_QWZ, HAM_DERIVATIVES_QWZ, &
                            chi_fix_x(ilambda), chi_fix_y(ilambda), &
                            chi_fs_x(ilambda), chi_fs_y(ilambda), &
                            chi_tot_x(ilambda), chi_tot_y(ilambda) &
                        )
    enddo

    ! ==========
    ! 计算reference OM
    ! ==========
    ! ----- 提取原点参数 ----- 
    ham_param_ref = [lambda_list(1), m, t_perp]

    ! ----------
    ! Reconstruct M
    ! ----------
    call M_FROM_CHI( N, nk, b1, b2, delta_h, mu, prefactor, zdiag, &
                    ham_param_ref, HAM_QWZ, &
                    lambda_list, chi_tot_x, chi_tot_y, &
                    Mx_from_chi, My_from_chi &
                )


    ! ==========
    ! Output
    ! ==========
    ! ----- character string -----
    write(s_m, '(F8.3)')m
    write(s_tperp, '(F8.3)')t_perp
    write(s_mu, '(F8.3)')mu 
    write(s_eta, '(F8.3)')eta
    write(s_nk, '(I0)')nk


    filename = '../data/qwz/fortran/susceptibility' &
                // '/qwz_chi_m_'//trim(adjustl(s_m)) &
                // '_tperp_'//trim(adjustl(s_tperp)) &
                // '_mu_'//trim(adjustl(s_mu)) &
                // '_eta_'//trim(adjustl(s_eta)) &
                // '_nk_'//trim(adjustl(s_nk))//'.dat'

    open(newunit=fileunit, file = trim(filename), status = 'replace')

    write(fileunit, '(A, ES32.16)') '# m = ', m
    write(fileunit, '(A, ES32.16)') '# t_perp = ', t_perp
    write(fileunit, '(A, ES32.16)') '# mu = ', mu
    write(fileunit, '(A, ES32.16)')  "# eta = ", eta
    write(fileunit, '(A, I0)') "# nk = ", nk
    write(fileunit, '(A)') '# columns: Lambda Mx_from_chi My_from_chi'  & 
                            // 'chi_fix_x chi_fix_y chi_fs_x chi_fs_y chi_tot_x chi_tot_y'

    ! ----------
    ! 写入数据
    ! ----------
    do ilambda = 1, nlambda

        write(fileunit, '(9(ES32.16E3,1X))') &
            lambda_list(ilambda), &
            Mx_from_chi(ilambda), &
            My_from_chi(ilambda), &
            chi_fix_x(ilambda), &
            chi_fix_y(ilambda), &
            chi_fs_x(ilambda), &
            chi_fs_y(ilambda), &
            chi_tot_x(ilambda), &
            chi_tot_y(ilambda)

    enddo

    close(fileunit)
    

end program QWZ_CHI



    





    







