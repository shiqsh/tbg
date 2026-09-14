module SUSCEPTIBILITY
    use CONSTANTS, only : PI
    use MODEL_INTERFACE, only : HAMILTONIAN
    use LINALG, only : myzheev
    use OM, only : CALCULATE_OM
    implicit none

    abstract interface
        ! ----- Hamiltonian derivative interface ----- 
        subroutine H_DERIVATIVES_INTERFACE(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)
            complex(8), intent(out) :: dH_dkx(:,:)
            complex(8), intent(out) :: dH_dky(:,:)
            complex(8), intent(out) :: dH_dlambda(:,:)
        end subroutine H_DERIVATIVES_INTERFACE

    end interface

contains

! ==========
! Gaussian function
! delta_eta(E - mu) = exp[- ((E - m u)/eta)^2]/(sqrt(pi)*eta)
! ==========
subroutine GAUSSIAN_DELTA(E, mu, eta, delta)
    implicit none

    real(8), intent(in) :: E
    real(8), intent(in) :: mu
    real(8), intent(in) :: eta

    real(8), intent(out) :: delta

    if (eta <= 0.d0) then
        error stop "DELTA_GAUSSIAN: eta must be positive"
    endif

    delta = exp(-((E-mu)/eta)**2)/(sqrt(PI)*eta)

end subroutine GAUSSIAN_DELTA


! ==========
! FS contribution at one k point
! [ H_beta_mm A_lambda(m) - Hlambda_mm A_beta(m) ]
! A_lambda(m) = sum_{n /= m} Re[Hlambda_mn z_nm]
! A_beta(m) = sum_{n /= m} Re[Hbeta_mn znm]
! ==========
subroutine CHI_FS_K(N, W, Hx_e, Hy_e, Hlambda_e, z_e, &
                    mu, eta, Gx, Gy)
    implicit none

    integer, intent(in) :: N
    real(8), intent(in) :: W(N)
    real(8), intent(in) :: mu
    real(8), intent(in) :: eta
    complex(8), intent(in) :: Hx_e(N,N)
    complex(8), intent(in) :: Hy_e(N,N)
    complex(8), intent(in) :: Hlambda_e(N,N)
    complex(8), intent(in) :: z_e(N,N)

    real(8), intent(out) :: Gx
    real(8), intent(out) :: Gy

    integer :: nn, mm
    real(8) :: delta_E
    real(8) :: A_x, A_y, A_lambda
    ! ----- error check ----- 
    real(8), parameter :: tol = 1.d-12

    complex(8) :: Hx_mm, Hy_mm, Hlambda_mm

    Gx = 0.d0
    Gy = 0.d0

    ! ----------
    ! 对所有band 指标 m 求和
    ! ----------
    do mm = 1, N

        call GAUSSIAN_DELTA(W(mm), mu, eta, delta_E)

        Hx_mm = Hx_e(mm, mm)
        Hy_mm = Hy_e(mm, mm)
        Hlambda_mm = Hlambda_e(mm, mm)

        ! -----
        ! Diagonal expection values of Hermitian operators should be real 
        ! -----
        if (abs(aimag(Hx_mm)) > tol) then
            error stop "CHI_FS_K: braket{m|dH_dkx|m} is not real"
        endif
        if (abs(aimag(Hy_mm)) > tol) then
            error stop "CHI_FS_K: braket{m|dH_dky|m} is not real"
        endif 
        if (abs(aimag(Hlambda_mm)) > tol) then
            error stop "CHI_FS_K: braket{m|dH_dlambda|m} is not real"
        endif 

        ! ----- origin value ----- 
        A_x = 0.d0
        A_y = 0.d0
        A_lambda = 0.d0

        ! -----
        ! 对n /= m进行求和
        ! ----- 
        do nn = 1, N 
            if ( nn == mm ) cycle 
            ! 下面的累加语句在 nn == mm 时不会执行
            A_x = A_x + real(Hx_e(mm, nn) * z_e(nn, mm), 8)
            A_y = A_y + real(Hy_e(mm, nn) * z_e(nn, mm), 8)
            A_lambda = A_lambda + real(Hlambda_e(mm, nn) * z_e(nn, mm), 8)
        enddo

        ! ----- 回到对mm的求和 ----- 
        Gx = Gx + delta_E * real(Hx_mm * A_lambda - Hlambda_mm * A_x, 8)
        Gy = Gy + delta_E * real(Hy_mm * A_lambda - Hlambda_mm * A_y, 8)

    enddo

end subroutine CHI_FS_K


! ==========
! FIX contribution at one K point
! n in occ, m in emp
! Fx = Re sum_{n, m} [Hx_nm Zlambda_mn - Hlambda_nm Zx_mn ]/(En-Em)
! Fy = Re sum_{n, m} [Hy_nm Zlambda_mn - Hlambda_nm Zy_mn ]/(En-Em)
! where Zlambda = {z,Hlambda}, Zx = {z,Hx}, Zy = {z,Hy}
! ==========
subroutine CHI_FIXED_K(N, W, nocc, Hx_e, Hy_e, Hlambda_e, z_e, &
                        Fx, Fy)
    implicit none

    integer, intent(in) :: N
    integer, intent(in) :: nocc
    real(8), intent(in) :: W(N)
    complex(8), intent(in) :: Hx_e(N, N)
    complex(8), intent(in) :: Hy_e(N, N)
    complex(8), intent(in) :: Hlambda_e(N, N)
    complex(8), intent(in) :: z_e(N, N)

    real(8), intent(out) :: Fx
    real(8), intent(out) :: Fy

    integer :: nn, mm
    real(8) :: denominator
    complex(8) :: term_x, term_y
    complex(8) :: Zx_e(N,N), Zy_e(N,N), Zlambda_e(N,N)

    ! ----- Anticommutator rep matrix ----- 
    Zx_e = matmul(z_e, Hx_e) + matmul(Hx_e, z_e)
    Zy_e = matmul(z_e, Hy_e) + matmul(Hy_e, z_e)
    Zlambda_e = matmul(z_e, Hlambda_e) + matmul(Hlambda_e, z_e)

    Fx = 0.d0
    Fy = 0.d0

    ! ----- occupied state check ----- 
    if (nocc == 0 .or. nocc == N) return

    ! ----------
    ! 循环
    ! n: occupied, m: empty
    ! ----------
    
    do nn = 1, nocc
        do mm = nocc + 1, N

            denominator = W(nn) - W(mm)

            ! ----- beta = x ----- 
            term_x = Hx_e(nn, mm) * Zlambda_e(mm, nn) - Hlambda_e(nn, mm) * Zx_e(mm, nn)
            Fx = Fx + real(term_x/denominator, 8)
            ! ----- beta = y ----- 
            term_y = Hy_e(nn, mm) * Zlambda_e(mm, nn) - Hlambda_e(nn, mm) * Zy_e(mm, nn)
            Fy = Fy + real(term_y/denominator, 8)
            
        enddo
    enddo

end subroutine CHI_FIXED_K


! ==========
! Integrate over the whole BZ
! epsilon_xy = 1
! chi_{FIX/FS}_x/y = \pm prefactor * int {F_x/y / G_x/y}
! ==========
subroutine CALCULATE_CHI(N, nk, b1, b2, mu, eta, prefactor, zdiag, ham_param, &
                        HAM, H_DERIVATIVES, &
                        chi_fix_x, chi_fix_y, chi_fs_x, chi_fs_y, &
                        chi_tot_x, chi_tot_y)
    implicit none

    integer, intent(in) :: N
    integer, intent(in) :: nk
    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: mu
    real(8), intent(in) :: eta
    real(8), intent(in) :: prefactor
    real(8), intent(in) :: ham_param(:)
    real(8), intent(in) :: zdiag(N)

    procedure(HAMILTONIAN) :: HAM
    procedure(H_DERIVATIVES_INTERFACE) :: H_DERIVATIVES

    real(8), intent(out) :: chi_fix_x
    real(8), intent(out) :: chi_fix_y
    real(8), intent(out) :: chi_fs_x
    real(8), intent(out) :: chi_fs_y
    real(8), intent(out) :: chi_tot_x
    real(8), intent(out) :: chi_tot_y
    
    ! ----- reciprocal lattice mesh coefficients ----- 
    integer :: iu, iv, nocc
    real(8) :: u, v
    real(8) :: k(2)
    ! ----- BZ area and integration weight ----- 
    real(8) :: area_bz
    real(8) :: weight
    ! ----- susceptibility ----- 
    real(8) :: sum_Fx, sum_Fy
    real(8) :: sum_Gx, sum_Gy

    real(8) :: Fx, Fy
    real(8) :: Gx, Gy
    ! ----- Hamiltonian and eigenvalues ----- 
    real(8) :: W(N)
    complex(8) :: H(N, N)
    complex(8) :: V_HAM(N, N)
    ! ----- z rep matrix ----- 
    integer :: ii 
    complex(8) :: zV(N, N)
    complex(8) :: z_e(N, N)
    ! ----- derivative matrix ----- 
    complex(8) :: dH_dkx(N, N), dH_dky(N, N), dH_dlambda(N, N)
    ! ----- derivative rep matrix ----- 
    complex(8) :: Hx_e(N, N), Hy_e(N, N), Hlambda_e(N, N)

    ! -----
    ! BZ area and integration weight
    ! -----
    area_bz = abs(b1(1)*b2(2) - b1(2)*b2(1))
    weight = area_bz/((2.d0*PI)*real(nk, 8))**2.d0
    ! -----
    ! original value
    ! -----
    sum_Fx = 0.d0
    sum_Fy = 0.d0

    sum_Gx = 0.d0
    sum_Gy = 0.d0

    ! =====
    ! BZ integration 
    ! =====
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP PRIVATE(ii,iu,iv,u,v,k,H,W,V_HAM,nocc,dH_dkx,dH_dky,dH_dlambda, & 
    !$OMP         Hx_e,Hy_e,Hlambda_e,zV,z_e,Fx,Fy,Gx,Gy) &
    !$OMP SHARED(N,nk,b1,b2,mu,eta,zdiag,ham_param) &
    !$OMP REDUCTION(+:sum_Fx,sum_Fy,sum_Gx,sum_Gy) &
    !$OMP COLLAPSE(2) SCHEDULE(DYNAMIC)
    do iu = 1, nk
        do iv = 1, nk
            ! -----Reduced coordinates-----
            u = -0.5d0 + real(iu-1, 8)/real(nk, 8)
            v = -0.5d0 + real(iv-1, 8)/real(nk, 8)
            ! -----Cartesian momentum-----
            k = u*b1 + v*b2
            ! -----Hamiltonian------
            call HAM(k, ham_param, H)
            ! ----- diagonalization ----- 
            call myzheev(N, W, H, V_HAM)
            ! ----- count occupied band ----- 
            nocc = count(W < mu)
            
            ! -----
            ! construct z rep matrix 
            ! -----
            do ii = 1, N
                zV(ii, :) = zdiag(ii) * V_HAM(ii, :)
            enddo
            z_e = matmul(transpose(conjg(V_HAM)), zV)

            ! ----------
            ! Hamiltonian derivatives operator
            ! ----------
            call H_DERIVATIVES(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
            ! -----
            ! operator to rep matrix   
            ! -----
            Hx_e = matmul(transpose(conjg(V_HAM)), matmul(dH_dkx, V_HAM))
            Hy_e = matmul(transpose(conjg(V_HAM)), matmul(dH_dky, V_HAM))
            Hlambda_e = matmul(transpose(conjg(V_HAM)), matmul(dH_dlambda, V_HAM))
            
            ! ----------
            ! single k point
            ! ----------
            ! ----- FS ----- 
            call CHI_FS_K(N, W, Hx_e, Hy_e, Hlambda_e, z_e, mu, eta, Gx, Gy)
            ! ----- FIX ----- 
            call CHI_FIXED_K(N, W, nocc, Hx_e, Hy_e, Hlambda_e, z_e, Fx, Fy)

            ! ----------
            ! accumulation
            ! ----------
            sum_Fx = sum_Fx + Fx
            sum_Fy = sum_Fy + Fy

            sum_Gx = sum_Gx + Gx
            sum_Gy = sum_Gy + Gy

        enddo
    enddo
    !$OMP END PARALLEL DO 
    
    ! =====
    ! susceptibility
    ! =====
    ! ----- FS ----- 
    chi_fs_x = prefactor*weight*sum_Gy
    chi_fs_y = -prefactor*weight*sum_Gx
    ! ----- FIX ----- 
    chi_fix_x = - prefactor*weight*sum_Fy
    chi_fix_y = prefactor*weight*sum_Fx
    ! ----- TOT ----- 
    chi_tot_x = chi_fix_x + chi_fs_x
    chi_tot_y = chi_fix_y + chi_fs_y

end subroutine CALCULATE_CHI


! ==========
! Integrated the susceptibility 
! ==========
subroutine M_FROM_CHI(N, nk, b1, b2, delta, mu, prefactor, zdiag, &
                        ham_param_ref, HAM, lambda_list, chi_x, chi_y, &
                        Mx_from_chi, My_from_chi)
    implicit none

    integer, intent(in) :: N, nk
    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: delta, mu, prefactor
    real(8), intent(in) :: zdiag(N)
    real(8), intent(in) :: ham_param_ref(:)

    procedure(HAMILTONIAN) :: HAM

    real(8), intent(in) :: lambda_list(:)
    real(8), intent(in) :: chi_x(:), chi_y(:)

    real(8), intent(out) :: Mx_from_chi(size(lambda_list))
    real(8), intent(out) :: My_from_chi(size(lambda_list))

    integer :: i
    integer :: nlambda
    real(8) :: dlambda
    real(8) :: Mx_ref, My_ref

    ! -----
    ! 维度检查
    ! -----
    nlambda = size(lambda_list)

    if (nlambda < 1) then
        error stop "M_FROM_CHI: lambda array is empty"
    endif

    if (size(chi_x) /= nlambda) then
        error stop "M_FROM_CHI: chi_x dimension mismatch"
    endif

    if (size(chi_y) /= nlambda) then
        error stop "M_FROM_CHI: chi_y dimension mismatch"
    endif

    ! ----------
    ! calculate the reference OM
    ! M_ref = M_{lambda_list(1)}
    ! ----------
    call CALCULATE_OM(N, nk, b1, b2, delta, mu, prefactor, zdiag, &
                        ham_param_ref, HAM, Mx_ref, My_ref)

    ! -----
    ! initial point
    ! M(lambda_1) = M_ref
    ! -----
    Mx_from_chi(1) = Mx_ref
    My_from_chi(1) = My_ref

    ! -----
    ! 积分
    ! -----
    do i = 2, nlambda

        dlambda = lambda_list(i) - lambda_list(i-1)

        Mx_from_chi(i) = Mx_from_chi(i-1) + 0.5d0*(chi_x(i) + chi_x(i-1))*dlambda
        My_from_chi(i) = My_from_chi(i-1) + 0.5d0*(chi_y(i) + chi_y(i-1))*dlambda

    enddo

end subroutine M_FROM_CHI


end module SUSCEPTIBILITY