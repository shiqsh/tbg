module SUSCEPTIBILITY_CUTOFF
    use CONSTANTS, only : PI, eV
    use MODEL_INTERFACE, only : HAMILTONIAN, H_DERIVATIVES_INTERFACE
    use LINALG, only : myzheev
    use stdlib_linalg, only : diag, hermitian
    use SUSCEPTIBILITY, only : GAUSSIAN_DELTA
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF
    implicit none
    
contains

! ==========
! FS contribution at one k point
! G_beta = 1/2 sum_m delta_eta(E_m-mu) Re[Hbeta_mm Zlambda_mm - Hlambda_mm Zbeta_mm]
! only the two active bands contribute to the Fermi surface
! ==========
subroutine CHI_FS_K_CUTOFF(N, ncut, W, Hx_e, Hy_e, Hlambda_e, Zx_e, Zy_e, Zlambda_e, &
                            mu, eta, Gx, Gy)
    implicit none

    integer, intent(in) :: N 
    integer, intent(in) :: ncut
    real(8), intent(in) :: W(N)
    real(8), intent(in) :: mu
    real(8), intent(in) :: eta

    complex(8), intent(in) :: Hx_e(:,:)
    complex(8), intent(in) :: Hy_e(:,:)
    complex(8), intent(in) :: Hlambda_e(:,:)
    complex(8), intent(in) :: Zx_e(:,:) ! {Hx, z}
    complex(8), intent(in) :: Zy_e(:,:)
    complex(8), intent(in) :: Zlambda_e(:,:) ! {Hlambda, z}

    real(8), intent(out) :: Gx
    real(8), intent(out) :: Gy

    integer :: iv, ic
    integer :: ilow 
    integer :: nb ! number of retained band: 2*ncut +2
    integer :: mm ! full active band index
    integer :: ii ! corresponding local index in cutoff matrices

    real(8) :: delta_E
    real(8), parameter :: tol = 1.d-12  ! ! Tolerance for imaginary diagonal elements

    complex(8) :: Hx_mm, Hy_mm, Hlambda_mm ! braket{m| dH_dxx|m}
    complex(8) :: Hx_diag(2*ncut + 2), Hy_diag(2*ncut + 2), Hlambda_diag(2*ncut + 2) ! diagonal elements of cutoff Hxx_e
    complex(8) :: Zx_diag(2*ncut + 2), Zy_diag(2*ncut + 2), Zlambda_diag(2*ncut + 2)

    ! ----- Dimension check  ----- 
    if (N < 2 .or. mod(N,2) /= 0) error stop 'CHI_FS_K_CUTOFF: invalid N'
    if (ncut < 0 .or. ncut > N/2 -1) error stop 'CHI_FS_K_CUTOFF: invalid ncut'
    if (eta <= 0.d0) error stop 'CHI_FS_K_CUTOFF: eta must be positive'

    nb = 2*ncut + 2

    if (any(shape(Hx_e) /= [nb, nb]) .or.  any(shape(Hy_e) /= [nb,nb]) .or. &
        any(shape(Hlambda_e) /= [nb,nb])) error stop 'CHI_FS_K_CUTOFF: H_e dimension mismatch'
    if (any(shape(Zx_e) /= [nb,nb]) .or. any(shape(Zy_e) /= [nb,nb]) .or. &
        any(shape(Zlambda_e) /= [nb,nb])) error stop 'CHI_FS_K_CUTOFF: Z_e dimension mismatch'

    ! ----- active bands and cutoff window ----- 
    iv = N/2
    ic = N/2 + 1
    ilow = iv - ncut 

    ! ----- diagonal matrix elements from stdlib_linalg ----- 
    Hx_diag = diag(Hx_e)
    Hy_diag = diag(Hy_e)
    Hlambda_diag = diag(Hlambda_e)
    Zx_diag = diag(Zx_e)
    Zy_diag = diag(Zy_e)
    Zlambda_diag = diag(Zlambda_e)

    ! ----- 初始化 ----- 
    Gx = 0.d0
    Gy = 0.d0

    ! ----------
    ! 对两个active band进行求和
    ! ----------
    do mm = iv, ic

        call GAUSSIAN_DELTA(W(mm), mu, eta, delta_E)
        ! -----
        ! 进入ncutmatrix
        ! if mm = N/2, ii = ncut + 1， if mm = N/2 + 1, ii = ncut + 2
        ! -----
        ii = mm - ilow + 1 

        Hx_mm = Hx_diag(ii)
        Hy_mm = Hy_diag(ii) 
        Hlambda_mm = Hlambda_diag(ii)

        ! -----
        ! Diagonal expection values of Hermitian operators should be real 
        ! -----
        if (abs(aimag(Hx_mm)) > tol*maxval(abs(Hx_e))) then
            error stop "CHI_FS_K: braket{m|dH_dkx|m} is not real"
        endif
        if (abs(aimag(Hy_mm)) > tol*maxval(abs(Hy_e))) then
            error stop "CHI_FS_K: braket{m|dH_dky|m} is not real"
        endif 
        if (abs(aimag(Hlambda_mm)) > tol*maxval(abs(Hlambda_e))) then
            error stop "CHI_FS_K: braket{m|dH_dlambda|m} is not real"
        endif 

        ! ----- 对m进行求和 ----- 
        Gx = Gx + 0.5d0 * delta_E * real(Hx_mm*Zlambda_diag(ii)-Hlambda_mm*Zx_diag(ii), 8)
        Gy = Gy + 0.5d0 * delta_E * real(Hy_mm*Zlambda_diag(ii)-Hlambda_mm*Zy_diag(ii), 8)

    enddo

end subroutine CHI_FS_K_CUTOFF
    

! ==========
! FIX contribution at one K point
! Z_beta is computed in the full BM space but projected into the cutoff space
! 不针对简并点的情况，或者说mu正好穿过Dirac点的情况
! ==========
subroutine CHI_FIXED_K_CUTOFF(N, ncut, noccA, W, Hx_e, Hy_e, Hlambda_e, Zx_e, Zy_e, Zlambda_e, & 
                            Fx, Fy)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: ncut
    integer, intent(in) :: noccA
    real(8), intent(in) :: W(N)

    complex(8), intent(in) :: Hx_e(:,:)
    complex(8), intent(in) :: Hy_e(:,:)
    complex(8), intent(in) :: Hlambda_e(:,:)
    complex(8), intent(in) :: Zx_e(:,:) ! {Hx, z}
    complex(8), intent(in) :: Zy_e(:,:)
    complex(8), intent(in) :: Zlambda_e(:,:)

    real(8), intent(out) :: Fx
    real(8), intent(out) :: Fy

    integer :: iv, ic
    integer :: ilow ! first retained full band index
    integer :: ihigh ! last retained full band index
    integer :: ilast ! last occupied full band index
    integer :: nb
    integer :: nn, mm
    integer :: ii, jj ! local occ/empty index in cutoff matrices

    real(8) :: denominator 
    complex(8) :: term_x, term_y

    ! ----- Dimension check -----
    if (N < 2 .or. mod(N,2) /= 0) error stop 'CHI_FIXED_K_CUTOFF: invalid N'
    if (ncut < 0 .or. ncut > N/2-1) error stop 'CHI_FIXED_K_CUTOFF: invalid ncut'
    if (noccA < 0 .or. noccA > 2) error stop 'CHI_FIXED_K_CUTOFF: invalid noccA'

    nb = 2*ncut + 2

    if (any(shape(Hx_e) /= [nb, nb]) .or.  any(shape(Hy_e) /= [nb,nb]) .or. &
        any(shape(Hlambda_e) /= [nb,nb])) error stop 'CHI_FIXED_K_CUTOFF: H_e dimension mismatch'
    if (any(shape(Zx_e) /= [nb,nb]) .or. any(shape(Zy_e) /= [nb,nb]) .or. &
        any(shape(Zlambda_e) /= [nb,nb])) error stop 'CHI_FIXED_K_CUTOFF: Z_e dimension mismatch'

    ! ----- Active bands and symmetric cutoff window -----
    iv = N/2
    ic = iv+1

    ilow = iv-ncut
    ihigh = ic+ncut
    
    ilast = iv+noccA-1

    ! ----- 初始化 ----- 
    Fx = 0.d0
    Fy = 0.d0

    ! ----------
    ! 求和
    ! n:occupied, m:empty
    ! ----------
    do nn = ilow, ilast
        do mm = ilast + 1, ihigh

            ii = nn - ilow + 1
            jj = mm - ilow + 1

            denominator = W(nn) - W(mm)
            ! ----- 简并检查 ----- 
            if (abs(denominator) < 1.d-12*eV) then
                !$OMP CRITICAL(FIX_DIAGNOSTIC)
                write(*,'(A,2I8)') 'FIX degeneracy: nn, mm = ', nn, mm
                write(*,'(A,I8)') 'noccA = ', noccA
                write(*,'(A,ES24.14)') 'E_n (eV) = ', W(nn)/eV
                write(*,'(A,ES24.14)') 'E_m (eV) = ', W(mm)/eV
                write(*,'(A,ES24.14)') 'Delta E (eV) = ', denominator/eV
                !$OMP END CRITICAL(FIX_DIAGNOSTIC)
                error stop 'CHI_FIXED_K_CUTOFF: occupied-empty degeneracy'
            endif
            ! ----- beta = x ----- 
            term_x = Hx_e(ii,jj)*Zlambda_e(jj,ii) - Hlambda_e(ii,jj)*Zx_e(jj,ii)
            Fx = Fx + real(term_x/denominator,8)
            ! ----- beta = y ----- 
            term_y = Hy_e(ii,jj)*Zlambda_e(jj,ii) - Hlambda_e(ii,jj)*Zy_e(jj,ii)
            Fy = Fy + real(term_y/denominator,8)
        
        enddo
    enddo

end subroutine CHI_FIXED_K_CUTOFF


! ==========
! Integrate over the whole BZ
! chi_{FIX/FS}_x/y = \pm prefactor * int {F_x/y / G_x/y}        
! ==========
subroutine CALCULATE_CHI_CUTOFF(N, nk, ncut, b1, b2, mu, eta, prefactor, zdiag, ham_param, &
                                HAM, H_DERIVATIVES, &
                                chi_fix_x, chi_fix_y, chi_fs_x, chi_fs_y, &
                                chi_tot_x, chi_tot_y)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: nk
    integer, intent(in) :: ncut
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

    ! ----- cutoff incdices ----- 
    integer :: iv, ic
    integer :: ilow, ihigh
    integer :: nb
    integer :: jj ! matrix column index in the original BM basis
    integer :: noccA 
    ! ----- reciprocal lattice mesh coefficients ----- 
    integer :: ik1, ik2
    real(8) :: u, v
    real(8) :: k(2)
    ! ----- weight ----- 
    real(8) :: weight
    ! ----- susceptibility ----- 
    real(8) :: sum_Fx, sum_Fy
    real(8) :: sum_Gx, sum_Gy
    real(8) :: Fx, Fy
    real(8) :: Gx, Gy
    ! ----- Hamiltonian  ----- 
    real(8), allocatable :: W(:)
    complex(8), allocatable :: H(:,:)
    complex(8), allocatable :: V_HAM(:,:)
    ! ----- 提取截断部分对应的本征向量 ----- 
    complex(8), allocatable :: Vs(:,:)         ! N-by-nb retained eigenvector columns
    complex(8), allocatable :: Vs_dagger(:,:)
    ! ----- derivative matrix ----- 
    complex(8), allocatable :: dH_dkx(:,:), dH_dky(:,:), dH_dlambda(:,:)
    ! ----- derivative rep matrix ----- 
    complex(8), allocatable :: Hx_e(:, :), Hy_e(:, :), Hlambda_e(:, :)
    ! ----- anti-commutator matrix ----- 
    complex(8), allocatable :: Zx(:,:), Zy(:,:), Zlambda(:,:)
    ! ----- anti-commutator rep matrix ----- 
    complex(8), allocatable :: Zx_e(:,:), Zy_e(:,:), Zlambda_e(:,:)

    ! ----- Dimension check -----
    if (N < 2 .or. mod(N,2) /= 0) error stop 'CALCULATE_CHI_CUTOFF: invalid N'
    if (nk <= 0) error stop 'CALCULATE_CHI_CUTOFF: nk must be positive'
    if (eta <= 0.d0) error stop 'CALCULATE_CHI_CUTOFF: eta must be positive'
    if (ncut < 0 .or. ncut > N/2-1) error stop 'CALCULATE_CHI_CUTOFF: invalid ncut' 

    ! ----- Full band indices and retained window ----- 
    iv = N/2
    ic = iv + 1

    ilow = iv - ncut
    ihigh = ic + ncut

    nb = 2*ncut + 2 ! ihigh - ilow + 1

    ! ----- weight ----- 
    weight = 1.d0/real(nk*nk, 8)

    ! ----- 初始化 ----- 
    sum_Fx = 0.d0
    sum_Fy = 0.d0

    sum_Gx = 0.d0
    sum_Gy = 0.d0

    ! ==========
    ! 开始积分
    ! ==========
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(N,nk,ncut,iv,ic,ilow,ihigh,nb,b1,b2,mu,eta,zdiag,ham_param) &
    !$OMP PRIVATE(ik1,ik2,jj,noccA,u,v,k,Fx,Fy,Gx,Gy,W,H,V_HAM,dH_dkx,dH_dky,dH_dlambda, &
    !$OMP         Zx,Zy,Zlambda,Vs,Vs_dagger,Hx_e,Hy_e,Hlambda_e,Zx_e,Zy_e,Zlambda_e) &
    !$OMP REDUCTION(+:sum_Fx,sum_Fy,sum_Gx,sum_Gy) 
        
        allocate(W(N),H(N,N),V_HAM(N,N))
        allocate(dH_dkx(N,N),dH_dky(N,N),dH_dlambda(N,N))
        allocate(Zx(N,N),Zy(N,N),Zlambda(N,N))
        allocate(Vs(N,nb),Vs_dagger(nb,N))
        allocate(Hx_e(nb,nb),Hy_e(nb,nb),Hlambda_e(nb,nb))
        allocate(Zx_e(nb,nb),Zy_e(nb,nb),Zlambda_e(nb,nb)) 

        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
        do ik1 = 1, nk
            do ik2 = 1, nk
                ! -----Reduced coordinates-----
                u = -0.5d0 + real(ik1-1, 8)/real(nk, 8)
                v = -0.5d0 + real(ik2-1, 8)/real(nk, 8)
                ! ----- momentum ----- 
                k = u*b1 + v*b2
                ! ----- Hamiltonian and relevant ----- 
                call HAM(k,ham_param,H)
                call myzheev(N,W,H,V_HAM)
                call H_DERIVATIVES(k,ham_param,dH_dkx,dH_dky,dH_dlambda)

                ! ----------
                ! anticommutators
                ! {z, H_beta}_{ij} = (z_i + z_j)*H_beta_jj
                ! ----------
                do jj = 1, N
                    Zx(:,jj) = (zdiag+zdiag(jj))*dH_dkx(:,jj)
                    Zy(:,jj) = (zdiag+zdiag(jj))*dH_dky(:,jj)
                    Zlambda(:,jj) = (zdiag+zdiag(jj))*dH_dlambda(:,jj)
                enddo

                ! ----------
                ! 提取截断部分对应的本征向量
                ! 得到表示矩阵
                ! ----------
                Vs = V_HAM(:,ilow:ihigh)
                Vs_dagger = hermitian(Vs)

                Hx_e = matmul(Vs_dagger,matmul(dH_dkx,Vs))
                Hy_e = matmul(Vs_dagger,matmul(dH_dky,Vs))
                Hlambda_e = matmul(Vs_dagger,matmul(dH_dlambda,Vs))
                Zx_e = matmul(Vs_dagger,matmul(Zx,Vs))
                Zy_e = matmul(Vs_dagger,matmul(Zy,Vs))
                Zlambda_e = matmul(Vs_dagger,matmul(Zlambda,Vs))

                ! ----------
                ! 对化学势进行检查
                ! ----------
                if (ncut > 0) then
                    if (any(W(ilow:iv-1) >= mu)) &
                        error stop 'CALCULATE_CHI_CUTOFF: remote valence crosses mu'
                    if (any(W(ic+1:ihigh) < mu)) &
                        error stop 'CALCULATE_CHI_CUTOFF: remote conduction crosses mu'
                endif
                ! ----- 计算占据的active band的数目 ----- 
                noccA = count(W(iv:ic) < mu)

                ! ----------
                ! Single k point 
                ! ----------
                ! ----- FS ----- 
                call CHI_FS_K_CUTOFF(N,ncut,W,Hx_e,Hy_e,Hlambda_e,&
                                     Zx_e,Zy_e,Zlambda_e,mu,eta,Gx,Gy)
                ! ----- FIX ----- 
                call CHI_FIXED_K_CUTOFF(N,ncut,noccA,W,Hx_e,Hy_e,Hlambda_e, &
                                        Zx_e,Zy_e,Zlambda_e,Fx,Fy)

                ! ----------
                ! accumulation
                ! ----------
                sum_Fx = sum_Fx + Fx
                sum_Fy = sum_Fy + Fy
    
                sum_Gx = sum_Gx + Gx
                sum_Gy = sum_Gy + Gy
            
            enddo
        enddo
        !$OMP END DO

        ! -----
        ! 释放数组
        ! -----
        deallocate(W,H,V_HAM,dH_dkx,dH_dky,dH_dlambda,Zx,Zy,Zlambda,Vs,Vs_dagger)
        deallocate(Hx_e,Hy_e,Hlambda_e,Zx_e,Zy_e,Zlambda_e)

    !$OMP END PARALLEL

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

end subroutine CALCULATE_CHI_CUTOFF


! ==========
! Integrated the susceptibility(CUTOFF)
! ==========
subroutine M_FROM_CHI_CUTOFF(N, nk, ncut, b1, b2, delta, mu, prefactor, zdiag, &
                            ham_param_ref, HAM, lambda_list, chi_x, chi_y, &
                            Mx_from_chi, My_from_chi)
    implicit none
    
    integer, intent(in) :: N, nk, ncut
    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: delta, mu, prefactor
    real(8), intent(in) :: zdiag(N)
    real(8), intent(in) :: ham_param_ref(:)

    procedure(HAMILTONIAN) :: HAM

    real(8), intent(in) :: lambda_list(:)
    real(8), intent(in) :: chi_x(:), chi_y(:)

    real(8), intent(out) :: Mx_from_chi(:)
    real(8), intent(out) :: My_from_chi(:)

    integer :: i
    integer :: nlambda
    real(8) :: dlambda
    real(8) :: Mx_ref, My_ref
    
    ! ----- Array dimension checks -----
    nlambda = size(lambda_list)

    if (nlambda < 1) error stop 'M_FROM_CHI_CUTOFF: empty lambda_list'
    if (size(chi_x) /= nlambda .or. size(chi_y) /= nlambda) &
        error stop 'M_FROM_CHI_CUTOFF: input chi dimensions mismatch'
    if (size(Mx_from_chi) /= nlambda .or. size(My_from_chi) /= nlambda) &
        error stop 'M_FROM_CHI_CUTOFF: output M dimensions mismatch'

    ! ----------
    ! calculate the reference OM
    ! M_ref = M_{lambda_list(1)}
    ! ----------
    call CALCULATE_OM_CUTOFF(N,nk,ncut,b1,b2,delta,mu,prefactor,zdiag, &
    ham_param_ref,HAM,Mx_ref,My_ref)

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

end subroutine M_FROM_CHI_CUTOFF


end module SUSCEPTIBILITY_CUTOFF