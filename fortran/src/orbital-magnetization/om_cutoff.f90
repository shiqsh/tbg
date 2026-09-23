module OM_CUTOFF
    use CONSTANTS, only : PI
    use LINALG, only : myzheev
    use stdlib_linalg, only : diag
    use MODEL_INTERFACE, only : HAMILTONIAN
    use OM, only : BUILD_Z
    implicit none

contains

! ==========
! Symmetric truncated projectors
! 单个K点对应的P
! ==========
subroutine PROJECTOR_CUTOFF(N, H, mu, ncut, P, Q, W, noccA)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: ncut
    complex(8), intent(in) :: H(N, N)
    real(8), intent(in) :: mu
    complex(8), intent(out) :: P(N, N)
    complex(8), intent(out) :: Q(N, N)
    real(8), intent(out) :: W(N)
    integer, intent(out) :: noccA

    integer :: iv, ic
    integer :: ilow, ihigh
    integer :: iband
    complex(8) :: V(N, N)

    ! ----- dimension check ----- 
    if (mod(N,2) /=0 ) then
        error stop "PROJECTOR_CUTOFF: N must be even"
    endif

    ! ----- diagonalization ----- 
    call myzheev(N, W, H, V)

    ! ----- 标记两个active bands ----- 
    iv  = N/2
    ic = N/2 + 1
    ! ----- cutoff check ----- 
    if (ncut < 0) then
        error stop "PROJECTOR_CUTOFF: ncut must be positive"
    endif
    if (ncut > iv - 1 .or. ncut > N - ic) then
        error stop "PROJECTOR_CUTOFF: ncut is too large"
    endif
    ! ----- symmetric remote bands window ----- 
    ilow = iv - ncut
    ihigh = ic + ncut
    ! ----- 初始化 ----- 
    P = (0.d0, 0.d0)
    Q = (0.d0, 0.d0)

    ! =====
    ! 分开构建remote band 的 P 和 Q
    ! =====
    P = matmul(&
    V(:, ilow : iv - 1), transpose(conjg(V(:, ilow : iv - 1))) &
    )
    
    Q = matmul(&
    V(:, ic + 1 : ihigh), transpose(conjg(V(:, ic + 1 : ihigh))) &
    )

    ! =====
    ! Active bands
    ! 占据取决于每一个k点的化学势
    ! =====
    noccA = 0

    do iband = iv, ic
        if (W(iband) < mu) then
            P = P +  matmul(&
                V(:,iband:iband), transpose(conjg(V(:,iband:iband)))&
                )
            noccA = noccA + 1
        else
            Q = Q +  matmul(&
                V(:,iband:iband), transpose(conjg(V(:,iband:iband)))&
                )
        endif
    enddo 

end subroutine PROJECTOR_CUTOFF


! ==========
! 中心差分， under symmetric cutoff
! 分别计算dP/dx(y),dQ/dx(y)
! ==========
subroutine CENTRAL_DIFF_CUTOFF(N, k, h, mu, ncut, ham_param, HAM, &
                                dP_dx, dP_dy, dQ_dx, dQ_dy)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: ncut
    real(8), intent(in) :: k(2)
    real(8), intent(in) :: h, mu ! h 是差分步长
    real(8), intent(in) :: ham_param(:)
    complex(8), intent(out) :: dP_dx(N,N), dP_dy(N, N)
    complex(8), intent(out) :: dQ_dx(N,N), dQ_dy(N, N)

    integer :: noccA
    real(8) :: kp(2), km(2)
    real(8) :: W(N)
    complex(8) :: HAM_K(N,N)
    complex(8) :: P_plus(N,N), P_minus(N,N)
    complex(8) :: Q_plus(N,N), Q_minus(N,N)

    if (h <= 0.d0) error stop "CENTRAL_DIFF: h must be positive"

    ! ----------
    ! d/dx
    ! ----------
    kp = k
    km = k
    kp(1) = k(1) + h
    km(1) = k(1) - h
    ! ----- kx + h ----- 
    call HAM(kp, ham_param, HAM_K)
    call PROJECTOR_CUTOFF(N, HAM_K, mu, ncut, P_plus, Q_plus, W, noccA)
    ! ----- kx - h ----- 
    call HAM(km, ham_param, HAM_K)
    call PROJECTOR_CUTOFF(N, HAM_K, mu, ncut, P_minus, Q_minus, W, noccA)
    ! ----- 中心差分 ----- 
    dP_dx = (P_plus - P_minus)/(2.d0*h)
    dQ_dx = (Q_plus - Q_minus)/(2.d0*h)

    ! ----------
    ! d/dy
    ! ----------
    kp = k
    km = k
    kp(2) = k(2) + h
    km(2) = k(2) - h
    ! ----- kx + h ----- 
    call HAM(kp, ham_param, HAM_K)
    call PROJECTOR_CUTOFF(N, HAM_K, mu, ncut, P_plus, Q_plus, W, noccA)
    ! ----- kx - h ----- 
    call HAM(km, ham_param, HAM_K)
    call PROJECTOR_CUTOFF(N, HAM_K, mu, ncut, P_minus, Q_minus, W, noccA)
    ! ----- 中心差分 ----- 
    dP_dy = (P_plus - P_minus)/(2.d0*h)
    dQ_dy = (Q_plus - Q_minus)/(2.d0*h)

end subroutine CENTRAL_DIFF_CUTOFF









end module OM_CUTOFF
