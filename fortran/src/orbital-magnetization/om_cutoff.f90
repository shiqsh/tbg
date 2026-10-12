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
    complex(8), allocatable :: V(:,:)

    ! ----- dimension check ----- 
    if (mod(N,2) /=0 ) then
        error stop "PROJECTOR_CUTOFF: N must be even"
    endif

    allocate(V(N,N))
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

    deallocate(V)

end subroutine PROJECTOR_CUTOFF


! ==========
! 中心差分， under symmetric cutoff
! Geometric central difference
! ==========
subroutine CENTRAL_DIFF_CUTOFF(N, k, h, ncut, noccA, ham_param, HAM, &
                                dP_dx, dP_dy)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: ncut
    integer, intent(in) :: noccA ! 这里本质上是通过化学势决定的
    real(8), intent(in) :: k(2)
    real(8), intent(in) :: h ! h 是差分步长
    real(8), intent(in) :: ham_param(:)
    procedure(HAMILTONIAN) :: HAM
    complex(8), intent(out) :: dP_dx(N,N), dP_dy(N, N)

    integer :: iv, ilow, ilast, idir
    real(8) :: kp(2), km(2)

    real(8), allocatable :: W(:)
    complex(8), allocatable :: HAM_K(:,:), V(:,:)
    complex(8), allocatable :: P_plus(:,:), P_minus(:,:)

    if (h <= 0.d0) error stop "CENTRAL_DIFF: h must be positive"

    ! ----------
    ! band index
    ! ----------
    iv  = N/2
    
    if (ncut < 0 .or. ncut > iv - 1) then
        error stop "CENTRAL_DIFF: invalid ncut"
    endif
    if (noccA < 0 .or. noccA > 2) then
        error stop "CENTRAL_DIFF: invalid noccA"
    endif
    ! ----- 占据态截断最低能带指标 ----- 
    ilow = iv - ncut
    ! ----- 占据态最高能带指标 ----- 
    ilast = iv + noccA - 1

    ! ----- 确定动态数组的大小 ----- 
    allocate(W(N))
    allocate(HAM_K(N,N), V(N,N), P_plus(N,N), P_minus(N,N))

    ! ----------
    ! 在两个方向上进行中心差分
    ! 同时要构造新的投影算符
    ! ----------
    do idir = 1, 2

        kp = k
        km = k
        kp(idir) = k(idir) + h
        km(idir) = k(idir) - h
        ! ----- k + h  ----- 
        call HAM(kp, ham_param, HAM_K)
        call myzheev(N, W, HAM_K, V)
        P_plus = (0.d0, 0.d0)
        if (ilast >= ilow) then
            P_plus = matmul(&
                V(:, ilow : ilast), transpose(conjg(V(:, ilow : ilast))) &
                )
        endif
        ! ----- k - h ----- 
        call HAM(km, ham_param, HAM_K)
        call myzheev(N, W, HAM_K, V)
        P_minus = (0.d0, 0.d0)
        if (ilast >= ilow) then
            P_minus = matmul(&
                V(:, ilow : ilast), transpose(conjg(V(:, ilow : ilast))) &
                )
        endif

        ! ----- 几何微分 ----- 
        if (idir == 1) then
            dP_dx = (P_plus - P_minus)/(2.d0*h)
        else
            dP_dy = (P_plus - P_minus)/(2.d0*h)
        endif   

    enddo

    deallocate(W, HAM_K, V, P_plus, P_minus)

end subroutine CENTRAL_DIFF_CUTOFF


! ==========
! 计算轨道磁化
! 输出结果为单位原胞内的轨道磁化强度，单位是波尔磁子（mu_B = e hbar/2 m_2)
! ==========
subroutine CALCULATE_OM_CUTOFF(N, nk, ncut, b1, b2, delta, mu, prefactor, zdiag, ham_param, HAM, &
                                Mx, My)
    implicit none

    integer, intent(in) :: N, nk, ncut
    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: delta, mu, prefactor ! delta 微分步长
    real(8), intent(in) :: zdiag(N)
    real(8), intent(in) :: ham_param(:)
    procedure(HAMILTONIAN) :: HAM

    real(8), intent(out) :: Mx, My

    integer :: iu, iv
    integer :: noccA

    real(8) :: u, v
    real(8) :: k(2)
    real(8) :: weight
    real(8) :: Fx, Fy
    real(8) :: sumFx, sumFy

    real(8), allocatable :: W(:)

    complex(8), allocatable :: H(:,:)
    complex(8), allocatable :: P(:,:), Q(:,:)
    complex(8), allocatable :: dP_dx(:,:), dP_dy(:,:)
    complex(8), allocatable :: Z(:,:)
    complex(8), allocatable :: A(:,:), B(:,:)


    ! ----- input check ----- 
    if (nk <= 0) then 
            error stop "CALCULATE_OM_CUTOFF: nk must be positive"
    endif
    if (delta <= 0.d0) then
            error stop "CALCULATE_OM_CUTOFF: delta must be positive"
    endif

    ! -----
    ! integration weight
    ! -----
    weight = 1/real(nk*nk,8)

    ! ----- initialize ----- 
    sumFx = 0.d0
    sumFy = 0.d0    

    ! ----------
    ! INTEGRATION
    ! ----------
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(iu,iv,u,v,k,H,P,Q,dP_dx,dP_dy,Z,W,A,B,noccA,Fx,Fy) &
    !$OMP SHARED(N, nk, ncut, b1, b2, delta, mu, prefactor, zdiag, ham_param, weight) &
    !$OMP REDUCTION(+:sumFx, sumFy) 

        ! ----- 每个线程分别分配自己的工作数组 ----- 
        allocate(W(N))
        allocate(H(N,N), P(N,N), Q(N,N))
        allocate(dP_dx(N,N), dP_dy(N,N), Z(N,N))
        allocate(A(N,N), B(N,N))
    
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
        do iu = 1, nk
            do iv = 1, nk
                ! ----- reduced coordinates ----- 
                u = -0.5d0 + real(iu-1, 8)/real(nk, 8)
                v = -0.5d0 + real(iv-1, 8)/real(nk, 8)
                ! ----- momentum ----- 
                k = u*b1 + v*b2
                ! ----- Hamiltonian ----- 
                call HAM(k, ham_param, H)
                ! ----- projector ----- 
                call PROJECTOR_CUTOFF(N, H, mu, ncut, P, Q, W, noccA)
                ! ----- central difference -----
                call CENTRAL_DIFF_CUTOFF(N, k, delta, ncut, noccA, ham_param, HAM, &
                                        dP_dx, dP_dy)   
                ! ----- anti-commutator Z  ----- 
                call BUILD_Z(H, zdiag, mu, Z)

                ! -----
                ! Fx = Re Tr [P (dP/dx) Q Z]
                ! Fy = Re Tr [P (dP/dy) Q Z]
                ! -----
                A = matmul(matmul(matmul(P, dP_dx), Q), Z)
                Fx = sum(real(diag(A)))

                B = matmul(matmul(matmul(P, dP_dy), Q), Z)
                Fy = sum(real(diag(B)))

                ! ----- accumulation ----- 
                sumFx = sumFx + Fx
                sumFy = sumFy + Fy
            enddo
        enddo
        !$OMP END DO

        ! 每个线程分别释放自己的工作数组
        deallocate(W, H, P, Q, dP_dx, dP_dy, Z, A, B)

    !$OMP END PARALLEL

    ! -----
    ! epsilon_xy = +1
    ! prefactor = 2 me/hbar^2
    ! -----
    Mx = -prefactor*weight*sumFy
    My = prefactor*weight*sumFx 

end subroutine CALCULATE_OM_CUTOFF

end module OM_CUTOFF
