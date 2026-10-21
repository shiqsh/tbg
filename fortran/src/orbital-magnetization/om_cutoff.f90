module OM_CUTOFF
    use CONSTANTS, only : PI
    use LINALG, only : myzheev
    use stdlib_linalg, only : diag
    use MODEL_INTERFACE, only : HAMILTONIAN， H_DERIVATIVES_INTERFACE
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

end subroutine 

! ==========
! 通过本征态来直接构造投影算符
! 这适用于哈密顿量不变的情况
! ==========
subroutine PROJECTORS_CUTOFF_FROM_V(N, V, ncut, noccA, P, Q)
    implicit none
    
    integer, intent(in) :: N
    complex(8), intent(in) :: V(:,:)
    integer, intent(in) :: ncut, noccA

    complex(8), intent(out) :: P(:,:), Q(:,:)

    integer :: iv, ic                       ! Two active-band indices
    integer :: ilow, ihigh, ilast

    ! ----- diemension check ----- 
    if (N < 2 .or. mod(N,2) /= 0) error stop 'PROJECTORS_CUTOFF_FROM_V: invalid N'

    ! -----
    ! band indices
    ! -----
    iv = N/2
    ic = iv + 1
    
    if (ncut < 0 .or. ncut > iv-1) error stop 'PROJECTORS_CUTOFF_FROM_V: invalid ncut'
    if (noccA < 0 .or. noccA > 2) error stop 'PROJECTORS_CUTOFF_FROM_V: invalid noccA'

    ! ----- Symmetric cutoff band indices -----
    ilow = iv-ncut
    ihigh = ic+ncut
    ilast = iv+noccA-1

    ! ----- initialize ----- 
    P = 0.d0
    Q = 0.d0

    ! -----
    ! occupied projector
    ! -----
    if (ilow <= ilast) then
        P = matmul(V(:,ilow:ilast), conjg(transpose(V(:,ilow:ilast))))
    endif
    ! -----
    ! empty projector
    ! -----
    if (ilast+1 <= ihigh) then
        Q = matmul(V(:,ilast+1:ihigh), conjg(transpose(V(:,ilast+1:ihigh))))
    endif

end subroutine PROJECTORS_CUTOFF_FROM_V


! ==========
! 计算哈密顿量固定的轨道磁化，只对化学势mu进行扫描
! 只做一次对角化，不要每个化学势变化一下就对角化
! ==========
subroutine CALCULATE_OM_CUTOFF_MU(N, nk, ncut, b1, b2, delta, mu_list, prefactor, zdiag, ham_param, HAM, &
                                Mx_list, My_list)
    implicit none
    
    integer, intent(in) :: N, nk, ncut
    real(8), intent(in) :: b1(2), b2(2), delta, mu_list(:), prefactor
    real(8), intent(in) :: zdiag(N),  ham_param(:)
    procedure(HAMILTONIAN) :: HAM   
    
    real(8), intent(out) :: Mx_list(size(mu_list)), My_list(size(mu_list))

    integer :: ik1, ik2
    integer :: imu
    integer :: iv, ic
    integer :: noccA

    real(8) :: u, v
    real(8) :: k(2), kp(2), km(2)
    real(8) :: weight
    real(8) :: mu
    real(8) :: Fx, Fy

    real(8), allocatable :: W(:), Wshift(:)
    real(8), allocatable :: sumFx(:), sumFy(:)

    complex(8), allocatable :: H0(:,:)
    complex(8), allocatable :: Hshift(:,:)

    complex(8), allocatable :: V0(:,:)  ! eigenvectors at central k
    complex(8), allocatable :: Vxp(:,:), Vxm(:,:)
    complex(8), allocatable :: Vyp(:,:), Vym(:,:)
    
    complex(8), allocatable :: P0(:,:), Pxp(:,:), Pxm(:,:), Pyp(:,:), Pym(:,:)
    complex(8), allocatable :: Q0(:,:), Qxp(:,:), Qxm(:,:), Qyp(:,:), Qym(:,:)

    complex(8), allocatable :: dP_dx(:,:), dP_dy(:,:)
    complex(8), allocatable :: Z(:,:)
    complex(8), allocatable :: A(:,:), B(:,:) ! P(dP)Q

    ! -----
    ! Dimension check
    ! -----
    if (N < 2 .or. mod(N,2) /= 0) error stop 'OM_MU: invalid N'
    if (nk <= 0 .or. delta <= 0.d0) error stop 'OM_MU: invalid nk or delta'
    if (ncut < 0 .or. ncut > N/2-1) error stop 'OM_MU: invalid ncut'
    if (size(mu_list) == 0) error stop 'OM_MU: empty mu_list'

    ! ----------
    ! Two active bands
    ! Symmetric remote bands 
    ! ----------
    iv = N/2
    ic = N/2 + 1

    ! -----
    ! 积分权重
    ! -----
    weight = 1.d0/real(nk*nk, 8)

    ! ----------
    ! allocate the final arrays
    ! ----------
    allocate(sumFx(size(mu_list)), sumFy(size(mu_list)))

    ! ----- 初始化 ----- 
    sumFx = 0.d0
    sumFy = 0.d0

    ! ==========
    ! INTEGRATION
    ! ==========
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(N,nk,ncut,iv,ic,b1,b2,delta,mu_list,zdiag,ham_param) &
    !$OMP PRIVATE(ik1,ik2,imu,noccA,u,v,k,kp,km,mu,Fx,Fy) &
    !$OMP PRIVATE(W,Wshift,H0,Hshift,V0,Vxp,Vxm,Vyp,Vym) &
    !$OMP PRIVATE(P0,Pxp,Pxm,Pyp,Pym,Q0,Qxp,Qxm,Qyp,Qym) &
    !$OMP PRIVATE(dP_dx,dP_dy,Z,A,B) &
    !$OMP REDUCTION(+:sumFx,sumFy)

        ! ----- Allocate workspace once for each thread -----
        allocate(W(N), Wshift(N))
        allocate(H0(N,N), Hshift(N,N))

        allocate(V0(N,N), Vxp(N,N), Vxm(N,N), Vyp(N,N), Vym(N,N))

        allocate(P0(N,N), Q0(N,N))
        allocate(Pxp(N,N), Qxp(N,N), Pxm(N,N), Qxm(N,N))
        allocate(Pyp(N,N), Qyp(N,N), Pym(N,N), Qym(N,N))

        allocate(dP_dx(N,N), dP_dy(N,N), Z(N,N), A(N,N), B(N,N))

        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
        do ik1 = 1, nk
            do ik2 = 1, nk
                ! ----- momentum ----- 
                u = -0.5d0 + real(ik1-1,8)/real(nk,8)
                v = -0.5d0 + real(ik2-1,8)/real(nk,8)
                k = u*b1 + v*b2  
                ! ----- Central Hamiltonian and diagonalization ----- 
                call HAM(k, ham_param, H0)
                call myzheev(N, W, H0, V0)
                ! ----- Diagonalization at kx pm delta ----- 
                kp = k
                km = k
                kp(1) = kp(1)+delta
                km(1) = km(1)-delta

                call HAM(kp, ham_param, Hshift)
                call myzheev(N, Wshift, Hshift, Vxp)

                call HAM(km, ham_param, Hshift)
                call myzheev(N, Wshift, Hshift, Vxm)

                ! ----- Diagonalization at ky pm delta ----- 
                kp = k
                km = k
                kp(2) = kp(2)+delta
                km(2) = km(2)-delta

                call HAM(kp, ham_param, Hshift)
                call myzheev(N, Wshift, Hshift, Vyp)

                call HAM(km, ham_param, Hshift)
                call myzheev(N, Wshift, Hshift, Vym)

                ! -----
                ! 扫描化学势
                ! -----
                do imu = 1, size(mu_list)

                    mu = mu_list(imu)

                    ! ----- 确定占据能带数量 ----- 
                    noccA = count(W(iv:ic) < mu)
                    ! ----- 构造五个点的投影算符 ----- 
                    call PROJECTORS_CUTOFF_FROM_V(N, V0, ncut, noccA, P0, Q0)
                    call PROJECTORS_CUTOFF_FROM_V(N, Vxp, ncut, noccA, Pxp, Qxp)
                    call PROJECTORS_CUTOFF_FROM_V(N, Vxm, ncut, noccA, Pxm, Qxm)
                    call PROJECTORS_CUTOFF_FROM_V(N, Vyp, ncut, noccA, Pyp, Qyp)
                    call PROJECTORS_CUTOFF_FROM_V(N, Vym, ncut, noccA, Pym, Qym)

                    ! ----- Geometric derivatives of P ----- 
                    dP_dx = (Pxp-Pxm)/(2.d0*delta)
                    dP_dy = (Pyp-Pym)/(2.d0*delta)

                    ! ----- Anticommutator  ----- 
                    call BUILD_Z(H0, zdiag, mu, Z)

                    ! -----
                    ! Fx = Re Tr [P (dP/dx) Q Z]
                    ! Fy = Re Tr [P (dP/dy) Q Z]
                    ! -----
                    A = matmul(matmul(matmul(P0, dP_dx), Q0), Z)
                    Fx = sum(real(diag(A)))

                    B = matmul(matmul(matmul(P0, dP_dy), Q0), Z)
                    Fy = sum(real(diag(B)))

                    sumFx(imu) = sumFx(imu) + Fx
                    sumFy(imu) = sumFy(imu) + Fy
                enddo

            enddo
        enddo
        !$OMP END DO 

        ! ----- 释放数组 ----- 
        deallocate(W,Wshift,H0,Hshift,V0,Vxp,Vxm,Vyp,Vym)
        deallocate(P0,Q0,Pxp,Qxp,Pxm,Qxm,Pyp,Qyp,Pym,Qym)
        deallocate(dP_dx,dP_dy,Z,A,B)

    !$OMP END PARALLEL

    Mx_list = -prefactor*weight*sumFy
    My_list =  prefactor*weight*sumFx

end subroutine CALCULATE_OM_CUTOFF_MU


! ==========
! 直接通过本征态求轨道磁化
! ==========
subroutine CALCULATE_OM_CUTOFF_KET(N, nk, ncut, b1, b2, mu, prefactor, zdiag, ham_param, HAM, &
                                    Mx, My)
    implicit none
    
    integer, intent(in) :: N
    integer, intent(in) :: nk
    integer, intent(in) :: ncut

    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: mu 
    real(8), intent(in) :: prefactor
    real(8), intent(in) :: zdiag(N)
    real(8), intent(in) :: ham_param(:)

    procedure(HAMILTONIAN) :: HAM
    procedure(H_DERIVATIVES_INTERFACE): H_DERIVATIVES

    real(8), intent(out) :: Mx, My

    integer :: iv, ic
    integer :: ilow, ihigh, ilast
    integer :: nb, noccA

    ! ----- loop indices  ----- 
    integer :: ik1, ik2
    integer :: iband, jband ! Full band indices
    integer :: i, j, icol ! local band indices

    real(8) :: u, v
    real(8) :: k(2), kref(2)
    real(8) :: weight
    real(8) :: denominator
    real(8) :: gap_tol ! degeneracy threshold
    real(8) :: factor ! En + Em - 2 mu / (En - Em)
    real(8) :: Fx, Fy
    real(8) :: sumFx, sumFy  

    real(8), allocatable :: W(:)
    complex(8), allocatable :: dH_dkx(:,:), dH_dky(:,:)
    complex(8), allocatable :: H0(:,:)
    complex(8), allocatable :: V0(:,:), Vcut(:,:)
    complex(8), allocatable :: Temp(:,:) ! matrix  workspace 
    
    ! ----- rep matrix ----- 
    complex(8), allocatable :: Hx_e(:,:)   ! <n|dH/dkx|m>
    complex(8), allocatable :: Hy_e(:,:)   ! <n|dH/dky|m>
    complex(8), allocatable :: z(:,:)







end module OM_CUTOFF
