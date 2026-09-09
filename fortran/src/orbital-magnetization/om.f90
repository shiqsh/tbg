module OM
    use CONSTANTS
    use LINALG, only : myzheev
    use stdlib_linalg, only : diag
    implicit none

    abstract interface
        subroutine HAMILTONIAN(k, ham_param, Ham)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)
            complex(8), intent(out) :: Ham(:,:)
        end subroutine HAMILTONIAN

    end interface
    
contains

! ==========
! Projector
! 单个k点对应的P
! ==========
subroutine PROJECTOR(N, H, mu, P, W, nocc)
    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: H(N, N)
    real(8), intent(in) :: mu
    complex(8), intent(out) :: P(N, N)
    real(8), intent(out) :: W(N)
    integer, intent(out) :: nocc

    complex(8) :: V(N, N)

    call myzheev(N, W, H, V)
    
    nocc = count(W < mu)

    P = (0.d0, 0.d0)

    if (nocc > 0) then
        P = matmul( & 
        V(:, 1: nocc), transpose(conjg(V(:, 1: nocc))) &
        ) 
    endif

end subroutine PROJECTOR 


! ==========
! 中心差分
! 在P_grid的基础上构造dP/dx, dP/dy
! =========
subroutine CENTRAL_DIFF(N, k, h, mu, ham_param, HAM, dP_dx, dP_dy)
    implicit none

    integer, intent(in) :: N
    real(8), intent(in) :: k(2)
    real(8), intent(in) :: h, mu ! h 差分步长
    real(8), intent(in) :: ham_param(:)
    procedure(HAMILTONIAN) :: HAM
    complex(8), intent(out) :: dP_dx(N,N)
    complex(8), intent(out) :: dP_dy(N,N)

    integer :: nocc
    real(8) :: kp(2), km(2)
    real(8) :: W(N)
    complex(8) :: HAM_K(N,N)
    complex(8) :: P_plus(N,N), P_minus(N,N)

    if (h <= 0.d0) error stop "CENTRAL_DIFF: h must be positive"


    ! ----------
    ! dP/dkx
    ! ----------
    kp = k
    km = k
    kp(1) = k(1) + h
    km(1) = k(1) - h
    ! P(kx+h, ky)
    call HAM(kp, ham_param, Ham_k)
    call PROJECTOR(N, Ham_k, mu, P_plus, W, nocc)
    ! P(kx-h, ky)
    call HAM(km, ham_param, Ham_k)
    call PROJECTOR(N, Ham_k, mu, P_minus, W, nocc)
    ! 中心差分
    dP_dx = (P_plus - P_minus)/(2.d0*h)
    ! ----------
    ! dP/dky
    ! ----------
    kp = k
    km = k
    kp(2) = k(2) + h
    km(2) = k(2) - h
    ! P(kx+h, ky)
    call HAM(kp, ham_param, Ham_k)
    call PROJECTOR(N, Ham_k, mu, P_plus, W, nocc)
    ! P(kx-h, ky)
    call HAM(km, ham_param, Ham_k)
    call PROJECTOR(N, Ham_k, mu, P_minus, W, nocc)
    ! 中心差分
    dP_dy = (P_plus - P_minus)/(2.d0*h)

end subroutine CENTRAL_DIFF


! ==========
! 构造mathcal{Z} = {H-mu, z}， 适用于z是对角阵
! ==========
subroutine BUILD_Z(H, zdiag, mu, Z)
    implicit none

    complex(8), intent(in) :: H(:,:)
    real(8), intent(in) :: zdiag(:)
    real(8), intent(in) :: mu
    complex(8), intent(out) :: Z(:,:)

    integer :: N, i, j

    N = size(H, 1)

     if (size(zdiag) /= N) then
        error stop "BUILD_Z: dimension error"
    endif

    do i = 1, N
        do j = 1, N
            Z(i,j) = (zdiag(i) + zdiag(j))*H(i,j)
        enddo
    enddo

    do i = 1, N
        Z(i,i) = Z(i,i) - 2.d0*mu*zdiag(i)
    enddo
    ! Z = matmul(z_+tbg, H) + matmul(H, z_tbg) - 2.d0*mu*z_tbg

end subroutine BUILD_Z


! ==========
! 计算轨道磁化
! ==========
subroutine CALCULATE_OM(N, nk, b1, b2, delta, mu, prefactor, zdiag, ham_param, Ham, &
                        Mx, My)
    implicit none

    integer, intent(in) :: N, nk
    real(8), intent(in) :: b1(2), b2(2)
    real(8), intent(in) :: delta, mu, prefactor
    real(8), intent(in) :: zdiag(N)
    real(8), intent(in) :: ham_param(:)
    procedure(HAMILTONIAN) :: HAM

    real(8), intent(out) :: Mx, My

    integer :: iu, iv, i
    integer :: nocc
    real(8) :: u, v
    real(8) :: k(2)
    real(8) :: detB, weight
    real(8) :: W(N)
    real(8) :: Fx, Fy
    real(8) :: sumFx, sumFy

    complex(8) :: H(N, N)
    complex(8) :: P(N,N), Q(N, N)
    complex(8) :: dP_dx(N, N), dP_dy(N, N)
    complex(8) :: Z(N, N)
    complex(8) :: A(N, N), B(N, N)

    ! -----
    ! Reciprocal-cell area，|b1 x b2|
    ! -----
    detB = b1(1)*b2(2) - b1(2)*b2(1)

    if (abs(detB) < 1.d-14) then
        error stop "CALCULATE_OM: invalid reciprocal lattice"
    endif

    ! -----
    ! integration weight
    ! -----
    weight = abs(detB)/(2.d0*PI)**2/real(nk*nk,8)

    sumFx = 0.d0
    sumFy = 0.d0

    ! ----------
    ! INTEGRATION
    ! ----------
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(iu,iv,i,u,v,k,H,P,Q,dP_dx,dP_dy,Z,A,B,W,nocc,Fx,Fy) &
    !$OMP SHARED(N,nk,b1,b2,delta,mu,prefactor,zdiag,ham_param,weight) &
    !$OMP REDUCTION(+:sumFx,sumFy) COLLAPSE(2) SCHEDULE(DYNAMIC)
    do iu = 1, nk
        do iv = 1, nk
            ! -----Reduced coordinates-----
            u = -0.5d0 + real(iu-1, 8)/real(nk, 8)
            v = -0.5d0 + real(iv-1, 8)/real(nk, 8)
            ! -----Cartesian momentum-----
            k = u*b1 + v*b2
            ! -----Hamiltonian------
            call HAM(k, ham_param, H)
            ! -----projector-----
            call PROJECTOR(N, H, mu, P, W, nocc)
            Q = -P
            do i = 1, N
                Q(i,i) = Q(i,i) + 1.d0
            enddo
            ! -----dP/dx, dP/dy-----
            call CENTRAL_DIFF(N, k, delta, mu, ham_param, HAM, dP_dx, dP_dy)

            call BUILD_Z(H, zdiag, mu, Z)

            ! -----
            ! Fx = Re Tr [P(dP/dx) Q Z]
            ! Fy = Re Tr [P(dP/dy) Q Z]
            ! -----
            A = matmul(matmul(matmul(P, dP_dx), Q), Z)
            Fx = sum(real(diag(A)))

            B = matmul(matmul(matmul(P, dP_dy), Q), Z)
            Fy = sum(real(diag(B)))
            
            sumFx = sumFx + Fx
            sumFy = sumFy + Fy

        enddo
    enddo
    !$OMP END PARALLEL DO

    ! -----
    ! epsilon_xy = +1
    ! -----
    Mx = -prefactor*weight*sumFy
    My = prefactor*weight*sumFx

end subroutine CALCULATE_OM

    
end module OM 