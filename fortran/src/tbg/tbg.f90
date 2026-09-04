module TBG
    use CONSTANTS
    implicit none
                        
contains

! Singel Layer Hamiltonian of Graphene, We set the K^t as the new origin
subroutine HSLG(p, phi, M)
    implicit none

    real(8), intent(in) :: p(2), phi
    complex(8), intent(out) :: M(2,2)

    real(8) :: theta_k

    theta_k = atan2(p(2), p(1))

    M = vF*sqrt(p(1)**2 + p(2)**2)*reshape([(0.d0, 0.d0), exp( CI*(theta_k + phi)), &
        exp(-CI*(theta_k + phi)), (0.d0, 0.d0)], [2,2])

end subroutine HSLG


subroutine DIAGONAL(a, D)
    implicit none

    real(8), intent(in) :: a(2)
    complex(8), intent(out) :: D(2, 2)

    D =vF*(a(1)*sigma_1 + a(2)*sigma_2)

end subroutine DIAGONAL


subroutine OFFDIAGONAL(w0, w1, T1, T2, T3)
    implicit none 

    real(8), intent(in) :: w0, w1
    complex(8), intent(out) :: T1(2, 2), T2(2, 2), T3(2, 2)

    T1 = w0*sigma_0 + w1*sigma_1
    T2 = w0*sigma_0 + w1*(cos(2.d0*PI/3.d0)*sigma_1 + sin(2.d0*PI/3.d0)*sigma_2)
    T3 = w0*sigma_0 + w1*(cos(2.d0*PI/3.d0)*sigma_1 - sin(2.d0*PI/3.d0)*sigma_2)

end subroutine OFFDIAGONAL


subroutine VECTORS(theta, q1, q2, q3, b1m, b2m)
    implicit none

    real(8), intent(in) :: theta
    real(8), intent(out) :: q1(2), q2(2), q3(2)
    real(8), intent(out) :: b1m(2), b2m(2)

    real(8) :: q

    q = 8.d0*PI*sin(theta/2.d0)/(3.d0*a)

    q1 = q*[0.d0, -1.d0]
    q2 = q*[sqrt(3.d0)/2.d0, 1.d0/2.d0]
    q3 = q*[-sqrt(3.d0)/2.d0, 1.d0/2.d0]
    b1m = sqrt(3.d0)*q*[1.d0/2.d0, -sqrt(3.d0)/2.d0]
    b2m = sqrt(3.d0)*q*[1.d0/2.d0,  sqrt(3.d0)/2.d0]

end subroutine VECTORS


subroutine ELEMENTS(tH, alpha, n1, n2, beta, m1, m2, M)
    implicit none

    complex(8), intent(inout) :: tH(N,N,2,2)

    integer, intent(in) :: alpha, n1, n2
    integer, intent(in) :: beta, m1, m2

    complex(8), intent(in) :: M(2,2)

    integer :: i, j

    i = (alpha - 1)*((2*tr+1)**2) + (n1+tr)*(2*tr+1) + (n2+tr) + 1
    j = (beta  - 1)*((2*tr+1)**2) + (m1+tr)*(2*tr+1) + (m2+tr) + 1
    
    tH(i,j,:,:) = M

    if (i /= j) tH(j,i,:,:) = transpose(conjg(M))

end subroutine ELEMENTS


subroutine ASMBL(N, H, tH)
    implicit none
    integer, intent(in) :: N
    complex(8), intent(in) :: tH(N, N, 2, 2)
    complex(8), intent(out) :: H(2*N, 2*N)
    

    integer :: i, j

    do i = 1, N
        do j = 1, N
            H(2*i-1:2*i, 2*j-1:2*j) = tH(i, j, : ,:)
        enddo
    enddo

end subroutine ASMBL


subroutine HTBG(k, te, theta, H)
    implicit none
    
    real(8), intent(in) :: k(2), theta, te
    complex(8), intent(out) :: H(:,:) 
    complex(8), allocatable :: tH(:,:,:,:)

    real(8) :: p(2)
    real(8) :: b1m(2), b2m(2)
    real(8) :: q1(2), q2(2), q3(2) !qb, qtr, qtl
    integer :: layer1, layer2, n1, n2 ! n1, n2  ket index

    complex(8) :: T1(2, 2), T2(2, 2), T3(2, 2) ! off-diagonal term
    complex(8) :: M(2, 2)

    call VECTORS(theta, q1, q2, q3, b1m, b2m)
    call OFFDIAGONAL(w0, w1, T1, T2, T3)

    allocate(tH(N, N, 2, 2)) 
    tH = 0.d0

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(p, layer1,layer2,n1,n2,M) &
    !$OMP SHARED(k,te, q1, b1m, b2m,T1,T2,T3,tH) COLLAPSE(4) SCHEDULE(dynamic)
    do layer1 = 1, 2
        do layer2 = 1, 2
            do n1 = -tr, tr
                do n2 = -tr, tr
                ! diagonal term
                if  (layer1 ==1 .and. layer2 == 1) then ! bottom layer
                    p = k - q1 - n1*b1m - n2*b2m
                    call HSLG(p, te/2.d0, M)    
                    call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2, M)
                elseif (layer1 ==2 .and. layer2 == 2) then ! top layer
                    p = k - n1*b1m - n2*b2m
                    call HSLG(p, -te/2.d0, M)
                    call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2, M) 
                ! off-diagonal term
                elseif (layer1 == 1 .and. layer2 == 2) then
                    ! 三个最近耦合的非对角块（T1, T2, T3）
                    call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2,   T1)        ! Δ(0,0)
                    if (abs(n2 - 1) <= tr) call ELEMENTS(tH, layer1, n1, n2, layer2, n1,   n2-1, T2) ! Δ(0,-1)
                    if (abs(n1 + 1) <= tr) call ELEMENTS(tH, layer1, n1, n2, layer2, n1+1, n2,   T3) ! Δ(+1,0)
                end if
                end do
            end do
        end do                                                                                              
    end do
    !$OMP END PARALLEL DO 

    H = 0.d0 
    call ASMBL(N, H, tH) 
    deallocate(tH)

end subroutine HTBG 


subroutine KPATH(numk, theta, k, dk)
    implicit none

    integer, intent(in) :: numk
    real(8), intent(in) :: theta
    real(8), intent(out) :: k(2,3*numk+1)
    real(8), intent(out) :: dk(3*numk+1)

    integer :: i, j, n
    real(8) :: l, t
    real(8) :: P(2,4)

    l = 8.d0*PI*sin(theta/2.d0)/(3.d0*a)

    P(:,1) = [0.d0, -l]
    P(:,2) = [sqrt(3.d0)*l/2.d0, -l/2.d0]
    P(:,3) = [0.d0, -l/2.d0]
    P(:,4) = [0.d0, 0.d0]

    n = 1
    k(:,1) = P(:,1)

    do i = 1, 3
        do j = 1, numk

            t = real(j,8)/real(numk,8)

            n = n + 1
            k(:,n) = (1.d0 - t)*P(:,i) + t*P(:,i+1)

        enddo
    enddo

    dk(1) = 0.d0

    do i = 2, 3*numk+1
        dk(i) = dk(i-1) + sqrt(sum((k(:,i) - k(:,i-1))**2))
    enddo

end subroutine KPATH

!====================================================================
! 计算 Moiré BZ 线性尺度， l = |b1m| = (8π sin(θ/2)) / (sqrt{3}*a)
!====================================================================
subroutine MOIRE_BZ_SCALE(theta, l)
    implicit none

    real(8), intent(in) :: theta
    real(8), intent(out) :: l

    l = (8.d0*PI*sin(theta/2.d0))/(sqrt(3.d0)*a)

end subroutine MOIRE_BZ_SCALE


end module TBG
