module BILAYER_QWZ
    use CONSTANTS
    use TBG, only : ASMBL
    use stdlib_linalg, only: kronecker_product, diag
    implicit none

contains

subroutine H_BILAYER_QWZ(k, lambda, m, t_perp, H)
    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: lambda, m
    real(8), intent(in) :: t_perp
    complex(8), intent(out) :: H(4, 4)

    real(8) :: kt(2), kb(2)
    complex(8) :: tH(2,2,2,2)

    kt = [k(1) + lambda, k(2)]
    kb = [k(1) - lambda, k(2)]

    tH(1,1,:,:) = sin(kt(1))*sigma_1 + sin(kt(2))*sigma_2 + (m + cos(kt(1)) + cos(kt(2)))*sigma_3
    tH(2,2,:,:) = sin(kb(1))*sigma_1 + sin(kb(2))*sigma_2 + (m + cos(kb(1)) + cos(kb(2)))*sigma_3
    tH(1,2,:,:) = t_perp*sigma_0
    tH(2,1,:,:) = t_perp*sigma_0

    call ASMBL(2, H, tH)

end subroutine H_BILAYER_QWZ


! ==========
! Bilayer QWZ Hamiltonian interface
!
! ham_param(1) = lambda
! ham_param(2) = m
! ham_param(3) = t_perp
! ==========
subroutine HAM_QWZ(k, ham_param, Ham)

    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: ham_param(:)

    complex(8), intent(out) :: Ham(:,:)

    real(8) :: lambda, m, t_perp


    if (size(ham_param) /= 3) then
        error stop "HAM_QWZ: ham_param dimension error"
    endif


    lambda = ham_param(1)
    m      = ham_param(2)
    t_perp = ham_param(3)


    call H_BILAYER_QWZ(k, lambda, m, t_perp, Ham)


end subroutine HAM_QWZ


! ==========
! Bilayer QWZ model derivatives 
! ==========
subroutine H_DERIVATIVES_QWZ(k, lambda, dH_dkx, dH_dky, dH_dlambda)
    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: lambda

    complex(8), intent(out) :: dH_dkx(:,:)
    complex(8), intent(out) :: dH_dky(:,:)
    complex(8), intent(out) :: dH_dlambda(:,:)

    real(8) :: kt(2), kb(2)
    complex(8) :: tHx(2,2,2,2), tHy(2,2,2,2), tHlambda(2,2,2,2)

    kt = [k(1) + lambda, k(2)]
    kb = [k(1) - lambda, k(2)]

    tHx = (0.d0, 0.d0)
    tHy = (0.d0, 0.d0)
    tHlambda = (0.d0, 0.d0)

    ! -----
    ! dH/dkx
    ! -----
    tHx(1,1,:,:) = cos(kt(1))*sigma_1 - sin(kt(1))*sigma_3
    tHx(2,2,:,:) = cos(kb(1))*sigma_1 - sin(kb(1))*sigma_3
    ! -----
    ! dH/dky
    ! -----
    tHy(1,1,:,:) = cos(kt(2))*sigma_2 - sin(kt(2))*sigma_3
    tHy(2,2,:,:) = cos(kb(2))*sigma_2 - sin(kb(2))*sigma_3
    ! -----
    ! dH/dlambda
    ! -----
    tHlambda(1,1,:,:) = cos(kt(1))*sigma_1 - sin(kt(1))*sigma_3
    tHlambda(2,2,:,:) = -cos(kb(1))*sigma_1 + sin(kb(1))*sigma_3

    call ASMBL(2, dH_dkx, tHx)
    call ASMBL(2, dH_dky, tHy)
    call ASMBL(2, dH_dlambda, tHlambda)
    
end subroutine H_DERIVATIVES_QWZ


! ==========
! Bilayer QWZ derivative Hamiltonian interface 
! ==========
subroutine HAM_DERIVATIVES_QWZ(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: ham_param(:)

    complex(8), intent(out) :: dH_dkx(:,:), dH_dky(:,:), dH_dlambda(:,:)

    real(8) :: lambda, m, t_perp

    if (size(ham_param) /= 3) then
        error stop "HAM_QWZ: ham_param dimension error"
    endif


    lambda = ham_param(1)
    m      = ham_param(2)
    t_perp = ham_param(3)

    call H_DERIVATIVES_QWZ(k, lambda, dH_dkx, dH_dky, dH_dlambda)

end subroutine HAM_DERIVATIVES_QWZ


subroutine KPATH_QWZ(numk, k, dk)
    implicit none

    integer, intent(in) :: numk
    real(8), intent(out) :: k(2, 3*numk+1)
    real(8), intent(out) :: dk(3*numk+1)

    integer :: i, j, ik
    real(8) :: t
    real(8) :: P(2,4)

    ! Gamma -> X -> M -> Gamma
    P(:,1) = [0.d0, 0.d0]
    P(:,2) = [PI,    0.d0]
    P(:,3) = [PI,    PI]
    P(:,4) = [0.d0, 0.d0]

    ik = 1
    k(:,1) = P(:,1)

    do i = 1, 3
        do j = 1, numk

            t = real(j,8)/real(numk,8)

            ik = ik + 1
            k(:,ik) = P(:,i) + t*(P(:,i+1) - P(:,i))
        enddo
    enddo

    dk(1) = 0.d0

    do i = 2, 3*numk+1
        dk(i) = dk(i-1) + sqrt(sum((k(:,i) - k(:,i-1))**2))
    enddo

end subroutine KPATH_QWZ


! ==========
! z方向上投影算符
! ==========
subroutine Z_BQWZ(d, z)
    implicit none

    real(8), intent(in) :: d
    real(8), intent(out) :: z(:)

    z(1) = d/2.d0;    z(2) = d/2.d0
    z(3) = -d/2.d0;    z(4) = -d/2.d0

    ! z = diag((d/2.d0) * kronecker_product(sigma_3, sigma_0))

end subroutine Z_BQWZ
    

end module BILAYER_QWZ
