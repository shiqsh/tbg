module BILAYER_QWZ
    use CONSTANTS
    use TBG, only : ASMBL
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
        tH(2,2,:,:) = sin(kb(1))*sigma_1 + sin(kb(2))*sigma_2 + (m + cos(kt(1)) + cos(kb(2)))*sigma_3
        tH(1,2,:,:) = t_perp*sigma_0
        tH(2,1,:,:) = t_perp*sigma_0

        call ASMBL(2, H, tH)

    end subroutine H_BILAYER_QWZ


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
    

end module BILAYER_QWZ
