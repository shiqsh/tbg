module MODEL_INTERFACE

    implicit none

    abstract interface

        subroutine HAMILTONIAN(k, ham_param, HAM)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)

            complex(8), intent(out) :: HAM(:,:)

        end subroutine HAMILTONIAN

        ! ----- Hamiltonian derivative interface ----- 
        subroutine H_DERIVATIVES_INTERFACE(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)
            complex(8), intent(out) :: dH_dkx(:,:)
            complex(8), intent(out) :: dH_dky(:,:)
            complex(8), intent(out), optional :: dH_dlambda(:,:)
        end subroutine H_DERIVATIVES_INTERFACE

    end interface

end module MODEL_INTERFACE