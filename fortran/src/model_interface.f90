module MODEL_INTERFACE

    implicit none

    abstract interface

        subroutine HAMILTONIAN(k, ham_param, Ham)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)

            complex(8), intent(out) :: Ham(:,:)

        end subroutine HAMILTONIAN

    end interface

end module MODEL_INTERFACE