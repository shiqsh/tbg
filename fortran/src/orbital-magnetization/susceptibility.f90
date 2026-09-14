module SUSCEPTIBILITY
    use CONSTANTS, only : PI
    use LINALG, only : myzheev
    use stdlib_linalg, only : diag
    implicit none

    abstract interface
        ! ----- Hamiltonian interface ----- 
         subroutine HAMILTONIAN(k, ham_param, Ham)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)
            complex(8), intent(out) :: Ham(:,:)
        end subroutine HAMILTONIAN
        ! ----- Hamiltonian derivative interface ----- 
        subroutine H_DERIVATIVE_INTERFACE(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
            implicit none

            real(8), intent(in) :: k(2)
            real(8), intent(in) :: ham_param(:)
            complex(8), intent(out) :: dH_dkx(:,:)
            complex(8), intent(out) :: dH_dky(:,:)
            complex(8), intent(out) :: dH_dlambda(:,:)
        end subroutine H_DERIVATIVE_INTERFACE

    end interface

contains

! ==========
! Gaussian function
! delta_eta(E - mu) = exp[- ((E - m u)/eta)^2]/(sqrt(pi)*eta)
! ==========
subroutine GAUSSIAN_DELTA(E, mu, eta, delta)
    implicit none

    real(8), intent(in) :: E
    real(8), intent(in) :: mu
    real(8), intent(in) :: eta

    if (eta <= 0.d0) then
        error stop "DELTA_GAUSSIAN: eta must be positive"
    endif

    delta = exp(-((E-mu)/eta)**2)/(sqrt(PI)*eta)

end subroutine GAUSSIAN_DELTA


! ==========
! FIX contribution at one k point
! ==========



            






end module SUSCEPTIBILITY
