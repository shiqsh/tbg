module CONSTANTS
    implicit none

! ==================
! Mathematical constants
! ==================
    complex(8), parameter :: CI = (0.d0, 1.d0)
    complex(8), parameter :: RI = (1.d0, 0.d0)
    real(8), parameter :: PI   = acos(-1.d0)

! ==================
! Physical constants
! ==================
    real(8), parameter :: eV   = 1.602176634d-19
    real(8), parameter :: hbar = 1.054572663d-34

! ==================
! Pauli matrices
! ==================
    complex(8), parameter :: sigma_0(2,2) = reshape([(1.d0, 0.d0), (0.d0, 0.d0), (0.d0, 0.d0), (1.d0, 0.d0)], [2, 2])
    complex(8), parameter :: sigma_1(2,2) = reshape([(0.d0, 0.d0), (1.d0, 0.d0), (1.d0, 0.d0), (0.d0, 0.d0)], [2, 2])
    complex(8), parameter :: sigma_2(2,2) = reshape([(0.d0, 0.d0), (0.d0, 1.d0), (0.d0, -1.d0), (0.d0, 0.d0)], [2, 2])
    complex(8), parameter :: sigma_3(2,2) = reshape([(1.d0, 0.d0), (0.d0, 0.d0), (0.d0, 0.d0), (-1.d0, 0.d0)], [2, 2])

! ==================
! Graphene constants
! ==================
    real(8), parameter :: d_graphene   = 1.42d-10
    real(8), parameter :: a_graphene    = sqrt(3.d0)*d_graphene
    real(8), parameter :: vF   = 5.944d-10*eV ! 注意，这里的vF实际上是vF*hbbar, 也就是说最初的vF = 5.944d-10*eV/1.054572663d-34

! ------
! grapphene geometry
! -----
    real(8), parameter :: a1(2) = a_graphene*[1.d0/2.d0, sqrt(3.d0)/2.d0 ]
    real(8), parameter :: a2(2) = a_graphene*[-1.d0/2.d0, sqrt(3.d0)/2.d0 ]
    real(8), parameter :: b1(2) = (4.d0*PI/(sqrt(3.d0)*a_graphene))*[sqrt(3.d0)/2.d0, 1.d0/2.d0 ]
    real(8), parameter :: b2(2) = (4.d0*PI/(sqrt(3.d0)*a_graphene))*[-sqrt(3.d0)/2.d0, 1.d0/2.d0 ]
    real(8), parameter :: K0(2) = (4.d0*PI/(3.d0*a_graphene))*[1.d0, 0.d0]

! ============================================================
! BM model parameters
! ============================================================
    real(8), parameter :: w0 = 0.0544d0*eV
    real(8), parameter :: w1 = 0.1249d0*eV

! ----------
! off-diagonal term
! ----------
    complex(8), parameter :: T1(2, 2) = w0*sigma_0 + w1*sigma_1
    complex(8), parameter :: T2(2, 2) = w0*sigma_0 + w1*(cos(2.d0*PI/3.d0)*sigma_1 + sin(2.d0*PI/3.d0)*sigma_2)
    complex(8), parameter :: T3(2, 2) = w0*sigma_0 + w1*(cos(2.d0*PI/3.d0)*sigma_1 - sin(2.d0*PI/3.d0)*sigma_2)

! ------
! wavefunction truncation
! -----
    integer, parameter :: tr = 4
    integer, parameter :: NTBG  = 2*((2*tr + 1)**2)


end module CONSTANTS