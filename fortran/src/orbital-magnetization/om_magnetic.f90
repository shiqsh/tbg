module OM_MAGNETIC
    use CONSTANTS, only : hbar, m_e, NTBG
    use STRAIN, only : HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: CALCULATE_OMEGA_B, CALCULATE_OM_MAGNETIC
    real(8), parameter :: e_charge = 1.602176634d-19

contains

! Raw, unregularized grand potential of the finite BM matrix, J / cell.
! One spin and one valley. All occupied eigenvalues are included (no ncut).
! A=(z By,-z Bx,0), electron charge -e, k -> k+e A/hbar.
subroutine CALCULATE_OMEGA_B(N, nk, b1, b2, mu, zdiag, ham_param, Bfield, Omega)
    integer, intent(in) :: N, nk
    real(8), intent(in) :: b1(2), b2(2), mu, zdiag(N), ham_param(:), Bfield(2)
    real(8), intent(out) :: Omega
    integer :: iu, iv, i, j, info
    real(8) :: k(2), sumOmega
    real(8), allocatable :: W(:), rwork(:)
    complex(8), allocatable :: H(:,:), Href(:,:), hx(:,:), hy(:,:), HB(:,:), work(:)

    if (N /= 2*NTBG .or. nk <= 0) error stop 'OMEGA_B: invalid N or nk'
    if (.not. all(ieee_is_finite([mu,Bfield]))) error stop 'OMEGA_B: nonfinite input'
    allocate(Href(N,N), hx(N,N), hy(N,N), HB(N,N))
    call HAM_STRAIN([0.d0,0.d0], ham_param, Href)
    call HAM_DERIVATIVES_STRAIN_PHI([0.d0,0.d0], ham_param, hx, hy)
    ! This cache is exact for the current k-linear, locally coupled BM model.
    ! The anticommutator guarantees Hermiticity; here [z,hx]=[z,hy]=0,
    ! so it is exactly the row multiplication in the proposed implementation.
    do j=1,N
        do i=1,N
            HB(i,j) = (e_charge/(2.d0*hbar))*(zdiag(i)+zdiag(j))* &
                (Bfield(2)*hx(i,j)-Bfield(1)*hy(i,j))
        enddo
    enddo
    sumOmega = 0.d0
    !$omp parallel default(none) shared(N,nk,b1,b2,mu,Href,hx,hy,HB) &
    !$omp private(iu,iv,k,H,W,work,rwork,info) reduction(+:sumOmega)
    allocate(H(N,N), W(N), work( max(1, 64*N) ), rwork(max(1,3*N-2)))
    !$omp do collapse(2) schedule(static)
    do iu=1,nk
        do iv=1,nk
            k = (-0.5d0+real(iu-1,8)/nk)*b1 + (-0.5d0+real(iv-1,8)/nk)*b2
            H = Href + k(1)*hx + k(2)*hy + HB
            ! Eigenvalues only; no eigenvectors are needed for Omega.
            call zheev('N','U',N,H,N,W,work,size(work),rwork,info)
            if (info /= 0) error stop 'OMEGA_B: zheev failed'
            sumOmega = sumOmega + sum(min(W-mu,0.d0))
        enddo
    enddo
    !$omp end do
    deallocate(H,W,work,rwork)
    !$omp end parallel
    Omega = sumOmega/real(nk,8)**2
end subroutine CALCULATE_OMEGA_B

! Raw moment, mu_B / cell. Check B0 convergence; this is NOT a UV subtraction.
subroutine CALCULATE_OM_MAGNETIC(N,nk,b1,b2,mu,zdiag,ham_param,B0,Mx,My,omegas)
    integer, intent(in) :: N,nk
    real(8), intent(in) :: b1(2),b2(2),mu,zdiag(N),ham_param(:),B0
    real(8), intent(out) :: Mx,My
    real(8), intent(out), optional :: omegas(4)
    real(8) :: o(4),mu_B
    if (.not. ieee_is_finite(B0)) error stop 'OM_MAGNETIC: nonfinite B0'
    if (B0 <= 0.d0) error stop 'OM_MAGNETIC: B0 must be positive'
    mu_B=e_charge*hbar/(2.d0*m_e)
    call CALCULATE_OMEGA_B(N,nk,b1,b2,mu,zdiag,ham_param,[ B0,0.d0],o(1))
    call CALCULATE_OMEGA_B(N,nk,b1,b2,mu,zdiag,ham_param,[-B0,0.d0],o(2))
    call CALCULATE_OMEGA_B(N,nk,b1,b2,mu,zdiag,ham_param,[0.d0, B0],o(3))
    call CALCULATE_OMEGA_B(N,nk,b1,b2,mu,zdiag,ham_param,[0.d0,-B0],o(4))
    Mx=-(o(1)-o(2))/(2.d0*B0*mu_B)
    My=-(o(3)-o(4))/(2.d0*B0*mu_B)
    if (present(omegas)) omegas=o
end subroutine CALCULATE_OM_MAGNETIC
end module OM_MAGNETIC
