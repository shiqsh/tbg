program check_dirac
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use CONSTANTS, only : PI, eV, NTBG
    use STRAIN, only : H_STRAIN, GEOMETRY_STRAIN
    use LINALG, only : myzheev
    use DIRAC, only : FIND_DIRAC_POINTS, CALCULATE_GAP, dirac_candidate
    implicit none

    real(8), parameter :: theta = 1.4d0*PI/180.d0
    real(8), parameter :: epsilon = 0.001d0
    real(8), parameter :: phi = 37.d0*PI/180.d0
    real(8), parameter :: spectral_tolerance = 1.d-10*eV
    real(8), parameter :: gap_tolerance = 1.d-8*eV
    real(8) :: bs(2,2), q1(2), k(2), kscale, test_points(2,3)
    real(8) :: gap, energy, reference_gap, reference_energy, seeds(2,2)
    type(dirac_candidate), allocatable :: unstrained(:), strained(:), periodic(:), limited(:)
    integer :: i, n_found, n_periodic, failures

    failures = 0

    ! Compare the new two-eigenvalue calculation with the existing full solver
    ! at momenta away from the optimized nodes, where both bands are separated.
    call GEOMETRY_STRAIN(theta, epsilon, phi, bs, q1)
    kscale = sqrt(sum(q1**2))
    test_points(:,1) = 0.07d0*bs(:,1) - 0.11d0*bs(:,2)
    test_points(:,2) = 0.5d0*q1 + 0.13d0*bs(:,1)
    test_points(:,3) = q1 - 0.09d0*bs(:,1) + 0.04d0*bs(:,2)
    do i = 1, size(test_points,2)
        k = test_points(:,i)
        call CALCULATE_GAP(k(1), k(2), theta, epsilon, phi, gap, energy)
        call full_spectrum(k, epsilon, phi, reference_gap, reference_energy)
        call check(abs(gap-reference_gap) < spectral_tolerance, 'arbitrary-k gap matches full spectrum')
        call check(abs(energy-reference_energy) < spectral_tolerance, 'arbitrary-k energy matches full spectrum')
    end do

    ! Zero strain recovers the two symmetry-related nodes and zero node energy.
    call GEOMETRY_STRAIN(theta, 0.d0, phi, bs, q1)
    kscale = sqrt(sum(q1**2))
    call FIND_DIRAC_POINTS(theta, 0.d0, phi, unstrained, n_found, &
                           gap_tolerance=gap_tolerance, k_tolerance=1.d-9)
    call check(size(unstrained) == 2, 'zero strain returns both candidate slots')
    call check(n_found == 2, 'zero strain finds two Dirac points')
    if (size(unstrained) == 2) then
        call check(sqrt(unstrained(1)%kx**2 + unstrained(1)%ky**2)/kscale < 1.d-5, &
                   'first zero-strain node is near the origin')
        call check(sqrt((unstrained(2)%kx-q1(1))**2 + (unstrained(2)%ky-q1(2))**2)/kscale < 1.d-5, &
                   'second zero-strain node is near q1')
        do i = 1, 2
            call check(unstrained(i)%valid .and. unstrained(i)%converged, 'zero-strain node is accepted and converged')
            call check(abs(unstrained(i)%energy) < 1.d-6*eV, 'zero-strain node energy is near zero')
            call check(unstrained(i)%gap < gap_tolerance, 'zero-strain node is gapless')
            k = [unstrained(i)%kx, unstrained(i)%ky]
            call verify_candidate(unstrained(i), k, 0.d0, phi)
        end do
    end if

    ! A uniaxial strain tensor is invariant under phi -> phi + pi. Do not
    ! require different orientations to have equal energies at finite strain.
    call FIND_DIRAC_POINTS(theta, epsilon, phi, strained, n_found, &
                           gap_tolerance=gap_tolerance, k_tolerance=1.d-9)
    call FIND_DIRAC_POINTS(theta, epsilon, phi+PI, periodic, n_periodic, &
                           gap_tolerance=gap_tolerance, k_tolerance=1.d-9)
    call check(size(strained) == 2 .and. size(periodic) == 2, 'strained calculations return both candidate slots')
    call check(n_found == 2 .and. n_periodic == 2, 'finite strain and periodic angle find both nodes')
    if (size(strained) == 2 .and. size(periodic) == 2) then
        do i = 1, 2
            call check(strained(i)%valid .and. strained(i)%converged, 'finite-strain node is accepted and converged')
            call check(periodic(i)%valid .and. periodic(i)%converged, 'periodic-angle node is accepted and converged')
            k = [strained(i)%kx-periodic(i)%kx, strained(i)%ky-periodic(i)%ky]
            call check(sqrt(sum(k**2))/kscale < 1.d-6, 'node position has pi angle periodicity')
            call check(abs(strained(i)%energy-periodic(i)%energy) < gap_tolerance, &
                       'node energy has pi angle periodicity')
            k = [strained(i)%kx, strained(i)%ky]
            call verify_candidate(strained(i), k, epsilon, phi)
        end do
    end if

    ! A deliberately inadequate iteration budget must remain visible to callers.
    ! Far-away starts prevent accidental convergence at a default exact seed.
    call GEOMETRY_STRAIN(theta, epsilon, phi, bs, q1)
    seeds(:,1) = 0.21d0*bs(:,1) + 0.13d0*bs(:,2)
    seeds(:,2) = q1 - 0.17d0*bs(:,1) + 0.23d0*bs(:,2)
    call FIND_DIRAC_POINTS(theta, epsilon, phi, limited, n_found, seeds=seeds, &
                           gap_tolerance=1.d-10*eV, k_tolerance=1.d-12, max_iterations=1)
    call check(size(limited) == 2, 'failed search preserves both candidate slots')
    if (size(limited) == 2) then
        call check(.not. any(limited%converged), 'one iteration does not claim convergence')
        call check(.not. any(limited%valid), 'unconverged candidates are rejected')
    end if
    call check(n_found == 0, 'failed search reports zero accepted nodes')

    if (failures /= 0) then
        print '(a,i0)', 'Dirac regression checks failed: ', failures
        error stop 1
    end if
    print '(a)', 'Dirac regression checks passed.'

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//message
        end if
    end subroutine check

    subroutine full_spectrum(momentum, eps, angle, central_gap, central_energy)
        real(8), intent(in) :: momentum(2), eps, angle
        real(8), intent(out) :: central_gap, central_energy
        complex(8), allocatable :: h(:,:), vectors(:,:)
        real(8), allocatable :: eigenvalues(:)
        integer, parameter :: n = 2*NTBG

        allocate(h(n,n), vectors(n,n), eigenvalues(n))
        call H_STRAIN(momentum, theta, eps, angle, h)
        call myzheev(n, eigenvalues, h, vectors)
        central_gap = eigenvalues(NTBG+1) - eigenvalues(NTBG)
        central_energy = 0.5d0*(eigenvalues(NTBG+1) + eigenvalues(NTBG))
    end subroutine full_spectrum

    subroutine verify_candidate(candidate, momentum, eps, angle)
        type(dirac_candidate), intent(in) :: candidate
        real(8), intent(in) :: momentum(2), eps, angle
        real(8) :: full_gap, full_energy

        call check(ieee_is_finite(candidate%kx) .and. ieee_is_finite(candidate%ky), 'candidate momentum is finite')
        call check(ieee_is_finite(candidate%energy) .and. ieee_is_finite(candidate%gap), 'candidate energies are finite')
        call check(candidate%evaluations > 0, 'candidate records function evaluations')
        call full_spectrum(momentum, eps, angle, full_gap, full_energy)
        call check(abs(candidate%gap-full_gap) < spectral_tolerance, 'node gap matches full spectrum')
        call check(abs(candidate%energy-full_energy) < spectral_tolerance, 'node energy matches full spectrum')
        call check(full_gap < gap_tolerance, 'full spectrum confirms a gapless node')
    end subroutine verify_candidate

end program check_dirac
