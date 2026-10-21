program check_dos
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use CONSTANTS, only : PI, eV, NTBG
    use STRAIN, only : H_STRAIN, GEOMETRY_STRAIN
    use LINALG, only : myzheev
    use DOS, only : CAL_DOS_TRIANGLE, DOS_FROM_BANDS
    implicit none

    real(8), parameter :: theta = 1.4d0*PI/180.d0
    real(8), parameter :: epsilon = 0.001d0
    real(8), parameter :: phi = 90.d0*PI/180.d0
    real(8), parameter :: analytic_tolerance = 2.d-11
    integer :: failures

    failures = 0
    call check_linear_band()
    call check_tent_band()
    call check_flat_bands()
    call check_empty_window()
    call check_narrow_triangle()
    call check_physical_sampling()

    if (failures /= 0) then
        print '(a,i0)', 'DOS regression checks failed: ', failures
        error stop 1
    end if
    print '(a)', 'DOS regression checks passed.'

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//message
        end if
    end subroutine check

    function integral(axis, density) result(value)
        real(8), intent(in) :: axis(:), density(:)
        real(8) :: value

        value = sum(0.5d0*(density(:size(density)-1)+density(2:))*(axis(2:)-axis(:size(axis)-1)))
    end function integral

    subroutine assert_regular_output(axis, density, lower, upper)
        real(8), intent(in) :: axis(:), density(:), lower, upper

        call check(all(ieee_is_finite(axis)) .and. all(ieee_is_finite(density)), 'output is finite')
        call check(all(density >= 0.d0), 'DOS is nonnegative')
        call check(abs(axis(1)-lower) < 1.d-14*eV, 'energy axis starts at requested energy in J')
        call check(abs(axis(size(axis))-upper) < 1.d-14*eV, 'energy axis ends at requested energy in J')
        call check(all(axis(2:) > axis(:size(axis)-1)), 'energy axis is strictly increasing')
    end subroutine assert_regular_output

    subroutine check_linear_band()
        ! E(x,y) = x eV has a uniform DOS of one state/eV on [0,1].
        ! This also exercises triangles whose two vertices have equal energy,
        ! including both sides of a shared grid line at E = 0.5 eV.
        real(8) :: bands(3,3,1), two_bands(3,3,2)
        real(8) :: axis(101), density(101), fine_axis(201), fine_density(201)
        real(8) :: clipped_axis(121), clipped_density(121)
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                bands(i,j,1) = real(i-1,8)*0.5d0*eV
            end do
        end do
        call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
        call assert_regular_output(axis, density, 0.d0, eV)
        call check(maxval(abs(density*eV-1.d0)) < analytic_tolerance, &
                   'linear band gives constant 1/eV DOS, including equal-energy vertices')
        call check(abs(integral(axis,density)-1.d0) < analytic_tolerance, 'one complete band integrates to one')

        call DOS_FROM_BANDS(bands, 0.d0, eV, 200, fine_axis, fine_density)
        call check(maxval(abs(fine_density*eV-1.d0)) < analytic_tolerance, &
                   'refining the energy axis does not multiply DOS by an extra 1/dE')
        call check(abs(density(31)-fine_density(61))*eV < analytic_tolerance, &
                   'same interior energy has the same DOS after energy-grid refinement')

        two_bands(:,:,1) = bands(:,:,1)
        two_bands(:,:,2) = bands(:,:,1)
        call DOS_FROM_BANDS(two_bands, 0.d0, eV, 100, axis, density)
        call check(maxval(abs(density*eV-2.d0)) < analytic_tolerance, 'identical bands add their DOS')
        call check(abs(integral(axis,density)-2.d0) < analytic_tolerance, 'two complete bands integrate to two')

        call DOS_FROM_BANDS(bands, 0.2d0*eV, 0.8d0*eV, 120, clipped_axis, clipped_density)
        call check(maxval(abs(clipped_density*eV-1.d0)) < analytic_tolerance, 'clipped energy window keeps local DOS')
        call check(abs(integral(clipped_axis,clipped_density)-0.6d0) < analytic_tolerance, &
                   'clipped energy window contains only its actual spectral weight')
    end subroutine check_linear_band

    subroutine check_tent_band()
        ! E(x,y) = (x+y) eV has the exact triangular density E below one
        ! and 2-E above one. Output nodes are means over midpoint bins,
        ! with half-width bins at the two ends of the energy window.
        real(8) :: bands(3,3,1), axis(201), density(201), expected(201)
        real(8) :: left, right, step
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                bands(i,j,1) = real(i+j-2,8)*0.5d0*eV
            end do
        end do
        call DOS_FROM_BANDS(bands, 0.d0, 2.d0*eV, 200, axis, density)
        call assert_regular_output(axis, density, 0.d0, 2.d0*eV)
        step = (axis(2)-axis(1))/eV
        do i = 1, size(axis)
            left = max(axis(1)/eV,axis(i)/eV-0.5d0*step)
            right = min(axis(size(axis))/eV,axis(i)/eV+0.5d0*step)
            expected(i) = (tent_cdf(right)-tent_cdf(left))/(right-left)
        end do
        call check(maxval(abs(density*eV-expected)) < analytic_tolerance, &
                   'E=x+y matches analytical triangle-bin probabilities')
        call check(abs(integral(axis,density)-1.d0) < analytic_tolerance, 'tent band preserves one state')
        call check(abs(density(101)*eV-(1.d0-step/4.d0)) < analytic_tolerance, &
                   'bin straddling the tent peak has its analytical mean density')
    end subroutine check_tent_band

    function tent_cdf(energy_ev) result(value)
        real(8), intent(in) :: energy_ev
        real(8) :: value

        if (energy_ev <= 0.d0) then
            value = 0.d0
        else if (energy_ev <= 1.d0) then
            value = 0.5d0*energy_ev**2
        else if (energy_ev < 2.d0) then
            value = 1.d0-0.5d0*(2.d0-energy_ev)**2
        else
            value = 1.d0
        end if
    end function tent_cdf

    subroutine check_flat_bands()
        real(8) :: bands(3,3,1), axis(101), density(101), expected(101)
        real(8), parameter :: centers(3) = [0.d0,0.373d0,1.d0]
        real(8) :: step
        integer :: i, node

        step = eV/100.d0
        do i = 1, size(centers)
            bands = centers(i)*eV
            call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
            call assert_regular_output(axis, density, 0.d0, eV)
            node = nint(centers(i)*100.d0)+1
            expected = 0.d0
            expected(node) = 1.d0/step
            if (node == 1 .or. node == size(axis)) expected(node) = 2.d0/step
            call check(maxval(abs(density-expected))*eV < analytic_tolerance, &
                       'flat band deposits delta weight on nearest energy node with endpoint correction')
            call check(abs(integral(axis,density)-1.d0) < analytic_tolerance, &
                       'flat band delta integrates to one, including at either endpoint')
        end do
    end subroutine check_flat_bands

    subroutine check_empty_window()
        real(8) :: bands(3,3,1), axis(101), density(101)
        integer :: i, j

        ! Nearby out-of-window flat bands must not be rounded onto an end node.
        bands = -0.001d0*eV
        call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
        call check(all(density == 0.d0), 'flat band below energy window contributes zero')
        bands = 1.001d0*eV
        call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
        call check(all(density == 0.d0), 'flat band above energy window contributes zero')

        do j = 1, 3
            do i = 1, 3
                bands(i,j,1) = (2.d0+real(i+j-2,8)*0.5d0)*eV
            end do
        end do
        call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
        call check(all(density == 0.d0), 'ordinary triangles outside energy window contribute zero')
    end subroutine check_empty_window

    subroutine check_narrow_triangle()
        ! A band far narrower than dE and crossing the energy-window edge
        ! must keep its fractional spectral weight, rather than becoming a delta.
        real(8) :: bands(2,2,1), axis(101), density(101)
        bands(1,:,1) = -0.25d-9*eV
        bands(2,:,1) = 0.75d-9*eV
        call DOS_FROM_BANDS(bands, 0.d0, eV, 100, axis, density)
        call check(abs(integral(axis,density)-0.75d0) < analytic_tolerance, &
                   'narrow triangles preserve partial state counts at the energy-window edge')
    end subroutine check_narrow_triangle

    subroutine check_physical_sampling()
        ! Build every k-point independently with the existing H_STRAIN and
        ! full eigenvectors. This checks the optimized H(0)+k construction,
        ! eigenvalues-only solver, geometry, all-band indexing, and J units.
        integer, parameter :: nx = 2, ny = 2, nh = 2*NTBG, ne = 200
        real(8) :: bands(nx+1,ny+1,nh), axis(ne+1), density(ne+1)
        real(8) :: reference_axis(ne+1), reference_density(ne+1), eigenvalues(nh)
        real(8) :: bs(2,2), k(2), scale
        complex(8), allocatable :: h(:,:), vectors(:,:)
        integer :: i, j

        allocate(h(nh,nh), vectors(nh,nh))
        call GEOMETRY_STRAIN(theta, epsilon, phi, bs)
        do j = 0, ny
            do i = 0, nx
                k = real(i,8)/real(nx,8)*bs(:,1)+real(j,8)/real(ny,8)*bs(:,2)
                call H_STRAIN(k, theta, epsilon, phi, h)
                call myzheev(nh, eigenvalues, h, vectors)
                bands(i+1,j+1,:) = eigenvalues
            end do
        end do
        call DOS_FROM_BANDS(bands, -0.05d0*eV, 0.05d0*eV, ne, reference_axis, reference_density)
        call CAL_DOS_TRIANGLE(theta, epsilon, phi, -0.05d0*eV, 0.05d0*eV, ne, &
                              axis, density, nkx=nx, nky=ny)
        call assert_regular_output(axis, density, -0.05d0*eV, 0.05d0*eV)
        scale = maxval(abs(reference_density))
        call check(scale > 0.d0, 'physical comparison window contains spectral weight')
        call check(maxval(abs(axis-reference_axis)) < 1.d-14*eV, 'physical and reference energy axes agree')
        call check(maxval(abs(density-reference_density)) < 1.d-7*max(scale,1.d0/eV), &
                   'optimized physical DOS agrees with independent full-spectrum sampling')
    end subroutine check_physical_sampling

end program check_dos
