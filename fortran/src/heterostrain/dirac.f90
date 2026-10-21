module DIRAC
    !! Shifted central-band crossings of H_STRAIN. Angles: rad; strain:
    !! fraction (0.001 = 0.1%); k: m^-1; returned energy/gap: J.
    use CONSTANTS, only : NTBG, eV, vF, K_graphene, CI
    use STRAIN, only : H_STRAIN, GEOMETRY_STRAIN
    use stdlib_linalg_lapack, only : heevr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter :: nh = 2*NTBG
    real(8), parameter :: default_gap_tolerance = 1.d-8*eV
    real(8), parameter :: default_k_tolerance = 1.d-8

    type, public :: dirac_candidate
        real(8) :: kx = 0.d0, ky = 0.d0
        real(8) :: energy = 0.d0, gap = huge(1.d0)
        logical :: converged = .false., valid = .false.
        integer :: evaluations = 0
    end type dirac_candidate

    ! One workspace per call; no saved mutable state or nested OpenMP search.
    type :: dirac_workspace
        complex(8), allocatable :: h0(:,:), h(:,:), work(:)
        real(8), allocatable :: rwork(:)
        integer, allocatable :: iwork(:)
        real(8) :: scale, w(nh)
    end type dirac_workspace

    public :: FIND_DIRAC_POINT, FIND_DIRAC_POINTS, CALCULATE_GAP

contains

    subroutine FIND_DIRAC_POINT(theta, epsilon, phi, candidate, seed, &
                               gap_tolerance, k_tolerance, max_iterations)
        !! Local search near seed (default: k=0). No global-BZ claim.
        real(8), intent(in) :: theta, epsilon, phi
        type(dirac_candidate), intent(out) :: candidate
        real(8), intent(in), optional :: seed(2), gap_tolerance, k_tolerance
        integer, intent(in), optional :: max_iterations
        type(dirac_workspace) :: ws
        real(8) :: start(2), step

        call INITIALIZE(theta, epsilon, phi, ws)
        start = 0.d0
        step = 0.03d0
        if (present(seed)) then
            if (.not. all(ieee_is_finite(seed))) error stop 'DIRAC: nonfinite seed'
            start = seed/ws%scale
            step = 0.005d0
        end if
        call MINIMIZE_GAP(ws, start, step, candidate, gap_tolerance, k_tolerance, max_iterations)
    end subroutine FIND_DIRAC_POINT

    subroutine FIND_DIRAC_POINTS(theta, epsilon, phi, candidates, n_found, seeds, &
                                gap_tolerance, k_tolerance, max_iterations)
        !! Two local searches initialized near 0 and q1s, or continuation seeds.
        !! Always return two diagnostics, including unsuccessful searches.
        !! n_found counts valid distinct crossings; array size is always two.
        real(8), intent(in) :: theta, epsilon, phi
        type(dirac_candidate), allocatable, intent(out) :: candidates(:)
        integer, intent(out) :: n_found
        real(8), intent(in), optional :: seeds(2,2), gap_tolerance, k_tolerance
        integer, intent(in), optional :: max_iterations
        type(dirac_workspace) :: ws
        type(dirac_candidate) :: trial
        real(8) :: bs(2,2), q1(2), start(2,2), midpoint(2), step, separation
        real(8) :: direct_distance, swapped_distance
        logical :: retry, duplicate
        real(8), parameter :: directions(2,4) = reshape([1.d0,0.d0, -1.d0,0.d0, &
                                                       0.d0,1.d0, 0.d0,-1.d0], [2,4])
        integer :: i, j, other, evaluations

        call INITIALIZE(theta, epsilon, phi, ws)
        call GEOMETRY_STRAIN(theta, epsilon, phi, bs, q1)
        start(:,1) = 0.d0
        start(:,2) = q1/ws%scale
        step = 0.03d0
        if (present(seeds)) then
            if (.not. all(ieee_is_finite(seeds))) error stop 'DIRAC: nonfinite seeds'
            start = seeds/ws%scale
            step = 0.005d0
        end if
        allocate(candidates(2))
        do i = 1, 2
            call MINIMIZE_GAP(ws, start(:,i), step, candidates(i), &
                              gap_tolerance, k_tolerance, max_iterations)
        end do
        separation = 10.d0*default_k_tolerance
        if (present(k_tolerance)) separation = max(10.d0*k_tolerance,1.d-8)
        retry = .true.
        if (present(max_iterations)) retry = max_iterations >= 20
        ! Recover failed searches or shared attraction basins with transverse
        ! starts. Branches are trajectories, not permanent layer identities.
        if (retry) then
            do j = 1, 2
                other = 3-j
                duplicate = candidates(j)%valid .and. candidates(other)%valid .and. &
                            SAME_POINT(candidates(j), candidates(other), ws%scale, separation)
                if (candidates(j)%valid .and. (.not. duplicate .or. j == 1)) cycle
                midpoint = start(:,j)
                if (candidates(other)%valid) midpoint = [candidates(other)%kx,candidates(other)%ky]/ws%scale
                evaluations = candidates(j)%evaluations
                do i = 1, 4
                    call MINIMIZE_GAP(ws, midpoint + 0.15d0*directions(:,i), 0.03d0, trial, &
                                      gap_tolerance, k_tolerance, max_iterations)
                    evaluations = evaluations + trial%evaluations
                    if (trial%valid) then
                        if (candidates(other)%valid) then
                            if (SAME_POINT(trial, candidates(other), ws%scale, separation)) cycle
                        end if
                        candidates(j) = trial
                        exit
                    else if (.not. candidates(j)%valid .and. trial%gap < candidates(j)%gap) then
                        candidates(j) = trial
                    end if
                end do
                candidates(j)%evaluations = evaluations
            end do
        end if
        if (candidates(1)%valid .and. candidates(2)%valid) then
            if (SAME_POINT(candidates(1), candidates(2), ws%scale, separation)) candidates(2)%valid = .false.
        end if
        if (present(seeds) .and. all(candidates%valid)) then
            direct_distance = norm2([candidates(1)%kx,candidates(1)%ky]-seeds(:,1)) + &
                              norm2([candidates(2)%kx,candidates(2)%ky]-seeds(:,2))
            swapped_distance = norm2([candidates(2)%kx,candidates(2)%ky]-seeds(:,1)) + &
                               norm2([candidates(1)%kx,candidates(1)%ky]-seeds(:,2))
            if (swapped_distance < direct_distance) then
                trial = candidates(1)
                candidates(1) = candidates(2)
                candidates(2) = trial
            end if
        end if
        n_found = count(candidates%valid)
    end subroutine FIND_DIRAC_POINTS

    subroutine CALCULATE_GAP(kx, ky, theta, epsilon, phi, gap, energy)
        !! Central eigenvalues NTBG and NTBG+1; also return their mean if requested.
        real(8), intent(in) :: kx, ky, theta, epsilon, phi
        real(8), intent(out) :: gap
        real(8), intent(out), optional :: energy
        type(dirac_workspace) :: ws
        real(8) :: center

        if (.not. all(ieee_is_finite([kx,ky]))) error stop 'DIRAC: nonfinite k'
        call INITIALIZE(theta, epsilon, phi, ws)
        call CENTRAL_BANDS(ws, [kx,ky]/ws%scale, gap, center)
        gap = gap*eV
        if (present(energy)) energy = center*eV
    end subroutine CALCULATE_GAP

    subroutine INITIALIZE(theta, epsilon, phi, ws)
        real(8), intent(in) :: theta, epsilon, phi
        type(dirac_workspace), intent(out) :: ws
        complex(8) :: work_query(1), z_dummy(1,1)
        real(8) :: rwork_query(1)
        integer :: iwork_query(1), support(4), m, info

        if (.not. all(ieee_is_finite([theta,epsilon,phi]))) error stop 'DIRAC: nonfinite parameters'
        if (theta <= 0.d0) error stop 'DIRAC: theta must be positive (radians)'
        ws%scale = theta*norm2(K_graphene)
        allocate(ws%h0(nh,nh), ws%h(nh,nh))
        call H_STRAIN([0.d0,0.d0], theta, epsilon, phi, ws%h0)
        ! Order-one eV units avoid an absolute LAPACK tolerance in joules.
        ws%h0 = ws%h0/eV
        ws%h = ws%h0
        call heevr('N', 'I', 'U', nh, ws%h, nh, 0.d0, 0.d0, NTBG, NTBG+1, &
                   0.d0, m, ws%w, z_dummy, 1, support, work_query, -1, &
                   rwork_query, -1, iwork_query, -1, info)
        if (info /= 0) error stop 'DIRAC: LAPACK workspace query failed'
        allocate(ws%work(max(2*nh,ceiling(real(work_query(1),8)))), &
                 ws%rwork(max(24*nh,ceiling(rwork_query(1)))), &
                 ws%iwork(max(10*nh,iwork_query(1))))
    end subroutine INITIALIZE

    subroutine CENTRAL_BANDS(ws, x, gap, energy)
        type(dirac_workspace), intent(inout) :: ws
        real(8), intent(in) :: x(2)
        real(8), intent(out) :: gap, energy
        complex(8) :: z_dummy(1,1), offdiagonal
        integer :: i, m, info, support(4)

        ws%h = ws%h0
        ! Exact k dependence of H_STRAIN: add the same 2x2 Dirac block
        ! to every site. Only the upper triangle is used by LAPACK.
        offdiagonal = (vF/eV)*ws%scale*(x(1) - CI*x(2))
        do i = 1, nh, 2
            ws%h(i,i+1) = ws%h(i,i+1) + offdiagonal
        end do
        call heevr('N', 'I', 'U', nh, ws%h, nh, 0.d0, 0.d0, NTBG, NTBG+1, &
                   0.d0, m, ws%w, z_dummy, 1, support, ws%work, size(ws%work), &
                   ws%rwork, size(ws%rwork), ws%iwork, size(ws%iwork), info)
        if (info /= 0 .or. m /= 2) error stop 'DIRAC: central-band eigensolver failed'
        gap = max(0.d0,ws%w(2)-ws%w(1))
        energy = 0.5d0*(ws%w(1)+ws%w(2))
    end subroutine CENTRAL_BANDS

    subroutine MINIMIZE_GAP(ws, start, step, result, gap_tolerance, k_tolerance, max_iterations)
        !! Two-dimensional Nelder-Mead, with cached values at simplex vertices.
        type(dirac_workspace), intent(inout) :: ws
        real(8), intent(in) :: start(2), step
        type(dirac_candidate), intent(out) :: result
        real(8), intent(in), optional :: gap_tolerance, k_tolerance
        integer, intent(in), optional :: max_iterations
        real(8) :: x(2,3), g(3), energy(3), center(2), xr(2), xt(2)
        real(8) :: gr, er, gt, et, tol_gap, tol_k, diameter
        integer :: iteration, limit, i

        tol_gap = default_gap_tolerance/eV
        tol_k = default_k_tolerance
        limit = 200
        if (present(gap_tolerance)) tol_gap = gap_tolerance/eV
        if (present(k_tolerance)) tol_k = k_tolerance
        if (present(max_iterations)) limit = max_iterations
        if (.not. all(ieee_is_finite([tol_gap,tol_k]))) error stop 'DIRAC: nonfinite tolerances'
        if (tol_gap <= 0.d0 .or. tol_k <= 0.d0 .or. limit < 1) error stop 'DIRAC: invalid search options'
        result = dirac_candidate()
        x(:,1) = start
        x(:,2) = start + [step,0.d0]
        x(:,3) = start + [0.d0,step]
        do i = 1, 3
            call EVALUATE(x(:,i), g(i), energy(i))
        end do
        do iteration = 1, limit
            call SORT_SIMPLEX(x, g, energy)
            diameter = max(norm2(x(:,2)-x(:,1)),norm2(x(:,3)-x(:,1)))
            if (diameter <= tol_k) exit
            center = 0.5d0*(x(:,1)+x(:,2))
            xr = 2.d0*center - x(:,3)
            call EVALUATE(xr, gr, er)
            if (gr < g(1)) then
                xt = center + 2.d0*(xr-center)
                call EVALUATE(xt, gt, et)
                if (gt < gr) then
                    call REPLACE_WORST(xt, gt, et)
                else
                    call REPLACE_WORST(xr, gr, er)
                end if
            else if (gr < g(2)) then
                call REPLACE_WORST(xr, gr, er)
            else
                if (gr < g(3)) then
                    xt = center + 0.5d0*(xr-center)
                else
                    xt = center + 0.5d0*(x(:,3)-center)
                end if
                call EVALUATE(xt, gt, et)
                if (gt < min(gr,g(3))) then
                    call REPLACE_WORST(xt, gt, et)
                else
                    do i = 2, 3
                        x(:,i) = x(:,1) + 0.5d0*(x(:,i)-x(:,1))
                        call EVALUATE(x(:,i), g(i), energy(i))
                    end do
                end if
            end if
        end do
        call SORT_SIMPLEX(x, g, energy)
        diameter = max(norm2(x(:,2)-x(:,1)),norm2(x(:,3)-x(:,1)))
        result%kx = x(1,1)*ws%scale
        result%ky = x(2,1)*ws%scale
        result%gap = g(1)*eV
        result%energy = energy(1)*eV
        result%converged = diameter <= tol_k
        result%valid = result%converged .and. g(1) <= tol_gap
    contains
        subroutine EVALUATE(point, gap, mean_energy)
            real(8), intent(in) :: point(2)
            real(8), intent(out) :: gap, mean_energy
            call CENTRAL_BANDS(ws, point, gap, mean_energy)
            result%evaluations = result%evaluations + 1
        end subroutine EVALUATE

        subroutine REPLACE_WORST(point, gap, mean_energy)
            real(8), intent(in) :: point(2), gap, mean_energy
            x(:,3) = point
            g(3) = gap
            energy(3) = mean_energy
        end subroutine REPLACE_WORST
    end subroutine MINIMIZE_GAP

    subroutine SORT_SIMPLEX(x, g, energy)
        real(8), intent(inout) :: x(2,3), g(3), energy(3)
        real(8) :: point(2), value
        integer :: i, j
        do i = 1, 2
            do j = i+1, 3
                if (g(j) >= g(i)) cycle
                point = x(:,i); x(:,i) = x(:,j); x(:,j) = point
                value = g(i); g(i) = g(j); g(j) = value
                value = energy(i); energy(i) = energy(j); energy(j) = value
            end do
        end do
    end subroutine SORT_SIMPLEX

    logical function SAME_POINT(a, b, scale, tolerance)
        type(dirac_candidate), intent(in) :: a, b
        real(8), intent(in) :: scale, tolerance
        SAME_POINT = norm2([a%kx-b%kx,a%ky-b%ky]) <= tolerance*scale
    end function SAME_POINT
end module DIRAC
