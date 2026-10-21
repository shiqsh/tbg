module DOS
    !! Legacy triangle integration on k = u*b1s + v*b2s, 0 <= u,v <= 1.
    !! Angles: radians; epsilon: fraction; energies: J; DOS: states/J/cell
    !! for one spin and one valley. Each complete band integrates to one.
    use CONSTANTS, only : NTBG, eV, vF, CI
    use STRAIN, only : H_STRAIN, GEOMETRY_STRAIN
    use stdlib_linalg_lapack, only : heev
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    integer, parameter :: nh = 2*NTBG
    public :: CAL_DOS_TRIANGLE, DOS_FROM_BANDS

contains

    subroutine CAL_DOS_TRIANGLE(theta, epsilon, phi, E_min, E_max, N_intervals, &
                               dos_energy, dos_value, nkx, nky)
        real(8), intent(in) :: theta, epsilon, phi, E_min, E_max
        integer, intent(in) :: N_intervals
        real(8), intent(out) :: dos_energy(:), dos_value(:)
        integer, intent(in), optional :: nkx, nky
        real(8), allocatable :: bands(:,:,:)
        integer :: nx, ny

        if (.not. all(ieee_is_finite([theta,epsilon,phi]))) error stop 'DOS: nonfinite physical parameters'
        if (theta <= 0.d0) error stop 'DOS: theta must be positive (radians)'
        call CHECK_ENERGY_AXIS(E_min, E_max, N_intervals, size(dos_energy), size(dos_value))
        nx = 100
        ny = 100
        if (present(nkx)) nx = nkx
        if (present(nky)) ny = nky
        if (nx < 1 .or. ny < 1) error stop 'DOS: nkx and nky must be positive'

        call SAMPLE_BANDS(theta, epsilon, phi, nx, ny, bands)
        call DOS_FROM_BANDS(bands, E_min, E_max, N_intervals, dos_energy, dos_value)
    end subroutine CAL_DOS_TRIANGLE

    subroutine SAMPLE_BANDS(theta, epsilon, phi, nx, ny, bands)
        real(8), intent(in) :: theta, epsilon, phi
        integer, intent(in) :: nx, ny
        real(8), allocatable, intent(out) :: bands(:,:,:)
        complex(8), allocatable :: h0(:,:), h(:,:), work(:)
        real(8), allocatable :: eigenvalues(:), rwork(:)
        complex(8) :: work_query(1), offdiagonal
        real(8) :: bs(2,2), k(2)
        integer :: i, j, block, info, lwork

        allocate(h0(nh,nh), h(nh,nh), eigenvalues(nh), rwork(3*nh-2))
        call GEOMETRY_STRAIN(theta, epsilon, phi, bs)
        call H_STRAIN([0.d0,0.d0], theta, epsilon, phi, h0)
        h0 = h0/eV
        h = h0
        call heev('N', 'U', nh, h, nh, eigenvalues, work_query, -1, rwork, info)
        if (info /= 0) error stop 'DOS: LAPACK workspace query failed'
        lwork = max(2*nh,ceiling(real(work_query(1),8)))
        deallocate(h, eigenvalues, rwork)
        allocate(bands(nx+1,ny+1,nh))

        ! One matrix and queried LAPACK workspace per thread. H(0) is shared;
        ! no eigenvectors, critical counters, or nested H_STRAIN calls at k.
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(nx,ny,bs,h0,bands,lwork) &
        !$OMP PRIVATE(h,work,eigenvalues,rwork,i,j,k,offdiagonal,block,info)
        allocate(h(nh,nh), work(lwork), eigenvalues(nh), rwork(3*nh-2))
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        do j = 0, ny
            do i = 0, nx
                k = real(i,8)/real(nx,8)*bs(:,1) + real(j,8)/real(ny,8)*bs(:,2)
                h = h0
                offdiagonal = (vF/eV)*(k(1)-CI*k(2))
                do block = 1, nh, 2
                    h(block,block+1) = h(block,block+1) + offdiagonal
                end do
                call heev('N', 'U', nh, h, nh, eigenvalues, work, lwork, rwork, info)
                if (info /= 0) error stop 'DOS: eigenvalue calculation failed'
                bands(i+1,j+1,:) = eigenvalues*eV
            end do
        end do
        !$OMP END DO
        deallocate(h, work, eigenvalues, rwork)
        !$OMP END PARALLEL
    end subroutine SAMPLE_BANDS

    subroutine DOS_FROM_BANDS(bands, E_min, E_max, N_intervals, dos_energy, dos_value)
        !! Integrate vertex energies in J using the original two triangles per
        !! cell. Return bin averages of the original tent, with midpoint bin
        !! boundaries and half-width endpoint bins. Trapezoidal integration
        !! then counts states exactly, even for triangles narrower than dE.
        real(8), intent(in) :: bands(:,:,:), E_min, E_max
        integer, intent(in) :: N_intervals
        real(8), intent(out) :: dos_energy(:), dos_value(:)
        real(8) :: dE, weight, e1, e2, e3, e4
        integer :: nx, ny, i, j, ib, ie

        call CHECK_ENERGY_AXIS(E_min, E_max, N_intervals, size(dos_energy), size(dos_value))
        nx = size(bands,1)-1
        ny = size(bands,2)-1
        if (nx < 1 .or. ny < 1 .or. size(bands,3) < 1) error stop 'DOS: invalid band-grid shape'
        if (.not. all(ieee_is_finite(bands))) error stop 'DOS: nonfinite band energy'
        dE = (E_max-E_min)/real(N_intervals,8)
        do ie = 1, N_intervals+1
            dos_energy(ie) = E_min + real(ie-1,8)*dE
        end do
        dos_energy(N_intervals+1) = E_max
        dos_value = 0.d0
        weight = 0.5d0/(real(nx,8)*real(ny,8))

        do ib = 1, size(bands,3)
            ! Piecewise-linear triangles cannot leave their vertex-energy range.
            if (minval(bands(:,:,ib)) > E_max .or. maxval(bands(:,:,ib)) < E_min) cycle
            do j = 1, ny
                do i = 1, nx
                    e1 = bands(i,j,ib)
                    e2 = bands(i+1,j,ib)
                    e3 = bands(i,j+1,ib)
                    e4 = bands(i+1,j+1,ib)
                    call ADD_TRIANGLE(e1, e2, e3, weight, dE, dos_energy, dos_value)
                    call ADD_TRIANGLE(e2, e4, e3, weight, dE, dos_energy, dos_value)
                end do
            end do
        end do
        ! The tent height already has units 1/J. Do not multiply by 1/dE.
    end subroutine DOS_FROM_BANDS

    subroutine ADD_TRIANGLE(e1, e2, e3, weight, dE, axis, density)
        real(8), intent(in) :: e1, e2, e3, weight, dE, axis(:)
        real(8), intent(inout) :: density(:)
        real(8) :: ea, eb, ec, tmp, height, lower, upper, center, mass
        integer :: first, last, ie, n

        ea = e1; eb = e2; ec = e3
        if (ea > eb) then
            tmp = ea; ea = eb; eb = tmp
        end if
        if (eb > ec) then
            tmp = eb; eb = ec; ec = tmp
        end if
        if (ea > eb) then
            tmp = ea; ea = eb; eb = tmp
        end if
        n = size(axis)
        if (ec < axis(1) .or. ea > axis(n)) return

        if (ec == ea) then
            ! A flat triangle contributes weight*delta(E-center). Deposit on
            ! the nearest energy node; endpoint trapezoidal weights are 1/2.
            center = (ea+eb+ec)/3.d0
            if (center < axis(1) .or. center > axis(n)) return
            ie = max(1,min(n,nint((center-axis(1))/dE)+1))
            height = weight/dE
            if (ie == 1 .or. ie == n) height = 2.d0*height
            density(ie) = density(ie) + height
            return
        end if

        ! Average the same linear-interpolation tent over each energy bin,
        ! rather than sampling its height and potentially missing a narrow peak.
        ! Clamping before index conversion also avoids remote-band overflow.
        first = max(1,floor((max(ea,axis(1))-axis(1))/dE+0.5d0)+1)
        last = min(n,floor((min(ec,axis(n))-axis(1))/dE+0.5d0)+1)
        do ie = first, last
            lower = max(axis(1),axis(ie)-0.5d0*dE)
            upper = min(axis(n),axis(ie)+0.5d0*dE)
            mass = max(0.d0,TRIANGLE_CDF(upper,ea,eb,ec)-TRIANGLE_CDF(lower,ea,eb,ec))
            density(ie) = density(ie) + weight*mass/(upper-lower)
        end do
    end subroutine ADD_TRIANGLE

    real(8) function TRIANGLE_CDF(energy, ea, eb, ec) result(value)
        real(8), intent(in) :: energy, ea, eb, ec
        ! Integral of the legacy tent of height 2/(ec-ea); its area is one.
        ! The branch conditions avoid division by zero for ea=eb or eb=ec.
        if (energy <= ea) then
            value = 0.d0
        else if (energy >= ec) then
            value = 1.d0
        else if (energy < eb) then
            value = ((energy-ea)/(eb-ea))*((energy-ea)/(ec-ea))
        else
            value = 1.d0-((ec-energy)/(ec-ea))*((ec-energy)/(ec-eb))
        end if
    end function TRIANGLE_CDF

    subroutine CHECK_ENERGY_AXIS(E_min, E_max, N_intervals, n_energy, n_dos)
        real(8), intent(in) :: E_min, E_max
        integer, intent(in) :: N_intervals, n_energy, n_dos
        if (.not. all(ieee_is_finite([E_min,E_max]))) error stop 'DOS: nonfinite energy bounds'
        if (E_max <= E_min .or. N_intervals < 1) error stop 'DOS: invalid energy range or interval count'
        if (n_energy /= N_intervals+1 .or. n_dos /= N_intervals+1) error stop 'DOS: output size must be N_intervals+1'
    end subroutine CHECK_ENERGY_AXIS
end module DOS
