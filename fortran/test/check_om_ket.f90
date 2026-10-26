program CHECK_OM_KET
    use CONSTANTS, only : NTBG, PI, vF
    use STRAIN, only : HAM_STRAIN, HAM_DERIVATIVES_STRAIN_PHI, HAM_DERIVATIVES_STRAIN_EPSILON
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF, CALCULATE_OM_CUTOFF_KET
    use omp_lib, only : omp_set_num_threads, omp_set_dynamic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none

    integer, parameter :: nstrain = 2*NTBG
    integer :: idir
    real(8) :: k(2), kp(2), km(2), param(3), pp(3), pm(3), step
    real(8) :: mx, my, mx_fd, my_fd, mx_parallel, my_parallel
    complex(8), allocatable :: hx(:,:), hy(:,:), hl(:,:), hx_only(:,:), hy_only(:,:)
    complex(8), allocatable :: hp(:,:), hm(:,:), fd(:,:)

    call omp_set_dynamic(.false.)
    call omp_set_num_threads(4)
    allocate(hx(nstrain,nstrain), hy(nstrain,nstrain), hl(nstrain,nstrain))
    allocate(hx_only(nstrain,nstrain), hy_only(nstrain,nstrain))
    allocate(hp(nstrain,nstrain), hm(nstrain,nstrain), fd(nstrain,nstrain))
    k = [1.d7, -2.d7]
    param = [1.4d0*PI/180.d0, 0.001d0, PI/6.d0]

    ! Omitting the optional lambda output must give the same velocity matrices.
    call HAM_DERIVATIVES_STRAIN_PHI(k, param, hx, hy, hl)
    call HAM_DERIVATIVES_STRAIN_PHI(k, param, hx_only, hy_only)
    if (maxval(abs(hx-hx_only)) > 1.d-14*vF) error stop 'phi: optional output changes hx'
    if (maxval(abs(hy-hy_only)) > 1.d-14*vF) error stop 'phi: optional output changes hy'

    step = 1.d-5
    pp = param
    pm = param
    pp(3) = pp(3) + step
    pm(3) = pm(3) - step
    call HAM_STRAIN(k, pp, hp)
    call HAM_STRAIN(k, pm, hm)
    fd = (hp-hm)/(2.d0*step)
    if (maxval(abs(fd-hl)) > 1.d-7*maxval(abs(hl))) error stop 'phi derivative mismatch'

    call HAM_DERIVATIVES_STRAIN_EPSILON(k, param, hx, hy, hl)
    call HAM_DERIVATIVES_STRAIN_EPSILON(k, param, hx_only, hy_only)
    if (maxval(abs(hx-hx_only)) > 1.d-14*vF) error stop 'epsilon: optional output changes hx'
    if (maxval(abs(hy-hy_only)) > 1.d-14*vF) error stop 'epsilon: optional output changes hy'
    step = 1.d-7
    pp = param
    pm = param
    pp(2) = pp(2) + step
    pm(2) = pm(2) - step
    call HAM_STRAIN(k, pp, hp)
    call HAM_STRAIN(k, pm, hm)
    fd = (hp-hm)/(2.d0*step)
    if (maxval(abs(fd-hl)) > 1.d-7*maxval(abs(hl))) error stop 'epsilon derivative mismatch'

    step = 1.d4
    do idir = 1, 2
        kp = k
        km = k
        kp(idir) = kp(idir) + step
        km(idir) = km(idir) - step
        call HAM_STRAIN(kp, param, hp)
        call HAM_STRAIN(km, param, hm)
        fd = (hp-hm)/(2.d0*step)
        if (idir == 1) then
            if (maxval(abs(fd-hx)) > 1.d-8*vF) error stop 'kx derivative mismatch'
        else
            if (maxval(abs(fd-hy)) > 1.d-8*vF) error stop 'ky derivative mismatch'
        endif
    enddo
    print *, 'PASS: optional and full strain derivatives, checked by finite differences'

    ! Full band window: compare the ket sum with the existing projector formula.
    ! More than one k point is required to detect a missing per-point reset.
    call omp_set_num_threads(1)
    call CALCULATE_OM_CUTOFF_KET(4, 3, 1, [1.d0,0.d0], [0.d0,1.d0], 0.d0, 1.d0, &
        [-1.d0,-1.d0,1.d0,1.d0], [0.d0], LINEAR_HAM, LINEAR_DERIVATIVES, mx, my)
    call CALCULATE_OM_CUTOFF(4, 3, 1, [1.d0,0.d0], [0.d0,1.d0], 1.d-4, 0.d0, 1.d0, &
        [-1.d0,-1.d0,1.d0,1.d0], [0.d0], LINEAR_HAM, mx_fd, my_fd)
    if (.not. all(ieee_is_finite([mx,my]))) error stop 'nonfinite ket magnetization'
    if (max(abs(mx-mx_fd),abs(my-my_fd)) > 1.d-8) error stop 'ket/projector mismatch'

    call omp_set_num_threads(4)
    call CALCULATE_OM_CUTOFF_KET(4, 3, 1, [1.d0,0.d0], [0.d0,1.d0], 0.d0, 1.d0, &
        [-1.d0,-1.d0,1.d0,1.d0], [0.d0], LINEAR_HAM, LINEAR_DERIVATIVES, mx_parallel, my_parallel)
    if (.not. all(ieee_is_finite([mx_parallel,my_parallel]))) error stop 'nonfinite parallel result'
    if (max(abs(mx-mx_parallel),abs(my-my_parallel)) > 1.d-13) error stop 'thread-dependent ket result'
    print *, 'PASS: ket/projector agreement and consistent results with 1 and 4 threads'

contains

    subroutine LINEAR_HAM(k, ham_param, ham)
        real(8), intent(in) :: k(2), ham_param(:)
        complex(8), intent(out) :: ham(:,:)
        complex(8) :: dx(4,4), dy(4,4)
        call LINEAR_DERIVATIVES(k, ham_param, dx, dy)
        ham = k(1)*dx + k(2)*dy
        ham(1,1) = -2.d0
        ham(2,2) = -0.4d0
        ham(3,3) = 0.7d0
        ham(4,4) = 3.d0
        ham(1,3) = ham(1,3) + (0.15d0,0.1d0)
        ham(3,1) = conjg(ham(1,3))
        ham(2,4) = ham(2,4) + (0.12d0,-0.08d0)
        ham(4,2) = conjg(ham(2,4))
    end subroutine LINEAR_HAM

    subroutine LINEAR_DERIVATIVES(k, ham_param, dx, dy, dl)
        real(8), intent(in) :: k(2), ham_param(:)
        complex(8), intent(out) :: dx(:,:), dy(:,:)
        complex(8), intent(out), optional :: dl(:,:)
        dx = (0.d0,0.d0)
        dy = (0.d0,0.d0)
        dx(1,3) = (0.2d0,0.1d0)
        dx(2,4) = (-0.3d0,0.05d0)
        dx(3,1) = conjg(dx(1,3))
        dx(4,2) = conjg(dx(2,4))
        dy(1,3) = (0.1d0,-0.2d0)
        dy(2,4) = (0.25d0,0.1d0)
        dy(3,1) = conjg(dy(1,3))
        dy(4,2) = conjg(dy(2,4))
        if (present(dl)) dl = (0.d0,0.d0)
    end subroutine LINEAR_DERIVATIVES

end program CHECK_OM_KET
