program test_environment

    use stdlib_kinds, only: dp
    use omp_lib, only: omp_get_num_threads

    implicit none

    integer :: nthreads
    integer :: info
    integer, parameter :: n = 2
    integer, parameter :: lwork = 4

    complex(dp) :: H(n,n)
    complex(dp) :: work(lwork)

    real(dp) :: eigenvalues(n)
    real(dp) :: rwork(3*n-2)

    interface
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            import :: dp

            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork

            complex(dp), intent(inout) :: a(lda,*)
            real(dp), intent(out) :: w(*)
            complex(dp), intent(out) :: work(*)
            real(dp), intent(out) :: rwork(*)

            integer, intent(out) :: info
        end subroutine zheev
    end interface


    print *, "======================================"
    print *, "TBG Fortran environment test"
    print *, "======================================"

    ! --------------------------------------------------
    ! Test 1: stdlib
    ! --------------------------------------------------

    print *
    print *, "[1] stdlib test"
    print *, "dp kind =", dp


    ! --------------------------------------------------
    ! Test 2: OpenMP
    ! --------------------------------------------------

    print *
    print *, "[2] OpenMP test"

    nthreads = 0

    !$omp parallel
    !$omp single
        nthreads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    print *, "OpenMP threads =", nthreads


    ! --------------------------------------------------
    ! Test 3: LAPACK / Apple Accelerate
    ! --------------------------------------------------

    print *
    print *, "[3] LAPACK ZHEEV test"

    H = cmplx(0.0_dp, 0.0_dp, kind=dp)

    H(1,1) = 1.0_dp
    H(2,2) = 2.0_dp

    H(1,2) = cmplx(0.0_dp, 1.0_dp, kind=dp)
    H(2,1) = conjg(H(1,2))

    call zheev( &
        'V', &
        'U', &
        n, &
        H, &
        n, &
        eigenvalues, &
        work, &
        lwork, &
        rwork, &
        info &
    )

    if (info /= 0) then
        print *, "ZHEEV failed."
        print *, "INFO =", info
        error stop
    endif

    print *, "Eigenvalues:"
    print *, eigenvalues

    print *
    print *, "All tests passed."

end program test_environment