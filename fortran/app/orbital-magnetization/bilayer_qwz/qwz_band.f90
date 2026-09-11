program QWZ_BAND
    use BILAYER_QWZ, only : H_BILAYER_QWZ, KPATH_QWZ
    use LINALG, only : myzheev
    implicit none 

    integer, parameter :: numk = 100
    integer, parameter :: nk = 3*numk + 1

    real(8), parameter :: m = 2.0d0
    real(8), parameter :: lambda = 0.0d0
    real(8), parameter :: t_perp = 0.5d0
    
    integer :: i, fileunit
    real(8) :: k(2, nk), dk(nk)
    real(8) :: W(4)
    complex(8) :: H(4,4), V(4,4)
    character(len=256) :: filename
    character(len=16) :: s_lambda, s_m, s_tperp

    write(s_lambda, '(F8.3)') lambda
    write(s_m, '(F8.3)') m
    write(s_tperp, '(F8.3)') t_perp

    filename = '../data/orbital-magnetization/qwz/fortran/band/qwz_band_lambda_'//trim(adjustl(s_lambda))//'_m_'//trim(adjustl(s_m))//'_tperp_'//trim(adjustl(s_tperp))//'.dat'

    call KPATH_QWZ(numk, k, dk)

    open(newunit=fileunit, file=trim(filename), status='replace')

    write(fileunit, '(A, ES32.16)') '# lambda = ', lambda
    write(fileunit, '(A, ES32.16)') '# m = ', m
    write(fileunit, '(A, ES32.16)') '# t_perp = ', t_perp

    do i = 1, nk
        call H_BILAYER_QWZ(k(:,i), lambda, m, t_perp, H)
        call myzheev(4, W, H, V)

        write(fileunit, '(7ES32.16)')dk(i), k(1,i), k(2,i), W
    enddo

    close(fileunit)


end program QWZ_BAND


        