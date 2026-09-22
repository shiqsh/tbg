program CHECK_STRAIN_BAND
    use CONSTANTS, only : PI, eV, NTBG
    use STRAIN, only : H_STRAIN, KPATH_STRAIN
    use LINALG, only : myzheev
    implicit none
    
    integer, parameter :: numk = 100
    integer, parameter :: nk = 3*numk + 1
    integer, parameter :: NH = 2*NTBG

    real(8) :: theta
    real(8) :: epsilon1, strain_phi
    real(8) :: k(2, nk), dk(nk)
    real(8) :: W(NH)
    complex(8) :: H(NH, NH), V(NH, NH)

    integer :: ik, iband
    integer :: unit_band

    ! ----- output file ----- 
    character(len=*), parameter :: filename = "/Users/shiqsh/Workspace/tbg/data/tbg/fortran/band/check_strain_band.dat"
    
    ! ----------
    ! strain parameters
    ! ----------
    theta = 1.05d0*PI/180.d0
    epsilon1 = 0.d0
    strain_phi = 0.d0*PI/180.d0

    ! ==========
    ! K-PATH
    ! ==========
    call KPATH_STRAIN( numk, theta, epsilon1, strain_phi, k, dk)

    ! ==========
    ! Output the data
    ! ==========
    open(newunit = unit_band, file = filename, status = "replace", action = "write")

    do ik = 1, nk

        call H_STRAIN(k(:,ik), theta, epsilon1, strain_phi, H)
        call myzheev(NH, W, H, V)

        write(unit_band, '(*(ES32.16, 1X))')dk(ik), (W(iband)/eV, iband=1,NH)
    
    enddo

    close(unit_band)

    print *, "Band data written to:"
    print *, trim(filename)

end program CHECK_STRAIN_BAND

