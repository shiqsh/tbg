program CHECK_LINALG

    use LINALG, only : mydot, myzdotc

    implicit none

    integer, parameter :: N = 3
    
    complex(8) :: X(N)
    complex(8) :: Y(N)

    complex(8) :: DOT_MYDOT
    complex(8) :: DOT_ZDOTC

! -----------
! Define test vectors
! X = (1+i, 2, -i)
! Y = (1, i, 2 + i)
! Exact the inner product is 3i 
! -----------
    X(1) = (1.d0, 1.d0);    X(2) = (2.d0, 0.d0);    X(3) = (0.d0, -1.d0)
    Y(1) = (1.d0, 0.d0);    Y(2) = (0.d0, 1.d0);    Y(3) = (2.d0, 1.d0)

    call mydot(N, X, Y, DOT_MYDOT)
    call myzdotc(N, X, Y, DOT_ZDOTC)

    print *, "mydot result:"
    print *, DOT_MYDOT
    print *
    print *, "myzdotc result:"
    print *, DOT_ZDOTC

end program CHECK_LINALG
    