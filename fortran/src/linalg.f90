module LINALG

    use stdlib_linalg_blas, only : stdlib_zdotc
    implicit none

    contains

!-------------------------------------------------
! Hermitian 矩阵特征值分解 (ZHEEV)
!-------------------------------------------------
subroutine myzheev(N, W, H, V)
    implicit none
    
    integer,intent(in) :: N
    real(8),intent(out) :: W(:)
    complex(8),intent(in) :: H(:,:)
    complex(8), intent(out) :: V(:,:) ! 特征向量矩阵（列向量为特征向量）
    integer :: INFO
    real(8) :: RWORK(3*N-2)
    complex(8) :: WORK(3*N)
    complex(8) :: H_TEMP(N,N) 
    
    H_TEMP = H ! 复制输入矩阵 H 到 H_TEMP（避免破坏原矩阵）
    call ZHEEV('V', 'U', N, H_TEMP, N, W, WORK, 3*N, RWORK, INFO) 
    V = H_TEMP
    
    ! 检查 INFO
    if (INFO /= 0) then
        print *, "ZHEEV 错误: INFO =", INFO
    end if
    
endsubroutine myzheev
    
!-------------------------------------------------
! 非对称复数矩阵特征值分解 (ZGEEV)
!-------------------------------------------------
subroutine myzgeev(N, W, U, VR)
    implicit none
    
    integer, intent(in) :: N
    complex(8), intent(out) :: W(N)       
    complex(8), intent(in) :: U(N,N)      
    complex(8), intent(out) :: VR(N,N)    ! 右特征向量矩阵（列向量为特征向量）
    integer :: INFO
    complex(8) :: A(N,N)                  
    complex(8) :: VL(N,N)                 
    real(8) :: RWORK(2*N)
    complex(8) :: WORK(4*N)               
    
    ! 复制输入矩阵 H 到 A（避免破坏原矩阵）
    A = U
    call ZGEEV( 'N', 'V', N, A, N, W, VL, N, VR, N, WORK, 4*N, RWORK, INFO) 
    
    ! 检查 INFO 是否成功
    if (INFO /= 0) then
        print *, "ZGEEV 错误: INFO =", INFO
    end if
end subroutine myzgeev

!-------------------------------------------------
! 实数矩阵奇异值分解 (DGESVD)
!-------------------------------------------------
subroutine mydgesvd(M, N, A, S, U, VT)
    implicit none
    integer, intent(in) :: M, N
    real(8), intent(in) :: A(M, N)      ! 输入矩阵
    real(8), intent(out) :: S(min(M,N)) ! 奇异值数组
    real(8), intent(out) :: U(M, M)     ! 左奇异向量矩阵（全尺寸）
    real(8), intent(out) :: VT(N, N)    ! 右奇异向量矩阵的转置（全尺寸）
    
    real(8) :: A_copy(M, N)
    integer :: INFO, LWORK, min_mn, max_mn
    real(8), allocatable :: WORK(:)
    
    ! 计算最小工作空间大小
    min_mn = min(M, N)
    max_mn = max(M, N)
    LWORK = max(3*min_mn + max_mn, 5*min_mn)  ! LAPACK建议的最小值
    allocate(WORK(LWORK))
    
    ! 复制矩阵以避免覆盖输入
    A_copy = A
    
    ! 调用 DGESVD
    call DGESVD('A', 'A', M, N, A_copy, M, S, U, M, VT, N, WORK, LWORK, INFO)
    
    ! 错误处理
    if (INFO /= 0) then
        print *, "DGESVD 错误: INFO =", INFO
    end if
    
    deallocate(WORK)
end subroutine mydgesvd

!-------------------------------------------------
! 复数矩阵奇异值分解 (ZGESVD)
!-------------------------------------------------
subroutine myzgesvd(M, N, A, S, U, VT)
    implicit none
    integer, intent(in) :: M, N
    complex(8), intent(in) :: A(M, N)   ! 输入矩阵
    real(8), intent(out) :: S(min(M,N)) ! 奇异值数组
    complex(8), intent(out) :: U(M, M)  ! 左奇异向量矩阵（全尺寸）
    complex(8), intent(out) :: VT(N, N) ! 右奇异向量矩阵的转置（全尺寸）
    
    complex(8) :: A_copy(M, N)
    integer :: INFO, LWORK, min_mn, max_mn
    complex(8), allocatable :: WORK(:)
    real(8) :: RWORK(5*min(M,N))         ! RWORK 大小为至少 5*min(M,N)
    
    ! 计算最小工作空间大小
    min_mn = min(M, N)
    max_mn = max(M, N)
    LWORK = max(2*min_mn + max_mn, 2*min_mn + M, 2*min_mn + N)  ! LAPACK建议的最小值
    allocate(WORK(LWORK))
    
    ! 复制矩阵以避免覆盖输入
    A_copy = A
    
    ! 调用 ZGESVD
    call ZGESVD('A', 'A', M, N, A_copy, M, S, U, M, VT, N, WORK, LWORK, RWORK, INFO)
    
    ! 错误处理
    if (INFO /= 0) then
        print *, "ZGESVD 错误: INFO =", INFO
    end if
    
    deallocate(WORK)

end subroutine myzgesvd


! ---------------------
! 复数向量共轭内积 
! dot_product(版本)
! ---------------------
subroutine mydot(N, X, Y, DOTC)
    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: X(:)
    complex(8), intent(in) :: Y(:)

    complex(8), intent(out) :: DOTC

    ! Dimension check
    if (size(X) /= N) then
        error stop 'mydot: X dimension error'
    end if

    if (size(Y) /= N) then
        error stop 'mydot: Y dimension error'
    end if

    ! -----
    ! DOTC = sum_i conjg(X(i)) * Y(i)
    ! -----
    DOTC = dot_product(X, Y)

end subroutine mydot


! -------------------
! 复数向量共轭内积
! zdotc版本
! -------------------
subroutine myzdotc(N, X, Y, DOTC)
    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: X(:)
    complex(8), intent(in) :: Y(:)

    complex(8), intent(out) :: DOTC

    if (size(X) /= N .or. size(Y) /= N) then
        error stop "myzdotc: dimension error"
    end if

    DOTC = stdlib_zdotc(N, X, 1, Y, 1)

end subroutine myzdotc




    
    



    
end module LINALG