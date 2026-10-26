module STRAIN
    use CONSTANTS
    use TBG, only : DIAGONAL, ELEMENTS, ASMBL
    use MODEL_INTERFACE, only : HAMILTONIAN 
    implicit none

    ! ----- STRAIN PARAMETERS ----- 
    real(8), parameter :: nu = 0.16d0
    real(8), parameter :: alpha = -1.0542d0*eV
    real(8), parameter :: gamma_vF = -3.3644d0*eV ! 请注意，这里是gamma*vF
    
contains

! ==========
! Pseudovector fields and Pseudoscalar fields
! ==========
subroutine PSEUDOFIELDS(epsilon1, strain_phi, S, vph, A)
    implicit none
    
    real(8), intent(in) :: epsilon1
    real(8), intent(in) :: strain_phi
    real(8), intent(out) :: vph
    real(8), intent(out) :: A(2)
    real(8), intent(out) :: S(2,2)

    real(8) :: epsilon2

    epsilon2 = -nu * epsilon1

    ! ----- pesudoscalar fields ----- 
    vph = epsilon1 + epsilon2
    
    ! ----- pesudovector fields ----- 
    A(1) = (epsilon1 - epsilon2)*cos(2.d0*strain_phi)
    A(2) = (epsilon1 - epsilon2)*sin(2.d0*strain_phi)

    ! ----- strain tensor ----- 
    S(1,1) = epsilon1*(cos(strain_phi)**2) + epsilon2*(sin(strain_phi)**2)
    S(1,2) = (-epsilon1 + epsilon2)*cos(strain_phi)*sin(strain_phi)
    S(2,1) = S(1,2)
    S(2,2) = epsilon1*(sin(strain_phi)**2) + epsilon2*(cos(strain_phi)**2)

end subroutine PSEUDOFIELDS


! ==========
! Heterostrain下的形变（包含应力与转角）
! ==========
subroutine DEFORMATION(theta, S, E)
    implicit none
    
    real(8), intent(in) :: theta    
    real(8), intent(in) :: S(2,2)
    real(8), intent(out) :: E(2,2)

    ! ----- 加入转角 ----- 
    E(1, 1) = S(1, 1)
    E(1, 2) = S(1, 2) - theta
    E(2, 1) = S(2, 1) + theta
    E(2, 2) = S(2, 2)

end subroutine DEFORMATION


! ==========
! Moire reciprocal lattice vectors under heterostrain
! bh is the 2 times 2 array of b1h and b2h
! ==========
subroutine GEOMETRY_STRAIN(theta, strain_epsilon, strain_phi, bs, q1s)
    implicit none 

    real(8), intent(in) :: theta, strain_epsilon, strain_phi
    real(8), intent(out) :: bs(2,2)
    real(8), intent(out), optional :: q1s(2)

    real(8) :: S(2,2), E(2,2)
    real(8) :: vph, A(2)

    call PSEUDOFIELDS(strain_epsilon, strain_phi, S, vph, A)
    call DEFORMATION(theta, S, E)

    bs = matmul(transpose(E), b_graphene)
    
    if (present(q1s)) then
        q1s = matmul(transpose(E), K_graphene)
    endif

end subroutine GEOMETRY_STRAIN

! ==========
! 构造异质应力下的哈密顿量
! ==========
subroutine H_STRAIN(k, theta, epsilon1, strain_phi, H)
    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: theta
    real(8), intent(in) :: epsilon1
    real(8), intent(in) :: strain_phi
    complex(8), intent(out) :: H(:, :)

    complex(8), allocatable :: tH(:,:,:,:)

    real(8) :: p(2)
    real(8) :: S(2,2)
    real(8) :: q1s(2), bs(2, 2) ! geometry under heterostrain
    real(8) :: vph, A(2)
    integer :: layer1, layer2, n1, n2
    complex(8) :: M(2, 2)

    ! ----- 赝磁场 -----    
    call PSEUDOFIELDS(epsilon1, strain_phi, S, vph, A)
    ! ----- moire reciprocal lattice vectors -----  
    call GEOMETRY_STRAIN(theta, epsilon1, strain_phi, bs, q1s)

    allocate(tH(NTBG, NTBG, 2, 2)) 
    tH = (0.d0, 0.d0)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(p, layer1, layer2, n1, n2, M) &
    !$OMP SHARED(k, q1s, bs, vph, A, tH) COLLAPSE(4) SCHEDULE(dynamic)
    do layer1 = 1, 2
        do layer2 = 1, 2
            do n1 = -tr, tr
                do n2 = -tr, tr
                    ! diagonal term
                    if  (layer1 ==1 .and. layer2 == 1) then ! bottom layer
                        p = k - q1s - n1*bs(:,1) - n2*bs(:,2)
                        call DIAGONAL(p*vF - gamma_vF*A, M)
                        M = M - alpha * vph * sigma_0
                        call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2, M)
                    elseif (layer1 ==2 .and. layer2 == 2) then ! top layer
                        p = k - n1*bs(:,1) - n2*bs(:,2)
                        call DIAGONAL(p*vF + gamma_vF*A, M)
                        M = M + alpha * vph * sigma_0
                        call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2, M)
                    elseif (layer1 == 1 .and. layer2 == 2) then
                        call ELEMENTS(tH, layer1, n1, n2, layer2, n1, n2,   T1)        ! Δ(0,0)
                        if (abs(n2 - 1) <= tr) call ELEMENTS(tH, layer1, n1, n2, layer2, n1,   n2-1, T2) ! Δ(0,-1)
                        if (abs(n1 + 1) <= tr) call ELEMENTS(tH, layer1, n1, n2, layer2, n1+1, n2,   T3) ! Δ(+1,0)
                    end if
                end do
            end do
        end do
    end do
    !$OMP END PARALLEL DO 

    H = 0.d0
    call ASMBL(NTBG, H, tH)
    deallocate(tH)

end subroutine H_STRAIN

    

! ==========
! Strain Hamiltonian interface
! ==========
subroutine HAM_STRAIN(k, ham_param, HAM)
    implicit none

    real(8), intent(in) :: k(2)
    real(8), intent(in) :: ham_param(:)

    complex(8), intent(out) :: HAM(:,:)
    
    real(8) :: theta, epsilon1, strain_phi

    if (size(ham_param) /= 3) then
        error stop "HAM_STRAIN: ham_param dimension error"
    endif

    theta = ham_param(1)
    epsilon1 = ham_param(2)
    strain_phi = ham_param(3)

    call H_STRAIN(k, theta, epsilon1, strain_phi, HAM)

end subroutine HAM_STRAIN


! ==========
! KPATH STRAIN : K_b -> Gamma -> M -> K_t
! ==========
subroutine KPATH_STRAIN(numk, theta, epsilon1, strain_phi, k, dk)
    implicit none
    
    integer, intent(in) :: numk
    real(8), intent(in) :: theta
    real(8), intent(in) :: epsilon1
    real(8), intent(in) :: strain_phi
    real(8), intent(out) :: k(2,3*numk+1)
    real(8), intent(out) :: dk(3*numk+1)

    integer :: i, j, nk
    real(8) :: t
    real(8) :: P(2, 4)
    real(8) :: q1s(2), bs(2, 2)

    call GEOMETRY_STRAIN(theta, epsilon1, strain_phi, bs, q1s)

    P(:,1) = q1s
    P(:,2) = (2.d0/3.d0)*bs(:,1) + (1.d0/3.d0)*bs(:,2)
    P(:,3) = 0.5d0*q1s
    P(:,4) = 0.d0

    nk = 1
    k(:, 1) = P(:, 1)
    do i = 1, 3
        do j = 1, numk
            t = real(j,8)/real(numk,8)

            nk = nk + 1
            k(:,nk) = (1.d0 - t)*P(:,i) + t*P(:,i+1)
        enddo
    enddo

    dk(1) = 0.d0

    do i = 2, 3*numk+1
        dk(i) = dk(i-1) + sqrt(sum((k(:,i) - k(:,i-1))**2))
    enddo

end subroutine KPATH_STRAIN


! ==========
! Derivatives of pesudofields and strain tensor
! ==========
subroutine PSEUDOFIELDS_DERIVATIVES(strain_epsilon, strain_phi, lambda_name, &
                                    dS_dlambda, dvph_dlambda, dA_dlambda)
    implicit none
    
    real(8), intent(in) :: strain_epsilon
    real(8), intent(in) :: strain_phi
    character(len=*), intent(in) :: lambda_name ! epsilon or phi

    real(8), intent(out) :: dS_dlambda(2, 2)
    real(8), intent(out) :: dvph_dlambda
    real(8), intent(out) :: dA_dlambda(2)

    select case (lambda_name) !相当于if ... else ...

    ! ----- Derivatives with respect to epsilon  ----- 
    case ('strain_epsilon')
        dvph_dlambda = 1.d0 - nu

        dA_dlambda(1) = (1.d0+nu)*cos(2.d0*strain_phi)
        dA_dlambda(2) = (1.d0+nu)*sin(2.d0*strain_phi)

        dS_dlambda(1,1) = cos(strain_phi)**2-nu*sin(strain_phi)**2
        dS_dlambda(1,2) = -(1.d0+nu)*cos(strain_phi)*sin(strain_phi)
        dS_dlambda(2,1) = dS_dlambda(1,2)
        dS_dlambda(2,2) = sin(strain_phi)**2-nu*cos(strain_phi)**2

    ! ----- Derivatives with respect to phi ----- 
    case ("strain_phi")
        dvph_dlambda = 0.d0

        dA_dlambda(1) = -2.d0*(1.d0+nu)*strain_epsilon*sin(2.d0*strain_phi)
        dA_dlambda(2) =  2.d0*(1.d0+nu)*strain_epsilon*cos(2.d0*strain_phi)

        dS_dlambda(1,1) = -(1.d0+nu)*strain_epsilon*sin(2.d0*strain_phi)
        dS_dlambda(1,2) = -(1.d0+nu)*strain_epsilon*cos(2.d0*strain_phi)
        dS_dlambda(2,1) = dS_dlambda(1,2)
        dS_dlambda(2,2) =  (1.d0+nu)*strain_epsilon*sin(2.d0*strain_phi)

    case DEFAULT

        error stop "PSEUDOFIELDS_DERIVATIVES: invalid lambda_name"

    end select 
    
end subroutine PSEUDOFIELDS_DERIVATIVES


! ==========
! Derivatives under heterostrain
! Hx = dH/dkx, Hy = dH/dky, Hlambda = dH/dlambda
! lambda = epsilon or phi 
! ==========
subroutine H_STRAIN_DERIVATIVES(strain_epsilon, strain_phi, lambda_name,& 
                                dH_dkx, dH_dky, dH_dlambda)
    implicit none
    
    real(8), intent(in) :: strain_epsilon
    real(8), intent(in) :: strain_phi
    character(len = *), intent(in) :: lambda_name   

    complex(8), intent(out) :: dH_dkx(:,:)
    complex(8), intent(out) :: dH_dky(:,:)
    complex(8), intent(out), optional :: dH_dlambda(:,:)

    complex(8), allocatable :: tHx(:,:,:,:) ! Block tensor for dH/dkx
    complex(8), allocatable :: tHy(:,:,:,:)
    complex(8), allocatable :: tHlambda(:,:,:,:)

    real(8) :: p(2)
    real(8) :: dS_dlambda(2,2)
    real(8) :: dvph_dlambda
    real(8) :: dA_dlambda(2)

    integer :: layer1, layer2
    integer :: n1, n2

    complex(8) :: M(2,2)
    complex(8) :: Mx(2,2), My(2,2)
    logical :: need_lambda

    
    ! ----------
    ! 判断是否需要计算dH_dlambda
    ! ----------
    need_lambda = present(dH_dlambda)
    if (need_lambda) then
        ! ----- derivatives of pesudofields and strain tensor ----- 
        call PSEUDOFIELDS_DERIVATIVES(strain_epsilon, strain_phi, lambda_name, &
                                    dS_dlambda, dvph_dlambda, dA_dlambda)
        allocate(tHlambda(NTBG, NTBG, 2, 2))
        tHlambda = 0.d0
    end if


    ! ----------
    ! dH/dkx  = hbar vF sigma_1
    ! dH/dky = hbar vF sigma_2
    ! ----------
    call DIAGONAL([vF, 0.d0], Mx)
    call DIAGONAL([0.d0, vF], My)

    ! ----- allocate ----- 
    allocate(tHx(NTBG, NTBG, 2, 2))
    allocate(tHy(NTBG, NTBG, 2, 2))
    

    ! ----- 初始化 ----- 
    tHx = 0.d0
    tHy = 0.d0

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP PRIVATE(layer1,layer2,n1,n2,p,M) &
    !$OMP SHARED(dS_dlambda,dA_dlambda,dvph_dlambda,Mx,My,tHx,tHy,tHlambda,need_lambda) &
    !$OMP COLLAPSE(4) SCHEDULE(dynamic)
    do layer1 = 1, 2
        do layer2 = 1, 2
            do n1 = -tr, tr
                do n2 = -tr, tr
                    ! bottom layer
                    if (layer1 == 1 .and. layer2 == 1) then
                        if (need_lambda) then
                            p = K_graphene + n1*b_graphene(:,1) + n2*b_graphene(:,2)
                            call DIAGONAL(-vF*matmul(dS_dlambda, p) - gamma_vF*dA_dlambda, M)
                            M = M - alpha * dvph_dlambda*sigma_0
                            call ELEMENTS(tHlambda, layer1, n1, n2, layer2, n1, n2, M)
                        endif
                        ! ----- dH_bot/dkx, dH_bot/dky ----- 
                        call ELEMENTS(tHx, layer1, n1, n2, layer2, n1, n2, Mx)
                        call ELEMENTS(tHy, layer1, n1, n2, layer2, n1, n2, My)
                    ! top layer
                    else if (layer1 == 2 .and. layer2 == 2) then 
                        if (need_lambda) then
                            p = n1*b_graphene(:,1) + n2*b_graphene(:,2)
                            call DIAGONAL(-vF*matmul(dS_dlambda, p) + gamma_vF*dA_dlambda, M)
                            M = M + alpha * dvph_dlambda*sigma_0
                            call ELEMENTS(tHlambda, layer1, n1, n2, layer2, n1, n2, M)
                        endif
                        ! ----- dH_top/dkx, dH_top/dky ----- 
                        call ELEMENTS(tHx, layer1, n1, n2, layer2, n1, n2, Mx)
                        call ELEMENTS(tHy, layer1, n1, n2, layer2, n1, n2, My)
                    end if

                enddo
            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    ! ----- Assemble full derivative Hamiltonians -----
    dH_dkx = 0.d0
    dH_dky = 0.d0

    call ASMBL(NTBG, dH_dkx, tHx)
    call ASMBL(NTBG, dH_dky, tHy)
    
    if (present(dH_dlambda)) then
        dH_dlambda = (0.d0,0.d0)
        call ASMBL(NTBG, dH_dlambda, tHlambda)
        deallocate(tHlambda)
    endif

    ! ----- 释放数组 ----- 
    deallocate(tHx, tHy)

end subroutine H_STRAIN_DERIVATIVES


! ==========
! Strain Hamiltonian derivatives interface
! lambda = strain_epsilon
! ==========
subroutine HAM_DERIVATIVES_STRAIN_EPSILON(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
    implicit none
    
    real(8), intent(in) :: k(2) 
    real(8), intent(in) :: ham_param(:)

    complex(8), intent(out) :: dH_dkx(:,:), dH_dky(:,:)
    complex(8), intent(out), optional ::  dH_dlambda(:,:)

    call H_STRAIN_DERIVATIVES(ham_param(2), ham_param(3), 'strain_epsilon', &
                              dH_dkx, dH_dky, dH_dlambda)

end subroutine HAM_DERIVATIVES_STRAIN_EPSILON


! ==========
! Strain Hamiltonian derivatives interface
! lambda = strain_phi   
! ==========
subroutine HAM_DERIVATIVES_STRAIN_PHI(k, ham_param, dH_dkx, dH_dky, dH_dlambda)
    implicit none
    
    real(8), intent(in) :: k(2) 
    real(8), intent(in) :: ham_param(:)

    complex(8), intent(out) :: dH_dkx(:,:), dH_dky(:,:)
    complex(8), intent(out), optional :: dH_dlambda(:,:)

    call H_STRAIN_DERIVATIVES(ham_param(2), ham_param(3), 'strain_phi', &
                              dH_dkx, dH_dky, dH_dlambda)

end subroutine HAM_DERIVATIVES_STRAIN_PHI


end module STRAIN
