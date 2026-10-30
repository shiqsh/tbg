program CHECK_OM_MAGNETIC
    use CONSTANTS, only : NTBG,PI,hbar,m_e,eV
    use STRAIN, only : HAM_STRAIN,HAM_DERIVATIVES_STRAIN_PHI,GEOMETRY_STRAIN
    use TBG, only : Z_TBG
    use LINALG, only : myzheev
    use OM_MAGNETIC, only : CALCULATE_OM_MAGNETIC,CALCULATE_OMEGA_B
    use omp_lib, only : omp_set_num_threads
    implicit none
    integer, parameter :: N=2*NTBG,nk=2
    real(8) :: param(3),bs(2,2),z(N),mu,k(2),hf(2),fd(2),tmp(2),B(2),w(N)
    real(8) :: omega,omega_direct,pref
    complex(8), allocatable :: ham(:,:),v(:,:),hx(:,:),hy(:,:)
    integer :: iu,iv,a,i
    call omp_set_num_threads(2)
    allocate(ham(N,N),v(N,N),hx(N,N),hy(N,N))
    param=[1.4d0*PI/180.d0,0.001d0,PI/4.d0]
    mu=0.5d-3*eV
    call GEOMETRY_STRAIN(param(1),param(2),param(3),bs)
    call Z_TBG(3.35d-10,z)
    call HAM_DERIVATIVES_STRAIN_PHI([0.d0,0.d0],param,hx,hy)
    hf=0.d0; omega_direct=0.d0; B=[0.013d0,-0.007d0]
    do iu=1,nk
        do iv=1,nk
            k=(-0.5d0+real(iu-1,8)/nk)*bs(:,1)+(-0.5d0+real(iv-1,8)/nk)*bs(:,2)
            call HAM_STRAIN(k,param,ham)
            call myzheev(N,w,ham,v)
            do a=1,count(w<mu)
                tmp=[real(dot_product(v(:,a),z*matmul(hy,v(:,a))),8), &
                    -real(dot_product(v(:,a),z*matmul(hx,v(:,a))),8)]
                hf=hf+tmp
            enddo
            ! Reproduce the user's direct H(k)+z*(By Hx-Bx Hy) prescription.
            do i=1,N
                ham(i,:)=ham(i,:)+(eV/hbar)*z(i)*(B(2)*hx(i,:)-B(1)*hy(i,:))
            enddo
            if(maxval(abs(ham-conjg(transpose(ham))))>1.d-12*maxval(abs(ham))) &
                error stop 'field Hamiltonian is not Hermitian'
            call myzheev(N,w,ham,v)
            omega_direct=omega_direct+sum(min(w-mu,0.d0))
        enddo
    enddo
    pref=2.d0*m_e/hbar**2
    hf=pref*hf/real(nk,8)**2
    omega_direct=omega_direct/real(nk,8)**2
    call CALCULATE_OMEGA_B(N,nk,bs(:,1),bs(:,2),mu,z,param,B,omega)
    if(abs(omega-omega_direct)>1.d-12*abs(omega_direct)) error stop 'cached/direct Omega mismatch'
    call CALCULATE_OM_MAGNETIC(N,nk,bs(:,1),bs(:,2),mu,z,param,0.01d0,fd(1),fd(2))
    if(maxval(abs(fd-hf))>5.d-6) error stop 'Omega/Hellmann-Feynman mismatch'
    write(*,'(A,2ES16.7)') 'PASS: magnetic FD minus HF (mu_B/cell): ',fd-hf
    print *, 'PASS: direct finite-field Hamiltonian and cached Omega agree'
end program CHECK_OM_MAGNETIC
