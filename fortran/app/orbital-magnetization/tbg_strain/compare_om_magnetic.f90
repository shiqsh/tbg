! Independent audit: raw grand potential, Hellmann-Feynman, band-window formula,
! and the BZ boundary term. Does not change production OM parameters or data.
program COMPARE_OM_MAGNETIC
    use CONSTANTS, only : NTBG,tr,PI,hbar,m_e,eV
    use TBG, only : Z_TBG
    use STRAIN, only : HAM_STRAIN,GEOMETRY_STRAIN,HAM_DERIVATIVES_STRAIN_PHI
    use LINALG, only : myzheev
    use OM_CUTOFF, only : CALCULATE_OM_CUTOFF_KET, CALCULATE_OM_CUTOFF
    use OM_MAGNETIC, only : CALCULATE_OM_MAGNETIC
    use omp_lib, only : omp_get_wtime
    implicit none
    integer, parameter :: N=2*NTBG, nc=7
    integer :: nk, cuts(nc), iu,iv,a,b,c,no,lo,hi, idir, side, nedge, verify, j
    real(8) :: param(3),bs(2,2),z(N),mu,B0,pref,k(2),start,fd(2),o(4)
    real(8) :: sumf(2,nc), ff(2,nc), hf(2), dd(2), active(2), aa(2), factor, gap
    real(8) :: sumg(2), gdiff(2), g, edge(2), invbs(2,2),det,err, h, gp,gm, numerical(2)
    real(8) :: mx,my, mu_arg, eps_arg,phi_arg
    real(8), allocatable :: W(:)
    complex(8), allocatable :: H0(:,:),Href(:,:),V(:,:),hx(:,:),hy(:,:),xe(:,:),ye(:,:),ze(:,:),tmp(:,:)
    character(100) :: arg
    nk=12; B0=0.01d0; phi_arg=0.d0; eps_arg=0.1d0; mu_arg=0.5d0; verify=0
    do j=1,command_argument_count()
        call get_command_argument(j,arg)
        select case(j)
        case(1); read(arg,*) nk
        case(2); read(arg,*) B0
        case(3); read(arg,*) phi_arg
        case(4); read(arg,*) eps_arg
        case(5); read(arg,*) mu_arg
        case(6); read(arg,*) verify
        case default; error stop 'usage: compare_om_magnetic [nk B0_T phi_deg epsilon_percent mu_meV verify]'
        end select
    enddo
    if(nk<=0 .or. B0<0.d0) error stop 'invalid nk or B0 (0 skips magnetic FD)'
    start=omp_get_wtime()
    param=[1.4d0*PI/180.d0,eps_arg/100.d0,phi_arg*PI/180.d0]
    mu=mu_arg*1.d-3*eV
    pref=2.d0*m_e/hbar**2
    cuts=min([0,20,40,80,100,140,N/2-1],N/2-1)
    call Z_TBG(3.35d-10,z)
    call GEOMETRY_STRAIN(param(1),param(2),param(3),bs)
    allocate(Href(N,N),hx(N,N),hy(N,N))
    call HAM_STRAIN([0.d0,0.d0],param,Href)
    call HAM_DERIVATIVES_STRAIN_PHI([0.d0,0.d0],param,hx,hy)
    allocate(H0(N,N))
    k=0.13d0*bs(:,1)-0.27d0*bs(:,2)
    call HAM_STRAIN(k,param,H0)
    err=maxval(abs(H0-Href-k(1)*hx-k(2)*hy))/maxval(abs(H0))
    if(err>1.d-12) error stop 'H is not linear in k: cache invalid'
    do a=1,N
        if(maxval(abs((z(a)-z)*hx(a,:)))>1.d-12*maxval(abs(hx))*maxval(abs(z))) &
            error stop '[z,hx] nonzero'
        if(maxval(abs((z(a)-z)*hy(a,:)))>1.d-12*maxval(abs(hy))*maxval(abs(z))) &
            error stop '[z,hy] nonzero'
    enddo
    deallocate(H0)
    write(*,'(A,4I8)') '# tr N nk verify: ',tr,N,nk,verify
    write(*,'(A,4ES23.14)') '# B0_T phi_deg epsilon_percent mu_meV: ',B0,phi_arg,eps_arg,mu_arg
    write(*,'(A,ES23.14)') '# Hamiltonian cache relative error: ',err
    hf=0.d0; active=0.d0; sumf=0.d0
    !$omp parallel default(none) shared(Nk,bs,mu,z,Href,hx,hy,cuts) &
    !$omp private(iu,iv,k,W,H0,V,xe,ye,ze,tmp,a,b,c,no,lo,hi,dd,aa,ff,factor,gap) &
    !$omp reduction(+:hf,active,sumf)
    allocate(W(N),H0(N,N),V(N,N),xe(N,N),ye(N,N),ze(N,N),tmp(N,N))
    !$omp do collapse(2) schedule(static)
    do iu=1,nk
        do iv=1,nk
            k=(-0.5d0+real(iu-1,8)/nk)*bs(:,1)+(-0.5d0+real(iv-1,8)/nk)*bs(:,2)
            H0=Href+k(1)*hx+k(2)*hy
            call myzheev(N,W,H0,V)
            no=count(W<mu)
            if(no<N/2-1 .or. no>N/2+1) error stop 'mu outside active bands: original cutoff occupation invalid'
            do a=1,N
                tmp(:,a)=z*V(:,a)
            enddo
            ze=matmul(conjg(transpose(V)),tmp)
            ! HF expectation uses z*H_i directly, independently of the interband sum.
            dd=0.d0; aa=0.d0
            tmp=matmul(hx,V)
            do a=1,no
                dd(1)=dd(1)+real(dot_product(V(:,a),z*tmp(:,a)),8)
                if(a>=N/2) aa(1)=aa(1)+real(dot_product(V(:,a),z*tmp(:,a)),8)
            enddo
            xe=matmul(conjg(transpose(V)),tmp)
            tmp=matmul(hy,V)
            do a=1,no
                dd(2)=dd(2)+real(dot_product(V(:,a),z*tmp(:,a)),8)
                if(a>=N/2) aa(2)=aa(2)+real(dot_product(V(:,a),z*tmp(:,a)),8)
            enddo
            ye=matmul(conjg(transpose(V)),tmp)
            hf=hf+[dd(2),-dd(1)]
            active=active+[aa(2),-aa(1)]
            ff=0.d0
            gap=100.d0*epsilon(1.d0)*maxval(abs(W))
            do c=1,nc
                lo=N/2-cuts(c); hi=N/2+1+cuts(c)
                do a=lo,no
                    do b=no+1,hi
                        if(abs(W(a)-W(b))<=gap) then
                            write(*,*) 'degeneracy: iu iv no a b E_a E_b mu gap ',iu,iv,no,a,b,W(a),W(b),mu,gap
                            error stop 'occupied-empty degeneracy'
                        endif
                        factor=(W(a)+W(b)-2.d0*mu)/(W(a)-W(b))
                        ff(:,c)=ff(:,c)+factor*[real(xe(a,b)*ze(b,a),8),real(ye(a,b)*ze(b,a),8)]
                    enddo
                enddo
            enddo
            sumf=sumf+ff
        enddo
    enddo
    !$omp end do
    deallocate(W,H0,V,xe,ye,ze,tmp)
    !$omp end parallel
    hf=pref*hf/real(nk,8)**2
    active=pref*active/real(nk,8)**2
    do c=1,nc
        sumf(:,c)=pref*[-sumf(2,c),sumf(1,c)]/real(nk,8)**2
        write(*,'(A,I5,2ES25.16)') 'formula ',cuts(c),sumf(:,c)
    enddo
    write(*,'(A,2ES25.16)') 'HF_raw ',hf
    write(*,'(A,2ES25.16)') 'HF_active_only ',active
    flush(6)
    if(B0>0.d0) then
        call CALCULATE_OM_MAGNETIC(N,nk,bs(:,1),bs(:,2),mu,z,param,B0,fd(1),fd(2),o)
        write(*,'(A,4ES25.16)') 'Omega_J_xp_xm_yp_ym ',o
        write(*,'(A,2ES25.16)') 'Omega_FD ',fd
        write(*,'(A,2ES25.16)') 'Omega_FD_minus_HF ',fd-hf
    endif
    ! Independent boundary integral of G=sum_occ (E-mu)<z>.
    ! k = bs*(u,v): grad_k G = transpose(inv(bs))*grad_(u,v) G.
    nedge=max(40,2*nk); sumg=0.d0
    !$omp parallel do default(none) shared(nedge,bs,Href,hx,hy,z,mu) &
    !$omp private(iu,idir,side,k,g,gdiff) reduction(+:sumg) schedule(static)
    do iu=1,nedge
        gdiff=0.d0
        do idir=1,2
            do side=-1,1,2
                if(idir==1) then
                    k=0.5d0*side*bs(:,1)+(-0.5d0+(iu-0.5d0)/nedge)*bs(:,2)
                else
                    k=(-0.5d0+(iu-0.5d0)/nedge)*bs(:,1)+0.5d0*side*bs(:,2)
                endif
                call energy_z(k,g)
                gdiff(idir)=gdiff(idir)+side*g
            enddo
        enddo
        sumg=sumg+gdiff
    enddo
    !$omp end parallel do
    det=bs(1,1)*bs(2,2)-bs(1,2)*bs(2,1)
    invbs(1,:)=[bs(2,2),-bs(1,2)]/det
    invbs(2,:)=[-bs(2,1),bs(1,1)]/det
    edge=matmul(transpose(invbs),sumg/nedge)
    edge=pref*[edge(2),-edge(1)]
    write(*,'(A,2ES25.16)') 'BZ_boundary ',edge
    write(*,'(A,2ES25.16)') 'HF_minus_boundary ',hf-edge
    write(*,'(A,2ES25.16)') 'closure_error ',hf-edge-sumf(:,nc)
    if(verify/=0) then
        c=min(100,N/2-1)
        call CALCULATE_OM_CUTOFF_KET(N,nk,c,bs(:,1),bs(:,2),mu,pref,z,param, &
            HAM_STRAIN,HAM_DERIVATIVES_STRAIN_PHI,mx,my)
        write(*,'(A,2ES25.16)') 'original_ket ',mx,my
        if(maxval(abs([mx,my]-sumf(:,5)))>1.d-9) error stop 'original ket mismatch'
        ! Small-grid projector cross-check: use a smaller k step than production.
        if(verify==2) then
            h=1.d-5*norm2(bs(:,1))
            call CALCULATE_OM_CUTOFF(N,nk,c,bs(:,1),bs(:,2),h,mu,pref,z,param,HAM_STRAIN,mx,my)
            write(*,'(A,2ES25.16)') 'original_projector ',mx,my
        endif
        ! Local independent derivative check of G at a generic momentum.
        k=0.13d0*bs(:,1)-0.27d0*bs(:,2); h=1.d-5*norm2(bs(:,1))
        do idir=1,2
            k(idir)=k(idir)+h; call energy_z(k,gp)
            k(idir)=k(idir)-2*h; call energy_z(k,gm)
            k(idir)=k(idir)+h
            numerical(idir)=(gp-gm)/(2*h)
        enddo
        write(*,'(A,2ES25.16)') 'local_gradG_moment ',pref*[numerical(2),-numerical(1)]
    endif
    write(*,'(A,F12.3)') 'wall_seconds ',omp_get_wtime()-start
contains
    subroutine energy_z(kval,gval)
        real(8), intent(in) :: kval(2)
        real(8), intent(out) :: gval
        real(8), allocatable :: ew(:)
        complex(8), allocatable :: hh(:,:),vv(:,:)
        integer :: band
        allocate(ew(N),hh(N,N),vv(N,N))
        hh=Href+kval(1)*hx+kval(2)*hy
        call myzheev(N,ew,hh,vv)
        gval=0.d0
        do band=1,count(ew<mu)
            gval=gval+(ew(band)-mu)*sum(z*abs(vv(:,band))**2)
        enddo
    end subroutine energy_z
end program COMPARE_OM_MAGNETIC
