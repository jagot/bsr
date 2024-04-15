!======================================================================
     Module DBS_gauss
!======================================================================
!    defines the values of splines at the gaussian points for each
!    interval of a grid:  (1:nv+1,1:ks,1:ks)
!
!    bsp  (i,m,ith) - values of the i+ith-1 B-spline in interval i
!                     at gausian point m
!    bspd (i,m,ith) - corresponding values of first derivative
!
!    bspdd(i,m,ith) - corresponding values of second derivative
!
!    bsp(nv+1,1,.) and bspd(nv+1,1,.,.) - corresponding values at
!                                         last knot point (rmax)
!----------------------------------------------------------------------
     Implicit none

!... arrays for Gaussian-quadrature data (1:nv;1:ks)

     Real(8), allocatable :: gx(:)    ! gaussian points (1:ks)
     Real(8), allocatable :: gw(:)    ! gaussian weights (1:ks)

     Real(8), allocatable :: gr (:,:) ! gr (i,m) - gaussian points m in the interval i
     Real(8), allocatable :: grm(:,:) ! grm(i,m) - reciprocal value of gr(i,m)
     Real(8), allocatable :: grw(:,:) ! grw(i,m) - corr. gaussian weights
     Real(8), allocatable :: ygw(:,:) ! working array

     Real(8), allocatable :: dbiatx(:,:,:) ! working array (used for calling zbsplvd routine)

!... standard B-splines: (not used, keep for possible applications)

     Real(8), allocatable :: bsp(:,:,:)
     Real(8), allocatable :: bsq(:,:,:)
     Real(8), allocatable :: bspd(:,:,:)
     Real(8), allocatable :: bspdd(:,:,:)

!... basis for P-functions (ks=ksp):

     Real(8), allocatable, target :: pbsp(:,:,:)
     Real(8), allocatable :: pbsd(:,:,:)
     Real(8), allocatable :: dbip(:,:,:)
     Real(8), allocatable :: tp(:)
     Real(8), allocatable :: fpbs(:,:)
     Real(8), allocatable :: fpb_nucl(:,:)

!... basis for Q-functions (ks=ksq):

     Real(8), allocatable, target :: qbsp(:,:,:)
     Real(8), allocatable :: qbsd(:,:,:)
     Real(8), allocatable :: dbiq(:,:,:)
     Real(8), allocatable :: tq(:)
     Real(8), allocatable :: fqbs(:,:)
     Real(8), allocatable :: fqb_nucl(:,:)

!... mixed pq arrays:

     Real(8), allocatable :: fpqbs(:,:)
     Real(8), allocatable :: fqpbs(:,:)
     Real(8), allocatable :: fpqbsd(:,:)
     Real(8), allocatable :: fqpbsd(:,:)
     Real(8), allocatable :: fppqq(:,:)

     Real(8) :: memory_DBS_gauss = 0.d0

     End Module DBS_gauss


!====================================================================
     Subroutine alloc_DBS_gauss
!====================================================================
! ... Allocate space and define the arrays in MODULE DBS_gauss
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Implicit none
     Integer :: i,m,nw

     if(allocated(bsp)) &
      Deallocate(bsp,bsq,bspd,bspdd,gr,grm,grw,dbiatx)
     nw = nv+1
     Allocate(bsp (nw,ks,ks), bsq  (nw,ks,ks), &
              bspd(nw,ks,ks), bspdd(nw,ks,ks), &
              gr(nv,ks), grm(nv,ks), grw(nv,ks), ygw(nv,ks), &
              dbiatx(nv,ks,ks), gx(ks),gw(ks) )

     Call gauleg(0.d0,1.d0,gx,gw,ks)

     dbiatx = 0.d0
     Do m=1,ks
      gr (1:nv,m)= (t(1+ks:nv+ks)-t(ks:nv+ks-1))*gx(m)+t(ks:nv+ks-1)
      grm(1:nv,m)= 1.d0/gr(1:nv,m)
      grw(1:nv,m)= (t(1+ks:nv+ks)-t(ks:nv+ks-1))*gw(m)

      Call zbsplvd(ns,ks,nv,t,ks,nv,gr(1,m),3,dbiatx)

      bsp  (1:nv,m,1:ks) = dbiatx(1:nv,1:ks,1)
      bspd (1:nv,m,1:ks) = dbiatx(1:nv,1:ks,2)
      bspdd(1:nv,m,1:ks) = dbiatx(1:nv,1:ks,3)

      Do i = 1,ks
       bsq(1:nv,m,i) = bspd(1:nv,m,i) - grm(1:nv,m)*bsp(1:nv,m,i)
      End do

     End do

!... store values at the last knot

     Call zbsplvd(ns,ks,nv,t,ns,1,t(ns+1),3,dbiatx)

     bsp  (nw,:,:) = 0.d0
     bspd (nw,:,:) = 0.d0
     bspdd(nw,:,:) = 0.d0
     bsq  (nw,:,:) = 0.d0

     bsp  (nw,1,1:ks) = dbiatx(1,1:ks,1)
     bspd (nw,1,1:ks) = dbiatx(1,1:ks,2)
     bspdd(nw,1,1:ks) = dbiatx(1,1:ks,3)
     bsq  (nw,1,1:ks) = bspd(nw,1,1:ks) - bsp(nw,1,1:ks)/t(ns+1)

!... store also the values at the first knot

     Call zbsplvd(ns,ks,nv,t,ks,1,t(1),3,dbiatx)

     bsp  (nw,2,1:ks) = dbiatx(1,1:ks,1)
     bspd (nw,2,1:ks) = dbiatx(1,1:ks,2)
     bspdd(nw,2,1:ks) = dbiatx(1,1:ks,3)

     if(t(1).ne.0.d0) then
      bsq  (nw,2,1:ks) = bspd(nw,2,1:ks) - bsp(nw,2,1:ks)/t(1)
     else
      bsq  (nw,2,1:ks) = bspd(nw,2,1:ks) - bsp(nw,2,1:ks)/1.d-20
     end if

! ... case of different ks for p- and q-functions:

     Call Def_pbsp(ksp)
     Call Def_qbsp(ksq)
     Call Def_pq_arrays

     memory_DBS_gauss = (nw*ks*ks*4*8 + nv*ks*4*8 + nv*ks*ks*8 +  &
                         nw*ks*ks*3*8 + ns*ns*2*8 +               &
                         nw*ks*ks*3*8 + ns*ns*2*8 +               &
                         ns*ns*4*8 + ms*ms*8 )/(1024d0*1024d0)

     End Subroutine alloc_DBS_gauss


!====================================================================
     Subroutine Def_pbsp(k)
!====================================================================
! ... Allocate space and define the arrays specific for P-functions
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Integer :: k, m, nw

     ksp=k; if(k.eq.0) ksp=ks; nsp=nv+ksp-1
     if(allocated(pbsp)) DeAllocate(pbsp); Allocate(pbsp(nv+1,ks,ksp))
     if(allocated(pbsd)) DeAllocate(pbsd); Allocate(pbsd(nv+1,ks,ksp))
     if(allocated(dbip)) DeAllocate(dbip); Allocate(dbip(nv,ksp,ksp))
     if(allocated(tp  )) DeAllocate(tp  ); Allocate(tp(nsp+ksp))
     if(allocated(fpbs)) DeAllocate(fpbs); Allocate(fpbs(ns,ns))
     if(allocated(fpb_nucl)) DeAllocate(fpb_nucl); Allocate(fpb_nucl(ns,ns))

     tp(1:ksp) = t(1)
     tp(ksp+1:nsp) = t(ks+1:ns)
     tp(nsp+1:nsp+ksp) = t(ns+1:ns+ksp)

     dbip = 0.d0
     Do m=1,ks
      Call zbsplvd(nsp,ksp,nv,tp,ksp,nv,gr(1,m),2,dbip)
      pbsp (1:nv,m,1:ksp) = dbip(1:nv,1:ksp,1)
      pbsd (1:nv,m,1:ksp) = dbip(1:nv,1:ksp,2)
     End do

     nw = nv + 1
     pbsp (nw,:,:) = 0.d0
     pbsd (nw,:,:) = 0.d0

!... store values at the last knot:

     Call zbsplvd(nsp,ksp,nv,tp,nsp,1,tp(nsp+1),2,dbip)

     pbsp (nw,1,1:ksp) = dbip (1,1:ksp,1)
     pbsd (nw,1,1:ksp) = dbip (1,1:ksp,2)


!... store also the values at the first knot

     Call zbsplvd(nsp,ksp,nv,tp,ksp,1,tp(1),2,dbip)

     pbsp (nw,2,1:ksp) = dbip (1,1:ksp,1)
     pbsd (nw,2,1:ksp) = dbip (1,1:ksp,2)

     Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,grw,ns,fpbs)

     End Subroutine Def_pbsp


!====================================================================
     Subroutine Def_qbsp(k)
!====================================================================
! ... Allocate space and define the arrays specific for Q-functions
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Integer :: k, m, nw

     ksq=k; if(k.eq.0) ksq=ks; nsq=nv+ksq-1

     if(allocated(qbsp)) DeAllocate(qbsp); Allocate(qbsp(nv+1,ks,ksq))
     if(allocated(qbsd)) DeAllocate(qbsd); Allocate(qbsd(nv+1,ks,ksq))
     if(allocated(dbiq)) DeAllocate(dbiq); Allocate(dbiq(nv,ksq,ksq))
     if(allocated(tq))   DeAllocate(tq);   Allocate(tq(nsq+ksq))
     if(allocated(fqbs)) DeAllocate(fqbs); Allocate(fqbs(ns,ns))
     if(allocated(fqb_nucl)) DeAllocate(fqb_nucl); Allocate(fqb_nucl(ns,ns))

     tq(1:ksq) = t(1)
     tq(ksq+1:nsq) = t(ks+1:ns)
     tq(nsq+1:nsq+ksq) = t(ns+1:ns+ksq)

     dbiq = 0.d0
     Do m=1,ks
      Call zbsplvd(nsq,ksq,nv,tq,ksq,nv,gr(1,m),2,dbiq)
      qbsp (1:nv,m,1:ksq) = dbiq(1:nv,1:ksq,1)
      qbsd (1:nv,m,1:ksq) = dbiq(1:nv,1:ksq,2)
     End do

     nw = nv + 1
     qbsp (nw,:,:) = 0.d0
     qbsd (nw,:,:) = 0.d0

!... store values at the last knot

     Call zbsplvd(nsq,ksq,nv,tq,nsq,1,tq(nsq+1),2,dbiq)

     qbsp (nw,1,1:ksq) = dbiq (1,1:ksq,1)
     qbsd (nw,1,1:ksq) = dbiq (1,1:ksq,2)

!... store also the values at the first knot

     Call zbsplvd(nsq,ksq,nv,tq,ksq,1,tq(1),2,dbiq)

     qbsp (nw,2,1:ksp) = dbiq (1,1:ksp,1)
     qbsd (nw,2,1:ksp) = dbiq (1,1:ksp,2)

     Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,grw,ns,fqbs)

     End Subroutine Def_qbsp


!====================================================================
     Subroutine Def_pq_arrays
!====================================================================
! ... Allocate space and define the B-splines arrays
! ... connected  P- and Q-functions:
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     ygw = grm * grw
     if(Allocated(fpqbs)) DeAllocate(fpqbs); Allocate(fpqbs(ns,ns))
     Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsp,ygw,ns,fpqbs)
     if(Allocated(fqpbs)) DeAllocate(fqpbs); Allocate(fqpbs(ns,ns))
     Call ZINTYm (nv,ks,ksq,ksp,qbsp,pbsp,ygw,ns,fqpbs)

     if(Allocated(fpqbsd)) DeAllocate(fpqbsd); Allocate(fpqbsd(ns,ns))
     Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsd,grw,ns,fpqbsd)
     if(Allocated(fqpbsd)) DeAllocate(fqpbsd); Allocate(fqpbsd(ns,ns))
     Call ZINTYm (nv,ks,ksq,ksp,qbsp,pbsd,grw,ns,fqpbsd)

     if(Allocated(fppqq)) DeAllocate(fppqq); Allocate(fppqq(ms,ms))
     fppqq = 0.d0
     fppqq(1:ns,1:ns)=fpbs(1:ns,1:ns)
     fppqq(ns+1:ms,ns+1:ms)=fqbs(1:ns,1:ns)

     End Subroutine Def_pq_arrays
