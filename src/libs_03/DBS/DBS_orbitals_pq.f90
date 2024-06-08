!=======================================================================
Module DBS_orbitals_pq
  !=======================================================================
  !     contains description of atomic orbitals in B-spline representation
  !-----------------------------------------------------------------------
  Implicit none

  Integer :: mbf = 0                  ! max. number of orbitals
  Integer :: nbf = 0                  ! current number of orbitals
  Integer :: ibf = 2**10              ! initial prediction for mbf
  Integer, allocatable :: nbs(:)      ! n-values
  Integer, allocatable :: kbs(:)      ! kappa-values
  Integer, allocatable :: lbs(:)      ! l-value
  Integer, allocatable :: jbs(:)      ! j-values
  Integer, allocatable :: ibs(:)      ! set indexes
  Integer, allocatable :: mbs(:)      ! number of splines
  Integer, allocatable :: ipbs(:)     ! additional pointer

  Character(5), allocatable :: ebs(:) ! spectroscopic notation

  Real(8), allocatable :: pq(:,:,:)   ! radial functions

  Real(8), allocatable :: bpq(:,:,:)  ! B x pq
  Real(8), allocatable :: dbs(:,:,:)  ! L- (or d-) integrals
  Real(8), allocatable :: vbs(:,:,:)  ! dipV-integrals

  Real(8), allocatable :: obs1(:)     ! auxiliary array
  Real(8), allocatable :: obs2(:)     ! auxiliary array

  Real(8) :: memory_DBS_orbitals = 0.d0
  Integer :: ns_bf = 0

  ! ... overlaps between bound and continuum orbitals:
  ! ... (with <.|i> vectors if any):

  Integer :: nv_ch = 0
  Real(8), allocatable :: V_ch(:,:)
  Integer, allocatable :: i_ch(:),j_ch(:)
  Integer, parameter   :: ib = 2**15
  Real(8) :: eps_v = 1.d-10

contains

  !======================================================================
  Subroutine read_pqbs(nu, e_orb)
    !======================================================================
    !     read B-spline w.f. from bsw-file (unit nu) only those orbitals
    !     which labels are already in the list
    !----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss

    Implicit none

    Integer, intent(in) :: nu
    real(8), optional, intent(out) :: e_orb(:)
    Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw,kp,kq
    Character(5) :: elw
    Integer, external :: Ifind_bsorb
    Real(8) :: tt(ns+ks), S

    rewind(nu)
    read(nu) itype,nsw,ksw,tt,kp,kq
    if(itype.ne.grid_type) then
       write(*,'(" Read_pqbs:  another grid_type: itype = ",i15," <> grid_type = ",i15)') itype, grid_type
       Stop 1
    end if
    if(ksw.ne.ks) then
       write(*,'(" Read_pqbs:  ksw =",i15," <> ks = ",i15)') ksw, ks
       Stop 1
    end if
    if(nsw.ne.ns) then
       write(*,'(" Read_pqbs:  nsw =",i15," <> ns = ",i15)') nsw, ns
       Stop 1
    end if
    if(ksp.ne.kp) then
       write(*,'(" Read_pqbs:  ksp =",i15," <> kp = ",i15)') ksp, kp
       Stop 1
    end if
    if(ksq.ne.kq) then
       write(*,'(" Read_pqbs:  ksq =",i15," <> kq = ",i15)') ksq, kq
       Stop 1
    end if
    Do i=1,ns+ks
       if(abs(t(i)-tt(i)) > 1.d-12) then
          write(*,'(" Read_pqbs:  t =",1e20.13," <> tt =",1e20.13)') t(i), tt(i)
          Stop 1
       end if
    End do

1   read(nu,end=2) elw,mw,S
    Call EL_NLJK(elw,n,k,l,j,i)
    m = Ifind_bsorb(n,k,i,2)
    if(present(e_orb)) e_orb(m) = S
    mbs(m)=mw
    pq(1:ns,1,m)=0.d0; read(nu) pq(1:mw,1,m)
    pq(1:ns,2,m)=0.d0; read(nu) pq(1:mw,2,m)
    bpq(:,1,m) = MATMUL(fpbs,pq(:,1,m))
    bpq(:,2,m) = MATMUL(fqbs,pq(:,2,m))
    go to 1
2   Close(nu)

  End subroutine read_pqbs


  !======================================================================
  Subroutine Read_dbsw(nu,mode,ishift)
    !======================================================================
    !     read B-spline r.w.f. from bsw-file (unit nu)
    !     mode = 0  - clear previuous w.f. in module DBS_orbitals_pq
    !            1  - only read radial functions for existing orbitals
    !            2  - read and add if needed
    !     ishift -  shift the set indexes for new orbitals
    !               (-1 -> nulify them)
    !----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss

    Implicit none
    Integer, intent(in) :: nu, mode, ishift
    Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw
    Integer, External :: Ifind_bsorb
    Character(5) :: elw
    Real(8) :: x

    if(mode.eq.0)  Call alloc_DBS_orbitals_pq(0,ns)

    rewind(nu)
    read(nu) itype,nsw,ksw
    if(ksw.ne.ks) Stop ' Read_dbsw:  ksw <> ks'
    if(nsw.ne.ns) Stop ' Read_dbsw:  nsw <> ns'
    if(grid_type.gt.0.and.itype.ne.grid_type) &
         Stop 'Stop in read_dbsw: another grid_type ?'

1   read(nu,end=2) elw,mw
    Call EL_NLJK(elw,n,k,l,j,i); i=i+ishift; if(ishift.eq.-1) i=0
    m = Ifind_bsorb(n,k,i,0)
    if(mode.eq.1.and.m.eq.0) then
       read(nu) x
       read(nu) x
       go to 1
    end if
    m = Ifind_bsorb(n,k,i,2)
    pq(1:ns,1:2,m)=0.d0
    mbs(m)=mw
    read(nu) pq(1:mw,1,m)
    read(nu) pq(1:mw,2,m)

    bpq(1:ns,1,m) = MATMUL(fpbs, pq(:,1,m))
    bpq(1:ns,2,m) = MATMUL(fqbs, pq(:,2,m))

    go to 1
2   Close(nu)
  End subroutine Read_dbsw
End Module DBS_orbitals_pq


!=======================================================================
      Subroutine alloc_DBS_orbitals_pq(m,ns)
!=======================================================================
!     allocate, deallocate or reallocates arrays in module DBS_orbitals_pq
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: m,ns
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: parr(:,:,:), oarr(:)
      Character(5), allocatable :: carr(:)

      if(m.le.0) then
       if(allocated(nbs)) &
        Deallocate(nbs,kbs,lbs,jbs,ibs,mbs,ipbs,ebs,pq,bpq,obs1,obs2)
       nbf=0; mbf=0
      elseif(.not.allocated(nbs)) then
       mbf=m
       Allocate(nbs(m),kbs(m),lbs(m),jbs(m),ibs(m),mbs(m),ebs(m),ipbs(m),&
                pq(ns,2,m),bpq(ns,2,m),obs1(m),obs2(m))
       nbs=0; kbs=0; lbs=0; jbs=0; ibs=0; mbs=0; ipbs=0; obs1=0.d0; obs2=0.d0
       pq=0.d0; bpq=0.d0; nbf=0
      elseif(m.le.mbf) then
       Return
      elseif(nbf.eq.0) then
       Deallocate(nbs,kbs,lbs,jbs,ibs,mbs,ipbs,ebs,pq,bpq,obs1,obs2)
       mbf=m
       Allocate(nbs(m),kbs(m),lbs(m),jbs(m),ibs(m),mbs(m),ebs(m),ipbs(m),&
                pq(ns,2,m),bpq(ns,2,m),obs1(m),obs2(m))
       nbs=0; kbs=0; lbs=0; jbs=0; ibs=0; mbs=0; ipbs=0; obs1=0.d0; obs2=0.d0
       pq=0.d0; bpq=0.d0
      else
       mbf=m
       Allocate(iarr(nbf))
       iarr(1:nbf)=nbs(1:nbf); Deallocate(nbs); Allocate(nbs(m)); nbs=0
       nbs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=kbs(1:nbf); Deallocate(kbs); Allocate(kbs(m)); kbs=0
       kbs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=lbs(1:nbf); Deallocate(lbs); Allocate(lbs(m)); lbs=0
       lbs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=jbs(1:nbf); Deallocate(jbs); Allocate(jbs(m)); jbs=0
       jbs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=ibs(1:nbf); Deallocate(ibs); Allocate(ibs(m)); ibs=0
       ibs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=mbs(1:nbf); Deallocate(mbs); Allocate(mbs(m)); mbs=0
       mbs(1:nbf)=iarr(1:nbf)
       iarr(1:nbf)=ipbs(1:nbf); Deallocate(ipbs); Allocate(ipbs(m)); ipbs=0
       ipbs(1:nbf)=iarr(1:nbf)
       Deallocate(iarr)

       Allocate(parr(ns,2,nbf))
       parr(:,:,1:nbf)=pq(:,:,1:nbf); Deallocate(pq); Allocate(pq(ns,2,mbf))
       pq=0.d0; pq(:,:,1:nbf)=parr(:,:,1:nbf)
       parr(:,:,1:nbf)=bpq(:,:,1:nbf); Deallocate(bpq); Allocate(bpq(ns,2,mbf))
       bpq=0.d0; bpq(:,:,1:nbf)=parr(:,:,1:nbf)
       Deallocate(parr)

       Allocate(carr(nbf))
       carr(1:nbf)=ebs(1:nbf); Deallocate(ebs); Allocate(ebs(m))
       ebs(1:nbf)=carr(1:nbf)
       Deallocate(carr)

       Allocate(oarr(nbf))
       oarr(1:nbf)=obs1(1:nbf); Deallocate(obs1); Allocate(obs1(m))
       obs1(1:nbf)=oarr(1:nbf)
       oarr(1:nbf)=obs2(1:nbf); Deallocate(obs2); Allocate(obs2(m))
       obs2(1:nbf)=oarr(1:nbf)
       Deallocate(oarr)
      end if

      memory_DBS_orbitals = (4.0*12*mbf +  4.0*8*mbf*ns)/(1024d0*1024d0)
      ns_bf = ns

      End Subroutine alloc_DBS_orbitals_pq


!=======================================================================
      Integer Function Ifind_bsorb(n,k,iset,job)
!=======================================================================
!     find orbital (n,k,iset) in the list DBS_orbitals_pq:
!     job = 0  -  no further actions
!     job = 1  -  stop if fail to find
!     job = 2  -  add new orbital
!----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Use DBS_grid, only: ns

      Implicit none
      Integer, intent(in) :: n,k,iset,job
      Integer :: i
      Character(5), external :: ELi

      Ifind_bsorb=0

      Do i=1,nbf
       if(n.ne.nbs(i)) Cycle
       if(k.ne.kbs(i)) Cycle
       if(iset.ne.ibs(i)) Cycle
       Ifind_bsorb = i
       Return
      End do
      if(job.eq.0) Return

      if(job.eq.1) then
       Write(*,'(a,a,5i5)') 'Ifind_bsorb can not find the orbital:',&
                            ' n,k,iset = ',n,k,iset
       Stop ' '
      end if

      if(nbf.ge.mbf) Call Alloc_DBS_orbitals_pq(mbf+ibf,ns)
      nbf = nbf + 1
      nbs(nbf) = n
      kbs(nbf) = k
      lbs(nbf) = k; if(k.lt.0) lbs(nbf)=-k-1
      jbs(nbf) = k+k-1; if(k.lt.0) jbs(nbf)=-k-k-1
      ibs(nbf) = iset
      ebs(nbf) = ELi(n,k,iset)
      Ifind_bsorb = nbf

      End Function Ifind_bsorb

!=======================================================================
      Subroutine Get_pv(i,v,ns)
!=======================================================================
!     get two-component Function "i"  as one vector "v"
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Implicit none
      Integer, intent(in) :: i,ns
      Real(8), intent(out) :: v(ns+ns)
      v(1:ns)=pq(1:ns,1,i)
      v(ns+1:ns+ns)=pq(1:ns,2,i)
      End Subroutine Get_pv


!=======================================================================
      Subroutine Get_qv(i,v,ns)
!=======================================================================
!     get two-component Function "i"  as one vector "v"
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Implicit none
      Integer, intent(in) :: i,ns
      Real(8), intent(out) :: v(ns+ns)
      v(1:ns)=bpq(1:ns,1,i)
      v(ns+1:ns+ns)=bpq(1:ns,2,i)
      End Subroutine Get_qv


!=======================================================================
      Subroutine Put_pv(i,v,ns)
!=======================================================================
!     put two-component vector "v" in p-array, location "i"
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Implicit none
      Integer, intent(in) :: i,ns
      Real(8), intent(in) :: v(ns+ns)
      pq(1:ns,1,i)=v(1:ns)
      pq(1:ns,2,i)=v(ns+1:ns+ns)
      End Subroutine put_pv


!=======================================================================
      Subroutine Alloc_nv_ch(m)
!=======================================================================
!     alloc the continuum-bound overlaps
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: m
      Integer :: k, i,j

      if(allocated(i_ch)) Deallocate(i_ch)
      if(allocated(j_ch)) Deallocate(j_ch)
      nv_ch = 0
      if(m.le.0) Return

      k = 0
      Do i=1,nbf;  Do j=1,nbf
       if(kbs(i).ne.kbs(j)) Cycle
       if(ipbs(i).eq.0) Cycle
       if(ipbs(i).ne.0.and.ipbs(j).ne.0) Cycle
       if(ibs(j).eq.0) Cycle
       k = k + 1
      End do; End do
      nv_ch = k

      if(nv_ch.eq.0) Return

      Allocate(i_ch(nv_ch),j_ch(nv_ch))
      i_ch = 0; j_ch = 0
      k = 0
      Do i=1,nbf;  Do j=1,nbf
       if(kbs(i).ne.kbs(j)) Cycle
       if(ipbs(i).eq.0) Cycle
       if(ipbs(i).ne.0.and.ipbs(j).ne.0) Cycle
       if(ibs(j).eq.0) Cycle
       k = k + 1
       i_ch(k) = i; j_ch(k) = j
      End do; End do

      End Subroutine Alloc_nv_ch


!======================================================================
      Integer Function KBORT(i,j)
!======================================================================
!     find overlap between channel "i" and bound orbital "j"
!----------------------------------------------------------------------
      Use DBS_orbitals_pq,  only: nv_ch,i_ch,j_ch

      Implicit none
      Integer, intent(in) :: i,j
      Integer :: k,l,m

      KBORT = 0; if(nv_ch.eq.0) Return

! ... search position (m) for given overlap

      k=1; l=nv_ch
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (i.lt.i_ch(m)) then;       l = m - 1
      elseif(i.gt.i_ch(m)) then;       k = m + 1
      else
       if    (j.lt.j_ch(m)) then;      l = m - 1
       elseif(j.gt.j_ch(m)) then;      k = m + 1
       else
        KBORT=m;  Return
       end if
      end if
      go to 1
    2 KBORT = 0

      End Function KBORT

!======================================================================
      Integer Function IBORT(io,jo)
!======================================================================
!     recover old definition IBORT
!----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: ipbs,ib

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j
      Integer, external :: KBORT

      if(ipbs(io).gt.0.and.ipbs(jo).eq.0) then
        i = io; j = jo
      elseif(ipbs(io).eq.0.and.ipbs(jo).gt.0) then
        j = io; i = jo
      elseif(ipbs(io).gt.0.and.ipbs(jo).gt.0) then
        i = max(io,jo); j = min(io,jo)
        IBORT = i*ib+j; Return
      else
        Stop 'JBORT for bound orbitals?'
      end if

      IBORT = KBORT(i,j)
      if(IBORT.gt.0) IBORT = i*ib+j

      End Function IBORT


!======================================================================
      Subroutine Nulify_bort(i,j)
!======================================================================
!     find orthogonal condition for orbitals io and jo
!----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: nv_ch,i_ch,j_ch

      Implicit none
      Integer, intent(in) :: i,j
      Integer :: m,k
      Integer, external :: KBORT

      m = KBORT(i,j)
      if(m.eq.0) Return
      nv_ch = nv_ch-1
      Do k=m,nv_ch
       i_ch(k) = i_ch(k+1)
       j_ch(k) = j_ch(k+1)
      End do

      End Subroutine Nulify_bort


!======================================================================
      Subroutine read_bort_jj(nuc)
!======================================================================
!     read from file nuc the imposed orthogonal conditions given as
!     < n1k1| n2k2>=0  (or 1 | 2)
!----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: nbf,kbs,nbs

      Character(15) :: Aort

      rewind(nuc)
    1 read(nuc,'(a)',end=2) Aort
      if(Aort(1:1).ne.'<') go to 1
      Call EL_nljk(Aort(2: 6),n1,kap1,l1,j1,k1)
      i1 = Ifind_jjorb(n1,kap1,k1,0)
      if(i1.eq.0.and.k1.ne.0) go to 1
      Call EL_nljk(Aort(8:12),n2,kap2,l2,j2,k2)
      i2 = Ifind_jjorb(n2,kap2,k2,0)
      if(i2.eq.0.and.k2.ne.0) go to 1
      if(kap1.ne.kap2) go to 1
      read(Aort(15:15),'(i1)') ii

      if(ii.ne.0) go to 1

      if(k1.gt.0.and.k2.gt.0) then
       Call Nulify_bort(i1,i2)
      elseif(k1.gt.0.and.k2.eq.0) then
       Do i2=1,nbf
        if(kap2.ne.kbs(i2).or.n2.ne.nbs(i2)) Cycle
        Call Nulify_bort(i1,i2)
       End do
      elseif(k1.eq.0.and.k2.gt.0) then
       Do i1=1,nbf
        if(kap1.ne.kbs(i1).or.n1.ne.nbs(i1)) Cycle
        Call Nulify_bort(i1,i2)
       End do
      elseif(k1.eq.0.and.k2.eq.0) then
       Do i1=1,nbf
        if(kap1.ne.kbs(i1).or.n1.ne.nbs(i1)) Cycle
       Do i2=1,nbf
        if(kap2.ne.kbs(i2).or.n2.ne.nbs(i2)) Cycle
        Call Nulify_bort(i1,i2)
       End do; End do
      end if

      go to 1
    2 Continue

      End Subroutine read_bort_jj


!=======================================================================
      Subroutine Get_v_ch(m)
!=======================================================================
!     define v_ch arrays if any
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq, ns => ns_bf

      Integer, intent(in) :: m

      if(allocated(v_ch)) Deallocate(v_ch)
      if(m.le.0) Return

      if(nv_ch.gt.0) Allocate(V_ch(ns+ns,nv_ch))

      Do k = 1,nv_ch
       i = i_ch(k); j = j_ch(k)
       if(IBORT(i,j).eq.0) Cycle
       Call Get_v (i,j,V_ch(1,k))
      End do

      memory_DBS_orbitals = (48.0*mbf +  4.0*8*mbf*ns)/(1024d0*1024d0)  &
                         +  (8.0*nv_ch + 16.0*ns*nv_ch)/(1024d0*1024d0)

      End Subroutine Get_v_ch



!======================================================================
      Subroutine get_V (kl,nl,v)
!======================================================================
!     provides <.|nl>  for <kl|nl>
!     if additionally <kl|n'l'> = 0  then  nl -> nl - <nl|n'l'> n'l'
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: kl,nl
      Real(8) :: v(ms), w(ms)
      Integer :: i
      Real(8), external :: obs
      Integer, external :: IBORT

      if(IBORT(kl,nl).eq.0) then
       write(*,'(2a6,10i5)') &
        ebs(kl),ebs(nl),IBORT(kl,nl),IBORT(nl,kl),kl,nl,nbf
       write(*,'(10a6)') ebs(1:nbf)
       Stop 'get_V: <kl|nl> = 0 ?'
      end if

      Call Get_qv(nl,v,ns)

      Do i=1,nbf
       if(i.eq.nl) Cycle
       if(ipbs(i).ne.0) Cycle
       if(kbs(i).ne.kbs(nl)) Cycle
       if(abs(OBS(i,nl)).lt.eps_v) Cycle
       if(IBORT(i,kl).ne.0) Cycle
       Call Get_qv(i,w,ns)
       v(:) = v(:) - OBS(i,nl)*w(:)
      End do

      End Subroutine get_V
