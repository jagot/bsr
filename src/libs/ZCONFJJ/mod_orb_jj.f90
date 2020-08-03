!=====================================================================
      Module orb_jj
!=====================================================================
!     defines the one-electron orbitals in jj-coupling
!---------------------------------------------------------------------
      Implicit none

! ... list of one-electron orbitals

      Integer :: mwf = 0        ! max. number of orbitals
      Integer :: nwf = 0        ! current number of orbitals
      Integer :: iwf = 2**10    ! initial prediction for mwf

      Integer, allocatable :: nef(:)   ! n-values
      Integer, allocatable :: kef(:)   ! kappa-value
      Integer, allocatable :: lef(:)   ! l-value
      Integer, allocatable :: jef(:)   ! j-value
      Integer, allocatable :: ief(:)   ! subset index
      Integer, allocatable :: ipef(:)  ! additional pointer

      Character(5), allocatable :: ELF(:) ! spectroscopic notation

! ... orbital orthogonality and AFTER conditions

      Integer :: JORT = 1

! ... possible set indexes for orthogonal subsets:

      Integer, parameter :: kset = 61
      Character(kset) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
      Integer, parameter :: ksmax = kset*kset

! ... shift for the set indexes:

      Integer :: kshift = 0

! ... current allocated memory (in Mb)

      Real(8) :: memory_orb_jj   

      End Module orb_jj


!======================================================================
      Subroutine alloc_orb_jj(m)
!======================================================================
!     allocate, deallocate or reallocate orbital list in module orb_jj
!----------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer, intent(in) :: m
      Integer :: i,j
      Integer, allocatable :: iarr(:)
      Character(5), allocatable :: carr(:)

      if(m.le.0) then
       if(allocated(nef)) Deallocate (nef,kef,lef,jef,ief,ipef,ELF)
       mwf = 0; nwf = 0
      elseif(.not.allocated(nef)) then
       mwf = m
       Allocate(nef(mwf),kef(mwf),lef(mwf),jef(mwf),ief(mwf),ipef(mwf), &
                ELF(mwf))
       nef=0; kef=0; lef=0; jef=0; ief=0; ipef=0
      elseif(m.le.mwf) then
       Return
      elseif(nwf.eq.0) then
       Deallocate (nef,kef,lef,jef,ief,ipef,ELF)
       mwf = m
       Allocate(nef(mwf),kef(mwf),lef(mwf),jef(mwf),ief(mwf),ipef(mwf), &
                ELF(mwf))
       nef=0; kef=0; lef=0; jef=0; ief=0; ipef=0
      else
       mwf=m
       Allocate(iarr(nwf))
       iarr(1:nwf)=nef(1:nwf); Deallocate(nef); Allocate(nef(m)); nef=0
       nef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=kef(1:nwf); Deallocate(kef); Allocate(kef(m)); kef=0
       kef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=lef(1:nwf); Deallocate(lef); Allocate(lef(m)); lef=0
       lef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=jef(1:nwf); Deallocate(jef); Allocate(jef(m)); jef=0
       jef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=ief(1:nwf); Deallocate(ief); Allocate(ief(m)); ief=0
       ief(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=ipef(1:nwf); Deallocate(ipef); Allocate(ipef(m)); ipef=0
       ipef(1:nwf)=iarr(1:nwf)
       Deallocate(iarr)

       Allocate(carr(nwf))
       carr(1:nwf)=ELF(1:nwf); Deallocate(ELF); Allocate(ELF(m))
       ELF(1:nwf)=carr(1:nwf)
       Deallocate(carr)

      end if

      memory_orb_jj = (29*mwf + 4*mwf*mwf)/(1024d0*1024d0)

      End Subroutine alloc_orb_jj


!=======================================================================
      Integer Function Ifind_jjorb(n,k,ii,job)
!=======================================================================
!     find orbital (n,k,iset) in the list "orb_jj"
!     Options:
!     job = 0  -  no further actions
!     job = 1  -  stop if fail to find
!     job = 2  -  add new orbital
!------------------------------------------------------------------------
      USE orb_jj

      Implicit none
      Integer :: n,k,i,ii,j,job
      Character(5), external :: ELi

      Ifind_jjorb=0;  i = ii !   + kshift    !   ???

      Do j=1,nwf
       if(n.ne.nef(j)) Cycle
       if(k.ne.kef(j)) Cycle
       if(i.ge.0.and.i.ne.ief(j)) Cycle
       Ifind_jjorb=j
       Return
      End do
      if(job.eq.0) Return

      if(job.eq.1) then
       Write(*,'(a,a,5i5)') 'Ifind_jjorb can not find the orbital:',&
                            ' N,K,iset = ',n,k,i
       Stop
      end if

      if(nwf+1.gt.mwf) Call Alloc_orb_jj(mwf+iwf)
      nwf = nwf + 1
      nef(nwf)=n
      kef(nwf)=k
      lef(nwf) = k; if(k.lt.0) lef(nwf)=-k-1
      jef(nwf) = k+k-1; if(k.lt.0) jef(nwf)=-k-k-1
      ief(nwf)=i
      ELF(nwf)=ELi(n,k,i)
      Ifind_jjorb = nwf

      End Function Ifind_jjorb


!======================================================================
      Subroutine Read_bsw_orb_jj(nu)
!======================================================================
!     read only spectroscopic notation from bsw-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j,k,l,n
      Character(5) :: elw
      Integer, external :: Ifind_jjorb
      Real(8) :: p

      rewind(nu)
      read(nu) i ! ignore the rest of the record
    1 read(nu,end=2) elw
      Call EL_NLJK(elw,n,k,l,j,i)
      i = Ifind_jjorb(n,k,i,2)
      read(nu) p
      read(nu) p
      go to 1
    2 Close(nu)

      End Subroutine Read_bsw_orb_jj


!=======================================================================
      Subroutine Get_orb_jj(i,n,k,iset,ip)
!=======================================================================
!     get parameters for orbital "i" (for external calls)
!------------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer, intent(in) :: i
      Integer, intent(out) :: n,k,iset,ip

      n = 0; k = 0;  iset=0; ip =0
      if(i.le.0.or.i.gt.nwf) Return
      n = nef(i)
      k = kef(i)
      iset = ief(i)
      ip = ipef(i)

      End Subroutine Get_orb_jj

