!=====================================================================
      Module orb_LS
!=====================================================================
!     description of atomic one-electron orbitals
!
!     ROUTINES:    Alloc_orb_LS (m)
!                  Ifind_nlk    (n,l,k)
!                  Read_bsw_orb (nu)
!                  Get_nlki     (j,n,l,k,i)
!
!---------------------------------------------------------------------
      Implicit none

! ... LIST OF ONE-ELECTRON ORBITALS

      Integer :: mwf = 0                ! max. number of orbitals
      Integer :: nwf = 0                ! current number of orbitals
      Integer :: iwf = 512              ! initial prediction for mwf
      Integer :: nwf1, nwf2, kset1, kset2

      Integer, allocatable :: NEF(:)    ! n-values
      Integer, allocatable :: LEF(:)    ! l-values
      Integer, allocatable :: KEF(:)    ! set number
      Integer, allocatable :: IEF(:)    ! additional pointer

      Character(4), allocatable :: ELF(:) ! spectroscopic notation

! ... orbital orthogonality and AFTER conditions

      Integer :: JORT = 1
      Integer, allocatable :: IORT(:,:)

      End Module orb_LS


!======================================================================
      Subroutine Alloc_orb_LS(m)
!======================================================================
!     allocate arrays in the module "orb_LS"
!----------------------------------------------------------------------
      Use orb_LS

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:), ib(:,:)
      Character(4), allocatable :: aa(:)

      if(m.le.0) then
       if(allocated(NEF)) Deallocate (NEF,LEF,KEF,IEF,ELF,IORT)
       mwf = 0; nwf = 0
       if(m.lt.0) then
        mwf = iwf
        Allocate(NEF(mwf),LEF(mwf),KEF(mwf),IEF(mwf),ELF(mwf), &
                 IORT(mwf,mwf))
        NEF = 0; LEF = -1; KEF = 0; IEF = 0; ELF = '    '
        IORT = 0
       end if
      elseif(.not.allocated(NEF)) then
       mwf = m
       Allocate(NEF(mwf),LEF(mwf),KEF(mwf),IEF(mwf),ELF(mwf), &
                IORT(mwf,mwf))
       NEF = 0; LEF = -1; KEF = 0; IEF = 0; ELF = '    '
       IORT = 0
      elseif(m.le.mwf) then
       Return
      elseif(nwf.eq.0) then
       Deallocate (NEF,LEF,KEF,IEF,ELF,IORT)
       mwf = m
       Allocate(NEF(mwf),LEF(mwf),KEF(mwf),IEF(mwf),ELF(mwf), &
                IORT(mwf,mwf))
       NEF = 0; LEF = -1; KEF = 0; IEF = 0; ELF = '    '
       IORT = 0
      else
       mwf = m
       Allocate(ia(nwf))
       ia=NEF(1:nwf); Deallocate(NEF)
       Allocate(NEF(mwf)); NEF(1:nwf)=ia
       ia=LEF(1:nwf); Deallocate(LEF)
       Allocate(LEF(mwf)); LEF(1:nwf)=ia
       ia=KEF(1:nwf); Deallocate(KEF)
       Allocate(KEF(mwf)); KEF(1:nwf)=ia
       ia=IEF(1:nwf); Deallocate(IEF)
       Allocate(IEF(mwf)); IEF(1:nwf)=ia
       Deallocate(ia)
       Allocate(aa(nwf))
       aa=ELF(1:nwf); Deallocate(ELF)
       Allocate(ELF(mwf)); ELF(1:nwf)=aa
       Deallocate(aa)
       Allocate(ib(nwf,nwf))
       ib=IORT(1:nwf,1:nwf); Deallocate(IORT)
       Allocate(IORT(mwf,mwf)); IORT(1:nwf,1:nwf)=ib
       Deallocate(ib)
      end if

      End Subroutine Alloc_orb_LS


!=======================================================================
      Integer function Ifind_nlk(n,l,k,job)
!=======================================================================
!     find orbital "nlk" in the orbital list, = 0, if not
!
!     job = 0  -  no further actions
!     job = 1  -  stop if fail to find
!     job = 2  -  add new orbital
!------------------------------------------------------------------------
      Use orb_LS

      Implicit none
      Integer, intent(in) :: n,l,k,job
      Integer :: i
      Character(4), external :: ELF4

      Ifind_nlk=0

      Do i=1,nwf
       if(n.ne.NEF(i)) Cycle
       if(l.ne.LEF(i)) Cycle
       if(k.ne.KEF(i)) Cycle
       Ifind_nlk = i
       Return
      End do
      if(job.eq.0) Return

      if(job.eq.1) then
       write(*,'(a,a,3i5,a6)') 'Ifind_nlk can not find the orbital:',&
                               ' N,L,K = ',n,l,k,ELF4(n,l,k)
       Stop
      end if

      if(nwf.ge.mwf) Call Alloc_orb_LS(mwf+iwf)
      nwf = nwf + 1
      NEF(nwf) = n; LEF(nwf) = l; KEF(nwf) = k; IEF(nwf) = 0
      ELF(nwf) = ELF4(n,l,k)
      Ifind_nlk = nwf

      End Function Ifind_nlk


!=======================================================================
      Subroutine Get_nlki (j,n,l,k,i)
!=======================================================================
!     provide nlk-values for orbital j
!------------------------------------------------------------------------
      Use orb_LS
      Implicit none
      Integer :: i,j,n,l,k
      n = NEF(j); l=LEF(j); k=KEF(j); i=IEF(j)
      End Subroutine Get_nlki


!======================================================================
      Subroutine Read_bsw_orb_LS(nu)
!======================================================================
!     read only spectroscopic notation from bsw-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,n,l,k
      Character(4) :: el
      Integer, External :: Ifind_nlk
      Real(8) :: p

      rewind(nu)
    1 read(nu,end=2) el
      Call EL4_nlk(el,n,l,k)
      i = Ifind_nlk(n,l,k,2)
      read(nu) p
      go to 1
    2 Close(nu)

      End Subroutine Read_bsw_orb_LS
