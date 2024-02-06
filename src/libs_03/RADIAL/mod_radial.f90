!====================================================================
      MODULE RADIAL
!====================================================================
!     contains the radial grid, radial wavefinctions and other
!     common variables, used in the 'radial' routines
!--------------------------------------------------------------------

      IMPLICIT NONE

! ... radial grid parameters:

      Real(8) :: RHO = -4.D0   ! initial value of logarifmic grid
      Real(8) :: H = 1./16.D0  ! step of logarifmic grid
      Real(8) :: H1            ! = 2*H/3
      Real(8) :: EH            ! = DEXP(-H)

      Integer :: NR = 220   ! max. number of radial points

      Integer :: mrf =  0   ! max. number of radial functions
      Integer :: nrf =  0   ! current number of radial functions
      Integer :: irf = 90   ! initial number of radial functions
      Integer :: jrf = 90   ! increment for radial functions
      Integer :: iscr= 90   ! file unit for re-allocation

! ... radial arrays:

      Real(8), Allocatable :: R(:), RR(:), R2(:)

! ... radial wave functions:

      Real(8), Allocatable ::  P(:,:), OVERLAPS(:,:)

      Integer, Allocatable :: nro(:),lro(:),kro(:),mro(:)
      Character(4), Allocatable :: ero(:)

      Real(8), Allocatable :: AZ(:), aexp(:), bexp(:)
      Integer, Allocatable :: mexp(:)

! ... potentials arrays:

      Real(8), Allocatable :: ZK(:), YK(:)
      Integer :: mzk, myk
      Real(8) :: azk,bzk, ayk,byk

! ... other parameters:

      Real(8)    :: Z    = 0.d0    ! charge of nuclear

      Integer :: kclosd = 0     ! number of common closed shells

      Logical :: REL  = .FALSE. ! pointer on the rel. corrections

      Logical :: OO   = .FALSE. ! pointer on orbit-orbit interaction

      Integer :: MASS   = 0     ! pointer to the mass correction

      Real(8)    :: RMASS  = 0.d0  ! mass-correction parameter

! ... the fine-structure constant:

      Real(8) :: FINE = 0.25d0/(137.03599976d0)**2

! ... commonly used constants:

      Real(8) :: D0=0.D0, D1=1.D0, D2=2.D0, D3=3.D0, D4=4.D0,  &
                 D5=0.5D0, D6=6.D0, D8=8.D0, D10=10.D0,        &
                 D12=12.D0, D16=16.D0, D30=30.D0

      Character(6) :: atom = 'atom', term = 'term'

      END MODULE RADIAL


!======================================================================
      Subroutine Alloc_radial(m)
!======================================================================

      USE RADIAL

      IMPLICIT NONE
      Integer, Intent(in) :: m
      Integer :: i

      if(m.le.0) then

       if(allocated(P)) Deallocate(nro,lro,kro,ero,mro, &
                                   AZ,mexp,aexp,bexp,P,OVERLAPS)
       mrf = 0; nrf = 0
       if(m.lt.0) then
       mrf = irf
       Allocate(nro(mrf),lro(mrf),kro(mrf),ero(mrf),mro(mrf),  &
                AZ(mrf),aexp(mrf),bexp(mrf),mexp(mrf),P(NR,mrf), &
                OVERLAPS(mrf,mrf))
       end if
      elseif(m.lt.mrf) then

       Return

      elseif(m.gt.mrf.and.nrf.eq.0) then

       if(allocated(P)) Deallocate(nro,lro,kro,ero,mro, &
                                   AZ,mexp,aexp,bexp,P,OVERLAPS)
       mrf = m
       Allocate(nro(mrf),lro(mrf),kro(mrf),ero(mrf),mro(mrf),  &
                AZ(mrf),aexp(mrf),bexp(mrf),mexp(mrf),P(NR,mrf), &
                OVERLAPS(mrf,mrf))

      elseif(m.gt.mrf.and.nrf.gt.0) then

       Open(iscr,status='SCRATCH',form='UNFORMATTED')
       rewind(iscr)
       Do i = 1,nrf
        write(iscr) nro(i),lro(i),kro(i),ero(i),mro(i),AZ(i), &
                    aexp(i),bexp(i),mexp(i)
        write(iscr) P(1:NR,i)
        write(iscr) Overlaps(1:nrf,i)
       End do
       Deallocate(nro,lro,kro,ero,mro, AZ,mexp,aexp,bexp,P,OVERLAPS)

       mrf = m
       Allocate(nro(mrf),lro(mrf),kro(mrf),ero(mrf),mro(mrf),  &
                AZ(mrf),aexp(mrf),bexp(mrf),mexp(mrf),P(NR,mrf), &
                OVERLAPS(mrf,mrf))
       OVERLAPS=0.d0
       rewind(iscr)
       Do i = 1,nrf
        read(iscr) nro(i),lro(i),kro(i),ero(i),mro(i),AZ(i), &
                   aexp(i),bexp(i),mexp(i)
        read(iscr) P(1:NR,i)
        read(iscr) Overlaps(1:nrf,i)
       End do

      end if

      END Subroutine ALLOC_RADIAL


!======================================================================
      Subroutine Get_overlaps(eps_c,pri)
!======================================================================
!     define overlas for given set of orbitals
!----------------------------------------------------------------------

      USE RADIAL

      IMPLICIT NONE

      Real(8), intent(in) :: eps_c
      Integer, intent(in) :: pri
      Real(8) :: C
      Integer :: i,j
      Real(8), External :: QUADR

      Do i=1,nrf
       Do j=i,nrf
        C = 0.d0
        if(lro(i).eq.lro(j))  C = QUADR(i,j,0)
        if(abs(C).lt.Eps_C)   C = 0.d0
        OVERLAPS(i,j)=C; OVERLAPS(j,i)=C
       End do
      End do

      if(pri.eq.0) Return
      write(pri,'(/a/)') ' Non-trivial one-elctron overlaps:'
      Do i=1,nrf
       Do j=i,nrf
        C = OVERLAPS(i,j)
        if(i.ne.j.and.abs(C).gt.2.D-5.and.nro(i).ne.nro(j)) &
         write(pri,'(''<'',a4,''|'',a4,''>='',f10.5)') &
         ero(i),ero(j),C
         C = C - 1.d0
        if(lro(i).eq.lro(j).and.nro(i).eq.nro(j).and. &
           abs(C).gt.2.D-5) &
         write(pri,'(''<'',a4,''|'',a4,''>='',f10.5)') &
         ero(i),ero(j),C
       End do
      End do

      End Subroutine Get_overlaps
