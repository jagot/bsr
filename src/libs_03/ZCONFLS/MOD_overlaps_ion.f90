!===========================================================================
      Module overlaps_ion
!===========================================================================
!     define overlaps for ion scattering channels 
!---------------------------------------------------------------------------
      Implicit none

      Integer :: novl =  0                !  number of overlaps               

      Integer, allocatable :: iipar(:)    !  index of partial wave
      Integer, allocatable :: nopen(:)    !  open channels 
      Real(8), allocatable :: Eovl (:)    !  energy of state

      Real(8), allocatable :: ovlr(:,:)   !  real overlaps
      Real(8), allocatable :: ovli(:,:)   !  imaginary overlaps
      Real(8), allocatable :: ovlt(:,:)   !  target probabilities
      Real(8), allocatable :: ui  (:,:)   !  sct.phases

      Real(8), parameter   :: eps_e = 1.d-7 

      End Module overlaps_ion


!======================================================================
      Subroutine R_overlaps_ion(nut,mch,ntarg)
!======================================================================
!     read from file 'nut' overlaps information 
!----------------------------------------------------------------------
      Use overlaps_ion
      
      Implicit none
      Integer, intent(in) :: nut,mch,ntarg
      Character(80) :: line
      Integer :: i,j,m,n
 
      novl = 0
      rewind(nut)
    1 read(nut,'(a)',end=2) line
      if(line(1:4).ne.'klsp') go to 1
      novl = novl + 1
      go to 1
    2 if(novl.eq.0) Stop 'R_overlaps: novl = 0'

      if(allocated(iipar)) Deallocate(iipar,nopen,Eovl,ovlr,ovli,ovlt)
      Allocate(iipar(novl),nopen(novl),Eovl(novl), &
               ovlr(mch,novl),ovli(mch,novl),ovlt(ntarg,novl),ui(mch,novl))      

      i = 0
      rewind(nut)
      Do 
       read(nut,'(a)') line
       if(line(1:4).ne.'klsp') Cycle
       j = INDEX(line,'=')+1
       i = i + 1
       read(line(j:),*) iipar(i),m,n,Eovl(i) 
       nopen(i)=n
       read(nut,*) ovlr(1:n,i)
       read(nut,*) ovli(1:n,i)
       read(nut,*) ovlt(1:ntarg,i)
!       read(nut,*) ui(1:n,i)
       if(i.eq.novl) Exit
      End do

      End Subroutine R_overlaps_ion


!======================================================================
      Subroutine Find_overlap(ion_state,ii_ion,E,or,oi,phase)
!======================================================================
!     find required overlap
!----------------------------------------------------------------------
      Use overlaps_ion
      
      Implicit none
      Integer, intent(in) :: ion_state,ii_ion
      Real(8), intent(in) :: E
      Real(8), intent(out) :: or,oi,phase
      Integer :: i

      or = 0.d0; oi = 0.d0; phase = 0.d0
      Do i=1,novl
       if(iipar(i).ne.ii_ion) Cycle
       if(abs(Eovl(i)-E).gt.eps_e) Cycle
       if(ion_state.gt.nopen(i)) Exit
       or = ovlr(ion_state,i)
       oi = ovli(ion_state,i)
       phase = ui(ion_state,i)
       Exit
      End do 

      End Subroutine Find_overlap


!======================================================================
      Subroutine Find_ovl_ion(ilsp_ion,E,ntarg,ot)
!======================================================================
!     find <psedo|continuum>^2  for ion states (ntarg) 
!----------------------------------------------------------------------
      Use overlaps_ion
      
      Implicit none
      Integer, intent(in) :: ilsp_ion,ntarg
      Real(8), intent(in) :: E
      Real(8), intent(out) :: ot(ntarg)
      Integer :: i

      ot = 0.d0
      Do i=1,novl
       if(iipar(i).ne.ilsp_ion) Cycle
       if(abs(Eovl(i)-E).gt.eps_e) Cycle
       ot(1:ntarg) = ovlt(1:ntarg,i)
       Exit
      End do 

!     if(SUM(ot).eq.0.d0) Stop 'Find_ovl_ion - ?' 

      End Subroutine Find_ovl_ion
















