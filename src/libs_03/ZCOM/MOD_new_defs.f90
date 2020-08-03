!=====================================================================
      Module new_defs
!=====================================================================
!     contains the list of i,j,S constructions
!     (used to save possible determinant factors between two channels)
!     new entries (i,j,S) are added or summed for the same indeces
!--------------------------------------------------------------------
      Implicit none
      Integer :: mndef = 0      !  max.number of new overlaps
      Integer :: nndef = 0      !  curent number of new overlaps
      Integer :: indef = 2**10  !  initial suggestion for mndef

!     values of the bound-orbitals overlap deferminants:

      Real(8), allocatable :: adef(:)

!     pointer to the overlaps with continuum orbitals:

      Integer, allocatable :: iof(:),jof(:)

      End Module new_defs


!======================================================================
      Subroutine Allocate_ndefs(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module new_defs
!----------------------------------------------------------------------
      Use new_defs

      Implicit none
      Integer, intent(in) :: m
      Real(8), allocatable :: rarray(:)
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(allocated(adef)) Deallocate(adef,iof,jof) 
        mndef = indef; nndef = 0
       Allocate(adef(mndef),iof(mndef),jof(mndef))
      elseif(m.eq.0) then
       if(allocated(adef)) Deallocate(adef,iof,jof) 
        mndef = 0; nndef = 0
      elseif(m.gt.mndef) then
       if(nndef.eq.0) then
        if(allocated(adef)) Deallocate(adef,iof,jof) 
        mndef = m; nndef = 0
        Allocate(adef(mndef),iof(mndef),jof(mndef))
       else 
        Allocate(rarray(nndef)); rarray(1:nndef)=adef(1:nndef)
        Deallocate(adef); Allocate(adef(m))
        adef(1:nndef)=rarray(1:nndef); Deallocate(rarray)
        Allocate(iarray(nndef)); iarray(1:nndef)=iof(1:nndef)
        Deallocate(iof); Allocate(iof(m))
        iof(1:nndef)=iarray(1:nndef)
        iarray(1:nndef)=jof(1:nndef)
        Deallocate(jof); Allocate(jof(m))
        jof(1:nndef)=iarray(1:nndef); Deallocate(iarray)
       end if
      end if

      End Subroutine Allocate_ndefs


!======================================================================
      Subroutine Iadd_ndefs(io,jo,S)
!======================================================================
!     add new entry in module new_defs
!----------------------------------------------------------------------
      Use new_defs

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S
      Integer :: i,m

      m = 0
      Do i=1,nndef
       if(io.ne.iof(i)) Cycle
       if(jo.ne.jof(i)) Cycle
       adef(i) = adef(i) + S
       m = 1
       Exit
      End do
      if(m.ne.0) Return

      if(nndef+1.gt.mndef) Call Allocate_ndefs(mndef+indef)
      nndef = nndef + 1
      iof(nndef) = io
      jof(nndef) = jo
      adef(nndef) = S

      End Subroutine Iadd_ndefs


