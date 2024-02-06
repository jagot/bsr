!======================================================================
      Module new_dets
!======================================================================
!     contains the list of i,j,S constructions
!     (used to save possible overlaps <kl|nl>)
!     new entries (i,j,S) are just added to the list
!----------------------------------------------------------------------
      Implicit none

      Integer :: mndet = 0      !  max.number of new overlaps
      Integer :: nndet = 0      !  curent number of new overlaps
      Integer :: indet = 2**10  !  initial suggestion for mndet

      Real(8) :: eps_det = 1.d-12    ! ???

!     value of the bound-orbitals overlap determinants:

      Real(8), allocatable :: ADET(:)

!     pointer on the continuum orbitals if any:

      Integer, allocatable :: IZOD(:),JZOD(:)

      End Module new_dets


!======================================================================
      Subroutine Allocate_ndets(m)
!======================================================================
!     allocate, deallocate or Reallocate arrays in module new_dets
!----------------------------------------------------------------------
      Use new_dets

      Implicit none
      Integer, Intent(in) :: m
      Real(8), allocatable :: rarray(:)
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(Allocated(ADET)) Deallocate(ADET,IZOD,JZOD)
       mndet = indet; nndet = 0
       Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
      elseif(m.eq.0) then
       if(Allocated(ADET)) Deallocate(ADET,IZOD,JZOD)
       mndet = 0; nndet = 0
      elseif(m.gt.mndet) then
       if(nndet.le.0) then
        if(Allocated(ADET)) Deallocate(ADET,IZOD,JZOD)
        mndet = m; nndet = 0
        Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
       else
        Allocate(rarray(nndet)); rarray(1:nndet)=ADET(1:nndet)
        Deallocate(ADET); Allocate(ADET(m))
        ADET(1:nndet)=rarray(1:nndet); Deallocate(rarray)
        Allocate(iarray(nndet)); iarray(1:nndet)=IZOD(1:nndet)
        Deallocate(IZOD); Allocate(IZOD(m))
        IZOD(1:nndet)=iarray(1:nndet)
        iarray(1:nndet)=JZOD(1:nndet)
        Deallocate(JZOD); Allocate(JZOD(m))
        JZOD(1:nndet)=iarray(1:nndet); Deallocate(iarray)
       end if
      end if

      End Subroutine Allocate_ndets


!======================================================================
      Subroutine Iadd_ndets(io,jo,S)
!======================================================================
!     add new entry in the list new_dets
!----------------------------------------------------------------------
      Use new_dets

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S

      if(abs(S).lt.eps_det) Return

      if(nndet+1.gt.mndet) Call Allocate_ndets(mndet+indet)

      nndet = nndet + 1
      IZOD(nndet) = io
      JZOD(nndet) = jo
      ADET(nndet) = S

      End Subroutine Iadd_ndets


