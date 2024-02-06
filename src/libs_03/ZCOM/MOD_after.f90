!======================================================================
      Module after_conditions
!======================================================================
!     contains the desription of one-electron radial overlaps
!----------------------------------------------------------------------
      Implicit none

      Integer :: mafter = 0      !  max.number of  overlaps
      Integer :: nafter = 0      !  curent number of overlaps
      Integer :: kafter = 2**10  !  initial suggestion for mafter

!     pointer to orbitals (iafter > = jafter):

      Integer, allocatable :: iafter(:),jafter(:)

      End Module after_conditions


!======================================================================
      Subroutine Alloc_after_conditions(m)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use after_conditions

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(allocated(iafter)) Deallocate(iafter,jafter)
       mafter = kafter; nafter = 0
       Allocate(iafter(mafter),jafter(mafter))
      elseif(m.eq.0) then
       if(allocated(iafter)) Deallocate(iafter,jafter)
       mafter = 0; nafter = 0
      elseif(m.gt.mafter) then
       if(nafter.le.0) then
        if(Allocated(iafter)) Deallocate(iafter,jafter)
        mafter = m; nafter = 0
        Allocate(iafter(mafter),jafter(mafter))
       else
        Allocate(iarray(nafter))
        iarray(1:nafter)=iafter(1:nafter);  Deallocate(iafter)
        Allocate(iafter(m)); iafter(1:nafter)=iarray(1:nafter)
        iarray(1:nafter)=jafter(1:nafter);  Deallocate(jafter)
        Allocate(jafter(m)); jafter(1:nafter)=iarray(1:nafter)
        Deallocate(iarray)
        mafter = m
       end if
      end if

      End Subroutine Alloc_after_conditions


!======================================================================
      Integer Function after(io,jo)
!======================================================================
!     find after_conditions for the orbitals io,jo
!----------------------------------------------------------------------
      Use after_conditions

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j, k,l,m, kk

      if(io.eq.jo) then; after=0; Return; end if

      if(mafter.eq.0) Call Alloc_after_conditions(kafter)

! ... search record for given position:

      i = io; j = jo

      k=1; l=nafter
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (i.lt.iafter(m)) then;     l = m - 1
      elseif(i.gt.iafter(m)) then;     k = m + 1
      else
       if    (j.lt.jafter(m)) then;    l = m - 1
       elseif(j.gt.jafter(m)) then;    k = m + 1
       else
        after=1;  Return
       end if
      end if
      go to 1
    2 kk=k

! ... search record for opposite position:

      i = jo; j = io

      k=1; l=nafter
    3 if(k.gt.l) go to 4
      m=(k+l)/2
      if    (i.lt.iafter(m)) then;     l = m - 1
      elseif(i.gt.iafter(m)) then;     k = m + 1
      else
       if    (j.lt.jafter(m)) then;    l = m - 1
       elseif(j.gt.jafter(m)) then;    k = m + 1
       else
        after=-1;  Return
       end if
      end if
      go to 3
    4 k=kk

! ... shift the rest data up:

      Do l=nafter,k,-1; m = l + 1
       iafter(m)=iafter(l); jafter(m)=jafter(l)
      End do

! ... add new overlap:

      iafter(k)=io; jafter(k)=jo; nafter=nafter+1
      if(nafter.eq.mafter) Call alloc_after_conditions(mafter+kafter)

      after = 1

      End Function after



