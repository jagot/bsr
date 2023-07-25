!======================================================================
      Module orthogonality
!======================================================================
!     contains the desription of one-electron radial overlaps 
!----------------------------------------------------------------------
      Implicit none
   
      Integer :: morth = 0      !  max.number of  overlaps
      Integer :: north = 0      !  curent number of overlaps
      Integer :: korth = 2**10  !  initial suggestion for morth

!     pointer to orbitals (iorth > = jorth):

      Integer, allocatable :: iorth(:),jorth(:), lorth(:)

      End Module orthogonality

 
!======================================================================
      Subroutine Alloc_orthogonality(m)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use orthogonality

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(allocated(lorth)) Deallocate(lorth,iorth,jorth) 
       morth = korth; north = 0
       Allocate(lorth(morth),iorth(morth),jorth(morth))
      elseif(m.eq.0) then
       if(allocated(lorth)) Deallocate(lorth,iorth,jorth) 
       morth = 0; north = 0
      elseif(m.gt.morth) then
       if(north.le.0) then
        if(Allocated(lorth)) Deallocate(lorth,iorth,jorth) 
        morth = m; north = 0
        Allocate(lorth(morth),iorth(morth),jorth(morth))
       else 
        Allocate(iarray(north))
        iarray(1:north)=iorth(1:north);  Deallocate(iorth)
        Allocate(iorth(m)); iorth(1:north)=iarray(1:north)
        iarray(1:north)=jorth(1:north);  Deallocate(jorth)
        Allocate(jorth(m)); jorth(1:north)=iarray(1:north)
        iarray(1:north)=lorth(1:north);  Deallocate(lorth)
        Allocate(lorth(m)); lorth(1:north)=iarray(1:north)
        Deallocate(iarray)
        morth = m
       end if
      end if

      End Subroutine Alloc_orthogonality


!======================================================================
      Subroutine Iadd_orth(io,jo,orth)
!======================================================================
!     add new entry( or replace) in the ordered list of overlaps
!----------------------------------------------------------------------
      Use orthogonality

      Implicit none
      Integer, intent(in) :: io,jo,orth
      Integer :: i,j, k,l,m

      if(morth.eq.0) Call Alloc_orthogonality(korth)

      i = max(io,jo)
      j = min(io,jo)

! ... search position (m) for given overlap

      k=1; l=north
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (i.lt.iorth(m)) then;     l = m - 1
      elseif(i.gt.iorth(m)) then;     k = m + 1
      else
       if    (j.lt.jorth(m)) then;    l = m - 1
       elseif(j.gt.jorth(m)) then;    k = m + 1
       else
        lorth(m)=orth;  Return 
       end if
      end if
      go to 1
    2 Continue 
      
! ... shift the rest data up:

      Do l=north,k,-1; m = l + 1
       lorth(m)=lorth(l); iorth(m)=iorth(l); jorth(m)=jorth(l)
      End do

! ... add new overlap:

      lorth(k)=orth; iorth(k)=i; jorth(k)=j; north=north+1
      if(north.eq.morth) Call alloc_orthogonality(morth+korth) 

      End Subroutine Iadd_orth


!======================================================================
      Integer Function IORT(io,jo)
!======================================================================
!     find orthogonality for the orbitals io,jo
!----------------------------------------------------------------------
      Use orthogonality

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j, k,l,m

      i = max(io,jo)
      j = min(io,jo) 

! ... search position (m) for given overlap

      k=1; l=north
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (i.lt.iorth(m)) then;     l = m - 1
      elseif(i.gt.iorth(m)) then;     k = m + 1
      else
       if    (j.lt.jorth(m)) then;    l = m - 1
       elseif(j.gt.jorth(m)) then;    k = m + 1
       else
        IORT=lorth(m);  Return 
       end if
      end if
      go to 1
    2 IORT = 0 

      End Function IORT


