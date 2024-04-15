!======================================================================
      Module multipole_integrals
!======================================================================
!     contains the desription of one-electron radial overlaps
!----------------------------------------------------------------------
      Implicit none

      Integer :: mdip = 0      !  max.number of  overlaps
      Integer :: ndip = 0      !  curent number of overlaps
      Integer :: kdip = 2**10  !  initial suggestion for mdip

!     value of the bound-bound one-electron overlap:

      Real(8), allocatable :: CL(:), CV(:)

!     pointer to orbitals:

      Integer, allocatable :: idip(:),jdip(:)

      End Module multipole_integrals


!======================================================================
      Subroutine Alloc_multipole_integrals(m)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use multipole_integrals

      Implicit none
      Integer, intent(in) :: m
      Real(8), allocatable :: rarray(:)
      Integer, allocatable :: iarray(:)

      if(m.lt.0) then
       if(allocated(CL)) Deallocate(CL,CV,idip,jdip)
       mdip = kdip; ndip = 0
       Allocate(CL(mdip),CV(mdip),idip(mdip),jdip(mdip))
      elseif(m.eq.0) then
       if(allocated(CL)) Deallocate(CL,CV,idip,jdip)
       mdip = 0; ndip = 0
      elseif(m.gt.mdip) then
       if(ndip.le.0) then
        if(Allocated(CL)) Deallocate(CL,CV,idip,jdip)
        mdip = m; ndip = 0
        Allocate(CL(mdip),CV(mdip),idip(mdip),jdip(mdip))
       else
        Allocate(rarray(ndip))
        rarray(1:ndip)=CL(1:ndip); Deallocate(CL)
        Allocate(CL(m)); CL(1:ndip)=rarray(1:ndip)
        rarray(1:ndip)=CV(1:ndip); Deallocate(CV)
        Allocate(CV(m)); CV(1:ndip)=rarray(1:ndip)
        Deallocate(rarray)
        Allocate(iarray(ndip))
        iarray(1:ndip)=idip(1:ndip); Deallocate(idip)
        Allocate(idip(m));  idip(1:ndip)=iarray(1:ndip)
        iarray(1:ndip)=jdip(1:ndip); Deallocate(jdip)
        Allocate(jdip(m)); jdip(1:ndip)=iarray(1:ndip)
        Deallocate(iarray)
        mdip = m
       end if
      end if

      End Subroutine Alloc_multipole_integrals


!======================================================================
      Subroutine Iadd_dip(io,jo,SL,SV)
!======================================================================
!     add new entry( or replace) in the ordered list of overlaps
!----------------------------------------------------------------------
      Use multipole_integrals

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: SL,SV
      Integer :: i,j, k,l,m

      if(mdip.eq.0) Call Alloc_multipole_integrals(kdip)

      i=io   ! i = max(io,jo)
      j=jo   ! j = min(io,jo)

! ... search position (m) for given overlap

      k=1; l=ndip
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (i.lt.idip(m)) then;     l = m - 1
      elseif(i.gt.idip(m)) then;     k = m + 1
      else
       if    (j.lt.jdip(m)) then;    l = m - 1
       elseif(j.gt.jdip(m)) then;    k = m + 1
       else
        CL(m)=SL; CV(m)=SV; Return
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:

      Do l=ndip,k,-1; m = l + 1
       CL(m)=CL(l); CV(m)=CV(l); idip(m)=idip(l); jdip(m)=jdip(l)
      End do

! ... add new overlap:

      CL(k)=SL; CV(k)=SV; idip(k)=i; jdip(k)=j; ndip=ndip+1
      if(ndip.eq.mdip) Call alloc_multipole_integrals(mdip+kdip)

      End Subroutine Iadd_dip


!======================================================================
      Subroutine dip(io,jo,SL,SV)
!======================================================================
!     find overlaps <io|jo>
!----------------------------------------------------------------------
      Use multipole_integrals

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j, k,l,m
      Real(8), intent(out) :: SL,SV

      i=io   ! i = max(io,jo)
      j=jo   ! j = min(io,jo)

! ... search position (m) for given overlap

      k=1; l=ndip
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (i.lt.idip(m)) then;     l = m - 1
      elseif(i.gt.idip(m)) then;     k = m + 1
      else
       if    (j.lt.jdip(m)) then;    l = m - 1
       elseif(j.gt.jdip(m)) then;    k = m + 1
       else
        SL=CL(m); SV=CV(m); Return
       end if
      end if
      go to 1
    2 SL = 0.d0;  SV = 0.d0

      End Subroutine dip


