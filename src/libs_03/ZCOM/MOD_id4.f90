!====================================================================
      Module id4_data
!====================================================================
!     contains a set of coefficients ordering accordint to
!     four indexes; similar to rk4_data;  used in dbsr_mchf
!--------------------------------------------------------------------
      Implicit none 
    
      Integer :: nid = 0       ! current number of coefficients
      Integer :: mid = 0       ! maximum dimension
      Integer :: iid = 10000   ! initial dimension

! ... coefficients:

      Real(8), allocatable :: cid(:)   
    
! ... their attributes:

      Integer, allocatable :: id1(:),id2(:),id3(:),id4(:)

      End Module id4_data


!======================================================================
      Subroutine alloc_id4_data(m)
!======================================================================
! ... allocate, deallocate, or reallocate the data
!----------------------------------------------------------------------
      Use id4_data

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(cid)) Deallocate (cid,id1,id2,id3,id4)
       nid = 0; mid = 0 
      elseif(.not.allocated(cid)) then
       mid = m; nid = 0
       Allocate(cid(mid),id1(mid),id2(mid),id3(mid),id4(mid))
      elseif(m.le.mid) then
       Return
      elseif(nid.eq.0) then
       Deallocate (cid,id1,id2,id3,id4)
       mid = m
       Allocate(cid(mid),id1(mid),id2(mid),id3(mid),id4(mid))
      else
       Allocate(rar(nid))
       rar=cid; Deallocate(cid); Allocate(cid(m)); cid(1:nid)=rar
       Deallocate(rar)
       Allocate(iar(nid))
       iar=id1; Deallocate(id1); Allocate(id1(m)); id1(1:nid)=iar
       iar=id2; Deallocate(id2); Allocate(id2(m)); id2(1:nid)=iar
       iar=id3; Deallocate(id3); Allocate(id3(m)); id3(1:nid)=iar
       iar=id4; Deallocate(id4); Allocate(id4(m)); id4(1:nid)=iar
       Deallocate(iar)
       mid = m
      end if

      End Subroutine alloc_id4_data


!======================================================================
      Subroutine Add_id4_data(k1,k2,k3,k4,C)
!======================================================================
!     add new data to the list 
!----------------------------------------------------------------------
      Use id4_data

      Implicit none
      Integer, intent(in) ::  k1,k2,k3,k4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m

      if(mid.eq.0) Call alloc_id4_data(iid)

! ... search position (k) for new integral

      k=1; l=nid
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (k1.lt.id1(m)) then;       l = m - 1
      elseif(k1.gt.id1(m)) then;       k = m + 1
      else
       if    (k2.lt.id2(m)) then;      l = m - 1
       elseif(k2.gt.id2(m)) then;      k = m + 1
       else
        if    (k3.lt.id3(m)) then;     l = m - 1
        elseif(k3.gt.id3(m)) then;     k = m + 1
        else
         if    (k4.lt.id4(m)) then;    l = m - 1
         elseif(k4.gt.id4(m)) then;    k = m + 1
         else
          cid(m)=cid(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=nid,k,-1
       m = i + 1
       cid(m)=cid(i)
       id1(m)=id1(i); id2(m)=id2(i); id3(m)=id3(i); id4(m)=id4(i)
      End do

! ... add new integral:

      cid(k)=C; id1(k)=k1; id2(k)=k2; id3(k)=k3; id4(k)=k4; nid=nid+1
      if(nid.eq.mid) Call alloc_id4_data(mid+iid) 

      End Subroutine Add_id4_data

