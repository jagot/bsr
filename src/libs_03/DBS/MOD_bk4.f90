!====================================================================
      Module bk4_data
!====================================================================
!     contains a set of coefficients ordering accordint to
!     four indexes
!--------------------------------------------------------------------
      Implicit none
      Integer :: nbk = 0       ! current number of coefficients
      Integer :: mbk = 0       ! maximum dimension
      Integer :: ibk = 10000   ! initial dimension
! ... coefficients:
      Real(8), allocatable :: cbk(:)
! ... their attributes:
      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)
      End Module bk4_data


!======================================================================
      Subroutine alloc_bk4_data(m)
!======================================================================
! ... allocate, deallocate, or reallocate the arrays in module bk4_data
!----------------------------------------------------------------------
      Use bk4_data

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(cbk)) Deallocate (cbk,kr1,kr2,kr3,kr4)
       nbk = 0; mbk = 0
      elseif(.not.allocated(cbk)) then
       mbk = m; nbk = 0
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      elseif(m.le.mbk) then
       Return
      elseif(nbk.eq.0) then
       Deallocate (cbk,kr1,kr2,kr3,kr4)
       mbk = m
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      else
       Allocate(rar(nbk))
       rar=cbk; Deallocate(cbk); Allocate(cbk(m)); cbk(1:nbk)=rar
       Deallocate(rar)
       Allocate(iar(nbk))
       iar=kr1; Deallocate(kr1); Allocate(kr1(m)); kr1(1:nbk)=iar
       iar=kr2; Deallocate(kr2); Allocate(kr2(m)); kr2(1:nbk)=iar
       iar=kr3; Deallocate(kr3); Allocate(kr3(m)); kr3(1:nbk)=iar
       iar=kr4; Deallocate(kr4); Allocate(kr4(m)); kr4(1:nbk)=iar
       Deallocate(iar)
       mbk = m
      end if

      End Subroutine alloc_bk4_data

!======================================================================
      Subroutine add_bk4_data(k1,k2,k3,k4,C)
!======================================================================
!     add new data to the list bk4_data
!----------------------------------------------------------------------
      Use bk4_data

      Implicit none
      Integer, intent(in) ::  k1,k2,k3,k4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m

      if(mbk.eq.0) Call alloc_bk4_data(ibk)

! ... search position (k) for new integral

      k=1; l=nbk
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (k1.lt.kr1(m)) then;       l = m - 1
      elseif(k1.gt.kr1(m)) then;       k = m + 1
      else
       if    (k2.lt.kr2(m)) then;      l = m - 1
       elseif(k2.gt.kr2(m)) then;      k = m + 1
       else
        if    (k3.lt.kr3(m)) then;     l = m - 1
        elseif(k3.gt.kr3(m)) then;     k = m + 1
        else
         if    (k4.lt.kr4(m)) then;    l = m - 1
         elseif(k4.gt.kr4(m)) then;    k = m + 1
         else
          cbk(m)=cbk(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:
      
      if(k.le.nbk) then
      Do i=nbk,k,-1
       m = i + 1
       cbk(m)=cbk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do
      end if

! ... add new integral:

      cbk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nbk=nbk+1
      if(nbk.eq.mbk) Call alloc_bk4_data(mbk+ibk)

      End Subroutine add_bk4_data

