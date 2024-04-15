!====================================================================
      Module rk4_data
!====================================================================
!     contains a set of coefficients ordering accordint to
!     four indexes
!--------------------------------------------------------------------
      Implicit none

      Integer :: nrk = 0       ! current number of coefficients
      Integer :: mrk = 0       ! maximum dimension
      Integer :: irk = 2**20   ! initial dimension

! ... coefficients:

      Real(8), allocatable :: crk(:)

! ... their attributes:

      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)

      End Module rk4_data


!======================================================================
      Subroutine alloc_rk4_data(m)
!======================================================================
! ... allocate, deallocate, or reallocate the arrays in "rk4_data"
!----------------------------------------------------------------------
      Use rk4_data

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(crk)) Deallocate (crk,kr1,kr2,kr3,kr4)
       nrk = 0; mrk = 0
      elseif(.not.allocated(crk)) then
       mrk = m; nrk = 0
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      elseif(m.le.mrk) then
       Return
      elseif(nrk.eq.0) then
       Deallocate (crk,kr1,kr2,kr3,kr4)
       mrk = m
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      else
       Allocate(rar(nrk))
       rar=crk; Deallocate(crk); Allocate(crk(m)); crk(1:nrk)=rar
       Deallocate(rar)
       Allocate(iar(nrk))
       iar=kr1; Deallocate(kr1); Allocate(kr1(m)); kr1(1:nrk)=iar
       iar=kr2; Deallocate(kr2); Allocate(kr2(m)); kr2(1:nrk)=iar
       iar=kr3; Deallocate(kr3); Allocate(kr3(m)); kr3(1:nrk)=iar
       iar=kr4; Deallocate(kr4); Allocate(kr4(m)); kr4(1:nrk)=iar
       Deallocate(iar)
       mrk = m
       write(*,'(a,i10/a,f10.2,a)') ' Realloc_rk4_data: new dimension = ', m, &
        ' memory requred = ', 24.d0*m/(1024*1024),'  Mb'
      end if

      End Subroutine alloc_rk4_data


!======================================================================
      Subroutine Add_rk4_data(k1,k2,k3,k4,C)
!======================================================================
!     add new data to the list
!----------------------------------------------------------------------
      Use rk4_data

      Implicit none
      Integer, intent(in) ::  k1,k2,k3,k4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m

      if(mrk.eq.0) Call alloc_rk4_data(irk)

! ... search position (k) for new integral

      k=1; l=nrk
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
          crk(m)=crk(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:

      Do i=nrk,k,-1
       m = i + 1
       crk(m)=crk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do

! ... add new integral:

      crk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nrk=nrk+1
      if(nrk.eq.mrk) Call alloc_rk4_data(mrk+irk)

      End Subroutine Add_rk4_data


