!====================================================================
      Module rk_data
!====================================================================
!     contains a set of coefficients ordering accordint to
!     four pointers (packed with nine parameters)
!--------------------------------------------------------------------
      Implicit none

      Integer :: nrk = 0       ! current number of coefficients
      Integer :: mrk = 0       ! maximum dimension
      Integer :: irk = 2**10   ! initial dimension

! ... coefficients:

      Real(8), allocatable :: crk(:)

! ... their attributes:

      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)

! ... structure of the pointers:

      Integer ::  ibr=2**15, ibk=2**20, ibi=2**27, ibc=2**15

! ... kr1 = ic*(ic-1)/2+jc
! ... kr2 = int*ibi+k*ibk+idf
! ... kr3 = i1*ibr+i2
! ... kr4 = i3*ibr+i4

! ... ic,jc, i1,i2,i3,i4 < 2**15
! ... int <  8
! ... k < 2**7 = 128

      Real(8), parameter :: eps_C = 1.d-12

      End Module rk_data


!======================================================================
      Subroutine alloc_rk_data(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module "rk_data"
!----------------------------------------------------------------------
      Use rk_data

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

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
       Allocate(ia(nrk))
       ia(1:nrk)=kr1(1:nrk); Deallocate(kr1)
       Allocate(kr1(m)); kr1(1:nrk)=ia(1:nrk)
       ia(1:nrk)=kr2(1:nrk); Deallocate(kr2)
       Allocate(kr2(m)); kr2(1:nrk)=ia(1:nrk)
       ia(1:nrk)=kr3(1:nrk); Deallocate(kr3)
       Allocate(kr3(m)); kr3(1:nrk)=ia(1:nrk)
       ia(1:nrk)=kr4(1:nrk); Deallocate(kr4)
       Allocate(kr4(m)); kr4(1:nrk)=ia(1:nrk)
       Deallocate(ia)
       Allocate(ra(nrk))
       ra(1:nrk)=crk(1:nrk); Deallocate(crk)
       Allocate(crk(m)); crk(1:nrk)=ra(1:nrk)
       Deallocate(ra)
       mrk = m
       write(*,'(a,i810,a,f10.2,a)') ' Realloc_rk_data: new dimension = ', m, &
        '   memory requred = ', 24.d0*m/(1024*1024),'  Mb'
      end if

      End Subroutine alloc_rk_data


!======================================================================
      Subroutine Add_rk_data(ic,jc,int,kr,idf,i1,i2,i3,i4,C)
!======================================================================
!     add new data to the list
!----------------------------------------------------------------------
      USE rk_data

      Implicit none
      Integer, intent(in) :: int,ic,jc,kr,i1,i2,i3,i4,idf
      Real(8), intent(in) :: C
      Integer :: k1,k2,k3,k4, i,k,l,m

      if(mrk.eq.0) Call alloc_rk_data(irk)

      if(abs(C).lt.eps_C) Return

      if(ic.lt.jc.and.int.gt.1) Stop 'ic < jc in Add_rk_data'

      k1=ic*(ic-1)/2+jc; if(int.le.1) k1=ic*ibc+jc
                         if(jc .eq.0) k1=ic

      k2=int*ibi+kr*ibk+idf
      k3=i1*ibr+i3
      k4=i2*ibr+i4

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
      if(nrk.eq.mrk) Call alloc_rk_data(mrk+irk)

      End Subroutine Add_rk_data


