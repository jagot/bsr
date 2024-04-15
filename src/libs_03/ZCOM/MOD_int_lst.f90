!====================================================================
      Module int_list
!====================================================================
!     contains a set of integrals with 6 identifiers
!--------------------------------------------------------------------
      Implicit none

      Integer :: mint =  0      !  max. number of coefficients
      Integer :: nint =  0      !  current number of coefficients
      Integer :: kint = 2**10   !  initial space

      Integer, allocatable :: Jcase(:),Jpol(:),J1(:),J2(:),J3(:),J4(:)

      End Module int_list


!======================================================================
      Subroutine alloc_int(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module "int_list"
!----------------------------------------------------------------------
      Use int_list

      Implicit none
      Integer, Intent(in) :: m
      Integer, allocatable :: iar(:)

      if(m.le.0) then
       if(allocated(Jcase)) Deallocate (Jcase,Jpol,J1,J2,J3,J4)
       mint = 0; nint = 0; Return
      elseif(.not.allocated(Jcase)) then
       mint = m; Allocate(Jcase(m),Jpol(m),J1(m),J2(m),J3(m),J4(m))
      elseif(m.le.mint) then
       Return
      elseif(nint.eq.0) then
       Deallocate (Jcase,Jpol,J1,J2,J3,J4);  mint = m
       Allocate(Jcase(m),Jpol(m),J1(m),J2(m),J3(m),J4(m))
      else
       Allocate(iar(nint))
       iar(1:nint)=Jcase(1:nint); Deallocate(Jcase)
       Allocate(Jcase(m)); Jcase(1:nint)=iar(1:nint)
       iar(1:nint)=Jpol(1:nint); Deallocate(Jpol)
       Allocate(Jpol(m)); Jpol(1:nint)=iar(1:nint)
       iar(1:nint)=J1(1:nint); Deallocate(J1)
       Allocate(J1(m)); J1(1:nint)=iar(1:nint)
       iar(1:nint)=J2(1:nint); Deallocate(J2)
       Allocate(J2(m)); J2(1:nint)=iar(1:nint)
       iar(1:nint)=J3(1:nint); Deallocate(J3)
       Allocate(J3(m)); J3(1:nint)=iar(1:nint)
       iar(1:nint)=J4(1:nint); Deallocate(J4)
       Allocate(J4(m)); J4(1:nint)=iar(1:nint)
       Deallocate(iar)
       mint = m
       write(*,'(a,i10/a,f10.2,a)') ' Realloc_int: new dimension = ', mint, &
        ' memory requred = ', 24.d0*m/(1024*1024),'  Mb'
      end if

      End Subroutine alloc_INT


!=======================================================================
      Integer Function Iadd_int(icase,kpol,i1,i2,i3,i4)
!=======================================================================
!     add the integral to the list 'int_list'
!-----------------------------------------------------------------------
      Use int_list

      Implicit none
      Integer, Intent(in) :: icase,kpol,i1,i2,i3,i4
      Integer :: i

      if(mint.eq.0) Call Alloc_int(kint)

! ... check if the same integral is already in list:

      Do i=1,nint
       if(icase.ne.Jcase(i)) Cycle
       if(kpol.ne.Jpol(i)) Cycle
       if(i1.ne.J1(i)) Cycle
       if(i2.ne.J2(i)) Cycle
       if(i3.ne.J3(i)) Cycle
       if(i4.ne.J4(i)) Cycle
       Iadd_int=i; Return
      End do

! ... add new integral:

      if(nint.eq.mint) Call Alloc_int(mint+kint)
      nint=nint+1; Jcase(nint)=icase; Jpol(nint)=kpol
      J1(nint)=i1; J2(nint)=i2; J3(nint)=i3; J4(nint)=i4
      Iadd_int = nint

      End Function Iadd_int
