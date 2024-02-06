!======================================================================
      Subroutine pri_matr(iout,ii,jj,n,m,R,A)
!======================================================================
!     print matrix R(n,m) with name A
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: iout,ii,jj,n,m
      Real(8), intent(in) :: R(n,m)
      Character(*) A
      Integer :: i,j

      write (iout,'(/a/)') A
      Do i = 1, ii
        write (iout,'(10D13.5)') (R(i,j),j = 1,jj)
      End do

      End Subroutine pri_matr
