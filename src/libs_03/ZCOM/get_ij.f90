!======================================================================
      Subroutine Get_ij (k,i,j)
!======================================================================
!     defines i,j from k = (i-1)*i/2+j,  i >= j
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k
      Integer, intent(out) :: i,j

      i = sqrt(real(2*k))
      j = k - i*(i-1)/2
      if(j.gt.i) then
       i = i+1
       j = k - i*(i-1)/2
      end if

      End Subroutine Get_ij


!======================================================================
      Integer Function Def_ij (i,j)
!======================================================================
!     (i-1)*i/2+j
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i,j
      Integer :: ii,jj

      ii = max(i,j); jj=min(i,j); Def_ij=ii*(ii-1)/2+jj

      End Function Def_ij


!======================================================================
      Integer(8) Function Def_ij8 (i,j)
!======================================================================
!     (i-1)*i/2+j
!----------------------------------------------------------------------
      Implicit none
      Integer(4), intent(in) :: i,j
      Integer(8) :: ii,jj

      ii = max(i,j); jj = min(i,j); Def_ij8 = ii*(ii-1)/2+jj

      End Function Def_ij8



!======================================================================
      Subroutine Get_ijmatr (m,A,n,matr)
!======================================================================
      Implicit none
      Integer, intent(in) :: m,n
      Real(8), intent(in) :: A(m,*)
      Real(8), intent(out) :: matr(*)
      Integer :: i,j,ij

      Do i=1,n; Do j=1,i
       ij = i*(i-1)/2+j; matr(ij)=A(i,j)
      End do; End do

      End Subroutine Get_ijmatr


!======================================================================
      Subroutine Put_ijmatr (m,A,n,matr)
!======================================================================
      Implicit none
      Integer, intent(in) :: m,n
      Real(8), intent(in) :: matr(*)
      Real(8), intent(out) :: A(m,*)
      Integer :: i,j,ij

      Do i=1,n; Do j=1,i
       ij = i*(i-1)/2+j; A(i,j)=matr(ij); A(j,i)=matr(ij)
      End do; End do

      End Subroutine Put_ijmatr

