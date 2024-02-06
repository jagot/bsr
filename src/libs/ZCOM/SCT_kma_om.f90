!======================================================================
      Subroutine kma_om (mch,n,g,kma,om)
!======================================================================
!     produce OMEGA-matrix from K-matrix, stored in symmetric form
!     n - number of open channels
!     g - stat.weihgt
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: mch,n
      Real(8), intent(in)  :: g, kma(*)
      Real(8), intent(out) :: om(mch,mch)
      Integer :: i,j,ij
      Real(8) :: x,y, AK(n,n),AK1(n,n),AK2(n,n)

      Call Put_ijmatr(n,AK,n,kma)      !  kma -> AK, as a full matrix

      AK2 = MATMUL(AK,AK)
      AK1 = AK2
      Do i=1,n; AK1(i,i)=AK1(i,i)+1; End do

      Call Inv(n,n,AK1)                 !  A --> A^-1

      Do i=1,n
       Do j=i,n
        x = -2.0*SUM(AK2(i,:)*AK1(:,j))
        y =  2.0*SUM(AK (i,:)*AK1(:,j))
        om(i,j) = (x*x + y*y) * g
        om(j,i) = om(i,j)
       End do
      End do

      End Subroutine kma_om
