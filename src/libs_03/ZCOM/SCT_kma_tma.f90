!======================================================================
      Subroutine kma_tma (n,kma,tr,ti)
!======================================================================
!     T = 1 - S = -2i K/(1-iK) = (2K^2 - 2iK) / (1+K^2) = -2i [P + i Q]
!
!     P = K /(..);   Q = -K^2/(..) = -(1 - 1/(..))= 1/(..)-1
!
!     n - number of open channels
!----------------------------------------------------------------------
      Implicit real(8) (a-h,o-z)
      Integer, intent(in)  :: n
      Real(8), intent(in)  :: kma(*)
      Real(8), intent(out) :: tr(*),ti(*)
      Real(8) :: AK(n,n),AK1(n,n),AK2(n,n)

      Do i=1,n                               !   AK --> K-matrix
       Do j=1,i
        ij=(i-1)*i/2+j; AK(i,j)=kma(ij); AK(j,i)=kma(ij)
       End do
      End do

      Do i=1,n                               !   AK1 --> 1 + K^2
       Do j=i,n                              !   AK2 -->     K^2
        x=SUM(AK(i,:)*AK(:,j))
        AK2(i,j)=x; AK2(j,i)=x; AK1(i,j)=x; AK1(j,i)=x
        if(i.eq.j) AK1(i,i)=AK1(i,i)+1
       End do
      End do

      Call Inv(n,n,AK1)                      !  AK1 --> AK1^-1

      Do i=1,n
       Do j=i,n
        x=SUM(AK2(i,:)*AK1(:,j))
        y=SUM(AK (i,:)*AK1(:,j))
        ij=(j-1)*j/2+i; TR(ij)= +2.d0*x; TI(ij)= -2.d0*y
       End do
      End do

      End Subroutine kma_tma


!======================================================================
      Subroutine kma_tmb (n,kma,tr,ti,nb,kb)
!======================================================================
!
!     T = 1 - S = -2iK/(1-iK) = (2K^2 - 2iK) / (1+K^2)
!
!     n - number of open channels
!     nb - number of required initial channels
!     kb - resulting size of T-matrix arrays
!----------------------------------------------------------------------

      Implicit real(8) (a-h,o-z)

      Integer, intent(in)  :: n,nb
      Real(8), intent(in)  :: kma(*)
      Real(8), intent(out) :: tr(*),ti(*)

      Real(8), Allocatable :: AK(:,:),AK1(:,:),AK2(:,:)

      Allocate (AK(n,n),AK1(n,n),AK2(n,n))

      Do i=1,n                               !   AK --> K-matrix
       Do j=1,i
        ij=(i-1)*i/2+j; AK(i,j)=kma(ij); AK(j,i)=kma(ij)
       End do
      End do

      Do i=1,n                               !   AK1 --> 1 + K^2
       Do j=i,n                              !   AK2 -->     K^2
        x=SUM(AK(i,:)*AK(:,j))
        AK2(i,j)=x; AK2(j,i)=x; AK1(i,j)=x; AK1(j,i)=x
        if(i.eq.j) AK1(i,i)=AK1(i,i)+1.d0
       End do
      End do

      Call Inv(n,n,AK1)                      !  AK1 --> AK1^-1

      kb = 0
      Do i=1,nb
       Do j=i,n
        x=SUM(AK2(i,:)*AK1(:,j))
        y=SUM(AK (i,:)*AK1(:,j))
        kb = kb +1; TR(kb)= +2.d0*x; TI(kb)= -2.d0*y
       End do
      End do

      End Subroutine kma_tmb

