!====================================================================
      Integer Function Jparity(n,l,q)
!====================================================================
!     defines the parity for given configuration {l^q}
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n
      Integer, intent(in) :: l(n),q(n)
      Integer :: m,i
      m=0;   Do i=1,n; m=m+l(i)*q(i);  End do
      Jparity=(-1)**m
      End Function Jparity


