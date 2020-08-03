!================================================================================
      Subroutine Laguerre(x,a,m,P)
!================================================================================
!     provides Lagguerre functions L^a_n at point x for n=0:m in array P(0:m)
!     using the recurrence formular
!           i L(a,i) = (2i+a-1-x) L(a,i-1) - (i-1+a) L(a,i-2)
!     with  L(a,0) = 1,  L(a,1) = -x+a+1
!--------------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: m                 !  >= 2
      Real(8), intent(in) :: x, a              !  x > 0;  a > -1
      Real(8), intent(out) :: P(0:m)
      Integer :: i

! ... check the input parameters:
      
      if(m.lt.2) Stop 'Laguerre:  m < 2'
      if(a.le.-1.d0) Stop 'Laguerre:  a < -1.d0'
      if(x.lt. 0.d0) Stop 'Laguerre:  a < -1.d0'

! ... recurrence calculations:

      P(0) = 1.d0;  P(1) = -x+a-1.d0 
      Do i=2,m
       P(i) = ( (2*i-1+a-x)*P(i-1) - (i-1+a)*P(i-2)) / i
      End do

      End Subroutine Laguerre 