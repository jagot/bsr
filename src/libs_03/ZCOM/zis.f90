!=====================================================================
      Real(8) Function ZIS(H,N,A)
!=====================================================================
!     the simplest routine to evaluate the integral value 
!     according to the Simpson rule;
!     integral function A is supposed to given in N equal-spaced points
!     with interval H.
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N
      Real(8), intent(in) :: A(N),H
      Integer :: i
      Real(8) :: S1, S2
      S1=0.0; Do i=3,N,2; S1=S1+A(i); End do
      S2=0.0; Do i=2,N,2; S2=S2+A(i); End do
      ZIS=(2.d0*S1 + 4.d0*S2 + A(1)) * H/3.d0

      End Function ZIS
