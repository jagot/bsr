!======================================================================
      Subroutine LEGPOL(T,mk,PL)
!======================================================================
!     CALCULATION OF LEGENDRE POLYNOMIALS BY USING A STANDARD
!     RECURRENCE RELATION
!
!     T  - ANGLE IN RADIANS
!     mk - maximun needed index
!     PL - result
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: mk
      Real(8), intent(in) :: T
      Real(8), intent(out) :: PL(mk)
      Integer :: n
      Real(8) :: A,B,C,X

      X=COS(T); A=1.0D0; B=X; PL(1)=1.D0; PL(2)=X
      Do n=2,mk-1
       C=(B*X*(n+n-1)-(n-1)*A)/n; A=B; B=C; PL(n+1)=C
      End do

      End Subroutine LEGPOL

