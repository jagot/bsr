!======================================================================
      Real(8) FUNCTION QUADR(I,J,KK)
!======================================================================
!
!                                 kk
!     Evaluates the integral of  r   P (r) P (r) with respect to r
!                                     i     j
!----------------------------------------------------------------------

      USE RADIAL, MX => mro

      IMPLICIT NONE
      Integer(4), Intent(in) :: i,j,kk
      Integer(4) :: k,m,mm
      REAL(8) :: A,B,DEN,ZR

! ... region (0,r1) ...

      m   = mexp(i) + mexp(j) + kk
      DEN = m + 1
      A  = aexp(i) + aexp(j)
      B  = bexp(i) + bexp(j) + aexp(i)*aexp(j)
      B  = (DEN+D1)*A**2 - D2*B/(DEN+D2)
      A  = - A/(DEN + D1)

      K = KK + 2
      ZR = Z*R(1)
      A  = P(1,I)*P(1,J)*R(1)**K * (((B*ZR + A)*ZR + D1)/(DEN*H1) + D5)

! ... region (r1,inf) ...

      MM = MIN0(MX(I),MX(J))

      DO m = 3,MM,2
       A = A + P(m,I)*P(m,J)*R(m)**K
      END DO

      B = D0
      DO m = 2,MM,2
       B = B + P(m,I)*P(m,J)*R(m)**K
      END DO

      QUADR = H1*(A + B + B)

      END FUNCTION QUADR
