!====================================================================
      REAL(8) FUNCTION XK (i1,i2,i3,i4)
!====================================================================
!
!     Evaluates integral:
!
!     Int(r;0:inf)  P_i1(r)*P_i2(r)*P_i3(r)*P_i4(r) / r^2
!
!--------------------------------------------------------------------

      USE RADIAL, MX => mro

      IMPLICIT NONE
      Integer(4), Intent(in) :: i1,i2,i3,i4
      Integer(4) :: i,m
      Real(8) :: A,B,ZR,DEN

! ... (0:r1) interval:

      m = mexp(i1) + mexp(i2) + mexp(i3) + mexp(i4) - 2
      a = aexp(i1) + aexp(i2) + aexp(i3) + aexp(i4)
      b = bexp(i1) + bexp(i2) + bexp(i3) + bexp(i4) + &
          aexp(i1) *(aexp(i2) + aexp(i3) + aexp(i4)) + &
          aexp(i2) *(aexp(i3) + aexp(i4)) + aexp(i3)*aexp(i4)

      DEN = m + 1
      B = a*a/(DEN+1) - 2*b/(DEN+3)
      A = - a/(DEN+1)

      ZR = Z*R(1)
      A  = P(1,i1)*P(1,i2)*P(1,i3)*P(1,i4)*R(1) * &
           ((D1 + ZR*(A + B*ZR))/(DEN*H1) + D5)

! ... (r1:inf) interval:

      M = MIN0(MX(i1),MX(i2),MX(i3),MX(i4))

      DO i = 3,M,2
       A = A + P(i,i1)*P(i,i2)*P(i,i3)*P(i,i4)*R(i)
      END DO

      B = D0
      DO i = 2,M,2
       B = B + P(i,i1)*P(i,i2)*P(i,i3)*P(i,i4)*R(i)
      END DO

      XK = H1*(A + B + B)

      END FUNCTION XK

