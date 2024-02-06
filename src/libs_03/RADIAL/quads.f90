!=======================================================================
      REAL(8) FUNCTION QUADS (I,J,KK)
!=======================================================================
!
!                                       kk
!       Evaluates the integral of  (1/r)   ZK(r) P (r) P (r)  with
!                                                 i     j
!       respect to r.
!
!-----------------------------------------------------------------------

      USE RADIAL, MX => mro

      IMPLICIT NONE
      Integer(4), Intent(in) :: i,j,kk
      Integer(4) :: k,m,mm
      Real(8) :: A,B,DEN

! ... (0,r1) region ...

      A = azk + aexp(i) + aexp(j)
      B = bzk + bexp(i) + bexp(j) +  &
          azk*aexp(i) + azk*aexp(j) + aexp(i)*aexp(j)
      m = mzk + mexp(i) + mexp(j) - kk

      DEN = m + 1
      B =  A*A/(DEN+D1) - D2*B/(DEN+D2)
      A = -A/(DEN+D1)

      B = D1 + Z*R(1)*(A + Z*R(1)*B)
      K = 2 - KK
      A  = ZK(1)*P(1,I)*P(1,J)*R(1)**K * (B/(DEN*H1) + D5)

! ... (r1,inf) region ...

      MM = MIN0(MX(I),MX(J))

      Do M = 3,MM,2
       A  = A  + ZK(M)*P(M,I)*P(M,J)*R(M)**K
      End do

      B = D0
      Do M = 2,MM,2
       B = B + ZK(M)*P(M,I)*P(M,J)*R(M)**K
      End do

      QUADS = H1*(A + B + B)

      END FUNCTION QUADS
