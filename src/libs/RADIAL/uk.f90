!==================================================================
      REAL(8)  FUNCTION UK (I,J,II,JJ,K)
!==================================================================
!
!     Evaluates UK integral:
!
!     alfa^2/4/(2k+1)
!                                              k
!      (k-1) * Int(r;0-inf) P_i  P_ii  1/r^2  Z (P_j,P'_jj)
!                                               k-1
!     -(k+2) * Int(r;0-inf) P_j  P'_jj  1/r^2  Z   (P_i,P_ii)
!
!     where P'(r) = (d/dr - 1/r) P(r)
!------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL(8) (A-H,O-Z)

      Integer(4), Intent(in) :: I,II,J,JJ,K

      Real(8), External :: QUADS

      JP = nrf+2;     Call P_derive(JJ,JP)

      Call DZK(J,JP,k);   UK1 = QUADS(I,II,2)

      Call DZK(I,II,k-1); UK2 = QUADS(J,JP,2)

      UK = UK1 + UK2

      UK = UK / (k+k+1) * FINE

      END FUNCTION UK


