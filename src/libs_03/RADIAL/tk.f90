!==================================================================
      REAL(8) FUNCTION TK (I,J,II,JJ,K)
!==================================================================
!
!     Evaluate TK integral:
!
!     alfa^2/4/(2k+1) Int(r1,0,inf) Int(r2,0,inf)
!
!     P_i(r1)*P_ii(r2) [r<^k / r>^(k+1)] P'_j(r1)*P'_jj(r2)
!
!     P'(r) = (d/dr - 1/r) P(r)
!------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL (8) (A-H,O-Z)

      Integer(4), Intent(in) :: I,II,J,JJ,K

      Real(8), External :: QUADS

      IP = nrf+1
      Call P_derive(II,IP)

      JP = nrf+2
      Call P_derive(JJ,JP)

      Call DZK(I,IP,k)
      TK1 = QUADS(J,JP,1)

      Call DZK(J,JP,k)
      TK2 = QUADS(I,IP,1)

      TK = TK1 + TK2

      TK = TK / (k+k+1) * FINE

      END FUNCTION TK
