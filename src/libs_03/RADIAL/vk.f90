!=====================================================================
      REAL(8) FUNCTION VK (I,J,II,JJ,K)
!=====================================================================
!
!     Evaluate VK integral:
!
!     alfa^2/4
!                                       k+1
!       Int(r;0-inf) P_i  P'_ii  1/r^2  Z   (P_j,P_jj)
!                                       k
!     + Int(r;0-inf) P_j  P_jj  1/r^2  Z (P_i,P'_ii)
!
!     where P'(r) = (d/dr - 1/r) P(r)
!------------------------------------------------------------------

      USE RADIAL

      IMPLICIT NONE

      Integer(4), Intent(in) :: I,II,J,JJ,K
      Integer(4) :: ip
      Real(8) :: VK1,VK2
      Real(8), External :: QUADS

      IP = nrf+1;  Call P_derive(II,IP)

      Call DZK(J,JJ,k+1);  VK1 = QUADS(I,IP,2)

      Call DZK(I,IP,k);    VK2 = QUADS(J,JJ,2)

      VK = VK1 + VK2

      VK = VK * FINE

      END FUNCTION VK


