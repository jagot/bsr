!----------------------------------------------------------------------
      REAL(8) FUNCTION NK (I,J,II,JJ,K)
!----------------------------------------------------------------------
!
!     returns value of NK-integral:
!
!     alfa^2/4  Int(r1;0-inf) Int(r2;0-r1)
!
!     P_i(r1)*P_j(r2) [r<^k / r>^(k+3)] P_ii(r1)*P_jj(r2)
!
!----------------------------------------------------------------------

      USE RADIAL

      Implicit real(8) (A-H,O-Z)

      Integer(4), Intent(in) :: i,ii,j,jj,k

      Real(8), External :: QUADS

      CALL DZK(J,JJ,K)

      NK = QUADS(I,II,3)*FINE

      END FUNCTION NK



