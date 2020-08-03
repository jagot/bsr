!====================================================================
      REAL(8) FUNCTION RK (i1,i2,i3,i4,k)
!====================================================================
!
!     Evaluate RK integral:
!
!     Int(r1,0,inf) Int(r2,0,inf)
!
!     P_i1(r1)*P_i2(r2) [r<^k / r>^(k+1)] P_i3(r1)*P_i4(r2)
!
!     with re-definotion on the inclusion relativistic and
!     mass corrections
!
!     r< = min(r1,r2);   r> = max(r1,r2)
!==================================================================

      USE RADIAL

      IMPLICIT NONE
      Integer(4), Intent(in) :: i1,i2,i3,i4,k
      Real(8), External :: QUADS, XK
      Real(8) :: RK1, RK2

      Call DZK(i1,i3,k);      RK1 = QUADS(i2,i4,1)

      Call DZK(i2,i4,k);      RK2 = QUADS(i1,i3,1)

      RK = RK1 + RK2

      if(rel) RK = RK + (k+k+1) * FINE * XK(i1,i2,i3,i4)
      IF(MASS.GT.0) Call RK_MASS(i1,i2,i3,i4,k,RK)

      END function RK



