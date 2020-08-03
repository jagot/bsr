!====================================================================
      Subroutine RK_MASS (i1,i2,i3,i4,k,RK)
!====================================================================
!
!     Mass-correction of RK integral
!
!==================================================================

      USE RADIAL

      IMPLICIT NONE
      Integer(4), Intent(in) :: i1,i2,i3,i4,k 
      Real(8), Intent(inout) :: RK
      Real(8), External :: GRAD, QUADR
     
      IF (MASS .GT. 0) THEN
         IF (MASS .EQ. 1) THEN
            IF (K .EQ. 1) RK = RK - RMASS*GRAD(i1,i3)*GRAD(i2,i4)
         ELSE
            RK = RK*(D1 + RMASS/D2)
            IF (K .EQ. 1) RK = RK + Z*RMASS/D2 *            &
                          (QUADR(i1,i3, 1)*QUADR(i2,i4,-2)+ &
                           QUADR(i1,i3,-2)*QUADR(i2,i4, 1))
         END IF
      END IF

      End Subroutine RK_MASS