!=======================================================================
      Double precision function NKy (I1,J1,I2,J2,K)
!=======================================================================
!                 k
!     Evaluates  N (i1, j1; i2, j2)  through the gaussian points
!
!
!     Calling sequence:
!
!                  NKy
!                   |
!           ----------------
!          /           |    \
!   allocate_slater   BZK   YVAL
!                      |
!                  -------------
!                 /     |       \
!              FACDZK  YVAL   dgbsl
!                |
!              dgbfa
!-----------------------------------------------------------------------
      USE spline_orbitals, p => pbs
      USE spline_grid
      USE spline_slater
      USE spline_atomic

      IMPLICIT NONE
      INTEGER, INTENT(in) :: I1,J1,I2,J2,K

      if(ky.eq.-100) Call Allocate_slater

      Call BZK (j1,j2,k,0)

      if (i1 .ne. ic1) then
        Call YVAL (0,0,0,p(1,i1),fc1)
        ic1 = i1
      end if

      if (i2 .ne. ic2) then
        Call YVAL (0,0,0,p(1,i2),fc2)
        ic2 = i2
      end if

      fc = grm*grm*grw*fc1*fc2*fyk
      nky = SUM(fc)
      nky = nky * fine

      END FUNCTION NKy

