!-----------------------------------------------------------------------
      Double precision function RKy (I1,J1,I2,J2,K)
!-----------------------------------------------------------------------
!                 k
!     Evaluates  R (i1, j1; i2, j2)  through diff. equations
!
!
!     Calling sequence:
!
!                  RKy
!           ---------------------
!          /        |            \
! Allocate_slater  YKF          YVAL
!                   |
!                 ------------------------
!                /     |     |       |    \
!             FACDYK  YVAL  VINTY  QUADR  dgbsl(LINPACK)
!               |
!             dgbfa(LINPACK)
!-----------------------------------------------------------------------

      USE spline_grid
      USE spline_orbitals, p => pbs
      USE spline_slater
      USE spline_atomic

      INTEGER, INTENT(in) :: i1,j1,i2,j2,k

      if(ky.lt.-1) Call Allocate_slater

      Call YKF (j1,j2,k)

!      if (i1 .ne. ic1) then
        Call YVAL (0,0,0,p(1,i1),fc1)
!        ic1 = i1
!      end if

!      if (i2 .ne. ic2) then
        Call YVAL (0,0,0,p(1,i2),fc2)
!        ic2 = i2
!      end if

      fc = grw*fc1*fc2*fyk

      rky = SUM(fc)

    END function RKy
