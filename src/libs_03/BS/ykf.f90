!=====================================================================
      Subroutine ykf (j1,j2,k)
!=====================================================================
!     Computes the spline solution of the differential equation
!
!           (d^2-k(k+1)/r^2)yk(r)= -(2k+1)(1/r) P_1(r)P_2(r)
!
!     on exit   (in module spline_slater)
!     -------
!     yk    spline expansion of the solution
!     fyk   values at the gaussian points * (1/r)
!
!   Calling sequence:
!
!              YKF
!       ------------------------
!      /     |     |       |    \
!   FACDYK  YVAL  VINTY  QUADR  dgbtrs(LAPACK) or dgbsl(LINPACK)
!     |                         
!   dgbfa(LALINPACK)
!---------------------------------------------------------------------
      USE spline_param
      USE spline_grid
      USE spline_slater
      USE spline_orbitals, p => pbs
      USE spline_atomic

      IMPLICIT NONE

      INTEGER, INTENT(in) :: j1,j2,k

! ..  local variables

      REAL(8) :: c
      REAL(8), EXTERNAL :: QUADR
      Integer(4) :: info

!  ... create the numerical values of the orbitals on the grid

       if (iy1 .ne. j1) then
         Call YVAL (0,0,0,p(1,j1),fy1)
         iy1 = j1
       end if

       if (iy2 .ne. j2) then
         Call YVAL (0,0,0,p(1,j2),fy2)
         iy2 = j2
       end if

!  ... set up the array of yk(i) = INTEGRAL [(1/r)fc1(r)fc2(r)*B_i(r)]

       fc = grw*grm*fy1*fy2

       Call Vinty(fc,yk)

       c = - (2*k+1)
       yk = c*yk

!  ... boundary conditions:

       yk(1)  = 0.d0
       yk(ns) = 0.d0
       if( k.eq.0 ) yk(ns) = QUADR(j1,j2,0)

!  ... set up and factor the differential operator ...

       if(k.ne.ky) Call FACDYK (ktx,k,ipvtd,dyk)
       ky = k

!  ... solve the matrix equation

!      CALL dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)

       Call DGBTRS ('N',ns,ks-1,ks-1,1,dyk,ktx,ipvtd,yk,ns,info)
       if( info .ne. 0 ) Stop 'YKF: dgbtrs failed (LAPACK)'

!  ... evaluates function (1/r)yk(r) at the gaussian points

       Call YVAL (0,0,-1,yk,fyk)

!  ... add a relativistic correction if irel > 0

        if (rel) then
          c = (k+k+1)*fine
          fyk = fyk + c*grm*grm*fy1*fy2
        end if

       END SUBROUTINE ykf
