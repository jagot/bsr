!=======================================================================
      Double precision function TKy (I1,J1,I2,J2,K)
!=======================================================================
!                 k
!     Evaluates  T (i1, j1; i2, j2)  through diff. equations
!
!
!     Calling sequence:
!
!                  TKy
!           -------------------
!          /        ||        |  
!  Allocate_slater  YTK      YVAL
!                    |
!                 ------------------------
!                /     |     |       |    \
!             FACDYK  YVAL  VINTY  QUADR  dgbsl(LINPACK)
!               |
!             dgbfa(LINPACK)
!-----------------------------------------------------------------------
      USE spline_orbitals, p => pbs, L => lbs
      USE spline_grid
      USE spline_slater
      USE spline_atomic

      IMPLICIT NONE
      INTEGER, INTENT(in) :: I1,J1,I2,J2,K

      if(ky.eq.-100) Call Allocate_slater

      Call Ytk(j1,j2,k)

      if (i1 .ne. ic1) then
        Call YVAL (0,0,0,p(1,i1),fc1)
        ic1 = i1
      end if

      if (i2 .ne. ic2) then
        Call YVAL (0,0,0,p(1,i2),fc2)
        ic2 = i2
      end if

      Call YVAL (1,0,0,p(1,i2),fc)

      fc = fc - grm*fc2
      fc = fc*fc1*fyk*grw

      tky = SUM(fc)

      tky = tky * fine / (k+k+1)

      CONTAINS



!=====================================================================
      Subroutine ytk(j1,j2,k)
!=====================================================================
!
!     Computes the spline solution of the differential equation
!
!           (d^2-k(k+1)/r^2)yk(r)= -(2k+1)(1/r) P_1(r)P'_2(r)
!
!     Calls:  yval,  vinty,  dgbtrs, dgbtrf  (lapack)
!                                   
!     on exit
!     -------
!     yk    spline expansion of the solution
!     fyk   values at the gaussian points * (1/r)
!---------------------------------------------------------------------
       USE spline_param;  USE spline_galerkin

       IMPLICIT NONE
       INTEGER, INTENT(in) :: j1,j2,k
       INTEGER :: jj,jp,info
       REAL(8) :: c
       REAL(8), EXTERNAL :: QUADR

! ... create the numerical values of the orbitals on the grid

      if (iy1 .ne. j1) then
        Call YVAL (0,0,0,p(1,j1),fy1)
        iy1 = j1
      end if

      if (iy2 .ne. j2) then
        Call YVAL (0,0,0,p(1,j2),fy2)
        iy2 = j2
      end if

! ... set up the array of yk(i) = INTEGRAL 'r^m*fc1(r)*fc2'(r)*B_i(r)'

      Call YVAL (1,0,0,p(1,j2),fc)
      fc = fc - grm*fy2
      fc = fc*fy1*grw*grm

      Call Vinty(fc,yk)

      c = - (2*k+1)
      yk = c*yk

! ... boundary conditions:

      yk(1)   = 0.d0
      yk(ns)  = 0.d0   !  Assume that p(ns)=0.d0 for input w.f.'s

      if( k.eq.0 .and. l(j1).eq.l(j2) ) then  ! INT( p_j1(r)*p'_j2(r) dr)
        c =  0.d0
        do jp = 2,ks
          do jj = jp+1,ns-jp+1
            c =  c  +    db1(jj,ks-jp+1) *   &
               ( p(jj,j2)*p(jj-jp+1,j1) - p(jj,j1)*p(jj-jp+1,j2) )
          end do
        end do
        yk(ns) = c - QUADR(j1,j2,-1)
      end if

! ... set up and factor the differential operator ...

      if(k.ne.ky) Call FACDYK (ktx,k,ipvtd,dyk)
      ky = k

! ... solve the matrix equation

      Call DGBTRS ('N',ns,ks-1,ks-1,1,dyk,3*ks-2,ipvtd,yk,ns,info)
      if( info .ne. 0 ) Stop 'TKy: dgbtrs (LAPACK) failed'

! ... evaluates the function (1/r)yk(r) at all the gaussian points

      Call YVAL (0,0,-1,yk,fyk)

      END SUBROUTINE Ytk

      END FUNCTION TKy

