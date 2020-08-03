!=======================================================================
      Double precision function VKy (I1,J1,I2,J2,K) 
!=======================================================================
!                 k
!     Evaluates  V (i1, j1; i2, j2)  through the diff. equations
!
!
!     Calling sequence:
!
!                   VKy
!           --------------------
!          /       ||          |    
! Allocate_slater  YVK        YVAL   
!                   |
!                -----------------------
!               /     |     |      |    \
!             BZK  YVAL  VINTY   FACDZK  dgbsl
!              |                    |
!            -----------          dgbfa
!            |      |   \
!          FACDZK  YVAL  dgbsl
!            |
!          dgbfa
!
!-----------------------------------------------------------------------
      Use spline_orbitals, p => pbs
      Use spline_grid
      Use spline_slater
      Use spline_atomic

      Implicit none
      Integer, Intent(in) :: I1,J1,I2,J2,K

      if(ky.eq.-100) Call Allocate_slater

      Call YVK (j1,j2,k)

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
      fc = fc*fc1*fyk*grm*grm*grw

      vky = SUM(fc)

      vky = vky * fine

      CONTAINS


!=====================================================================
      Subroutine yvk(j1,j2,k)
!=====================================================================
!
!     Computes the spline solution of the differential equation
!
!           (d - (k+kk)/r) yk(r)= -(2k+kk)(1/r) zk(r)
!
!     with yk(ns) = zk(ns)
!
!     Calls:  bzk, yval, vinty, dgbtrs (lapack) or dgbfa (linpack)
!                               dgbtrf             dgbsl
!     on exit
!     -------
!     yk    spline expansion of the solution
!     fyk   values at the gaussian points * (1/r)
!---------------------------------------------------------------------

       USE spline_param; USE spline_galerkin

       IMPLICIT NONE

       INTEGER, INTENT(in) :: j1,j2,k

       REAL(8) :: fk, yns
       Integer(4) :: info 

!  ... calculation of Zk function ...

       Call BZK (j1,j2,k,1)
       yns = yk(ns)

!  ... set up the array of yk(i) = INTEGRAL '1/r zk(r) B_i(r)'

       fc = grw * fyk

       Call Vinty(fc,yk)

       fk = - (k+k+3)
       yk = fk * yk

!  ... set up and factor the differential operator ...

       if(kz.ne.-(k+3)) CALL FACDZK(ktx,-(k+3),dzk,ipvtz)
       kz = -(k+3)

!  ... boundary conditions:    yk(ns) = zk(ns)

        yk(ns) = yns

!  ... solve the matrix equation

!      CALL dgbsl(dzk,ktx,ns,ks-1,ks-1,ipvtz,yk,0) ! LINPACK

       Call DGBTRS ('N',ns,ks-1,ks-1,1,dzk,3*ks-2,ipvtz,yk,ns,info)
       if( info .ne. 0 ) Stop 'VKy: dgbtrs (LAPACK) failed'

! ... evaluates the function (1/r)yk(r) at all the gaussian points

       Call YVAL (0,0,-1,yk,fyk)

       END SUBROUTINE YVK

       END FUNCTION VKy
      
      
