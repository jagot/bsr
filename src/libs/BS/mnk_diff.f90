!======================================================================
      Subroutine MNK_diff (k)
!======================================================================
!
!     Sets up matrix nk which stores the integral
!
!           <B_isp, (1/r^3)zk{jsp,jsp+jth-1}, B_isp+ith-1>
!
!     where isp=1..n, ith=1..ks, and zk{...} is the solution
!     of the differential equation
!
!           (d + k/r)zk(r) = B_jsp(r)*B_jsp+jth-1(r).
!
!   Calling sequence:
!
!                 mnk_diff
!             ------------------
!            /     ||     |    |
!         facdzk bspzk  yval minty
!           |      |
!         dgbfa  dgbsl
!
!
!     on entry
!     --------
!     k    mupltipole index
!
!     on exit
!     -------
!     rbk    a four dimension array nk(isp,jsp,ith,jth)
!            in the symmetric upper-column storage mode
!---------------------------------------------------------------------

      USE spline_param
      USE spline_grid
      USE spline_galerkin
      USE spline_integrals
      USE spline_atomic

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: k

      ! .. local variables

      INTEGER :: jsp, jth
      REAL(KIND=8), DIMENSION(ns) :: yk
      REAL(KIND=8), DIMENSION(nv,ks) :: fyk
      REAL(KIND=8), DIMENSION(ns,ks) :: nki
      REAL(KIND=8), DIMENSION(3*ks-2,ns) :: dzk
      INTEGER, DIMENSION(ns) :: ipvtz

      ! .. check the need of calculations

      if(itype == 'aaa') Call allocate_integrals
      if(itype == 'nk ' .and. krk == k) Return

      ! .. factorize the Zk operator:

      Call FACDZK (3*ks-2,k,dzk,ipvtz)

      do jth = 1,ks
      do jsp = 1,ns

        ! ... consider element (jsp,jsp+jth-1) --> (j,j') at j<j'

        if (jsp+jth-1.gt.ns) Exit

        ! ... the coefficients of the B-spline basis for zk{jsp,jth}

        Call BSPZK (jsp,jsp+jth-1)

        ! ... zk{jsp,jth}/r^3 at all the gaussian points (weighted)

        Call YVAL (0,1,-3,yk,fyk)

        ! ... integrates yk with B_isp and B_isp+ith-1 --> superdiagonals

        Call MINTY (1,fyk,nki)

        rkb(1:ns,jsp,1:ks,jth) = nki(1:ns,1:ks)

      end do
      end do

      rkb = rkb * fine

      krk = k
      itype = 'nk '


      CONTAINS

!======================================================================
      Subroutine BSPZK (jf,js)
!======================================================================
!
!     Computes the B-spline expansion for the solution zk
!     of the equation using the spline Galerking method:
!
!       (d + k/r) zk(r)= B_jf*B_js(r),  at is>jf
!
!     Calls:  dgbtrs  (LAPACK) or  dgbsl (LINPACK)
!
!     on entry
!     --------
!     jf   index of the first B-spline
!     js   index of the second B-spline
!
!     on exit
!     -------
!     yk   the coefficients of the B-spline expansion of the solution
!---------------------------------------------------------------------

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: jf,js

       ! .. local variables

       INTEGER :: jv,jv1,jv2, ith, info
       REAL(KIND=8), DIMENSION(ks) :: v


!                                               ks     ks
! ...  the interval stretch of the two splines B     B
!                                               jf     js

       jv1 = max0( 1,js-ks+1)
       jv2 = min0(jf,ns-ks+1)

!  ... set up the array yk = integral of three B-splines  ...

       yk = 0.d0

       do jv = jv1,jv2
        v(:) =  grw(jv,:)*bsp(jv,:,jf-jv+1)*bsp(jv,:,js-jv+1)
        do ith = 1,ks
         yk(jv+ith-1) = yk(jv+ith-1) + SUM(v(:)*bsp(jv,:,ith))
        end do
       end do

! ... boundary condition ...

      yk(1) = 0.d0
      if(k.eq.-1) yk(2)=0.d0

! ... solve the matrix equation ...

!     CALL dgbsl(dzk,3*ks-2,ns,ks-1,ks-1,ipvtz,yk,0) ! LINPACK

      Call DGBTRS ('N',ns,ks-1,ks-1,1,dzk,3*ks-2,ipvtz,yk,ns,info)
      if(info.ne.0) Stop 'MNK_diff: dgbtrs (LAPACK) failed'

      End Subroutine BSPZK

    End Subroutine MNK_diff

