!======================================================================
      Subroutine MVK_diff(k)
!======================================================================
!
!     Sets up matrix vk which stores the integral
!
!           <B_isp, (1/r^3) yk{jsp,jsp+jth-1}, B'_isp+ith-1>
!
!     where isp=1..n, ith=1..ks, and yk{...} is the solution
!     of the differential equation
!
!           (d - (k+3)/r) yk(r) = -(2k+3)/r * Zk
!
!     and zk{...} is the solution  of the differential equation
!
!           (d + k/r)zk(r) = r * B_jsp(r)*B_jsp+jth-1(r).
!
!-----------------------------------------------------------------------
!
!     Calling sequence:
!
!                mvk
!                 |
!          -----------------
!         /     ||     |    \
!      facdzk bspvyk  yval dinty
!        |       |
!      dgbfa   ----------
!             /     |    \
!            yval  vinty  dgbsl
!-----------------------------------------------------------------------
!     on entry
!     --------
!     k    multipole index
!
!     on exit
!     -------
!     rkb   a four dimension array vk(isp,jsp,ith,jth)
!           in symmetric upper-column storage mode for j-indeces
!           and non-symmetric storage mode for i-indeces
!=====================================================================

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. local variables

    INTEGER :: jsp, jth
    REAL(KIND=8),DIMENSION(ns) :: yk
    REAL(KIND=8),DIMENSION(ns,2*ks-1) :: vk, vkd
    REAL(KIND=8),DIMENSION(nv,ks) :: fyk
    REAL(KIND=8),DIMENSION(3*ks-2,ns) :: dyk,dzk
    INTEGER, DIMENSION(ns) :: ipvtz,ipvty

    ! .. check the need of calculations

    if(itype == 'aaa') Call allocate_integrals
    if(itype == 'vk ' .and. krk == k) Return

    ! .. factorize the diff. operators

    Call FACDZK  (3*ks-2,  k   ,dzk,ipvtz)
    Call FACDZK  (3*ks-2,-(k+3),dyk,ipvty)

    do jth = 1,ks
    do jsp = 1,ns

      ! .. consider element (jsp,jsp+jth-1) --> (j,j') at j<=j'

      if (jsp+jth-1 .gt. ns) Exit

      ! .. B-spline represantation for yk{jsp,jth}

      Call BSPVYK (jsp,jsp+jth-1)

      ! .. yk{jsp,jth}/r^3 at all the gaussian points (weighted)

      Call YVAL (0,1,-3,yk,fyk)

      ! .. integrates yk with B_isp and B'_isp+ith-1 for all elements

      Call DINTY (1,fyk,vkd)

      ! .. yk{isp,ith}/r^4 at all the gaussian points (weighted)

      Call YVAL (0,1,-4,yk,fyk)

      ! .. integrates yk with B_isp and B_isp+ith-1 for all elements

      Call DINTY (0,fyk,vk)

      rkb(:,jsp,:,jth) = vkd(:,:) - vk(:,:)

    end do
    end do

    rkb = rkb * fine

    krk = k
    itype = 'vk '

    CONTAINS


!======================================================================
      Subroutine BSPVYK (jf,js)
!======================================================================
!
!     Computes the B-spline expansion for the solution yk
!
!     (d-(k+3)/r) yk(r)= -(2k+3)/r zk,  where
!
!     (d + k/r)zk(r) = r * B_is(r)*B_if(r), at is>jf
!
!     on entry
!     --------
!     jf   index of the first B-spline
!     js   index of the second B-spline
!
!     on exit
!     -------
!     yk   coefficients of the B-spline expansion of the solution
!
!     Calls: YVAL, VINTY,  dgbtrs (LAPPACK)  or  dgbsl (LINPACK)
!---------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jf,js

    ! .. local variables

    INTEGER :: jv,jv1,jv2, ith, info
    REAL(8), DIMENSION(ks) :: v
    REAL(8) :: yns

!----------------------------------------------------------------------
!                                                          zk-function:

!                                             ks     ks
! .. the interval stretch of the two splines B     B
!                                             jf     js

       jv1 = max0( 1,js-ks+1)
       jv2 = min0(jf,ns-ks+1)

!  ... set up the array yk = integral of three B-splines * r ...

       yk = 0.d0

       do jv = jv1,jv2
        v(:) =  grw(jv,:)*gr(jv,:)*bsp(jv,:,jf-jv+1)*bsp(jv,:,js-jv+1)
        do ith = 1,ks
         yk(jv+ith-1) = yk(jv+ith-1) + SUM(v(:)*bsp(jv,:,ith))
        end do
       end do

!  ... boundary condition ...

       yk(1) = 0.d0

!  ... solve the matrix equation ...

!      CALL dgbsl(dzk,3*ks-2,ns,ks-1,ks-1,ipvtz,yk,0) ! LINPACK

       Call DGBTRS ('N',ns,ks-1,ks-1,1,dzk,3*ks-2,ipvtz,yk,ns,info)
       if( info .ne. 0 ) Stop 'BSPVYK: dgbtrs failed (LAPACK)'

!----------------------------------------------------------------------
!                                                          yk-function:
       yns = yk(ns)

!  ... set up the array of yk(i) = INTEGRAL '1/r zk(r) B_i(r)'

!  ... (1/r)zk(r) at all the gaussian points (weighted)

       Call YVAL (0,1,-1,yk,fyk)

       Call Vinty(fyk,yk)

       yk = - (k+k+3) * yk

!  ... boundary conditions:

       yk(1)  = 0.d0
       yk(ns) = yns

!  ... solve the matrix equation

!      CALL dgbsl(dyk,3*ks-2,ns,ks-1,ks-1,ipvty,yk,0)  ! LINPACK

       Call DGBTRS ('N',ns,ks-1,ks-1,1,dyk,3*ks-2,ipvty,yk,ns,info)
       if(info.ne.0) Stop 'BSPVYK: dgbtrs (LAPACK) failed'

       End Subroutine BSPVYK

      End Subroutine MVK_diff

