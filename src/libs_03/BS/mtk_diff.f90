!=====================================================================
    SUBROUTINE mtk_diff(k)
!=====================================================================
!
!   Sets up matrix tk which stores the integral
!                 1
!       <B_isp  ---  yk{jsp,jth} B'_isp+ith-1>
!                 r
!   where isp=1..n, ith=1..ks, and yk{jsp,jth} is the solution
!   of the differential equation:
!
!        d^2    k(k+1)         (2k+1)
!       (--- - ------)yk(r)= - ------ B_jsp(r) B'_jsp+jth-1(r).
!        dr^2    r^2              r
!
!-----------------------------------------------------------------------
!
!     Calling sequence:
!
!                mtk
!                 |
!         ------------------
!        /     ||     |    |
!     facdyk bsptyk  yval dinty
!       |      |
!     dgbfa  dgbsl
!
!---------------------------------------------------------------------
!     on entry
!     --------
!     k    mupltipole index
!
!     on exit
!     -------
!     rbk    a four dimension array rk(isp,jsp,ith,jth)
!            in the symmetric upper-column storage mode
!---------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(in) :: k

    ! .. local variables

    INTEGER ::  isp, ith, jsp, jth, jv, jv1,jv2
    REAL(KIND=8), DIMENSION(ns) :: yk
    REAL(KIND=8), DIMENSION(nv,ks) :: fyk
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: tk, tkd
    REAL(KIND=8), DIMENSION(3*ks-2,ns) :: dyk
    REAL(KIND=8) :: c
    INTEGER, DIMENSION(ns) :: ipvtd

!---------------------------------------------------------------------
    ! .. check the need of calculations

    if(itype == 'aaa') Call allocate_integrals
    if(itype == 'tk ' .and. krk == k) Return

!---------------------------------------------------------------------
    ! .. set up and factorize the differential operator

    CALL facdyk(3*ks-2,k,ipvtd,dyk)

    DO ith = 1,2*ks-1
      DO isp = 1,ns

        jsp = isp+ith-ks

        ! ... consider element (isp,jsp) --> (j,j') at j<j'

        if (jsp > ns) EXIT

        ! .. B-spline represantation for yk{isp,ith}

        CALL bsptyk(isp,jsp)

        ! .. yk{jsp,jth}/r at all the gaussian points (weighted)

        Call YVAL (0,1,-1,yk,fyk)

        ! .. integrates yk with B_isp and B'_isp+ith-1 for all elements

        Call DINTY (1,fyk,tkd)

        ! .. yk{isp,ith}/r^2 at all the gaussian points (weighted)

        Call YVAL (0,1,-2,yk,fyk)

        ! .. integrates yk with B_isp and B_isp+ith-1 for all elements

        Call DINTY (0,fyk,tk)

        rkb(:,isp,:,ith) = tkd(:,:) - tk(:,:)

      END DO
    END DO

    c = fine / (k+k+1)
    rkb = rkb * c

    itype='tk '
    krk=k

    CONTAINS


!======================================================================
    SUBROUTINE bsptyk(jf,js)
!======================================================================
!
!   Computes the B-spline expansion for the solution yk
!   of the equation using the spline Galerking method:
!
!       (d^2-k(k+1)/r^2) yk(r) = -(2k+1) B_jf(r) B'_js(r)/r
!
!   with boundary conditions:
!       at r=0, yk(r) = 0;
!       at r=rmax, yk(r) = const, if k=0,
!                  dy/dr + (k/r)y = B_jf(r) B'_js(r) if k>0.
!
!   SUBROUTINES called:    dgbsl (LINPACK)
!
!---------------------------------------------------------------------
!
!   on entry
!   --------
!       jf    index of the first B-spline
!       js    index of the second B-spline
!
!   on exit
!   -------
!       yk    the coefficients of the B-spline expansion of yk(r)
!
!---------------------------------------------------------------------

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: jf,js

       ! .. local variables

       REAL(8), DIMENSION(ks) :: v
       INTEGER(4) :: info

!                                               ks     ks
! ...  the interval stretch of the two splines B     B
!                                               jf     js

       jv1 = max0( 1,jf-ks+1,js-ks+1)
       jv2 = min0(jf,js,ns-ks+1)

!  ... set up the array yk  ...

       yk = 0.d0

       do jv = jv1,jv2
        v(:) =  grw(jv,:)*grm(jv,:)*bsp(jv,:,jf-jv+1)*     &
                (bspd(jv,:,js-jv+1,1) - grm(jv,:)*bsp(jv,:,js-jv+1))
        do jth = 1,ks
         yk(jv+jth-1) = yk(jv+jth-1) + SUM(v(:)*bsp(jv,:,jth))
        end do
       end do

       c = -(2*k+1)
       yk = c*yk

       ! .. the boundary condition at the origin

       yk(1) = 0.d0

       ! .. the boundary condition at rmax

       if (k == 0) then
         if(js.gt.jf) then
          yk(ns) = db1(js,jf-js+ks) - rm1(js,jf-js+ks)
         else
          yk(ns) = -db1(jf,js-jf+ks) - rm1(jf,js-jf+ks)
         end if
       else
         if(jf.ne.ns) then
          yk(ns) = 0.d0
         elseif(js.eq.ns) then
          yk(ns) = bspd(nv+1,1,ks,1) - 1.d0/t(ns+1)
         elseif(js.eq.ns-1) then
          yk(ns) = bspd(nv+1,1,ks-1,1)
         end if
       end if

       ! .. solve the matrix equation

!      CALL dgbsl(dyk,3*ks-2,ns,ks-1,ks-1,ipvtd,yk,0) ! LINPACK

       Call DGBTRS ('N',ns,ks-1,ks-1,1,dyk,3*ks-2,ipvtd,yk,ns,info)
       if(info.ne.0) Stop 'mtk_diff: dgbtrs (LAPACK) failed'

     END SUBROUTINE bsptyk

  END SUBROUTINE mtk_diff
