!=====================================================================
    SUBROUTINE mrk_diff(k)
!=====================================================================
!
!   Sets up matrix rk which stores the integral
!                1
!       <B_isp  ---  yk{jsp,jth} B_isp+ith-1>
!                r
!   where isp=1..n, ith=1..ks, and yk{jsp,jth} is the solution
!   of the differential equation:
!
!        d^2    k(k+1)         (2k+1)
!       (--- - ------)yk(r)= - ------ B_jsp(r) B_jsp+jth-1(r).
!        dr^2    r^2              r
!
!-----------------------------------------------------------------------
!
!     Calling sequence:
!
!                mrk
!                 |
!         ------------------
!        /     ||     |    |
!     facdyk bspyk  yval minty
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
    REAL(KIND=8), DIMENSION(ns,ks) :: rki
    REAL(KIND=8), DIMENSION(3*ks-2,ns) :: dyk
    REAL(KIND=8) :: c
    INTEGER, DIMENSION(ns) :: ipvtd

!---------------------------------------------------------------------
    ! .. check the need of calculations

    if(itype == 'aaa') Call allocate_integrals
    if(itype == 'rk ' .and. krk == k) Return

!---------------------------------------------------------------------
    ! .. set up and factorize the differential operator

    CALL facdyk(3*ks-2,k,ipvtd,dyk)

    DO ith = 1,ks
      DO isp = 1,ns

        jsp = isp+ith-1

        ! ... consider element (isp,jsp) --> (j,j') at j<j'

        if (jsp > ns) EXIT

        ! .. B-spline represantation for yk{isp,ith}

        CALL bspyk(isp,jsp)

        ! .. evaluates the function yk{isp,ith} at all the gaussian points

        CALL yval(0,1,-1,yk,fyk)

         IF (rel) THEN
           c = (2*k+1)*fine
           jv1 = max0(  1,jsp-ks+1)
           jv2 = min0(isp, ns-ks+1)
           do jv = jv1,jv2
              fyk(jv,:) = fyk(jv,:) + c*grw(jv,:)*grm(jv,:)*grm(jv,:)*  &
                          bsp(jv,:,isp-jv+1)*bsp(jv,:,jsp-jv+1)
           end do
         END IF

        ! .. integrates yk with B_jsp and B_jsp+jth-1; superdiagonals

        CALL minty(1,fyk,rki)

        rkb(1:ns,isp,1:ks,ith) = rki(1:ns,1:ks)

      END DO
    END DO

    itype='rk '
    krk=k

    CONTAINS


!======================================================================
    SUBROUTINE bspyk(jf,js)
!======================================================================
!
!   Computes the B-spline expansion for the solution yk
!   of the equation using the spline Galerking method:
!
!       (d^2-k(k+1)/r^2) yk(r) = -(2k+1) B_jf(r) B_js(r)/r
!
!   with boundary conditions:
!       at r=0, yk(r) = 0;
!       at r=rmax, yk(r) = const(=ykbc), if k=0,
!                  dy/dr + (k/r)y = B_jf(r) B_js(r) if k>0.
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
       Integer(4) :: info

!                                               ks     ks
! ...  the interval stretch of the two splines B     B
!                                               jf     js

       jv1 = max0( 1,js-ks+1)
       jv2 = min0(jf,ns-ks+1)

!  ... set up the array yk = integral of three B-splines / r ...

       yk = 0.d0

       do jv = jv1,jv2
        v(:) =  grw(jv,:)*grm(jv,:)*bsp(jv,:,jf-jv+1)*bsp(jv,:,js-jv+1)
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
         yk(ns) = sb(js,jf-js+ks)
       else
         yk(ns) = 0.0d0
         if(js==ns .and. jf==ns) yk(ns) = 1.d0
       end if

       ! .. solve the matrix equation

!      CALL dgbsl(dyk,3*ks-2,ns,ks-1,ks-1,ipvtd,yk,0) ! LINPACK

      Call DGBTRS ('N',ns,ks-1,ks-1,1,dyk,3*ks-2,ipvtd,yk,ns,info)
      if(info.ne.0) Stop 'mrk_diff: dgbtrs (LAPACK) failed'

     END SUBROUTINE bspyk

  END SUBROUTINE mrk_diff
