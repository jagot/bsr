!=========================================================================
    SUBROUTINE mk_moments(k)
!=========================================================================
!
!   Defines moments for Mk-integrals in the B-spline cells
!
!   Calling sequence:          mk_moments
!                              ----------
!                               /    \\
!                           moments mk_pdiag
!                                     ||
!                                   mk_triang
!                                    /   \
!                                 gauss  vbsplvb
!
!-------------------------------------------------------------------------
!
!   on entry    k  -  multipole index
!   --------
!
!   on exit     rkd1,rkd2,rkd - off-diagonal and diagonal moments
!   -------                     (in module spline_moments)
!
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_moments
    USE spline_integrals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_moments
    if(mtype == 'mk ' .and. kmk == k) Return

    CALL moments(  k   , rkd1,'s','b')
    CALL moments(-(k+3), rkd2,'s','b')
    CALL mk_pdiag

    mtype = 'mk '
    kmk = k

    CONTAINS

!===================================================================
    SUBROUTINE mk_pdiag
!===================================================================
!
!   Controls the scaling propeties for diagonal B-spline Nk-interals
!
!   Calls:  mk_triang
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv,ik
    REAL(KIND=8) :: hp

    ! .. non-exponential grid

    if(me.eq.0) then
     DO iv=1,nv
      CALL mk_triang(iv)
     END DO
     Return
    end if

    ! .. the first equal step region

    DO iv=1,ml+ks-1
      CALL mk_triang(iv)
    END DO

    ! .. the exponential region - using scaling law

    hp=h+1.d0
    ik = ks*(ks+1)/2
    DO iv=ml+ks,ml+me-ks+2
     rkd(1:ik,1:ik,iv) = rkd(1:ik,1:ik,iv-1) / hp
    END DO

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL mk_triang(iv)
    END DO

    END SUBROUTINE mk_pdiag


!========================================================================
    SUBROUTINE mk_triang(iv)
!========================================================================
!
!   Returns the two-dimensional array of B-splin Mk-integrals
!         <B_i B_j| r2^k/r1^(k+3) |B_i' B_j'>
!   over the given triangle diagonal cell
!
!   Calls:   gauss, vbsplvd
!
!   On entry   iv  -  the index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of B-spline Mk-integrals for given
!   --------                 interval iv in the reduced-dimension mode
!                            (in module spline_moments)
!----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

    ! .. local variables

    INTEGER :: i,j, ip,jp, m,left, ii,jj, ik
    REAL(KIND=8) :: xbase, c
    REAL(KIND=8), DIMENSION(ks) :: x,w, gx,gw
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp
    REAL(KIND=8), DIMENSION(ks,ks,ks) ::  Int
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
    REAL(KIND=8), DIMENSION(ks*(ks+1)/2,ks*(ks+1)/2) :: a

    left=iv+ks-1
    xbase=t(left)

    ! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
        Call vbsplvd(t,left,1,gx(i),1,dbiatx)
        bspTmp (i,:) = dbiatx(1,:,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:) * gx(:)**k

!            / r(iv,m)                             k
! .. Int =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      c = grw(iv,m) * grm(iv,m)**(k+3)
      DO j=1,ks
        gx(:) = gw(:)*bspTmp(:,j)
        DO jp=j,ks
          Int(j,jp,m) = SUM(gx(:)*bspTmp(:,jp)) * c
        END DO
      END DO

    END DO    ! over m

! .. second integration ..

    ii = 0
    DO i=1,ks
     DO ip=i,ks
      ii = ii + 1
      gx(:) = bsp(iv,:,i)*bsp(iv,:,ip)
      jj = 0
      DO j=1,ks
       DO jp=j,ks
         jj = jj + 1
         a(ii,jj) = SUM(gx(:)*INT(j,jp,:))
       END DO
      END DO
     END DO
    END DO

    ik = ks*(ks+1)/2
    rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)

    END SUBROUTINE mk_triang

    END SUBROUTINE mk_moments

