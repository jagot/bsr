!=========================================================================
    SUBROUTINE qk_moments(k)
!=========================================================================
!
!   Defines moments for Qk-integrals in the B-spline cells
!
!   Calling sequence:          qk_moments
!                              ----------
!                               /    \\
!                           moments qk_pdiag
!                                     ||
!                                   qk_triang
!                                    /   \
!                                 gauss  qbsplvb
!
!-------------------------------------------------------------------------
!
!   on entry
!   --------
!       k     the multipole index
!
!   on exit
!   -------
!      rkd    four-dimensional array of qk integrals in diagonal cells
!      rkd1,rkd2,rkd3,rkd4 - non-diagonal shells
!
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_moments

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_integrals
    if(mtype == 'qk ' .and. kmk == k) Return

    ! .. compute the moments in the spline basis

    if(mtype == 'aaa') Call allocate_moments

    CALL moments(  k   , rkd3, 'n','d')
    CALL moments(  k-1 , rkd2, 'n','b')
    rkd3 = 2 * rkd3 + (k+2)*rkd2
    CALL moments(-(k+2), rkd2, 'n','b')

    CALL moments(-(k+3), rkd4, 'n','d')
    CALL moments(-(k+4), rkd1, 'n','b')
    rkd4 = 2 * rkd4 - (k+1)*rkd1
    CALL moments(  k+1 , rkd1, 'n','b')

    CALL qk_pdiag

    mtype = 'qk '
    kmk = k

!-----------------------------------------------------------------------
    CONTAINS
!-----------------------------------------------------------------------


!======================================================================
    SUBROUTINE qk_pdiag
!======================================================================
!
!   Controls the scaling propeties for diagonal B-spline Qk-interals
!
!   Calls:  qk_triang
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER(4) :: iv,ik,jk
    REAL(8) :: hp

    ! .. non-exponential grid

    if(me.eq.0) then
     DO iv=1,nv
      CALL qk_triang(iv)
     END DO
     Return
    end if

    ! .. the first equal step region.

    Do iv=1,ml+ks-1
      CALL qk_triang(iv)
    End do

    ! .. the log region --- using scaling law.

    hp=h+1.d0
    Do iv=ml+ks,ml+me-ks+2
      rkd(:,:,iv) = rkd(:,:,iv-1) / hp
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL qk_triang(iv)
    END DO

    END SUBROUTINE qk_pdiag

!========================================================================
    SUBROUTINE qk_triang(iv)
!========================================================================
!
!    Returns the "qk matrix element" in the diagonal cell   iv
!
!------------------------------------------------------------------------
!
!   SUBROUTINES called:
!       gauss
!       bsplvd
!
!---------------------------------------------------------------------
!
!   On entry
!   --------
!       k:     the indices of of the bsplines
!       iv:    the index of the integration region
!
!---------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

! .. local variables

    INTEGER :: i,j, ip,jp, ii,jj, m, left
    REAL(KIND=8) :: xbase, c
    REAL(KIND=8), DIMENSION(ks) :: x,w,gx,gw,gv, bi
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp,bspdTmp
    REAL(KIND=8), DIMENSION(ks,ks,ks) ::Int1,Int2
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx

    left = iv+ks-1
    xbase = t(left)

! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks                  !  loop over old knots

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call vbsplvd(t,left,1,gx(i),2,dbiatx)
       bspTmp (i,1:ks) = dbiatx(1,1:ks,1)
       bspdTmp(i,1:ks) = 2*dbiatx(1,1:ks,2) + k*bspTmp(i,1:ks)/gx(i)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:) * gx(:)**k


!              / r(iv,m)         ___           k
! .. Int1  =  |      bsp(iv,:,i) bsp(iv,:,ip) r  dr
!             / r_iv


      c = grm(iv,m)**(k+2) * grw(iv,m)
      DO i=1,ks
       bi(:) = gw(:)*bspTmp(:,i)
       DO ip=1,ks
        Int1(i,ip,m)= SUM(bi(:)*bspdTmp(:,ip)) * c
       END DO
      END DO

!              / r(iv,m)                       k+1
! .. Int2  =  |      bsp(iv,:,j) bsp(iv,:,jp) r    dr
!             / r_iv

      gv = gw * gx
      c = grm(iv,m)**(k+3) * grw(iv,m)
      DO j=1,ks
       bi(:) = gv(:)*bspTmp(:,j)
       DO jp=1,ks
        Int2(j,jp,m)= SUM(bi(:)*bspTmp(:,jp)) * c
       END DO
      END DO

    END DO	!  over m

! ... second integration

      jj = 0
      DO j=1,ks
       DO jp=1,ks
        jj = jj + 1
        gx (:) = bsp(iv,:,j)*bsp(iv,:,jp)
        ii = 0
        DO i=1,ks
         DO ip=1,ks
          ii = ii + 1
          rkd(ii,jj,iv) = SUM( gx(:)*Int1(i,ip,:))
         END DO
        END DO
       END DO
      END DO

      Do ip=1,ks
       bspTmp(:,ip) = 2*bspd(iv,:,ip,1) - (k+3)*grm(iv,:)*bsp(iv,:,ip)
      End do

      ii = 0
      DO i=1,ks
       DO ip=1,ks
        ii = ii + 1
        gx(:) = bsp(iv,:,i)*bspTmp(:,ip)
        jj = 0
        DO j=1,ks
         DO jp=1,ks
          jj = jj + 1
          rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM( gx(:)*Int2(j,jp,:))
         END DO
        END DO
       END DO
      END DO

    END SUBROUTINE qk_triang

    END SUBROUTINE qk_moments
