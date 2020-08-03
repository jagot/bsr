!=========================================================================
    SUBROUTINE tk_moments(k)
!=========================================================================
!
!   Defines moments for Tk-integrals in the B-spline cells
!
!   Calling sequence:          tk_moments             
!                              ----------             
!                               /    \\           
!                           moments tk_pdiag      
!                                     ||          
!                                   tk_triang     
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
    USE spline_galerkin
    USE spline_moments

    IMPLICIT NONE
    INTEGER(4), INTENT(in) :: k  

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_moments
    if(mtype == 'tk ' .and. kmk == k) Return

    ! .. compute the moments in the spline basis

    CALL moments(  k   , rkd1,'n','d')
    CALL moments(-(k+1), rkd2,'n','d')

    CALL tk_pdiag

    mtype = 'tk '
    kmk = k

    CONTAINS


!======================================================================
    SUBROUTINE tk_pdiag
!======================================================================
!
!   Controls the scaling propeties for diagonal B-spline Tk-interals
!
!   Calls:  tk_triang
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv, jk
    REAL(KIND=8) :: hp

    ! .. non-exponential grid

    if(me.eq.0) then
     DO iv=1,nv
      CALL tk_triang(iv)
     END DO
     Return
    end if

    ! .. the first equal step region

    Do iv=1,ml+ks-1
      CALL tk_triang(iv)
    End do

    ! .. the exponential region - using scaling law

    hp=h+1.d0
    jk = ks*ks
    Do iv=ml+ks,ml+me-ks+2
      rkd(1:jk,1:jk,iv) = rkd(1:jk,1:jk,iv-1) / hp
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL tk_triang(iv)
    END DO

    END SUBROUTINE tk_pdiag



!========================================================================
    SUBROUTINE tk_triang(iv)
!========================================================================
!
!   Returns the two-dimensional array of B-splin Tk-integrals
!                                  _    _
!         <B_i B_j| r2^k/r1^(k+1) |B_i' B_j'>
!
!   over the given triangle diagonal cell 
!                                                           
!   Calls:   gauss, vbsplvd
!
!   On entry   iv  -  the index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of B-spline Tk-integrals for given 
!   --------                 interval iv in the reduced-dimension mode
!                            (in module spline_moments)      
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: iv

! .. local variables

      INTEGER :: i,j, ii,jj, ip,jp, m, left, jk
      REAL(KIND=8) :: xbase
      REAL(KIND=8), DIMENSION(ks) :: x,w, gx,gw
      REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp, bspdTmp
      REAL(KIND=8), DIMENSION(ks,ks,ks) :: INT
      REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
      REAL(KIND=8), DIMENSION(ks*ks,ks*ks) :: a

      left=iv+ks-1
      xbase=t(left)

! .. setup the gaussian points

      CALL gauss(ks,x,w)

! .. first integration 

      DO m=1,ks

! .. the absolute coordinate at the new gaussian point

       gx(1:ks) = (gr(iv,m)-xbase)*x(1:ks) + xbase

! .. the bspline values at the new gaussian points

       DO i=1,ks
        Call vbsplvd(t,left,1,gx(i),2,dbiatx)
        bspTmp (i,1:ks) = dbiatx(1,1:ks,1)
        bspdTmp(i,1:ks) = dbiatx(1,1:ks,2)-bspTmp(i,1:ks)/gx(i)
      END DO

! .. and the corresponding gaussian weights

       gw(1:ks) = (gr(iv,m)-xbase)*w(1:ks)
       gx(1:ks) = gx(1:ks)**k * gw(1:ks)

!            / r(iv,m)         ___           k
! .. INT =  |      bsp(iv,:,j) bsp(iv,:,jp) r  dr
!           / r_iv

       DO j=1,ks
        gw(1:ks) = gx(1:ks) * bspTmp(1:ks,j)
        DO jp=1,ks
          INT(j,jp,m)= SUM(gw(1:ks) * bspdTmp(1:ks,jp))
        END DO
       END DO

      END DO    ! over m

! .. second integration

      gx(1:ks) = grw(iv,1:ks)*grm(iv,1:ks)**(k+1)

      Do ip=1,ks
       bspTmp(1:ks,ip) = bsq(iv,1:ks,ip)*gx(1:ks)
      End do

      ii = 0
      DO i=1,ks
       DO ip=1,ks
        ii = ii + 1
        gx(1:ks) =  bsp(iv,1:ks,i)*bspTmp(1:ks,ip)
        jj = 0
        DO j=1,ks
         DO jp=1,ks
          jj = jj + 1
          a(ii,jj) = SUM(gx(1:ks)*INT(j,jp,1:ks))
         END DO
        END DO
 
       END DO
      END DO

      jk = ks*ks
      rkd(1:jk,1:jk,iv) = a + Transpose(a)

      END SUBROUTINE tk_triang

      END SUBROUTINE tk_moments
