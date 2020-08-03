!=========================================================================
    SUBROUTINE vk_moments(k)
!=========================================================================
!
!   Defines moments for Vk-integrals in the B-spline cells
!
!   Calling sequence:          vk_moments             
!                              ----------             
!                               /    \\           
!                           moments vk_pdiag      
!                                     ||          
!                                   vk_triang     
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
    INTEGER, INTENT(in) :: k

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_moments
    if(mtype == 'vk ' .and. kmk == k) Return

    ! .. compute the moments in the spline basis

    CALL moments(  k+1 , rkd1,'s','b')
    CALL moments(-(k+2), rkd2,'s','b')
    CALL moments(  k   , rkd3,'n','d')
    CALL moments(-(k+3), rkd4,'n','d')
    CALL vk_pdiag

    mtype = 'vk '
    kmk = k

    CONTAINS


!======================================================================
    SUBROUTINE vk_pdiag
!======================================================================
!
!   Controls the scaling propeties for diagonal B-spline Vk-interals
!
!   Calls:  vk_triang
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER(4) :: iv,ik,jk
    REAL(8) :: hp

    ! .. non-exponential grid

    if(me.eq.0) then
     DO iv=1,nv
      CALL vk_triang(iv)
     END DO
     Return
    end if

    ! .. the first equal step region.

    Do iv=1,ml+ks-1
      CALL vk_triang(iv)
    End do

    ! .. the log region --- using scaling law.

    jk = ks*ks
    ik = ks*(ks+1)/2
    hp=h+1.d0
    Do iv=ml+ks,ml+me-ks+2
      rkd(1:jk,1:ik,iv) = rkd(1:jk,1:ik,iv-1) / hp
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL vk_triang(iv)
    END DO

    END SUBROUTINE vk_pdiag


!========================================================================
    SUBROUTINE vk_triang(iv)
!========================================================================
!
!   Returns the two-dimensional array of B-splin Vk-integrals
!                                  _    
!         <B_i B_j| (r<)^k/(r>)^(k+1) |B_i' B_j' r2 >
!
!   over the given triangle diagonal cell 
!                                                           
!   Calls:   gauss, vbsplvd
!
!   On entry   iv  -  the index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of B-spline Vk-integrals for given 
!   --------                 interval iv in the reduced-dimension mode
!                            (in module spline_moments)      
!----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

! .. local variables

    INTEGER(4) :: i,j, ip,jp, m, left, ii,jj
    REAL(8) :: xbase, c
    REAL(8), DIMENSION(ks) :: x,w,gx,gw, bi
    REAL(8), DIMENSION(ks,ks) :: bspTmp,bspdTmp
    REAL(8), DIMENSION(ks,ks,ks) ::Int1,Int2
    REAL(8), DIMENSION(nv,ks,ks) :: dbiatx	          

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
      bspdTmp(i,1:ks) = dbiatx(1,1:ks,2) - bspTmp(i,1:ks)/gx(i)
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

      gw = gw * gx
      c = grm(iv,m)**(k+3) * grw(iv,m)
      DO j=1,ks
       bi(:) = gw(:)*bspTmp(:,j)
       DO jp=1,ks
        Int2(j,jp,m)= SUM(bi(:)*bspTmp(:,jp)) * c
       END DO
      END DO

     END DO	!  over m

! .. second integration

     jj = 0
     DO j=1,ks
      DO jp=j,ks
       jj = jj + 1
       gx(:) = bsp(iv,:,j)*bsp(iv,:,jp)
       ii = 0
       DO i=1,ks
        DO ip=1,ks
         ii = ii + 1
         rkd(ii,jj,iv) = SUM(gx(:)*Int1(i,ip,:))
        END DO
       END DO
      END DO
     END DO

     ii = 0
     DO i=1,ks
      DO ip=1,ks
       ii = ii + 1
       gx(:) = bsp(iv,:,i)*bsq(iv,:,ip)
       jj = 0
       DO j=1,ks
        DO jp=j,ks
         jj = jj + 1
         rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM(gx(:)*Int2(j,jp,:)) 
       END DO
      END DO
     END DO
    END DO

    END SUBROUTINE vk_triang

    END SUBROUTINE vk_moments
