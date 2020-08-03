!=========================================================================
      SUBROUTINE mvk_cell(k)
!=========================================================================
!
!     Defines matrix of Vk integrals in the B-spline basis
!     by cell algorithm
!
!     Calls: vk_moments
!----------------------------------------------------------------------
!
!     on entry      k        multipole index
!     --------
!       
!     on exit       rkb     four-dimensional array of Vk integrals 
!     -------               of power k in the B-spline basis
!                           (in module spline-integrals)
!-------------------------------------------------------------------------

      USE spline_param
      USE spline_integrals
      USE spline_moments
      USE spline_atomic
    
      IMPLICIT NONE
      INTEGER, INTENT(in) :: k
    
      ! .. local variables
    
      INTEGER(4) :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      REAL(8) :: c
    
      ! .. check the need of calculations
    
      if(itype == 'aaa') Call allocate_integrals
      if(itype == 'vk ' .and. krk == k) Return
    
      ! .. compute the moments in the spline basis
    
      CALL vk_moments(k)
    
      ! .. assemble the moments
    
      rkb = 0.d0
    
      DO jv=1,nv
       jj = 0
       DO jh=1,ks
        j = jv + jh - 1
        DO jhp=jh,ks
         jp = jhp - jh + 1
         jj = jj + 1
    
         DO iv=1,nv
          ii = 0
          DO ih=1,ks
           i = iv + ih - 1
           DO ihp=1,ks
            ip = ihp - ih + ks
            ii = ii + 1
    
            IF( iv < jv ) THEN
             c  =  rkd3(ii,iv)*rkd2(jj,jv)
            ELSE IF( iv > jv ) THEN
             c  =  rkd1(jj,jv)*rkd4(ii,iv)
            ELSE
             c  =  rkd(ii,jj,iv)
            END IF
    
            rkb(i,j,ip,jp) = rkb(i,j,ip,jp) + c
                
           END DO
          END DO
         END DO

        END DO
       END DO
      END DO
    
      rkb = rkb * fine
    
      itype = 'vk '
      krk = k
    
      END SUBROUTINE mvk_cell
