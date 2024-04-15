!=========================================================================
      SUBROUTINE mtk_cell(k)
!=========================================================================
!
!     Defines matrix of Tk integrals in the B-spline basis
!     by cell algorithm
!
!     Calls: tk_moments
!
!----------------------------------------------------------------------
!
!     on entry      k        multipole index
!     --------
!
!     on exit       rkb     four-dimensional array of Tk integrals
!     -------               of power k in the B-spline basis
!                           (in module spline-integrals)
!-------------------------------------------------------------------------

      USE spline_param
      USE spline_atomic
      USE spline_integrals
      USE spline_moments

      IMPLICIT NONE
      INTEGER, INTENT(in) :: k

      ! .. local variables

      INTEGER(4) :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      REAL(8) :: c

      ! .. check the need of calculations

      if(itype == 'aaa') Call allocate_integrals
      if(itype == 'tk ' .and. krk == k) Return

      ! .. compute the moments in the spline basis

      CALL tk_moments(k)

      ! .. assemble the moments

      rkb = 0.d0

      DO jv=1,nv
       jj = 0
       DO jh=1,ks
        j = jv + jh - 1
        DO jhp=1,ks
         jp = jhp - jh + ks
         jj = jj + 1

         DO iv=1,nv
          ii = 0
          DO ih=1,ks
           i = iv + ih -1
           DO ihp=1,ks
            ip = ihp - ih + ks
            ii = ii + 1

            IF( iv < jv ) THEN
             c = rkd1(ii,iv)*rkd2(jj,jv)
            ELSE IF ( iv > jv ) THEN
             c = rkd1(jj,jv)*rkd2(ii,iv)
            ELSE
             c = rkd(ii,jj,iv)
            END IF

            rkb(i,j,ip,jp) = rkb(i,j,ip,jp) + c

           END DO
          END DO
         END DO

        END DO
       END DO
      END DO

      c = fine / (k+k+1)
      rkb = rkb * c

      itype = 'tk '
      krk = k

      END SUBROUTINE mtk_cell
