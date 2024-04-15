!======================================================================
      REAL(8) FUNCTION qk(i1,j1,i2,j2,k,meth)
!======================================================================
!               k
!     Returns  Q (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mqk_diff or mqk_cell, depending of the parameter 'meth':
!     meth = 'd' - differential equation method
!          = 'c' - cell integration method
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals
      USE spline_integrals

      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: i1,j1,i2,j2,k
      CHARACTER(1), INTENT(in), OPTIONAL :: meth

      ! .. local variables

      INTEGER(4) :: i,ip, j,jp, imin,imax, jmin,jmax
      REAL(8), DIMENSION(ns,ks+ks-1) :: a,b
      REAL(8) :: qki

      ! .. check the B-spline integrals in module spline-integrals

      if(k.ne.krk.or.itype.ne.'qk ') then
!       if(meth.eq.'d') then
!         Call Mqk_diff(k)
!       else
         Call Mqk_cell(k)
!       end if
      end if

      ! .. form cross-products

      Call density (ns,ks,a,pbs(1,i1),pbs(1,i2),'n')
      Call density (ns,ks,b,pbs(1,j1),pbs(1,j2),'n')

      ! .. assembling the B-spline integrals

      qk = 0.d0

      do jp = 1,ks+ks-1
       jmin=max( 1, 1 + ks-jp)
       jmax=min(ns,ns + ks-jp)
       do j = jmin,jmax
        qki = 0.d0
        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
           qki = qki + a(i,ip)*rkb(i,j,ip,jp)
          end do
        end do
        qk = qk + b(j,jp)*qki
       end do
      end do

      END FUNCTION qk
