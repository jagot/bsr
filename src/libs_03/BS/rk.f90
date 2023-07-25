!======================================================================
      REAL(8) FUNCTION rk (i1,j1,i2,j2,k,meth)
!======================================================================
!               k
!     Returns  R (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mrk_diff or mrk_cell, depending of the parameter 'meth':
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
  
      INTEGER(4) :: i,ip, j,jp
      REAL(8), DIMENSION(ns,ks) :: a,b
      REAL(8) :: rkj
  
      ! .. check the B-spline integrals in module spline-integrals
  
      if(k.ne.krk.or.itype.ne.'rk ') then
       if(meth.eq.'d') then
         Call MRK_diff(k)
       else
         Call MRK_cell(k)
       end if
      end if
  
      ! .. form cross-products
  
      Call density (ns,ks,a,pbs(1,i1),pbs(1,i2),'s')
      Call density (ns,ks,b,pbs(1,j1),pbs(1,j2),'s')
  
      ! .. assembling the B-spline integrals

      rk = 0.d0
      do ip = 1,ks
        do i = 1,ns-ip+1
          rkj = 0.d0
          do jp = 1,ks
            do j = 1,ns-jp+1
              rkj = rkj+b(j,jp)*rkb(j,i,jp,ip)
            end do
          end do
          rk = rk + a(i,ip)*rkj
        end do
      end do
  
      END FUNCTION rk
  