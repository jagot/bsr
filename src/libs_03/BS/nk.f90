!======================================================================
      REAL(8) FUNCTION nk (i1,j1,i2,j2,k,meth)
!======================================================================
!               k
!     Returns  N (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mnk_diff or mnk_cell, depending of the parameter 'meth':
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
      REAL(8) :: nkj
  
      ! .. check the B-spline integrals in module spline-integrals
  
      if(k.ne.krk.or.itype.ne.'nk ') then
       if(meth.eq.'d') then
         Call MNK_diff(k)
       else
         Call MNK_cell(k)
       end if
      end if
  
      ! .. form cross-products
  
      Call density (ns,ks,a,pbs(1,i1),pbs(1,i2),'s')
      Call density (ns,ks,b,pbs(1,j1),pbs(1,j2),'s')
  
      ! .. assembling the B-spline integrals

      nk = 0.d0
      do ip = 1,ks
        do i = 1,ns-ip+1
          nkj = 0.d0
          do jp = 1,ks
            do j = 1,ns-jp+1
              nkj = nkj+b(j,jp)*rkb(i,j,ip,jp)
            end do
          end do
          nk = nk + a(i,ip)*nkj
        end do
      end do
  
      END FUNCTION nk
  