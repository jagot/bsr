!====================================================================
      SUBROUTINE facsbl(l)
!====================================================================
!
!     Sets up the overlap bs which is a transpose of sb, <B_i,B_j>,
!     with the correct boundary conditions, and then factorizes bs.
!
!     SUBROUTINES called:   DPBTRF (from LAPACK)
!
!--------------------------------------------------------------------
    
      USE spline_param; USE spline_galerkin

      IMPLICIT NONE

      INTEGER(4), intent(in) :: l

      INTEGER(4) :: i,j, ierr

      bs = TRANSPOSE(sb)

! ... apply zero boundary condition at r=0:

      Do i = 1,l+1
       bs(ks,i) = 1.d0
       Do j = 1,ks-1
        bs(j,ks-j+i) = 0.d0
       End do
      End do

! ... apply zero boundary condition at r=a:

      Do i = 1,ks-1
       bs(i,ns)=0.d0
       bs(i,ns-1)=0.d0
      End do
      bs(ks,ns  ) = 1.d0
      bs(ks,ns-1) = 1.d0
 
      Call DPBTRF('U',ns,ks-1,bs,ks,ierr)
      if (ierr .ne. 0 )  Stop 'facsbl: dpbtrf (LAPACK) failed'

      END SUBROUTINE facsbl



