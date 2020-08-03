!=========================================================================
      Subroutine coulom(z, l, a, b)
!=========================================================================
!     Builds full a and b matrices for the coulomb problem
!  
!     a(i,j) = integration of B_i*(DD+2*z/r-l(l+1)/r**2)*B_j
!     b(i,j) = integration of B_i*B_j
!  
!     from the banded symmetric form of these elementary operators
!-------------------------------------------------------------------------
!     on exit
!     -------
!         a      full matrix for coulomb operator not assumed
!                to be symmetric
!         b      full matrix of the overlap operator
!---------------------------------------------------------------------
      Use spline_param
      Use spline_galerkin
      Use spline_hl
   
      Implicit none
      Real(8), intent(in) :: z
      Integer, intent(in) :: l
      REAL(8), intent(out) :: a(ns,ns), b(ns,ns)
   
      ! ..  LOCAL variables
   
      Integer :: i,j
   
      ! .. initialize
   
      a = 0.d0
      b = 0.d0
   
      Call HLM(l)
   
      ! .. lower matrix
   
      Do j = 1,ks
        Do i = ks-j+1,ns
          a(i,i-ks+j) = hl(i,j)
          b(i,i-ks+j) = sb(i,j)
        End do
      End do
   
      ! .. upper matrix by symmetry
   
      Do j = ks+1,2*ks-1
        Do i = j+1-ks,ns
          a(i-j+ks,i) = hl(i,2*ks-j)
          b(i-j+ks,i) = sb(i,2*ks-j)
        End do
      End do
   
      ! .. correction for asymmetry in db2
   
      a(ns-1,ns) = hl(1,1)
   
      End Subroutine coulom
   