!=========================================================================
      Subroutine bxv(k,n,b,v,y)
!=========================================================================
!     Computes   y = b * v    where b is a symmetric, banded matrix,
!     in lower-band storage mode,  and v, y are vectors.
!-------------------------------------------------------------------------
!     on entry
!     --------
!         k       the number of diagonals
!         n       the order of the matrix
!         b       the symmetric, banded matrix in column storage mode
!         v       vector
!
!     on exit
!     -------
!         y       y = b*v
!-------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n, k
      Real(8), intent(in) :: v(n), b(n,k)
      Real(8), intent(out) :: y(n)

      ! .. Local variables

      Integer :: i, j, jp

! ... contribution from central diagonal (jp=k)

      Do i=1,n
        y(i) = b(i,k)*v(i)
      End do

! ... off diagonal

      Do jp = 1,k-1
       Do i = k-jp+1,n
        j = i-k+jp
        y(i) = y(i) + b(i,jp)*v(j)             ! sub_diagonals
        y(j) = y(j) + b(i,jp)*v(i)             ! super-diagonals
       End do
      End do

      End Subroutine bxv
