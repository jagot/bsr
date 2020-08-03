!=======================================================================
    SUBROUTINE dxv(k,n,b,v,y)
!=======================================================================
!   Computes   y = b * v    where b is a anti-symmetric, banded matrix,
!   in lower-band storage mode,  and v, y are vectors.
!-----------------------------------------------------------------------
!   on entry
!   --------
!       k       the number of diagonals
!       n       the order of the matrix
!       b       the symmetric, banded matrix in column storage mode
!       v       vector
!
!   on exit
!   -------
!       y       y = b*v
!-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, k
    REAL(8), INTENT(IN) :: v(n), b(n,k)
    REAL(8), INTENT(out) :: y(n)

! ... Local variables

    INTEGER :: i, j, jp

! ...   off diagonal elements only

        Do jp = 1,k-1;   Do i = k-jp+1,n;  j = i-k+jp

          y(i) = y(i) + b(i,jp)*v(j)             ! sub_diagonals
          y(j) = y(j) - b(i,jp)*v(i)             ! super-diagonals

        End do;  End do

     END SUBROUTINE dxv
