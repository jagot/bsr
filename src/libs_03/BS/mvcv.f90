!====================================================================
    SUBROUTINE mvcv(l,vc)
!====================================================================
!
!   Computes the matrix elements for the mass-velocity correction
!   in the B-spline basis.  The lack of symmetry in the d^2/dr^2
!   operator is ignored.
!
!     VC(i,j) = INT [  (d^2/dr^2 - l(l+1)/r) B_i(r) *
!                      (d^2/dr^2 - l(l+1)/r) B_j(r)  ] dr
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       l    the angular momentum
!
!   on exit
!   -------
!       vc   the mass velocity correction in symmetric storage mode
!
!--------------------------------------------------------------------

    USE spline_param; USE spline_atomic;  USE spline_grid

    IMPLICIT NONE

    INTEGER(4), INTENT(in) :: l
    REAL(8), INTENT(inout), Dimension(ns,ks) :: vc

    INTEGER(4) :: m, ith, jth, i, irow, jcol
    REAL(8) :: fll, y1, y2, S
    Real(8), External :: AZL

    ! .. initialize the vc array

    vc = 0.d0;  fll = l*(l+1);  nv = ns-ks+1

    ! .. compute the matrix elements

    do m = 1,ks
      do i = 1,nv
        S = fll*grm(i,m)*grm(i,m)

! ... cutoff correction

!        B = gr(i,m)/(gr(i,m)+2*fine*Z);  B = B*B*B

        do ith = 1,ks
          irow = i+ith-1
          do jth = 1,ith
          jcol = jth-ith+ks

            y1 = bspd(i,m,ith,2) - S*bsp(i,m,ith)
            y2 = bspd(i,m,jth,2) - S*bsp(i,m,jth)
            vc(irow,jcol) = vc(irow,jcol) + grw(i,m)*y1*y2 !* B

          end do
        end do
      end do
    end do

    vc = vc * fine

! ... one-electron Darwin correction:

    if(l.eq.0) then
      S = azl(z,h,ks,l+1);  vc(2,ks) = vc(2,ks) - z*S*S*fine
    end if


    END SUBROUTINE mvcv


