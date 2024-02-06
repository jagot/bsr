!======================================================================
      Double precision function  QUADR(i,j,m)
!======================================================================
!     Evaluates   <P_i | r^m | P_j>     with respect to r
!----------------------------------------------------------------------

      USE spline_param
      USE spline_galerkin
      USE spline_orbitals, p => pbs

      IMPLICIT NONE

      INTEGER, INTENT(in) :: i,j,m

      REAL(KIND=8), EXTERNAL :: BVMV
      REAL(KIND=8) :: rm(ns,ks)

      if     ( m .eq. 1 ) then
        quadr = BVMV (ns,ks, r1,'s',p(1,i),p(1,j))
      elseif ( m .eq. 0 ) then
        quadr = BVMV (ns,ks, sb,'s',p(1,i),p(1,j))
      elseif ( m .eq.-1 ) then
        quadr = BVMV (ns,ks,rm1,'s',p(1,i),p(1,j))
      elseif ( m .eq.-2 ) then
        quadr = BVMV (ns,ks,rm2,'s',p(1,i),p(1,j))
      else
        Call MRM(m,rm)
        quadr = BVMV(ns,ks,rm,'s',p(1,i),p(1,j))
      end if

      End function  QUADR
