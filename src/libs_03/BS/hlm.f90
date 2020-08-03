!====================================================================
    SUBROUTINE hlm(l)
!====================================================================
!   Sets up the matrix hl for the coulomb operator
!
!       hl(i,j) = DD(i,j) - l*(l+1)*rm2(i,j) + 2 * z*rm1(i,j)
!
!   in symmetric band storage mode, where i=1,ns, j=1,ks and
!   hl(1,1) stores <B_{ns-1},H B_{ns}>.
!--------------------------------------------------------------------
!   on entry
!   --------
!       z    nuclear charge
!       l    orbital angular momentum
!--------------------------------------------------------------------
    USE spline_param
    Use spline_atomic
    USE spline_grid
    Use spline_galerkin
    Use spline_hl

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l
    REAL(8) :: fl, zz, a
    Real(8), external :: AZL

    if(lh.eq.-1) Call Allocate_hl
    if(lh.eq.l) Return

    fl = l*(l+1.d0)
    zz = 2.d0*z

!.. set up hl matrix

    hl = db2 - fl*rm2 + zz*rm1

!.. store the (ns-1,ns) asymmetric value in hl(1,1)

    hl(1,1) = hl(ns,ks-1) + (db2(1,1)-db2(ns,ks-1))

!.. relativistic shift

    if(rel) then
     Call mvc(l)
     hl = hl + fine*vc
     if(l.eq.0) then
      a = azl(z,h,ks,l+1)
      hl(2,ks) = hl(2,ks) - z*a*a*fine
     end if
    end if

    lh = l

    END SUBROUTINE hlm
