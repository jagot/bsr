!====================================================================
    SUBROUTINE facsb
!====================================================================
!
!   Factorizes bs matrix which is a transpose of overlap matrix sb,
!   <B_i,B_j>,  with the correct boundary condition at r=0
!
!   SUBROUTINES called:  dpbfa (from LINPACK)
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin

    IMPLICIT NONE
    INTEGER :: m, ierr

    ! .. copy the array, converting to row oriented band storage mode

    bs = TRANSPOSE(sb)

    ! .. apply boundary condition at r=0

    do m = 1,ks-1
      bs(m,ks-m+1)=0.d0
    end do
    bs(ks,1) = 1.d0

!    CALL dpbfa(bs,ks,ns,ks-1,ierr)
!    if (ierr /= 0 ) STOP 'facsb: dpbfa (LINPACK) failed'

    Call DPBTRF('U',ns,ks-1,bs,ks,ierr)
    if (ierr.ne.0)  Stop 'facsb: dpbtrf (LAPACK) failed'

  END SUBROUTINE facsb
