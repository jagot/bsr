!====================================================================
      MODULE spline_hl
!====================================================================
!     contains the spline representation of L operator
!--------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: lh = -1    ! l-value for current L-operator
    
      REAL(8), ALLOCATABLE :: hl(:,:), vc(:,:)

! ... hl(1:ns,1:ks) - matrix of L-operator in the B-spline basis
!                    (in almost symmetric lower-column mode)
! ... vc(1:ns,1:ks) - matrix of mass-velocity correction

      END MODULE spline_hl

!====================================================================
      SUBROUTINE allocate_hl
!====================================================================
      USE spline_param
      USE spline_hl

      if(allocated(hl)) Deallocate(hl,vc)
      ALLOCATE( hl(ns,ks), vc(ns,ks))
	  
      hl = 0.d0; vc = 0.d0; lh = -2

      END SUBROUTINE allocate_hl






