!====================================================================
    MODULE spline_atomic
!====================================================================
!   atomic parameters related to splines
!--------------------------------------------------------------------
    IMPLICIT NONE
    REAL(8) :: z  = 1.d0        ! nuclear charge
    REAL(8) :: EC = 0.d0        ! core energy
    REAL(8) :: fine = 0.25D0/(137.036D0)**2
    LOGICAL :: rel = .FALSE.    ! relativistic corrections
    INTEGER :: irel  =  0       ! the same
    INTEGER :: kclosd = 0       ! closed shells
    INTEGER :: MASS =   0       ! mass-corrections
    INTEGER :: ioo =    0       ! orbit-orbit interaction
    END MODULE spline_atomic

